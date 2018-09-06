-- The neural network for WLNM
-- *author: Muhan Zhang, Washington University in St. Louis

require'paths'
require 'torch'
require 'nn'
require 'cunn'
require 'cutorch'
require 'optim'
require 'svm'
require 'nnsparse'


cmd = torch.CmdLine()
cmd:text()
cmd:text()
cmd:text('Train a neural network for link prediction.')
cmd:text()
cmd:text('Options')
cmd:option('-inputdim', 190, 'first layer input dimension, should be equal to K(K-1)/2')
cmd:option('-ith_experiment', 1, 'the ID of the current experiment')
cmd:text()


-- parse input params
params = cmd:parse(arg)
ith_experiment = params.ith_experiment


num_gpus = cutorch.getDeviceCount()
cutorch.setDevice(ith_experiment % num_gpus + 1)
--cutorch.setDevice((ith_experiment-1) % 3 + 1)

-- network configuration
l1 = params.inputdim
l2 = 32
l3 = 32
l4 = 16

net = nn.Sequential()
net:add(nn.Linear(l1, l2))
--net:add(nn.BatchNormalization(l2))
net:add(nn.ReLU())
--net:add(nn.Dropout(0.5))
net:add(nn.Linear(l2, l3))
--net:add(nn.BatchNormalization(l3))
net:add(nn.ReLU())
--net:add(nn.Dropout(0.5))
net:add(nn.Linear(l3, l4))
--net:add(nn.BatchNormalization(l4))
net:add(nn.ReLU())
--net:add(nn.Dropout(0.5))
net:add(nn.Linear(l4, 2))
net:add(nn.LogSoftMax())


-- load dataset
trainname = paths.concat(paths.cwd(), string.format('tempdata/traindata_%d', ith_experiment))
testname = paths.concat(paths.cwd(), string.format('tempdata/testdata_%d', ith_experiment))
trainvaldata = svm.ascread(trainname)
testdata = svm.ascread(testname)


-- function for converting libsvm format to table
function libsvm2tensor(data)
   local n = table.getn(data)
   local set = {
      label = torch.Tensor(n), 
      size = function() return n end,
      data = torch.Tensor(n, l1)
   }
   for i = 1, n do
      set.label[i] = data[i][1] + 1  -- positive = 2, negative = 1
      local tmp = torch.cat(data[i][2][1]:double(), data[i][2][2]:double(), 2):densify(0, l1)
      set.data[i] = tmp
   end
   return set
end

-- split traindata into train and validation
local ntv = table.getn(trainvaldata)  -- number of train and val
local shuffle = torch.randperm(ntv)
traindata = {}
valdata = {}
local valratio = 0.1
local valnum = torch.floor(valratio * ntv)
for i = 1, valnum do
   valdata[i] = trainvaldata[shuffle[i]]
end
for i = valnum + 1, ntv do
   traindata[i - valnum] = trainvaldata[shuffle[i]]
end


-- convert libsvm format to tensor
trainset = libsvm2tensor(traindata)
valset = libsvm2tensor(valdata)
testset = libsvm2tensor(testdata)


-- prepare trainset for function 'train'
setmetatable(trainset, 
  {__index = function(t, i) return {t.data[i], t.label[i]} end}
)


runthis = false
if runthis then
mean = {} -- store the mean, to normalize the test set in the future
stdv  = {} -- store the standard-deviation for the future
for i=1,3 do -- over each image channel
    mean[i] = trainset.data[{ {}, {i}, {}, {}  }]:mean() -- mean estimation
    print('Channel ' .. i .. ', Mean: ' .. mean[i])
    trainset.data[{ {}, {i}, {}, {}  }]:add(-mean[i]) -- mean subtraction
    
    stdv[i] = trainset.data[{ {}, {i}, {}, {}  }]:std() -- std estimation
    print('Channel ' .. i .. ', Standard Deviation: ' .. stdv[i])
    trainset.data[{ {}, {i}, {}, {}  }]:div(stdv[i]) -- std scaling
end
end

net = net:cuda()

criterion = nn.ClassNLLCriterion()
criterion = criterion:cuda()
trainset.data = trainset.data:cuda()
trainset.label = trainset.label:cuda()
valset.data = valset.data:cuda()
valset.label = valset.label:cuda()
testset.data = testset.data:cuda()
testset.label = testset.label:cuda()


-- retrieve parameters and gradients
parameters,gradParameters = net:getParameters()


-- training function
function train(dataset)
   net:training()
   -- epoch tracker
   epoch = epoch or 1

   -- local vars
   local time = sys.clock()
   local trainError = 0

   -- do one epoch
   print('<trainer> on training set:')
   print("<trainer> online epoch # " .. epoch .. ' [batchSize = ' .. opt.batchSize .. ']')
   shuffle = torch.randperm(dataset:size())
   batchNum = torch.floor(dataset:size() / opt.batchSize)
   for batch = 1, batchNum do
      -- disp progress
      xlua.progress((batch - 1) * opt.batchSize, dataset:size())

      -- create mini batch
      local inputs = torch.zeros(opt.batchSize, l1):cuda()
      local targets = torch.zeros(opt.batchSize):cuda()
      local batchCount = 0
      for i = (batch - 1) * opt.batchSize + 1, batch * opt.batchSize do
         batchCount = batchCount + 1
         local input = dataset.data[shuffle[i]]
         local target = dataset.label[shuffle[i]]
         inputs[batchCount] = input
         targets[batchCount] = target
      end

      -- create closure to evaluate f(X) and df/dX
      local feval = function(x)
         -- get new parameters
         if x ~= parameters then
            parameters:copy(x)
         end
         local output = net:forward(inputs)
         tmp, outputLabels = torch.max(output, 2)
         local f = criterion:forward(output, targets)
         -- reset gradients
         gradParameters:zero()
         local df_do = criterion:backward(output, targets)
         net:backward(inputs, df_do)
         confusion:batchAdd(output, targets)
         -- normalize gradients and f(X)
         trainError = trainError + f * batchCount
         --print(torch.norm(opt.learningRate*gradParameters))
         return f,gradParameters
      end

      -- optimize on current mini-batch
      if opt.optimization == 'CG' then
         Config = Config or {maxIter = opt.maxIter}
         optim.cg(feval, parameters, Config)

      elseif opt.optimization == 'LBFGS' then
         Config = Config or {learningRate = opt.learningRate,
                             maxIter = opt.maxIter,
                             nCorrection = 10}
         optim.lbfgs(feval, parameters, Config)

      elseif opt.optimization == 'SGD' then
         Config = Config or {learningRate = opt.learningRate,
                             weightDecay = opt.weightDecay,
                             momentum = opt.momentum,
                             learningRateDecay = opt.learningRateDecay
                             }
         optim.sgd(feval, parameters, Config)
      elseif opt.optimization == 'ASGD' then
         Config = Config or {eta0 = opt.learningRate,
                             t0 = nbTrainingPatches * opt.t0}
         _,_,average = optim.asgd(feval, parameters, Config)

      elseif opt.optimization == 'ADAM' then
         Config = Config or {learningRate = opt.learningRate,  
                             weightDecay = opt.weightDecay
       }
         optim.adam(feval, parameters, Config)
      else
         error('unknown optimization method')
      end
   end

   -- train error
   trainError = trainError / math.floor(dataset:size())
   print(trainError)

   -- time taken
   time = sys.clock() - time
   time = time / dataset:size()
   --print(time)
   --print("<trainer> time to learn 1 sample = " .. (time*1000) .. 'ms')

   -- print confusion matrix
   print(confusion)
   local trainAccuracy = confusion.totalValid * 100
   confusion:zero()

   -- next epoch
   epoch = epoch + 1

   return trainAccuracy, trainError
end


-- test function
function test(dataset, writeflag)
   net:evaluate()
   writeflag = writeflag or false
   -- local vars
   local testError = 0
   local time = sys.clock()
   local testBatch = 32

   -- averaged param use?
   if average then
      cachedparams = parameters:clone()
      parameters:copy(average)
   end

   -- test over given dataset
   print('<trainer> on testing Set:')
   -- disp progress

   scores = torch.Tensor(dataset:size()):cuda()
   for t = 1, torch.ceil(dataset:size()/testBatch) do
     -- get new sample
     local input = dataset.data[{ {(t-1)*testBatch+1, math.min(t*testBatch, dataset:size())}, {}}]
     local target = dataset.label[{ {(t-1)*testBatch+1, math.min(t*testBatch, dataset:size())} }]
     -- test sample
     local pred = net:forward(input)
     scores[{ {(t-1)*testBatch+1, math.min(t*testBatch, dataset:size())}}] = pred[{{}, {2}}]
     confusion:batchAdd(pred, target)
     -- compute error
     err = criterion:forward(pred, target)
     testError = testError + err
     -- timing
     time = sys.clock() - time
     time = time / dataset:size()
     --print("<trainer> time to test a batch = " .. (time/1000) .. 'ms')
   end
   -- testing error estimation
   testError = testError / dataset:size()
   -- print confusion matrix
   print(confusion)
   local testAccuracy = confusion.totalValid * 100
   confusion:zero()

   -- write predictions to file, also print auc
   if writeflag == true then
      file = torch.DiskFile(string.format('tempdata/test_log_scores_%d.asc', ith_experiment), 'w')
      for ex = 1, dataset:size() do
         file:writeDouble(scores[ex])
      end
      -- print auc
      metrics = require 'metrics'
      roc_points, thresholds = metrics.roc.points(scores, dataset.label, 1, 2)
      auc = metrics.roc.area(roc_points)
      print(auc)
   end
   -- averaged param use?
   if average then
      -- restore parameters
      parameters:copy(cachedparams)
   end

   return testAccuracy, testError
 end

dname, fname = sys.fpath()
opt = {batchSize = 128, optimization = 'ADAM', learningRate = 0.001, weightDecay = 1e-3, momentum = 0.9, learningRateDecay = 1e-5, save2 = fname:gsub('.lua', ''), save = string.format('tempdata/nDNN_%d', ith_experiment)}

-- this matrix records the current confusion across classes
confusion = optim.ConfusionMatrix({'neg link', 'pos link'})

-- log results to files
accLogger = optim.Logger(paths.concat(opt.save, 'accuracy.log'))
errLogger = optim.Logger(paths.concat(opt.save, 'error.log'   ))

bestTestAcc = 0
maxIter = 100
for iter = 1, maxIter do
   -- train/test
   trainAcc, trainErr = train(trainset)
   testAcc,  testErr  = test(valset)
   
   if testAcc > bestTestAcc then
      bestTestAcc = testAcc
      -- save/log current net
      filename = paths.concat(opt.save, string.format('bestNet_%d.t7', ith_experiment))
      os.execute('mkdir -p ' .. paths.dirname(filename))
      print('<trainer> saving network to '..filename)
      torch.save(filename, net)
   end
      
   -- update logger
   accLogger:add{['% train accuracy'] = trainAcc, ['% val accuracy'] = testAcc}
   errLogger:add{['% train error']    = trainErr, ['% val error']    = testErr}

   -- plot logger
   accLogger:style{['% train accuracy'] = '-', ['% val accuracy'] = '-'}
   errLogger:style{['% train error']    = '-', ['% val error']    = '-'}
   
end


-- load the best net
net =  torch.load(filename)
testAcc, testErr = test(testset, true)

-- save weights for visualization
runthis = false
if runthis then
local weights, gradients = net:parameters()
weights_file = torch.DiskFile(paths.concat(opt.save, string.format('l1_weights.asc')), 'w')
weights_file:noAutoSpacing()
for i = 1, weights[1]:size(1) do
   for j = 1, weights[1]:size(2) do
      weights_file:writeDouble(weights[1][i][j])
      weights_file:writeChar(32)
   end
   weights_file:writeChar(10)
end
end



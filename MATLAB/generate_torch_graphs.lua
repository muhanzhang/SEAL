-- Usage: for converting .mat graphs in tempdata/data_name into .dat format read by torch
-- require 'matio' installed
-- *author: Muhan Zhang, Washington University in St. Louis

matio = require 'matio'
require 'paths'
require 'os'

local cmd = torch.CmdLine()
cmd:option('-dataName',         'train_data',   'Specify the name of the data set to convert')
cmd:option('-ith_experiment',   1,              'Specify which experiment is running')
cmd:option('-multiLabel',       false,          'whether the given label is a vector of multiple labels')

local opt = cmd:parse(arg or {})
print('Converting dataset '..opt.dataName..' to Torch format...')

t0 = os.time()
dataname = opt.dataName
local instance = {}
local label = {}
local datapath = 'tempdata/'..dataname..'_'..opt.ith_experiment..'/'
tmp_label = {}
tmp_instance =  {}
local i = 1
while true do
   print('Torch is reading split '..i)
   if paths.filep(datapath..'split_'..i..'.mat') then
      local current_split = matio.load(datapath..'split_'..i..'.mat')
      local current_label = current_split[string.lower('l'..dataname)]
      local current_instance = current_split[dataname][1]
      for j = 1, current_label:size(1) do
         if opt.multiLabel then
            table.insert(tmp_label, current_label[j]) -- each label is a vector
         else
            table.insert(tmp_label, current_label[j][1]) -- each label is a number
         end
         table.insert(tmp_instance, current_instance[j])
      end
      i = i + 1
   else
      break
   end
end
t1 = os.time()
print(t1-t0)

if not opt.multiLabel then
   -- transform labels to standard 1, 2, ..., n classes
   local label_map = {}
   for i = 1, #tmp_label do
      label_map[tmp_label[i]] = tmp_label[i]
   end
   local unique_label = {}
   for k, v in pairs(label_map) do
      table.insert(unique_label, k)
   end
   table.sort(unique_label)
   local label_map = {}
   for k, v in pairs(unique_label) do
      label_map[v] = k
   end
   for i = 1, #tmp_label do
      local l = tmp_label[i]
      tmp_label[i] = label_map[l]
   end
end

t2 = os.time()
print(t2-t1)

-- save in torch serialization
for i = 1, #tmp_label do
   if i % 100 == 0 then print(i) end
   ins = tmp_instance[i]
   local node_information = ins['nl'].values
   local tmp1 = {ins['am'], node_information}
   local tmp2 = tmp_label[i]
   instance[i] = tmp1
   label[i] = tmp2
end
local dataset = {instance = instance, label = label}
torch.save(datapath..dataname..'_'..opt.ith_experiment..'.dat', dataset)

t3 = os.time()
print(t3-t2)

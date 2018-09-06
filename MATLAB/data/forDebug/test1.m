clear;
linklist = [ 1 2 1; 2 3 1; 3 1 1; 1 4 1; 4 5 1; 5 4 1; 5 5 1];
train = spconvert(linklist);
train(5,5) = 0;
test = train;
tempauc = SimRank(train, test);  

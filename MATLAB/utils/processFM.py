# Transform to libfm format for running FM

import sys

f = open(sys.argv[1],'r')
f1 = open(sys.argv[1]+'.libfm','w')
for line in f.readlines():
  #print line
  #if line == '':
  #	break
  l = line.split()
  tmp = l[0]
  for i in range(1,len(l)):
    tmp += ' '+l[i]+':1'
  tmp += '\n'
  f1.write(tmp)
f.close()
f1.close()


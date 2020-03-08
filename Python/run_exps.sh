#! /bin/bash

data=${1}
for n in 50 100 200
do
  python Main.py --num-paths ${n} --hop "auto" --data-name ${data}
done

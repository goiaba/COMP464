#!/bin/bash

for n in 100 1000 10000; do
  for p in 1 2 4 8 16 32 64 128; do
    echo "====== ${n}n ${p}p 4t start ======"
    folder=strong.${n}n.${p}p.4t
    echo -n "Creating ${folder}... "
    mkdir $folder
    if [ $? -eq 0 ]; then
      echo "done"
    else
      echo "fail"
      exit 1
    fi
    ibrun -n $p -o 0 ./parallel.exec ${n} > ${folder}.out
    mv *.out $folder
    echo -e "====== ${n}n ${p}p 4t done ======\n"
  done
done

mult=1
for p in 1 4 16 64 256; do
  for nloop in 100 500 1000; do
    n=$(($nloop*$mult))
    echo "====== ${n}n ${p}p 4t start ======"
    folder=weak.${n}n.${p}p.4t
    echo -n "Creating ${folder} folder... "
    mkdir $folder
    if [ $? -eq 0 ]; then
        echo "done"
    else
        echo "fail"
        exit 1
    fi
    ibrun -n $p -o 0 ./parallel.exec ${n} > ${folder}.out
    mv *.out $folder
    echo -e "====== ${n}n ${p}p 4t done ======\n"
  done
  mult=$(($mult*2))
done


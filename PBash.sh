#!/usr/bin/bash
echo $"Resultados para STD" >> Results.txt
#bash for loop
cmake --build ParallelAssembleBenchmark_build/ --config Debug --target pbench -j 14 --
for a in 2 4 8 12 16 24 32 64
do
echo $a > PTests.txt
echo Number of threads: $a
echo Number of threads: $a >>  Results.txt
cd ParallelAssembleBenchmark_build/;./pbench;cd ..
done



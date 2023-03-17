#!/usr/bin/bash
echo $"Resultados para TBB" >> Results.txt
#bash for loop
/usr/bin/cmake --build build --config Debug --target pbench -j 14 --
for a in 2 4 8 16 32 64
do
echo $a > PTests.txt
echo Number of threads: $a
echo Number of threads: $a >>  Results.txt
cd build;./pbench;cd ..
done



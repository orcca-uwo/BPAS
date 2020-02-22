"Maple-18 Multiplication" > maple.txt
make;
for (( d=1024; d<=1048576; d=2*d ))
do
        ./mul $d $d 6 0;
        maple -q mul.mm >> maple.txt
done

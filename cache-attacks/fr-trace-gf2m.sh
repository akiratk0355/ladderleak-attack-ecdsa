mkdir -p results-gf2m
rm results-gf2m/*.dat
./FR-trace -r 60000 -F results-gf2m/out.%05d.dat -s 5000 -c 10000 -l 500 -p 2 -H -f ../apps/openssl -m ec2_mult.c:153 -m ec2_mult.c:161 -t bn_gf2m.c:412 -t bn_gf2m.c:412+64 &

../../cache-attacks/FR-trace -r 60000 -F results-gf2m/out.%05d.dat -s 5000 -c 10000 -l 500 -p 4 -H -f dudect_sect163r1_-O2 -m ec2_mult.c:141 -t bn_gf2m.c:412 -t bn_gf2m.c:412+64 > /dev/null &

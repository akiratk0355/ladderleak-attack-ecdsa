mkdir -p results-gfp
rm results-gfp/*.dat
./FR-trace -r 10000 -F results-gfp/out.%05d.dat -s 5000 -c 10000 -l 500 -p 2 -H -f openssl-1.0.2u/apps/openssl -m ecp_smpl.c:927 -m ecp_smpl.c:934  -t bn_lib.c:514 -t bn_lib.c:543 &

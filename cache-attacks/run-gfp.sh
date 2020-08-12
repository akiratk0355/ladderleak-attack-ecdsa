#!/bin/sh
openssl-1.0.2u/apps/openssl ecparam -out secp192r1.pem -name secp192r1 -genkey
openssl-1.0.2u/apps/openssl ec -in secp192r1.pem -pubout -out secp192r1.pub
touch test.txt
./fr-trace-gfp.sh
sleep 1
./ecdsa_sign.sh 10000 secp192r1 results-gfp
killall -9 ./FR-trace
python analysis-gfp.py results-gfp
rm test.txt secp192r1.pem secp192r1.pub

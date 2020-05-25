#!/bin/sh
../apps/openssl ecparam -out sect163r1.pem -name sect163r1 -genkey
../apps/openssl ec -in sect163r1.pem -pubout -out sect163r1.pub
touch test.txt
./fr-trace-gf2m.sh
sleep 1
./ecdsa_sign.sh 10000 sect163r1 results-gf2m
killall -9 ./FR-trace-modified
python analysis-gf2m.py results-gf2m

wget -c https://www.openssl.org/source/old/1.0.2/openssl-1.0.2u.tar.gz
tar xzvf openssl-1.0.2u.tar.gz
cd openssl-1.0.2u; ./Configure -g -fno-inline-functions linux-x86_64 && make clean && make depend && make

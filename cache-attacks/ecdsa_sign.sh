#!/bin/sh
i=1
max=$1
rm $3/$2.log

while : 
do
        echo "----------------------$i"
	date "+%a %b %d %H:%M:%S %Y" >> $3/$2.log
        ../apps/openssl dgst -sha1 -sign $2.pem test.txt > /dev/null 2>>$3/$2.log
        sleep 1
        i=$((i+1))
        if [ $i -gt $max ]; then break; fi
done
echo "DONE"

#!/usr/bin/bash
# Drive the program make_drove.sh
for file in include lib; do
if [[ ! -d ${file} ]]; then
    mkdir ${file}
fi
done

for x in Class OptMap _fnc _fmt; do
    root -l<<EOF
        .L src/io${x}.cc+
EOF
    if [[ ! -f include/io${x}.h ]]; then
        cd include
        ln -s ../src/io${x}.h .
        cd ..
    fi
    if [[ ! -f lib/libio${x}.so ]]; then
        cd lib
        ln -s ../src/io${x}_cc.so libio${x}.so
        cd ..
    fi
done

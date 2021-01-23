#!/usr/bin/bash
# Drive the program make_drove.sh
for file in include lib; do
if [[ ! -d ${file} ]]; then
    mkdir ${file}
fi
done

while read line; do
# for x in Class OptMap _fnc _fmt; do
    root -l<<EOF
        .L src/io${line}.cc+
EOF
    if [[ ! -f include/io${line}.h ]]; then
        cd include
        ln -s ../src/io${line}.h .
        cd ..
    fi
    if [[ ! -f lib/libio${line}.so ]]; then
        cd lib
        ln -s ../src/io${line}_cc.so libio${line}.so
        cd ..
    fi
done < lib_list

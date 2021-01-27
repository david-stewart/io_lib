#!/usr/bin/bash
for file in include lib; do
if [[ ! -d ${file} ]]; then
    mkdir ${file}
fi
done

# make the src_oi files
while read line; do
    make -f src_oi/Makefile lib/liboi${line}.so
    if [ ! -L include/oi${line}.h ]; then
        cd include
        ln -s ../src_oi/oi${line}.h .
        cd ..
    fi
done < oi_lib_list

# Make the src_io files
while read line; do
    root -l<<EOF
        .L src_io/io${line}.cc+
EOF
    if [ ! -L include/io${line}.h ]; then
        cd include
        ln -s ../src_io/io${line}.h .
        cd ..
    fi
    if [ ! -L lib/libio${line}.so ]; then
        cd lib
        ln -s ../src_io/io${line}_cc.so libio${line}.so
        cd ..
    fi
    if [ ! -L lib/io${line}_cc_ACLiC_dict_rdict.pcm ]; then
        cd lib
        ln -s ../src_io/io${line}_cc_ACLiC_dict_rdict.pcm
        cd ..
    fi
done < io_lib_list


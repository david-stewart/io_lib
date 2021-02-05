#!/usr/bin/bash
for file in include lib obj; do
if [[ ! -d ${file} ]]; then
    mkdir ${file}
fi
done

# ------------------------------
# working on src_io files
# ------------------------------

# line the header files
cd include
while read line; do
    if [ ! -L io${line}.h ]; then
        ln -s ../src_io/io${line}.h .
    fi
done < ../io_lib_list
cd ..


# make a script that will compile all the the files
# it looks like:
#    root -l <<EOF
#    gSystem->Load("${ROOUNFOLD}/libRooUnfold.so");
#    .L src_io/ioOptMap.cc+
#    .L src_io/ioXsec_pAu2015.cc+
#    etc...
#    EOF
temp_script=root_L_script.sh
echo 'root -l <<EOF' > ${temp_script}
echo 'gSystem->Load("${ROOUNFOLD}/libRooUnfold.so");' >> ${temp_script}
while read line; do
    echo ".L src_io/io${line}.cc+" >> ${temp_script}
done < io_lib_list
echo 'EOF' >> ${temp_script}
chmod u+x ${temp_script}
./${temp_script}
rm ${temp_script}

cd lib
while read line; do
    if [ ! -L libio${line}.so ]; then
        ln -s ../src_io/io${line}_cc.so libio${line}.so
    fi
    if [ ! -L io${line}_cc_ACLiC_dict_rdict.pcm ]; then
        ln -s ../src_io/io${line}_cc_ACLiC_dict_rdict.pcm
    fi
done < ../io_lib_list
cd ..



# make the src_oi files
while read line; do
    make -f src_oi/Makefile lib/liboi${line}.so
    if [ ! -L include/oi${line}.h ]; then
        cd include
        ln -s ../src_oi/oi${line}.h .
        cd ..
    fi
done < oi_lib_list

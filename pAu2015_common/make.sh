#!/usr/bin/bash
# echo HERE NOW

temp_script=root_L_script.sh
echo "cd src" > ${temp_script}
echo 'root -l <<EOF' >> ${temp_script}

# load all the compilations macros in
while read line; do
    echo ".L io${line}.cc+" >> ${temp_script}
done < src_list

echo 'EOF' >> ${temp_script}
echo "cd .." >> ${temp_script}

chmod u+x ${temp_script}
./${temp_script}
rm ${temp_script}

cd ../include
while read line; do
# echo `pwd`
# echo ${line}
    if [ ! -L io${line}.h ]; then
        ln -s ../pAu2015_common/src/io${line}.h .
    fi
done < ../pAu2015_common/src_list

cd ../lib
while read line; do
# echo `pwd`
# echo ${line}
    if [ ! -L libio${line}.so ]; then
        ln -s ../pAu2015_common/src/io${line}_cc.so libio${line}.so
    fi
    if [ ! -L io${line}_cc_ACLiC_dict_rdict.pcm ]; then
        ln -s ../pAu2015_common/src/io${line}_cc_ACLiC_dict_rdict.pcm
    fi
done < ../pAu2015_common/src_list
cd ../pAu2015_common





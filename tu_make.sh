#!/usr/bin/env bash

# set up files
for file in include lib obj; do
if [[ ! -d ${file} ]]; then
    mkdir ${file}
fi
done

# ------------------------------
# working on src files
# ------------------------------

# link the header files
cd include
while read line; do
    if [[ $line = \#* ]] ; then continue; fi
    # echo $line
    cline=${line%$'\r'}
    if [ ! -L tu${cline}.h ]; then
        ln -s ../src/tu${cline}.h .
    fi
done < ../tu_lib_list
cd ..

# compile files
temp_script=root_L_script.sh
echo 'root -l <<EOF' > ${temp_script}
# echo 'gSystem->Load("${ROOUNFOLD}/libRooUnfold.so");' >> ${temp_script}
while read line; do
    if [[ $line = \#* ]] ; then continue; fi
    cline=${line%$'\r'}
    echo "cout << \".L src/tu${cline}.cc+\" << endl;" >> ${temp_script}
    echo ".L src/tu${cline}.cc+" >> ${temp_script}
done < tu_lib_list
echo 'EOF' >> ${temp_script}
chmod u+x ${temp_script}
./${temp_script}

cd lib
while read line; do
    if [[ $line = \#* ]] ; then continue; fi
    cline=${line%$'\r'}
    if [ ! -L libtu${cline}.so ]; then
        ln -s ../src/tu${cline}_cc.so libtu${cline}.so
    fi
    if [ ! -L tu${cline}_cc_ACLiC_dict_rdict.pcm ]; then
        ln -s ../src/tu${cline}_cc_ACLiC_dict_rdict.pcm
    fi
done < ../tu_lib_list
cd ..

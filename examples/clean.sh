find ./ -name 'tmp.*' -exec rm -r {} \;
find ./ -name 'qepy_input_tmp.in' -exec rm  {} \;
find ./jupyter/ -name 'pwscf.*' -exec rm -r {} \;

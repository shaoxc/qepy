find ./ -name 'tmp.*' -exec rm -r {} \;
find ./ -name 'qepy_input_tmp.in' -exec rm  {} \;
find ./jupyter/ -name 'pwscf.*' -exec rm -r {} \;
find ./jupyter/ -name '.virtual_documents' -exec rm -r {} \;
find ./opt/ase/ -name 'opt.*' -exec rm  {} \;

source="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/01.prep"

cd "$source"

for file in "$source"/*.txt; do

  # isolate .txt file name
  file_name=$(basename "$file")

  # assign modified file names
  new_file="cut_${file_name}"
  new_2_file="space_${file_name}"
  new_3_file="final_${file_name}"

  # delete vestigial text
  sed -e 's/"1"//g' -e 's/"x"//g' -e 's/\"//g' ${file_name} > ${new_file}

  # add 6 spaces before each entry
  awk '{for (i=1; i<=NF; i++) $i = "      " $i} 1' ${new_file} > ${new_2_file}

  # replace single quotes with double
  sed "s/'/\"/g" ${new_2_file} > ${new_3_file}

  # consider adding command to remove cut and space files after processing
done
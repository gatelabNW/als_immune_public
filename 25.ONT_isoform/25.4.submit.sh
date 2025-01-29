input_dir="/projects/b1042/Gate_Lab/projects/als-project/ont_seq/04.isoquant/02.ymal"
input_files=($input_dir/*)

num_files=$(ls ${input_files} | wc -l)
echo "Files in array:"
for file in "${input_files[@]}"; do
    echo "$file"
done

sbatch --array=1-${#input_files[@]} 25.4.batch_isoquant.sh

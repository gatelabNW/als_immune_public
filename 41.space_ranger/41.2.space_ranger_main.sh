


fastq_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.raw"
cytassit_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment"
output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004"
out_log="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/logs.txt"
sample_info="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.meta/ALS_spatial_meta_latest.csv"
sample_id_column="unique_idx"
transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
feature="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
shopt -s extglob
> "$out_log"
column_number=$(head -1 "$sample_info" | awk -v col_name="$sample_id_column" -F',' '{
    for(i=1; i<=NF; i++) {
        if ($i == col_name) {
            print i;
            exit;
        }
    }
}')

# Check if the column was found
if [ -z "$column_number" ]; then
    echo "Column '$column_name' not found."
    exit 1
fi

declare -A processed_values

awk -v col_num="$column_number" -F',' 'NR>1 && $col_num!="" {print $col_num}' "$sample_info" | while IFS= read -r value; do
    if [[ -n "${processed_values[$value]}" ]]; then
        echo "Skip searched value $value"
        continue
    fi
    dir_path="${output_dir}/${value}"
    # create sample's spaceranger output root dir
    # Check if the directory exists
    if [ ! -d "$dir_path" ]; then
        # Directory does not exist, so create it
        mkdir -p "$dir_path"
        echo "Directory created: $dir_path"
    fi

    # if there is output, do not submit job
    # TODO: Commented for rerun purpose; Uncomment for regular logic
    if [[ -d "${output_dir}/${value}/${value}/outs" ]]; then
        echo "${value} output exists. Skipped!" >> "$out_log"
        echo "<============================================================>" >> "$out_log"
        continue
    fi
    donor_id="${value%%___*}"
    temp="${value#*___}"
    slide_id="${temp%%___*}"
    index="${temp#*___}"

    # Mark this value as processed
    processed_values[$value]=1
    # search for matching fastq files
    found_files=($(find "$fastq_dir" -type f -name "*$donor_id*" -and -name "*$index*"))

    # Count the number of files found
    file_count=${#found_files[@]}

    # Check if no files are found
    if [[ $file_count -eq 0 ]]; then
        echo "Error: Sample $donor_id - Incorrect number of files found: $file_count"
        echo "$value has $file_count associated files" >> "$out_log"
        continue
    fi

    # check if more sample is ran more than once
    if [[ $file_count -gt 8 ]]; then
      echo "Warning: Sample $value - Incorrect number of files found: $file_count"
      echo "Check for multiple runs" >> "$out_log"
      continue
    fi
    # extract dir and file name for input file
    # Extract file names from paths and store in a new array
    file_names=()
    declare -A array_hash
    for file_path in "${found_files[@]}"; do
        file_name=$(basename "$file_path")
        extracted_string=${file_name%%_S+([0-9])*}
        if [[ -z ${array_hash[$extracted_string]} ]]; then
            # Add element to the hash map
            array_hash[$extracted_string]=1
            # Add element to the result array
            file_names+=("$extracted_string")
        fi
    done

    # Write library input file for each sample
    cur_sample_lib_file="${dir_path}/lib.csv"
    # clear it if there is content
    > "$cur_sample_lib_file"
    echo "fastqs,sample,library_type" >> "$cur_sample_lib_file"
    for element in "${file_names[@]}"; do
      if [[ $element == *"TS"* ]]; then
        echo "${fastq_dir},${element},Gene Expression" >> "$cur_sample_lib_file"
      elif [[ $element == *"NT"* ]]; then
        echo "${fastq_dir},${element},Antibody Capture" >> "$cur_sample_lib_file"
      fi
    done

    # find corresponding CytAssist stuff
    loupe_dir=$(find "$cytassit_dir" -type d -name "*$donor_id*" -and -name "*$slide_id*" -print -quit)

    # check if ready to rerun
    # TODO: comment out for now, uncomment when rerun for binarized IF intensity
#    rerun_folder="${loupe_dir}/old"
#    if [ ! -d "${rerun_folder}" ]; then
#        echo "${value} not ready to rerun. Skipped!" >> "$out_log"
#        echo "<============================================================>" >> "$out_log"
#        continue
#    fi

    # Check if the dir was found
    if [ -z "$loupe_dir" ]; then
        echo "$donor_id loupe dir not found."
        echo "${donor_id} has no cytassist files" >> "$out_log"
        continue
    fi
    cyt_img=$(find "$loupe_dir" -type f -name "*CAVG*" -print)

    # TODO: Change if running for IF or original
    # regular
#    dark_img=$(find "$loupe_dir" -type f -name "*tdp43-multipage.tif" -maxdepth 1 -print)
    # for IF
    dark_img=$(find "$loupe_dir" -type f -name "*small.tif" -maxdepth 1 -print)
    # running with original json file here, so not binarized cleaned signal that would be added at later stage
    json_file=$(find "$loupe_dir" -type f -name "*.IF.json" -maxdepth 1 -print)
#    json_file=$(find "$loupe_dir" -type f -name "*.tdp43.json" -maxdepth 1 -print)

    # if any is empty continue
    if [[ -z $cyt_img || -z $dark_img || -z $json_file ]]; then
        echo "Check input in cytassit dir. Missing at least a file!"
        echo "${donor_id} has no cytassist files" >> "$out_log"
        continue
    fi

    # TODO Change if running for IF or original
    json_based=$(echo "$json_file" | grep -oP '(?<=/)[^/]+(?=.IF.json)')
#    json_based=$(echo "$json_file" | grep -oP '(?<=/)[^/]+(?=.tdp43.json)')
    slide_id=${json_based%-[[:alnum:]]*}
    temp=${json_based##*-}
    area_id=${temp%.*}

    echo "<============================================================>" >> "$out_log"
    echo "" >> "$out_log"
    echo "ID: ${donor_id}" >> "$out_log"
    echo "Transcriptome ref: ${transcriptome_ref}" >> "$out_log"
    echo "Probe set: ${probe_set}" >> "$out_log"
    echo "Feature ref: ${feature}" >> "$out_log"
    echo "Library: ${cur_sample_lib_file}" >> "$out_log"
    echo "CytAssist image: ${cyt_img}" >> "$out_log"
    echo "Highres image: ${dark_img}" >> "$out_log"
    echo "Alignment file: ${json_file}" >> "$out_log"
    echo "Slide id: ${slide_id}" >> "$out_log"
    echo "Capture area: ${area_id}" >> "$out_log"
    echo "Output dir: ${dir_path}" >> "$out_log"
    echo "" >> "$out_log"
    echo "${value} submitted" >> "$out_log"
    sbatch /projects/p31535/zzhang/als/als_repo/41.space_ranger/41.1.batch_space_ranger.sh \
      "${value}" "${transcriptome_ref}" "${probe_set}" "${feature}" "${cur_sample_lib_file}" "${cyt_img}" "${dark_img}"\
      "${json_file}" "${slide_id}" "${area_id}" "${dir_path}"

done
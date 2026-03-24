#!/bin/bash

# Parse command-line arguments

while getopts s:d:p: flag; do

  case "${flag}" in

    s) source_dir=${OPTARG};;
    d) destination_dir=${OPTARG};;
    p) scriptdir=${OPTARG};;
  esac

done

# Pre-split files

$scriptdir/pre_splitfiles.sh -s "$source_dir" -n 16

output_prefix="splitted_folder"
folders=($(ls -d ${source_dir}/${output_prefix}_*))

# Calculate the number of folders

num_folders=${#folders[@]}

# Generate the Slurm submission script for the first job

cat <<EOL > job1.sh
#!/bin/bash

#SBATCH -p intel-sc3,amd-ep2,intel-sc3-32c
#SBATCH -q normal
#SBATCH -J splitted
#SBATCH -c 16
#SBATCH --output=/storage/caishangLab/fangxiunan/log/exec_%A_%a.log
#SBATCH --error=/storage/caishangLab/fangxiunan/log/exec_%A_%a.err
#SBATCH --array=0-$(($num_folders - 1))
#SBATCH --mem=32G

source_dir="${source_dir}"
destination_dir="${destination_dir}"
scriptdir="${scriptdir}"

output_prefix="splitted_folder"
folders=(\$(ls -d \${source_dir}/\${output_prefix}_*))

# Calculate the number of folders

num_folders=\${#folders[@]}

# Print debug information

echo "SLURM_ARRAY_TASK_ID: \${SLURM_ARRAY_TASK_ID}"
echo "Number of total splitted jobs: \${num_folders}"

if [ \${SLURM_ARRAY_TASK_ID} -lt \${num_folders} ]; then

    splitted_folder=\${folders[\${SLURM_ARRAY_TASK_ID}]}
    echo "Processing folder: \${splitted_folder}"
    splitted_jobs=\$(basename "\${splitted_folder}")
    splitted_resdir="\${destination_dir}/\${splitted_jobs}"

    mkdir -p "\${splitted_resdir}"

    # Run Step0

    echo "Running Step0 for folder \${splitted_folder}"
    \${scriptdir}/step0_test.sh -s "\${splitted_folder}" -d "\${splitted_resdir}"
    if [ \$? -ne 0 ]; then

        echo "Error: step0_test.sh failed for folder \${splitted_folder}"
        exit 1

    else

        echo "Step0 completed successfully for folder \${splitted_folder}"
    fi

    # Run Step1

    echo "Running Step1 for folder \${splitted_folder}"
    \${scriptdir}/step1_hostalignment.sh -s "\${splitted_folder}" -d "\${splitted_resdir}"
    if [ \$? -ne 0 ]; then

        echo "Error: step1_hostalignment.sh failed for folder \${splitted_folder}"
        exit 1

    else

        echo "Step1 completed successfully for folder \${splitted_folder}"
    fi

    # Run Step2

    echo "Running Step2 for folder \${splitted_folder}"
    \${scriptdir}/step2_dedupereads.sh -s "\${splitted_folder}" -d "\${splitted_resdir}"
    if [ \$? -ne 0 ]; then

        echo "Error: step2_dedupereads.sh failed for folder \${splitted_folder}"
        exit 1

    else

        echo "Step2 completed successfully for folder \${splitted_folder}"
    fi

else

    echo "Array index \${SLURM_ARRAY_TASK_ID} is out of range (0 to \${num_folders}-1)"
    exit 1

fi

EOL

# Make the generated script executable

chmod +x job1.sh

## find file differences from two folder and copy file to third folder

rsync -rvcm --compare-dest={parent_folder} {daughter_folder} {diff_folder}

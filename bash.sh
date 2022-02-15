## find file differences from two folder and copy file to third folder

diff -r /Volumes/G_DRIVE/Ujjwal_lab/multiplexing/Didem/Analysis/3_aligned_stacks/ /Volumes/G_DRIVE/Ujjwal_lab/multiplexing/Didem/Analysis/4_artifacts_output/ | grep /Volumes/G_DRIVE/Ujjwal_lab/multiplexing/Didem/Analysis/3_aligned_stacks/ | awk '{print $4}' > defin.txt
sed "/and/d" defin.txt > defin1.txt; rm defin.txt

while read FILENAME; 
do
find $PWD/Analysis/3_aligned_stacks -name "*${FILENAME}*" -exec cp '{}' $PWD/Analysis/5_artifacts_removal \; 
done < defin1.txt



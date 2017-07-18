hawkDir=/scratch/atif/hawk

caseCount=$(cat case_sorted_files.txt | wc -l);
controlCount=$(cat control_sorted_files.txt | wc -l);

$hawkDir/hawk $caseCount $controlCount

$hawkDir/bonf_fasta  
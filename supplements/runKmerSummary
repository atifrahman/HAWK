hawkDir=/scratch/atif/hawk
isCaseSequence=1	# set to 0 if sequence is associated with control


if [ $isCaseSequence -eq 1 ]
then
	file=case_out_wo_bonf.kmerDiff
else
	file=control_out_wo_bonf.kmerDiff
fi

caseCount=$(cat case_sorted_files.txt | wc -l);
controlCount=$(cat control_sorted_files.txt | wc -l);

echo -n $1 > seed.txt
${hawkDir}/kmersearch
grep -f kmers.txt $file > test.txt
${hawkDir}/kmersummary $caseCount $controlCount
mv test.txt kmerstats.txt
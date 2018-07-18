for file in Python-output/*_c1_len1_VSL2b_* 
do
	echo $file
	python printUrl.py -i $file &
done


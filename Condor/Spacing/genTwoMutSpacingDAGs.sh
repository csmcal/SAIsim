#!/bin/bash
#
# genTwoMutSpacingDAGs.sh
# A bash wrapper for running specific parameterizations of the writeTwoMutSpacingDAG.py python script
# 
# Run as:
# ./genTwoMutSpacingDAGs.sh [test]

#
# Intending to run with and without inversions separately:
# 



# Processing inputs, args[0] assumed to be a flag changing the DAG to testing mode
args=("$@")
if [ "${#args[@]}" -gt 0 ]; then
	if [ "${args[0]}" == 0 -o "${args[0]}" == test ]; then
		test=true
	else
		test=false
	fi
else
	test=false
fi

# The path to the DAG generating python script
DAGSCRIPTPATH='writeTwoMutSpacingDAG.py'
DAGFOLDERPATH='DAGs/'

# Combinable arg inputs
TESTINGFLAG='--test'
STDFLAGS=''
INVFLAGS='--inv'
LOWFREQFLAGS='--freq 0.05'

spacing_names=(
	"std_spacing.dag"
	"inv_spacing.dag"
	"rare_std_spacing.dag"
	"rare_inv_spacing.dag"
)

spacing_partials=(
	"--out ${DAGFOLDERPATH}${spacing_names[0]} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${spacing_names[1]} ${INVFLAGS}"
	"--out ${DAGFOLDERPATH}${spacing_names[2]} ${STDFLAGS} ${LOWFREQFLAGS}"
	"--out ${DAGFOLDERPATH}${spacing_names[3]} ${INVFLAGS} ${LOWFREQFLAGS}"
)

# Add the testing flag if given
if [ "$test" = true ]; then
	for i in ${!spacing_partials[@]}; do
		spacing_partials[$i]=${spacing_partials[$i]}' '${TESTINGFLAG}
	done
fi

# Write the DAG list
for i in ${!spacing_partials[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${spacing_partials[$i]}
	python ${DAGSCRIPTPATH} ${spacing_partials[$i]}
done

# Combine them to a single file
spacing_combName='comb_spacing.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${spacing_combName}
for i in ${!spacing_names[@]}; do
	cat "${DAGFOLDERPATH}${spacing_names[$i]}" | tail -r | tail -n +2 | tail -r >> "${DAGFOLDERPATH}${spacing_combName}"
done
tail -n 1 "${DAGFOLDERPATH}${spacing_names[0]}" >> "${DAGFOLDERPATH}${spacing_combName}"
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${spacing_combName}



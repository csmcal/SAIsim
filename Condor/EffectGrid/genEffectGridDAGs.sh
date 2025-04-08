#!/bin/bash
#
# genEffectGridDAGs.sh
# A bash wrapper for running specific parameterizations of the writeEffectGridDAG.py python script
# 
# Run as:
# ./genEffectGridDAGs.sh [test]

#
# Intending to run these parameters separately:
# 
# Encounter number = 10
# No male survival costs
# No female survival costs
# Population size = 10000



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
DAGSCRIPTPATH='writeEffectGridDAG.py'
DAGFOLDERPATH='DAGs/'

# Combinable arg inputs
NUMREPLICATES='50'
REPNFLAG='--rep_N '${NUMREPLICATES}
TESTINGFLAG='--test'
NMCFLAGS='--no_male_costs'
NFCFLAGS='--no_female_costs'
ENC10FLAGS='--enc_N 10'
N4POPFLAGS='--N 1e4 --g 2e5 --exec SAIsim_EffectGrid.sh'


effect_grid_names=(
	"std_effect_grid.dag"
	"nmc_effect_grid.dag"
	"nfc_effect_grid.dag"
	"E10_effect_grid.dag"
	"N4_effect_grid.dag"
)

effect_grid_job_names=(
	"std_effect_grid"
	"nmc_effect_grid"
	"nfc_effect_grid"
	"E10_effect_grid"
	"N4_effect_grid"
)

effect_grid_partials=(
	"--out ${DAGFOLDERPATH}${effect_grid_names[0]} --sub ${effect_grid_job_names[0]} ${REPNFLAG}"
	"--out ${DAGFOLDERPATH}${effect_grid_names[1]} --sub ${effect_grid_job_names[1]} ${REPNFLAG} ${NMCFLAGS}"
	"--out ${DAGFOLDERPATH}${effect_grid_names[2]} --sub ${effect_grid_job_names[2]} ${REPNFLAG} ${NFCFLAGS}"
	"--out ${DAGFOLDERPATH}${effect_grid_names[3]} --sub ${effect_grid_job_names[3]} ${REPNFLAG} ${ENC10FLAGS}"
	"--out ${DAGFOLDERPATH}${effect_grid_names[4]} --sub ${effect_grid_job_names[4]} ${REPNFLAG} ${N4POPFLAGS}"
)

if [ "$test" = true ]; then
	for i in ${!effect_grid_partials[@]}; do
		effect_grid_partials[$i]=${effect_grid_partials[$i]}' '${TESTINGFLAG}
	done
fi

# Write the DAG list
#	Useful default references: --rep 100

for i in ${!effect_grid_partials[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${effect_grid_partials[$i]}
	python ${DAGSCRIPTPATH} ${effect_grid_partials[$i]}
done

# Combine them to a single file
effect_grid_combName='comb_effect_grid.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${effect_grid_combName}
for i in ${!effect_grid_names[@]}; do
	cat "${DAGFOLDERPATH}${effect_grid_names[$i]}" | tail -r | tail -n +2 | tail -r >> "${DAGFOLDERPATH}${effect_grid_combName}"
done
tail -n 1 "${DAGFOLDERPATH}${effect_grid_names[0]}" >> "${DAGFOLDERPATH}${effect_grid_combName}"
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${effect_grid_combName}



n3_effect_grid_names=(
	"std_effect_grid.dag"
	"nmc_effect_grid.dag"
	"nfc_effect_grid.dag"
	"E10_effect_grid.dag"
)

# Combine them to a single file
n3_effect_grid_combName='n3_comb_effect_grid.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${n3_effect_grid_combName}
for i in ${!n3_effect_grid_names[@]}; do
	cat "${DAGFOLDERPATH}${n3_effect_grid_names[$i]}" | tail -r | tail -n +2 | tail -r >> "${DAGFOLDERPATH}${n3_effect_grid_combName}"
done
tail -n 1 "${DAGFOLDERPATH}${n3_effect_grid_names[0]}" >> "${DAGFOLDERPATH}${n3_effect_grid_combName}"
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${n3_effect_grid_combName}


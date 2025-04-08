#!/bin/bash
#
# A bash wrapper for running specific parameterizations of the writeEffectGridDAG.py python script
# 
# Run as:
# ./genParCombEffectGridDAGs.sh [test]

#
# Intending to generate  parameters separately:
# 
# Encounter numbers = [2,3,4,10,100]
# No male survival costs
# No female survival costs



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

NAMESUFFIX='_effect_grid'

# Combinable arg inputs
NUMREPLICATES='50'
REPNFLAG='--rep_N '${NUMREPLICATES}
TESTINGFLAG='--test'
NMCFLAGS='--no_male_costs'
NFCFLAGS='--no_female_costs'
ENC10FLAGS='--enc_N 10'
N4POPFLAGS='--N 1e4 --g 2e5 --exec SAIsim_EffectGrid.sh'
# GRIDFLAGS='--min_S 0 --max_S 0.975 --step_S 0.025 --min_R 0.025 --max_R 1 --step_R 0.025'

enc_nums=(
	'2'
	'3'
	'4'
	'10'
	'100'
	'200'
)

# cost_prefixes=(
# 	"std_"
# 	"nmc_"
# 	"nfc_"
# )

for i in ${!enc_nums[@]}; do
	enc_suffix="_e${enc_nums[i]}${NAMESUFFIX}"
	name="std${enc_suffix}"
	effect_grid_names+=("${name}.dag")
	effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG}")
	name="nmc${enc_suffix}"
	effect_grid_names+=("${name}.dag")
	effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG} ${NMCFLAGS}")
	name="nfc${enc_suffix}"
	effect_grid_names+=("${name}.dag")
	effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG} ${NFCFLAGS}")
done

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
effect_grid_combName='comb_enc_cost_effect_grid.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${effect_grid_combName}
for i in ${!effect_grid_names[@]}; do
	cat "${DAGFOLDERPATH}${effect_grid_names[$i]}" | tail -r | tail -n +2 | tail -r >> "${DAGFOLDERPATH}${effect_grid_combName}"
done
tail -n 1 "${DAGFOLDERPATH}${effect_grid_names[0]}" >> "${DAGFOLDERPATH}${effect_grid_combName}"
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${effect_grid_combName}




# Writing DAGs for the near-neutral range at higher resolution
FINEGRIDFLAGS='--min_S 0.8 --max_S 0.995 --step_S 0.005 --min_R 0.005 --max_R .2 --step_R 0.005'

for i in ${!enc_nums[@]}; do
	enc_suffix="_e${enc_nums[i]}_fine${NAMESUFFIX}"
	name="std${enc_suffix}"
	fine_effect_grid_names+=("${name}.dag")
	fine_effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG} ${FINEGRIDFLAGS}")
	name="nmc${enc_suffix}"
	fine_effect_grid_names+=("${name}.dag")
	fine_effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG} ${FINEGRIDFLAGS} ${NMCFLAGS}")
	name="nfc${enc_suffix}"
	fine_effect_grid_names+=("${name}.dag")
	fine_effect_grid_partials+=("--out ${DAGFOLDERPATH}${name}.dag --sub ${name} --enc_N ${enc_nums[i]} ${REPNFLAG} ${FINEGRIDFLAGS} ${NFCFLAGS}")
done

if [ "$test" = true ]; then
	for i in ${!fine_effect_grid_partials[@]}; do
		fine_effect_grid_partials[$i]=${fine_effect_grid_partials[$i]}' '${TESTINGFLAG}
	done
fi

# Write the DAG list
#	Useful default references: --rep 100

for i in ${!fine_effect_grid_partials[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${fine_effect_grid_partials[$i]}
	python ${DAGSCRIPTPATH} ${fine_effect_grid_partials[$i]}
done

# Combine them to a single file
fine_effect_grid_combName='comb_enc_cost_fine_effect_grid.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${fine_effect_grid_combName}
for i in ${!fine_effect_grid_names[@]}; do
	cat "${DAGFOLDERPATH}${fine_effect_grid_names[$i]}" | tail -r | tail -n +2 | tail -r >> "${DAGFOLDERPATH}${effect_grid_combName}"
done
tail -n 1 "${DAGFOLDERPATH}${fine_effect_grid_names[0]}" >> "${DAGFOLDERPATH}${fine_effect_grid_combName}"
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${fine_effect_grid_combName}


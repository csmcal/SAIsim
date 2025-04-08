#!/bin/bash
#
# dmelFullArrangementDAGs.sh
# A bash wrapper for running specific parameterizations of the writeFullDAG.py python script

#
# Based on Drosophila melanogaster parameters of:
#

# At autosomal c of 2.18e-8 (Comeron et al. 2012), expect 0.436 crossovers in the 20Mbp region 
# γ=1.25E-7 events/bp/female meiosis * 518 bp (Comeron et al. 2012) = 6.475e-5 conversions/bp/female meiosis (DO MALES CONVERT?)
# 5.21E-9 mutations/bp/gen on autosomes (Huang et al. 2016) Unknown rate of SA mutation? 1 in 100,000 for 5.21e-14? 1 in 1000 for 5.21e-12?
# Inversion rate maybe 100 times less common than SA variation?


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

# Script and output directory filepaths
DAGSCRIPTPATH='writeFullDAG.py'
DAGFOLDERPATH='DAGs/Arrangements/'

ARREXEC='SAIsim_Full_Arrangements.sh'
ARRSCRIPT='SAIsim_Full_Arrangements.py'


# Shared arg inputs
TESTINGFLAG='--test --rep 5'
STDFLAGS='--wf True --nmc False --nfc False'
NMCFLAGS='--wf True --nmc True --nfc False'
NFCFLAGS='--wf True --nmc False --nfc True'


N3POPARGS='--N 1e3 --g 4e4'

# For Scaling from Ne of 2e6 to 1e3; a factor of 2e3 for scaled /bp values, with more common mutation rates
#	r = 4.36e-5
#	g = 1.295e-1
#	µ = 1.042e-8
#	i = 1.042e-10
# For a track length of 1 Mbp:
#	r = 43.6
#	µ = 1.042e-2
#	i = 1.042e-4
# For a track length of 10 Kbp:
#	r = 4.36e-1
#	µ = 1.042e-4
#	i = 1.042e-6


N3L1MbpARGS='--m 1.042e-2 --i 1.042e-4 1.042e-2 --e 10 100 --r 4.36e1 --c 1.295e-1'
N3L10KbpARGS='--m 1.042e-4 --i 1.042e-6 1.042e-4 --e 10 100 --r 4.36e-1 --c 1.295e-1'

N3_mel_arrangement_names=(
	"N3_1Mbp_mel_Full_Arrangements.dag"
	"N3_1Mbp_mel_nmc_Full_Arrangements.dag"
	"N3_1Mbp_mel_nfc_Full_Arrangements.dag"
	"N3_10Kbp_mel_Full_Arrangements.dag"
	"N3_10Kbp_mel_nmc_Full_Arrangements.dag"
	"N3_10Kbp_mel_nfc_Full_Arrangements.dag"
)

N3_mel_arr=(
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[0]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpARGS} ${N3POPARGS} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[1]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpARGS} ${N3POPARGS} ${NMCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[2]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpARGS} ${N3POPARGS} ${NFCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[3]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpARGS} ${N3POPARGS} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[4]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpARGS} ${N3POPARGS} ${NMCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_arrangement_names[5]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpARGS} ${N3POPARGS} ${NFCFLAGS}"
)

if [ "$test" = true ]; then
	for i in ${!N3_mel_arr[@]}; do
		N3_mel_arr[$i]=${N3_mel_arr[$i]}' '${TESTINGFLAG}
	done
fi

# Write the Ne = 1e3 DAGs
#	Useful default references: --checkpointTime 7200 (seconds, =2hr)	--rep 1000

for i in ${!N3_mel_arr[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${N3_mel_arr[$i]}
	python ${DAGSCRIPTPATH} ${N3_mel_arr[$i]}
done

# Combine them to a single file
N3_mel_arr_combName='N3_mel_Full_Arrangements.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${N3_mel_arr_combName}
for i in ${!N3_mel_arr[@]}; do
	cat ${DAGFOLDERPATH}${N3_mel_arrangement_names[$i]} | tail -r | tail -n +2 | tail -r >> ${DAGFOLDERPATH}${N3_mel_arr_combName}
done
tail -n 1 ${DAGFOLDERPATH}${N3_mel_arrangement_names[0]} >> ${DAGFOLDERPATH}${N3_mel_arr_combName}
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${N3_mel_arr_combName}""




# For Scaling from Ne of 2e6 to 1e3; a factor of 2e3 for scaled /bp values, with rarer mutations
#	r = 4.36e-5
#	g = 1.295e-1
#	µ = 1.042e-10
#	i = 1.042e-12
# For a track length of 1 Mbp:
#	r = 43.6
#	µ = 1.042e-4
#	i = 1.042e-6
# For a track length of 10 Kbp:
#	r = 4.36e-1
#	µ = 1.042e-6
#	i = 1.042e-8

N3L1MbpRareMutARGS='--m 1.042e-4 --i 1.042e-6 1.042e-4 0 --e 10 100 --r 4.36e1 --c 1.295e-1'
N3L10KbpRareMutARGS='--m 1.042e-6 --i 1.042e-8 1.042e-6 0 --e 10 100 --r 4.36e-1 --c 1.295e-1'


N3_mel_rareMut_arr_names=(
	"N3_1Mbp_mel_rareMut_Full_Arrangements.dag"
	"N3_1Mbp_mel_rareMut_nmc_Full_Arrangements.dag"
	"N3_1Mbp_mel_rareMut_nfc_Full_Arrangements.dag"
	"N3_10Kbp_mel_rareMut_Full_Arrangements.dag"
	"N3_10Kbp_mel_rareMut_nmc_Full_Arrangements.dag"
	"N3_10Kbp_mel_rareMut_nfc_Full_Arrangements.dag"
)

N3_mel_rareMut_arr=(
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[0]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpRareMutARGS} ${N3POPARGS} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[1]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpRareMutARGS} ${N3POPARGS} ${NMCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[2]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpRareMutARGS} ${N3POPARGS} ${NFCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[3]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpRareMutARGS} ${N3POPARGS} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[4]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpRareMutARGS} ${N3POPARGS} ${NMCFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[5]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpRareMutARGS} ${N3POPARGS} ${NFCFLAGS}"
)


if [ "$test" = true ]; then
	for i in ${!N3_mel_rareMut_arr[@]}; do
		N3_mel_rareMut_arr[$i]=${N3_mel_rareMut_arr[$i]}' '${TESTINGFLAG}
	done
fi

# Write the Ne = 1e3 DAGs
#	Useful default references: --checkpointTime 7200 (seconds, =2hr)	--rep 1000

for i in ${!N3_mel_rareMut_arr[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${N3_mel_rareMut_arr[$i]}
	python ${DAGSCRIPTPATH} ${N3_mel_rareMut_arr[$i]}
done

# Combine them to a single file
N3_mel_rareMut_arr_combName='N3_mel_rareMut_Full_Arrangements.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${N3_mel_rareMut_arr_combName}
for i in ${!N3_mel_arr[@]}; do
	cat ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[$i]} | tail -r | tail -n +2 | tail -r >> ${DAGFOLDERPATH}${N3_mel_rareMut_arr_combName}
done
tail -n 1 ${DAGFOLDERPATH}${N3_mel_rareMut_arr_names[0]} >> ${DAGFOLDERPATH}${N3_mel_rareMut_arr_combName}
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${N3_mel_rareMut_arr_combName}""




# Std inv mutation rate, no SA mutation simulations


N3L1MbpNoMutARGS='--m 0 --i 1.042e-6 1.042e-4 1.042e-2 --e 10 100 --r 4.36e1 --c 1.295e-1'
N3L10KbpNoMutARGS='--m 0 --i 1.042e-8 1.042e-6 1.042e-4 --e 10 100 --r 4.36e-1 --c 1.295e-1'

N3_mel_noMut_arr_names=(
	"N3_1Mbp_mel_noMut_Full_Arrangements.dag"
	"N3_10Kbp_mel_noMut_Full_Arrangements.dag"
)

N3_mel_noMut_arr=(
	"--out ${DAGFOLDERPATH}${N3_mel_noMut_arr_names[0]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L1MbpNoMutARGS} ${N3POPARGS} ${STDFLAGS}"
	"--out ${DAGFOLDERPATH}${N3_mel_noMut_arr_names[1]} --exec ${ARREXEC} --pyscript ${ARRSCRIPT} ${N3L10KbpNoMutARGS} ${N3POPARGS} ${STDFLAGS}"
)

if [ "$test" = true ]; then
	for i in ${!N3_mel_noMut_arr[@]}; do
		N3_mel_noMut_arr[$i]=${N3_mel_noMut_arr[$i]}' '${TESTINGFLAG}
	done
fi

for i in ${!N3_mel_noMut_arr[@]}; do
	echo 'python '${DAGSCRIPTPATH}' '${N3_mel_noMut_arr[$i]}
	python ${DAGSCRIPTPATH} ${N3_mel_noMut_arr[$i]}
done


# Combine them to a single file
N3_mel_noMut_arr_combName='N3_mel_noMut_Full_Arrangements.dag'
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${N3_mel_noMut_arr_combName}
for i in ${!N3_mel_noMut_arr[@]}; do
	cat ${DAGFOLDERPATH}${N3_mel_noMut_arr_names[$i]} | tail -r | tail -n +2 | tail -r >> ${DAGFOLDERPATH}${N3_mel_noMut_arr_combName}
done
tail -n 1 ${DAGFOLDERPATH}${N3_mel_noMut_arr_names[0]} >> ${DAGFOLDERPATH}${N3_mel_noMut_arr_combName}
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${N3_mel_noMut_arr_combName}""




# Combine all N3 DAGs to a single file
N3_mel_arr_all_combName='N3_mel_Full_Arrangements_all.dag'
N3_mel_arr_all_names=("${N3_mel_arrangement_names[@]}" "${N3_mel_noMut_arr_names[@]}")
# Make or overwrite an empty file
:>| ${DAGFOLDERPATH}${N3_mel_arr_all_combName}
for i in ${!N3_mel_arr_all_names[@]}; do
	cat ${DAGFOLDERPATH}${N3_mel_arr_all_names[$i]} | tail -r | tail -n +2 | tail -r >> ${DAGFOLDERPATH}${N3_mel_arr_all_combName}
done
tail -n 1 ${DAGFOLDERPATH}${N3_mel_arr_all_names[0]} >> ${DAGFOLDERPATH}${N3_mel_arr_all_combName}
printf '%s\n\n' "Concatenated into "${DAGFOLDERPATH}${N3_mel_arr_all_combName}












# N4POPARGS='--N 1e4 --g 2e5'

# # For Scaling from Ne of 2e6 to 1e4; a factor of 2e2 for scaled /bp values, with rarer mutations
# #	r = 4.36e-6
# #	g = 1.295e-2
# #	µ = 1.042e-11
# #	i = 1.042e-13
# # For a track length of 1 Mbp:
# #	r = 4.36
# #	µ = 1.042e-5
# #	i = 1.042e-7
# # For a track length of 10 Kbp:
# #	r = 4.36e-2
# #	µ = 1.042e-7
# #	i = 1.042e-9


# N4L1MbpRareMutARGS='--m 1.042e-4 --i 1.042e-6 --e 2 10 100 --r 4.36e0 --c 1.295e-1'
# N4L10KbpRareMutARGS='--m 1.042e-6 --i 1.042e-8 --e 2 10 100 --r 4.36e-2 --c 1.295e-1'

# N4_mel_rareMut=(
# 	"--out N4_1Mbp_mel_rareMut_Full.dag ${N4L1MbpRareMutARGS} ${N4POPARGS} ${STDFLAGS}"
# 	"--out N4_1Mbp_mel_rareMut_nmc_Full.dag ${N4L1MbpRareMutARGS} ${N4POPARGS} ${NMCFLAGS}"
# 	"--out N4_1Mbp_mel_rareMut_nfc_Full.dag ${N4L1MbpRareMutARGS} ${N4POPARGS} ${NFCFLAGS}"
# 	"--out N4_10Kbp_mel_rareMut_Full.dag ${N4L10KbpRareMutARGS} ${N4POPARGS} ${STDFLAGS}"
# 	"--out N4_10Kbp_mel_rareMut_nmc_Full.dag ${N4L10KbpRareMutARGS} ${N4POPARGS} ${NMCFLAGS}"
# 	"--out N4_10Kbp_mel_rareMut_nfc_Full.dag ${N4L10KbpRareMutARGS} ${N4POPARGS} ${NFCFLAGS}"
# )

# if [ "$test" = true ]; then
# 	for i in ${!N4_mel_rareMut[@]}; do
# 		N4_mel_rareMut[$i]=${N4_mel_rareMut[$i]}' '${TESTINGFLAG}
# 	done
# fi

# # Write the Ne = 1e3 DAGs
# #	Useful default references: --checkpointTime 7200 (seconds, =2hr)	--rep 1000

# for i in ${!N4_mel_rareMut[@]}; do
# 	echo 'python ../writeFullDAG.py '${N4_mel_rareMut[$i]}
# 	python ../writeFullDAG.py ${N4_mel_rareMut[$i]}
# done

# # Combine them to a single file
# cat N4_1Mbp_mel_rareMut_Full.dag | tail -r | tail -n +2 | tail -r > N4_mel_rareMut_Full.dag
# cat N4_1Mbp_mel_rareMut_nmc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_rareMut_Full.dag
# cat N4_1Mbp_mel_rareMut_nfc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_rareMut_Full.dag
# cat N4_10Kbp_mel_rareMut_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_rareMut_Full.dag
# cat N4_10Kbp_mel_rareMut_nmc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_rareMut_Full.dag
# cat N4_10Kbp_mel_rareMut_nfc_Full.dag >> N4_mel_rareMut_Full.dag


# # Write the Ne = 1e4 DAGs
# #	Useful default references: --checkpointTime 7200 (seconds, =2hr)	--rep 1000
# python ../writeFullDAG.py --out N4_1Mbp_mel_Full.dag --m 1.042e-5 --i 1.042e-7 --e 10 100 --r 4.36e0 --c 1.295e-2 --N 1e4 --g 2e5 --wf True --nmc False --nfc False
# python ../writeFullDAG.py --out N4_10Kbp_mel_Full.dag --m 1.042e-7 --i 1.042e-9 --e 10 100 --r 4.36e-2 --c 1.295e-2 --N 1e4 --g 2e5 --wf True --nmc False --nfc False

# # Combine them to a single file
# cat N4_1Mbp_mel_Full.dag | tail -r | tail -n +2 | tail -r > N4_mel_Full.dag
# cat N4_10Kbp_mel_Full.dag >> N4_mel_Full.dag



# # For Scaling from Ne of 2e6 to 1e3; a factor of 2e3 for scaled /bp values
# #	r = 4.36e-5
# #	g = 1.295e-1
# #	µ = 1.042e-8
# #	i = 1.042e-10
# # For a track length of 1 Mbp:
# #	r = 4.36
# #	µ = 1.042e-3
# #	i = 1.042e-5
# # For a track length of 10 Kbp:
# #	r = 4.36e-2
# #	µ = 1.042e-5
# #	i = 1.042e-7



# N4L1MbpARGS='--m 1.042e-2 --i 1.042e-4 --e 10 100 --r 4.36e0 --c 1.295e-1'
# N4L10KbpARGS='--m 1.042e-4 --i 1.042e-6 --e 10 100 --r 4.36e-2 --c 1.295e-1'

# N4_mel=(
# 	"--out N4_1Mbp_mel_Full.dag ${N4L1MbpARGS} ${N4POPARGS} ${STDFLAGS}"
# 	"--out N4_1Mbp_mel_nmc_Full.dag ${N4L1MbpARGS} ${N4POPARGS} ${NMCFLAGS}"
# 	"--out N4_1Mbp_mel_nfc_Full.dag ${N4L1MbpARGS} ${N4POPARGS} ${NFCFLAGS}"
# 	"--out N4_10Kbp_mel_Full.dag ${N4L10KbpARGS} ${N4POPARGS} ${STDFLAGS}"
# 	"--out N4_10Kbp_mel_nmc_Full.dag ${N4L10KbpARGS} ${N4POPARGS} ${NMCFLAGS}"
# 	"--out N4_10Kbp_mel_nfc_Full.dag ${N4L10KbpARGS} ${N4POPARGS} ${NFCFLAGS}"
# )

# if [ "$test" = true ]; then
# 	for i in ${!N4_mel[@]}; do
# 		N4_mel[$i]=${N4_mel[$i]}' '${TESTINGFLAG}
# 	done
# fi

# # Write the Ne = 1e3 DAGs
# #	Useful default references: --checkpointTime 7200 (seconds, =2hr)	--rep 1000

# for i in ${!N4_mel[@]}; do
# 	echo 'python ../writeFullDAG.py '${N4_mel[$i]}
# 	python ../writeFullDAG.py ${N4_mel[$i]}
# done

# # Combine them to a single file
# cat N4_1Mbp_mel_Full.dag | tail -r | tail -n +2 | tail -r > N4_mel_Full.dag
# cat N4_1Mbp_mel_nmc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_Full.dag
# cat N4_1Mbp_mel_nfc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_Full.dag
# cat N4_10Kbp_mel_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_Full.dag
# cat N4_10Kbp_mel_nmc_Full.dag | tail -r | tail -n +2 | tail -r >> N4_mel_Full.dag
# cat N4_10Kbp_mel_nfc_Full.dag >> N4_mel_Full.dag



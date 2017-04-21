#!/usr/bin/env bash

################################################################################
# PRIMER REMOVAL
################################################################################
# requires:
# 1. cutadapt
# 2. seqtk
# 3. revcom function
revcom (){

  echo $1 |\
    rev   |\
    tr ACGTWSMKRYBDHVNacgtwsmkrybdhvn TGCASWKMYRVHDBNtgcaswkmyrvhdbn

}

PRIMER_REMOVAL_INPUT="${1}"
NO_PRIMERS="${2}"
PRIMER1="${3}"
PRIMER2="${4}"


# reverse complement primers
PRIMER1RC=$( revcom "${PRIMER1}" )
PRIMER2RC=$( revcom "${PRIMER2}" )

# count lines in primer removal input
echo $(date +%Y-%m-%d\ %H:%M) "Counting sequences in primer removal input..."
seq_N_demult_concat=$( grep -e '^>' --count "${PRIMER_REMOVAL_INPUT}" )
echo $(date +%Y-%m-%d\ %H:%M) "${seq_N_demult_concat}" "sequences found in primer removal input"
echo

echo $(date +%Y-%m-%d\ %H:%M) "Beginning primer removal..."
# remove primer 1 from left side of sequences
#
# #Alternative script to create 4 folders with sequences with each primer
#
#
(cutadapt\
	-a primerR1_removed="${PRIMER1RC}" \
	-a primerR2_removed="${PRIMER2RC}" --no-trim	-o "${PRIMER_REMOVAL_INPUT}"-{name}.fasta "${PRIMER_REMOVAL_INPUT}")
	#	\
		#-a primerR2_removed=^"${PRIMER2}"\
			#-a primerR1_removed="${PRIMER2RC}"$ -a primerR2_removed="${PRIMER1RC}"$\
			#
#Now reverse complement those reads in which we found the reversed  Fwd primer

seqtk seq -r "${PRIMER_REMOVAL_INPUT}"-primerR1_removed.fasta > "${PRIMER_REMOVAL_INPUT}"-R1_RC.fasta
#First change is to remove the file already used
rm "${PRIMER_REMOVAL_INPUT}"-primerR1_removed.fasta
#See how it goes - not much improvement. Now change the cat to add file 1 into 2, doesn't improve
cat "${PRIMER_REMOVAL_INPUT}"-R1_RC.fasta >> "${PRIMER_REMOVAL_INPUT}"-primerR2_removed.fasta
rm "${PRIMER_REMOVAL_INPUT}"-R1_RC.fasta
mv "${PRIMER_REMOVAL_INPUT}"-primerR2_removed.fasta "${PRIMER_REMOVAL_INPUT}"

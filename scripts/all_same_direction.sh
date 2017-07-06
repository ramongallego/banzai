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
echo $(date +%Y-%m-%d\ %H:%M) "Counting sequences in reorder input"
seq_N_demult_concat=$( grep -e '^>' --count "${PRIMER_REMOVAL_INPUT}" )
echo $(date +%Y-%m-%d\ %H:%M) "${seq_N_demult_concat}" "sequences found in reorder input"
echo

echo $(date +%Y-%m-%d\ %H:%M) "Beginning primer search..."
# remove primer 1 from left side of sequences
#
# #Alternative script to create 4 folders with sequences with each primer
#
#
(cutadapt\
	-g primerF1="${PRIMER1}" \
	-g primerR1="${PRIMER2}" --no-trim	-o "${PRIMER_REMOVAL_INPUT}"-{name}.fasta "${PRIMER_REMOVAL_INPUT}")
	#	\
		#-a primerR2_removed=^"${PRIMER2}"\
			#-a primerR1_removed="${PRIMER2RC}"$ -a primerR2_removed="${PRIMER1RC}"$\
			#
#Now reverse complement those reads in which we found the reversed  Fwd primer

seqtk seq -r "${PRIMER_REMOVAL_INPUT}"-primerR1.fasta >> "${PRIMER_REMOVAL_INPUT}"-primerF1.fasta
#First change is to remove the file already used
rm "${PRIMER_REMOVAL_INPUT}"-primerR1.fasta
#See how it goes - not much improvement. Now change the cat to add file 1 into 2, doesn't improve

mv "${PRIMER_REMOVAL_INPUT}"-primerF1.fasta "${PRIMER_REMOVAL_INPUT}"

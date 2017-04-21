#!/usr/bin/env bash

source "banzai_params.sh"

Rscript sr.r "${OUTPUT_DIRECTORY}" "${SCRIPT_DIR}" "r.output.txt"

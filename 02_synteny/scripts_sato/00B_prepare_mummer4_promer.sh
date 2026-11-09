#!/bin/bash

set -euo pipefail

BASE=/ru-auth/local/home/asato/scratch/synteny_4species/reanalysis_20260824
MUMMER_ENV=/ru-auth/local/home/asato/scratch/conda_envs/mummer4_env

SOURCE=${MUMMER_ENV}/bin/promer
OUTPUT=${BASE}/software_wrappers/promer_mummer4_4.0.1
CHECK=${BASE}/checks/software/mummer4_promer_configuration.tsv

BIN_DIR=${MUMMER_ENV}/bin
LIBEXEC_DIR=${MUMMER_ENV}/libexec/mummer
LIB_DIR=${MUMMER_ENV}/lib/mummer

for required in \
    "${SOURCE}" \
    "${BIN_DIR}/mummer" \
    "${BIN_DIR}/show-coords" \
    "${LIBEXEC_DIR}/mgaps" \
    "${LIBEXEC_DIR}/prepro" \
    "${LIBEXEC_DIR}/postpro"
do
    if [[ ! -x "${required}" ]]; then
        echo "ERROR: missing executable: ${required}" >&2
        exit 1
    fi
done

cp "${SOURCE}" "${OUTPUT}"

sed -i \
    -e "s|my \\\$BIN_DIR = '@BIN_DIR';|my \\\$BIN_DIR = '${BIN_DIR}';|" \
    -e "s|my \\\$AUX_BIN_DIR = '@LIBEXEC_DIR';|my \\\$AUX_BIN_DIR = '${LIBEXEC_DIR}';|" \
    -e "s|my \\\$LIB_DIR = '@LIB_DIR';|my \\\$LIB_DIR = '${LIB_DIR}';|" \
    "${OUTPUT}"

chmod +x "${OUTPUT}"

if grep -qE '@BIN_DIR|@LIBEXEC_DIR|@LIB_DIR' "${OUTPUT}"; then
    echo "ERROR: unresolved placeholders remain in ${OUTPUT}" >&2
    exit 1
fi

{
    printf "metric\tvalue\n"
    printf "conda_environment\t%s\n" "${MUMMER_ENV}"
    printf "conda_package\tmummer4=4.0.1\n"
    printf "source_promer\t%s\n" "${SOURCE}"
    printf "configured_promer\t%s\n" "${OUTPUT}"
    printf "BIN_DIR\t%s\n" "${BIN_DIR}"
    printf "LIBEXEC_DIR\t%s\n" "${LIBEXEC_DIR}"
    printf "LIB_DIR\t%s\n" "${LIB_DIR}"
    printf "status\tPASS\n"
} > "${CHECK}"

echo "===== Configured PROmer ====="
grep -nE \
    'my \$BIN_DIR|my \$AUX_BIN_DIR|my \$LIB_DIR' \
    "${OUTPUT}"

echo "===== Version ====="
"${OUTPUT}" --version

echo "===== Configuration record ====="
cat "${CHECK}"

echo "FINAL: PASS"

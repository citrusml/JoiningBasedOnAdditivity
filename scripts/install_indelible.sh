#!/bin/bash
# ==============================================================================
# Build and install INDELible v1.03 binary in <repo_root>/bin/indelible
# Safe for concurrent invocation across LSF job array tasks (uses atomic mv).
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${ROOT_DIR}/bin"
TARGET_BIN="${BIN_DIR}/indelible"

mkdir -p "${BIN_DIR}"

# 1. Check if TARGET_BIN already exists and runs on the current OS/architecture
if [ -x "${TARGET_BIN}" ]; then
    CHECK_OUT=$(echo "" | "${TARGET_BIN}" 2>&1 || true)
    if echo "${CHECK_OUT}" | grep -q "INDELible"; then
        echo "INDELible is already installed and functional at: ${TARGET_BIN}"
        exit 0
    fi
fi

# 2. Check if sibling PhyloTree_Methods_Comparison/bin/indelible exists and runs on this OS
SIBLING_BIN="${ROOT_DIR}/../PhyloTree_Methods_Comparison/bin/indelible"
if [ -x "${SIBLING_BIN}" ]; then
    CHECK_SIB=$(echo "" | "${SIBLING_BIN}" 2>&1 || true)
    if echo "${CHECK_SIB}" | grep -q "INDELible"; then
        TMP_CP="${TARGET_BIN}.tmp.$$"
        cp "${SIBLING_BIN}" "${TMP_CP}"
        chmod +x "${TMP_CP}"
        mv -f "${TMP_CP}" "${TARGET_BIN}"
        echo "Copied working INDELible binary from ${SIBLING_BIN} -> ${TARGET_BIN}"
        exit 0
    fi
fi

# 3. Compile INDELible v1.03 from source
echo "Building INDELible v1.03 from evolbioinfo/indelible..."
TMP_BUILD_DIR="$(mktemp -d -t indelible_build_XXXXXX)"
trap 'rm -rf "${TMP_BUILD_DIR}"' EXIT

git clone --depth 1 https://github.com/evolbioinfo/indelible.git "${TMP_BUILD_DIR}/indelible"
cd "${TMP_BUILD_DIR}/indelible/src"

CXX="${CXX:-g++}"
if ! command -v "${CXX}" >/dev/null 2>&1; then
    if command -v clang++ >/dev/null 2>&1; then
        CXX="clang++"
    else
        echo "Error: Neither g++ nor clang++ was found in PATH." >&2
        exit 1
    fi
fi

echo "Compiling with ${CXX} -O3..."
${CXX} -O3 -o indelible indelible.cpp -lm

TMP_BIN="${TARGET_BIN}.tmp.$$"
cp indelible "${TMP_BIN}"
chmod +x "${TMP_BIN}"
mv -f "${TMP_BIN}" "${TARGET_BIN}"

echo "INDELible v1.03 successfully installed at: ${TARGET_BIN}"

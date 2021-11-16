#!/usr/bin/env sh

set -e

C_FILES=$(find src -type f -name '*.c')
H_FILES=$(find src -type f -name '*.h')

clang-format -i ${C_FILES} ${H_FILES}

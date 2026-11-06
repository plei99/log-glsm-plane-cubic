#!/usr/bin/env bash
# Run every regression suite with per-suite timing.
#
#   ./run_tests.sh            run all suites
#   FAST=1 ./run_tests.sh     skip the slow localization/orchestration suites
#
# Every suite prints its own pass line; this script adds wall-clock timing,
# a summary, and a nonzero exit code if anything failed.  DOT_SAGE defaults
# to the user environment; export it first to redirect Sage's profile.

set -u
cd "$(dirname "${BASH_SOURCE[0]}")"
export PYTHONPATH="$PWD${PYTHONPATH:+:$PYTHONPATH}"

# Suites excluded by FAST=1, including the longer complete-intersection
# comparisons added in November 2026. The full suite is the release check.
SLOW_SUITES="
test_givental_teleman
test_hodge_integrals
test_o3_fixed_locus_graphs
test_genus_four_marked
test_grassmannian_laurent
test_projective_complete_intersections
test_toric_complete_intersections
"

is_slow() {
    for slow in $SLOW_SUITES; do
        if [ "$1" = "$slow" ]; then
            return 0
        fi
    done
    return 1
}

failures=""
skipped=0
total_start=$(date +%s)

for suite_file in tests/test_*.sage; do
    suite=$(basename "${suite_file%.sage}")
    if [ "${FAST:-}" = "1" ] && is_slow "$suite"; then
        skipped=$((skipped + 1))
        printf '%-46s skipped (FAST=1)\n' "$suite"
        continue
    fi
    start=$(date +%s)
    if output=$(sage "$suite_file" 2>&1); then
        elapsed=$(( $(date +%s) - start ))
        summary=$(printf '%s\n' "$output" | tail -1)
        printf '%-46s %4ss  %s\n' "$suite" "$elapsed" "$summary"
    else
        elapsed=$(( $(date +%s) - start ))
        printf '%-46s %4ss  FAILED\n' "$suite" "$elapsed"
        printf '%s\n' "$output" | tail -12 | sed 's/^/    /'
        failures="$failures $suite"
    fi
done

start=$(date +%s)
if output=$(python3 -m unittest discover -s tests -p 'test_*.py' -q 2>&1); then
    elapsed=$(( $(date +%s) - start ))
    printf '%-46s %4ss  %s\n' "python_unittests" "$elapsed" \
        "python unittest passed"
else
    elapsed=$(( $(date +%s) - start ))
    printf '%-46s %4ss  FAILED\n' "python_unittests" "$elapsed"
    printf '%s\n' "$output" | tail -12 | sed 's/^/    /'
    failures="$failures python_unittests"
fi

total=$(( $(date +%s) - total_start ))
echo "----------------------------------------------------------------"
if [ -n "$failures" ]; then
    echo "FAILED (${total}s):$failures"
    exit 1
fi
if [ "$skipped" -gt 0 ]; then
    echo "OK (${total}s, $skipped slow suite(s) skipped; run without FAST=1 before merging)"
else
    echo "OK (${total}s)"
fi

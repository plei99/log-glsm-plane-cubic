#!/bin/sh
# Run every regression suite with per-suite timing.
#
#   ./run_tests.sh            run all suites
#   FAST=1 ./run_tests.sh     skip the slow localization/orchestration suites
#
# Every suite prints its own pass line; this script adds wall-clock timing,
# a summary, and a nonzero exit code if anything failed.  DOT_SAGE defaults
# to the user environment; export it first to redirect Sage's profile.

set -u

# The suites measurably slower than ~10s (2026-08-25 full run: 168s total;
# givental_teleman 51s, hodge_integrals 34s, o3_fixed_locus_graphs 18s;
# everything else under 10s).
SLOW_SUITES="
test_givental_teleman
test_hodge_integrals
test_o3_fixed_locus_graphs
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

for suite_file in test_*.sage; do
    suite=${suite_file%.sage}
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
if output=$(python3 -m unittest -q test_elliptic_cubic_gw 2>&1); then
    elapsed=$(( $(date +%s) - start ))
    printf '%-46s %4ss  %s\n' "test_elliptic_cubic_gw" "$elapsed" \
        "python unittest passed"
else
    elapsed=$(( $(date +%s) - start ))
    printf '%-46s %4ss  FAILED\n' "test_elliptic_cubic_gw" "$elapsed"
    printf '%s\n' "$output" | tail -12 | sed 's/^/    /'
    failures="$failures test_elliptic_cubic_gw"
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

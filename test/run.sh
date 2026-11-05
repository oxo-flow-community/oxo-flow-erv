#!/usr/bin/env bash
# Acceptance test for the oxo-flow-erv port.
# Usage: bash test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run, paired-end (default)"
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-erv-dryrun-pe-$$.txt 2>&1
grep -qE "would (execute|run)" /tmp/oxo-erv-dryrun-pe-$$.txt

echo "==> dry-run, single-end"
"$OXO" dry-run main.oxoflow --samples first:1 layout=se > /tmp/oxo-erv-dryrun-se-$$.txt 2>&1
grep -qE "would (execute|run)" /tmp/oxo-erv-dryrun-se-$$.txt

echo "==> the two layouts select different alignment rules"
# every rule appears in the plan, but the inactive variant must be skipped.
# Assert on the run/skip marker, not on the reason text or the up-to-date
# state, so the check also passes when a results/ tree already exists.
grep -qE "fastp_pe.*\[(run|skip):" /tmp/oxo-erv-dryrun-pe-$$.txt
grep -q "fastp_se.*\[skip" /tmp/oxo-erv-dryrun-pe-$$.txt
! grep -q "fastp_se.*\[run:" /tmp/oxo-erv-dryrun-pe-$$.txt
grep -qE "fastp_se.*\[(run|skip):" /tmp/oxo-erv-dryrun-se-$$.txt
grep -q "fastp_pe.*\[skip" /tmp/oxo-erv-dryrun-se-$$.txt
! grep -q "fastp_pe.*\[run:" /tmp/oxo-erv-dryrun-se-$$.txt

echo "==> debug: expanded commands contain no literal placeholders"
OUT=$("$OXO" debug main.oxoflow 2>&1)
echo "$OUT" | grep -q '{sample}' && { echo "unexpanded {sample} in debug output"; exit 1; } || true
echo "$OUT" | grep -q '{config\.' && { echo "unexpanded {config.*} in debug output"; exit 1; } || true

echo "==> fixture scripts run"
python3 scripts/build_locus_from_rmsk.py test/fixtures/reference/rmsk.txt \
    "ERV1,ERVK,ERVL,ERVL-MaLR" /tmp/erv-locus-$$.gtf /tmp/erv-map-$$.tsv \
    /tmp/erv-exon-$$.gtf "^chr([0-9]+|X|Y)$" 2>/dev/null
test "$(($(wc -l < /tmp/erv-map-$$.tsv) - 1))" -eq 10
python3 scripts/make_squire_clean.py /tmp/erv-map-$$.tsv /tmp/erv-squire-$$ 2>/dev/null
python3 scripts/annotate_ervmap.py <(sort -k1,1V -k2,2n test/fixtures/reference/ERVmap.bed) \
    /tmp/erv-map-$$.tsv /tmp/erv-annot-$$.tsv 2>/dev/null
test "$(($(wc -l < /tmp/erv-annot-$$.tsv) - 1))" -eq 4
rm -rf /tmp/erv-locus-$$.gtf /tmp/erv-map-$$.tsv /tmp/erv-exon-$$.gtf \
    /tmp/erv-squire-$$ /tmp/erv-annot-$$.tsv /tmp/oxo-erv-dryrun-*-$$.txt

echo "PASS"

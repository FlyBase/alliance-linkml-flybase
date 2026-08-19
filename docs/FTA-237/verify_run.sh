#!/usr/bin/env bash
# FTA-237 run verification.
# Usage (on flysql26): verify_run.sh /data/alliance/PERSISTENT_main_FTA-237_balancer_rename
# Or from a laptop:    ssh flysql26 "bash -s" < docs/FTA-237/verify_run.sh <rundir>
set -euo pipefail
DIR="${1:?usage: verify_run.sh <run directory>}"
LOG=$(ls "$DIR"/allele_curation_*.log)
TSV=$(ls "$DIR"/allele_curation_*.tsv | grep -v -e _notes -e _renames -e _merges -e _failures)
RENAMES=$(ls "$DIR"/*_balancer_renames.tsv)

echo '--- rename tallies from the log (expect 24 flags, 24 + 1 hard-coded = 25) ---'
grep -E 'balancer rename flags|Renamed .* aberrations from curated flags|hard-coded rename of' "$LOG"

echo '--- FBba retrieval (the nested BalancerHandler must have run) ---'
grep -E 'obtained [0-9]+ FBba balancers' "$LOG"

echo '--- errors (expect 0) ---'
grep -c ' ERROR ' "$LOG" || true
grep ' ERROR ' "$LOG" || true

echo '--- name-collision warnings (expect 0 of each) ---'
echo -n 'many current symbols: '
grep -c 'many current symbols' "$LOG" || true
echo -n 'many current full_names: '
grep -c 'many current full_names' "$LOG" || true

echo '--- rename TSV row count (expect 25) ---'
tail -n +2 "$RENAMES" | wc -l

echo '--- exactly one hard-coded row, and it must be FBab0004818 ---'
awk -F'\t' 'NR>1 && $5=="hard-coded" {print}' "$RENAMES"

echo '--- worked example: In(1)Basc exports as Basc / Bar-apricot-scute ---'
grep -m1 'FB:FBab0004219' "$TSV"

echo '--- SM6 case: symbol SM6 in the primary TSV ---'
grep -m1 'FB:FBab0004818' "$TSV"

echo '--- every renamed aberration keeps its old symbol as a synonym ---'
missing=0
while IFS=$'\t' read -r fbab old; do
  line=$(grep -m1 "^FB:$fbab" "$TSV" || true)
  if [ -z "$line" ]; then
    echo "NOT IN PRIMARY TSV: $fbab"
    missing=$((missing + 1))
  elif ! printf '%s' "$line" | grep -qF "$old"; then
    echo "MISSING old symbol synonym: $fbab ($old)"
    missing=$((missing + 1))
  fi
done < <(awk -F'\t' 'NR>1 {print $1"\t"$7}' "$RENAMES" | sed 's/^FB://')
echo "old-symbol synonym problems: $missing"

echo '--- new symbol must be the exported symbol, not just a synonym ---'
while IFS=$'\t' read -r fbab new; do
  symbol=$(grep -m1 "^FB:$fbab" "$TSV" | cut -f2)
  [ "$symbol" = "$new" ] || echo "SYMBOL MISMATCH: $fbab exported as '$symbol', expected '$new'"
done < <(awk -F'\t' 'NR>1 {print $1"\t"$6}' "$RENAMES" | sed 's/^FB://')

echo '--- done ---'

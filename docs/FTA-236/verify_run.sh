#!/usr/bin/env bash
# FTA-236 run verification.
# Usage (on flysql26): verify_run.sh /data/alliance/PERSISTENT_main_FTA-236_balancer_merge
# Or from a laptop:    ssh flysql26 "bash -s" < docs/FTA-236/verify_run.sh <rundir>
#
# Expected counts come from docs/FTA-236/expected_counts.sql (verified 2026-08-14):
#   38 merge flags, 37 parents, 96 symbol synonyms, 35 fullname synonyms, 27 comments,
#   68 internal notes, 440 references, 326 carried_on rows of which 53 are held back.
# Until the AM1 (FBba0000688) note correction loads (week of 2026-08-17), expect 37 merges
# and exactly one ERROR naming FBba0000688 - that one is expected, not a regression.
set -euo pipefail
DIR="${1:?usage: verify_run.sh <run directory>}"
LOG=$(ls "$DIR"/allele_curation_*.log)
TSV=$(ls "$DIR"/allele_curation_*.tsv | grep -v -e _notes -e _renames -e _merges -e _failures)
MERGES=$(ls "$DIR"/*_balancer_merges.tsv)

echo '--- FBba retrieval and merge tallies ---'
grep -E 'obtained [0-9]+ FBba balancers|balancer merge flags|Merged .* balancers into|Moved [0-9]+ (synonyms|secondary_ids|comments|internal_notes|references)' "$LOG"

echo '--- carries move (only present when ADD_ALLELE_ALLELE_ASSOC=YES) ---'
grep -E 'balancer "carries" associations|Held back' "$LOG" || echo '(gate off: no carries lines)'

echo '--- errors: expect 0, or exactly 1 naming FBba0000688 before the note correction loads ---'
grep -c ' ERROR ' "$LOG" || true
grep ' ERROR ' "$LOG" || true

echo '--- name-collision warnings (expect 0 of each; this is what demoted_synonym_ids prevents) ---'
echo -n 'many current symbols: '
grep -c 'many current symbols' "$LOG" || true
echo -n 'many current full_names: '
grep -c 'many current full_names' "$LOG" || true

echo '--- merge TSV row count (expect 38, or 37 before the correction) ---'
tail -n +2 "$MERGES" | wc -l

echo '--- column totals moved, to compare with expected_counts.sql ---'
awk -F'\t' 'NR>1 {syn+=$6; sec+=$7; com+=$8; note+=$9; ref+=$10; car+=$11; exc+=$12}
            END {printf "synonyms=%d secondary_ids=%d comments=%d internal_notes=%d references=%d carries=%d excluded=%d\n",
                 syn, sec, com, note, ref, car, exc}' "$MERGES"

echo '--- In(2LR)SM6 must show two merge rows (SM6a and SM6b) ---'
grep -c 'FBab0004818' "$MERGES"

echo '--- the three excluded balancers must move 0 carries and hold back 15/20/18 ---'
awk -F'\t' 'NR>1 && ($3=="FBba0000011" || $3=="FBba0000039" || $3=="FBba0000040") {print $1"\t"$3"\tcarries="$11"\texcluded="$12}' "$MERGES"

echo '--- worked example: In(1)Basc must gain M5 as a synonym and FBba0000014 as a secondary ID ---'
grep -m1 'FB:FBab0004219' "$TSV"

echo '--- the N[opa33b] comment must have moved to the aberration ---'
grep -c 'FBab0004219.*opa33b' "$DIR"/allele_curation_*_notes.tsv || echo '(not found - check the notes TSV)'

echo '--- every merged parent must carry the balancer ID as a secondary ID ---'
missing=0
while IFS=$'\t' read -r fbab fbba; do
  line=$(grep -m1 "^FB:$fbab" "$TSV" || true)
  if [ -z "$line" ]; then
    echo "NOT IN PRIMARY TSV: $fbab"
    missing=$((missing + 1))
  elif ! printf '%s' "$line" | grep -qF "$fbba"; then
    echo "MISSING secondary ID: $fbab should list $fbba"
    missing=$((missing + 1))
  fi
done < <(awk -F'\t' 'NR>1 {print $1"\t"$3}' "$MERGES" | sed 's/^FB://')
echo "secondary-ID problems: $missing"

echo '--- done ---'

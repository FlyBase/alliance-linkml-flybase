-- FTA-236 expectation queries. Run on flysql26:
--   PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado -f expected_counts.sql
-- The numbers these produce are the table in docs/FTA-236/plan.md ("Verified input data").
-- NB psql here runs with standard_conforming_strings = off, so set it before using \y in a regex.
set standard_conforming_strings = on;

\echo === 1. FTA flag families on FBba (expect: merge 38, symbol/fullname 24) ===
select case when fp.value like 'FTA: Balancer - merge with parent%' then 'merge with parent (FTA-236)'
            when fp.value like 'FTA: Balancer - use balancer symbol and fullname%' then 'use symbol/fullname (FTA-237)'
            else 'other FTA: flag' end as flag,
       count(*) as props, count(distinct f.feature_id) as fbba
from feature f
join featureprop fp on fp.feature_id = f.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
where f.uniquename like 'FBba%' and cvt.name = 'internal_notes' and fp.value like 'FTA:%'
group by 1 order by 2 desc;

\echo === 2. merge flags vs parents vs the FTA-235 is_balancer set (expect 38 / 37 / 36 / 0) ===
with merges as (
  select f.uniquename as fbba, substring(fp.value from 'FBab[0-9]+') as parent_fbab
  from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
),
flagged as (
  select f.uniquename as fbab
  from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - mark this aberration%'
)
select count(*) as merge_flags,
       count(distinct parent_fbab) as distinct_parents,
       count(distinct parent_fbab) filter (where parent_fbab in (select fbab from flagged)) as parents_also_is_balancer,
       count(*) filter (where parent_fbab is null) as unparsed_parent
from merges;

\echo === 3. parents named by more than one balancer (expect only FBab0004818, SM6a + SM6b) ===
select substring(fp.value from 'FBab[0-9]+') as parent_fbab, count(*) as n_balancers,
       string_agg(f.uniquename, ', ' order by f.uniquename) as fbba
from feature f
join featureprop fp on fp.feature_id = f.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
group by 1 having count(*) > 1 order by 2 desc;

\echo === 4. parents named by a merge flag that are obsolete or absent (expect FBab0007127 until fixed) ===
select b.uniquename as fbba, substring(fp.value from 'FBab[0-9]+') as parent_fbab,
       p.is_obsolete as parent_is_obsolete
from feature b
join featureprop fp on fp.feature_id = b.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
left join feature p on p.uniquename = substring(fp.value from 'FBab[0-9]+')
where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
  and (p.feature_id is null or p.is_obsolete is true)
order by 1;

\echo === 5. volume to move per source field (expect 96 / 35 / 10 / 326 / 10 / 27 / 68 / 440) ===
with fbba as (
  select distinct f.feature_id, f.uniquename
  from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
)
select 'symbol synonyms (AB1b)' as source, count(*) as rows, count(distinct b.feature_id) as fbba_with_data
  from fbba b
  join feature_synonym fs on fs.feature_id = b.feature_id
  join synonym s on s.synonym_id = fs.synonym_id
  join cvterm t on t.cvterm_id = s.type_id
  where t.name = 'symbol'
union all
select 'fullname synonyms (AB2b)', count(*), count(distinct b.feature_id)
  from fbba b
  join feature_synonym fs on fs.feature_id = b.feature_id
  join synonym s on s.synonym_id = fs.synonym_id
  join cvterm t on t.cvterm_id = s.type_id
  where t.name = 'fullname'
union all
select 'secondary IDs (feature_dbxref FBba)', count(*), count(distinct b.feature_id)
  from fbba b
  join feature_dbxref fd on fd.feature_id = b.feature_id
  join dbxref dx on dx.dbxref_id = fd.dbxref_id
  join db on db.db_id = dx.db_id
  where db.name = 'FlyBase' and dx.accession like 'FBba%'
union all
select 'carries alleles (AB5b, carried_on, FBba as object)', count(*), count(distinct b.feature_id)
  from fbba b
  join feature_relationship fr on fr.object_id = b.feature_id
  join cvterm t on t.cvterm_id = fr.type_id
  where t.name = 'carried_on'
union all
select 'transposon insertions (AB5a, associated_with, FBba as subject)', count(*), count(distinct b.feature_id)
  from fbba b
  join feature_relationship fr on fr.subject_id = b.feature_id
  join cvterm t on t.cvterm_id = fr.type_id
  where t.name = 'associated_with'
union all
select 'other info (AB6, misc prop)', count(*), count(distinct b.feature_id)
  from fbba b
  join featureprop fp on fp.feature_id = b.feature_id
  join cvterm t on t.cvterm_id = fp.type_id
  where t.name = 'misc'
union all
select 'internal_notes props (incl. the merge flag itself)', count(*), count(distinct b.feature_id)
  from fbba b
  join featureprop fp on fp.feature_id = b.feature_id
  join cvterm t on t.cvterm_id = fp.type_id
  where t.name = 'internal_notes'
union all
select 'references (feature_pub)', count(*), count(distinct b.feature_id)
  from fbba b
  join feature_pub fpb on fpb.feature_id = b.feature_id
order by 1;

\echo === 6. carried_on rows total vs those held back by the 3 exclusions (expect 326 / 54 / 53) ===
with merges as (
  select f.feature_id as fbba_id, f.uniquename as fbba
  from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
)
select count(distinct fr.feature_relationship_id) as carried_on_rows,
       count(distinct fr.subject_id) as distinct_alleles,
       count(distinct fr.feature_relationship_id) filter (where m.fbba in ('FBba0000011', 'FBba0000039', 'FBba0000040')) as excluded_rows
from merges m
join feature_relationship fr on fr.object_id = m.fbba_id
join cvterm t on t.cvterm_id = fr.type_id and t.name = 'carried_on';

\echo === 7. the worked example from the ticket: Basc (FBba0000014) -> In(1)Basc (FBab0004219) ===
select 'symbol synonym' as datum, s.name as value
  from feature f
  join feature_synonym fs on fs.feature_id = f.feature_id
  join synonym s on s.synonym_id = fs.synonym_id
  join cvterm t on t.cvterm_id = s.type_id
  where f.uniquename = 'FBba0000014' and t.name = 'symbol'
union all
select 'carried allele', a.name
  from feature b
  join feature_relationship fr on fr.object_id = b.feature_id
  join cvterm t on t.cvterm_id = fr.type_id
  join feature a on a.feature_id = fr.subject_id
  where b.uniquename = 'FBba0000014' and t.name = 'carried_on'
union all
select 'comment (AB6)', fp.value
  from feature b
  join featureprop fp on fp.feature_id = b.feature_id
  join cvterm t on t.cvterm_id = fp.type_id
  where b.uniquename = 'FBba0000014' and t.name = 'misc'
order by 1, 2;

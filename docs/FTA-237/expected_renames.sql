-- FTA-237 expectation queries. Run on flysql26:
--   PGPASSWORD=ilongdenpgsql psql -h flysql24 -U ilongden -d production_chado -f expected_renames.sql
-- Query 1 is the oracle for the *_balancer_renames.tsv rows with source=flag; the curator spreadsheet
-- linked from the ticket needs auth, so this is the machine-checkable equivalent.
set standard_conforming_strings = on;

\echo === 1. the 24 renames: new names (from the balancer) and old names (from the parent) ===
with renames as (
  select f.feature_id as fbba_id, f.uniquename as fbba,
         substring(fp.value from 'FBab[0-9]+') as parent_fbab
  from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - use balancer symbol and fullname%'
),
curr as (
  select fs.feature_id, t.name as syn_type, s.name as value
  from feature_synonym fs
  join synonym s on s.synonym_id = fs.synonym_id
  join cvterm t on t.cvterm_id = s.type_id
  where fs.is_current is true and t.name in ('symbol', 'fullname')
)
select r.parent_fbab,
       r.fbba,
       (select string_agg(distinct value, ' / ') from curr where feature_id = r.fbba_id and syn_type = 'symbol') as new_symbol,
       (select string_agg(distinct value, ' / ') from curr where feature_id = r.fbba_id and syn_type = 'fullname') as new_full_name,
       (select string_agg(distinct value, ' / ') from curr where feature_id = p.feature_id and syn_type = 'symbol') as old_symbol,
       coalesce((select string_agg(distinct value, ' / ') from curr where feature_id = p.feature_id and syn_type = 'fullname'), '') as old_full_name
from renames r
left join feature p on p.uniquename = r.parent_fbab
order by r.parent_fbab;

\echo === 2. every rename-flagged balancer must have both a current symbol and a current fullname (expect 24 / 24 / 24 / 0) ===
with renames as (
  select f.feature_id as fbba_id from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - use balancer symbol and fullname%'
)
select count(*) as balancers,
       count(*) filter (where has_symbol) as with_current_symbol,
       count(*) filter (where has_fullname) as with_current_fullname,
       count(*) filter (where not has_fullname) as missing_current_fullname
from (
  select r.fbba_id,
         exists (select 1 from feature_synonym fs
                 join synonym s on s.synonym_id = fs.synonym_id
                 join cvterm t on t.cvterm_id = s.type_id
                 where fs.feature_id = r.fbba_id and fs.is_current and t.name = 'symbol') as has_symbol,
         exists (select 1 from feature_synonym fs
                 join synonym s on s.synonym_id = fs.synonym_id
                 join cvterm t on t.cvterm_id = s.type_id
                 where fs.feature_id = r.fbba_id and fs.is_current and t.name = 'fullname') as has_fullname
  from renames r
) x;

\echo === 3. no parent may be named by two rename flags, and none may be obsolete/absent (both expect 0 rows) ===
select substring(fp.value from 'FBab[0-9]+') as parent_fbab, count(*) as n_flags
from feature f
join featureprop fp on fp.feature_id = f.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - use balancer symbol and fullname%'
group by 1 having count(*) > 1;

select b.uniquename as fbba, substring(fp.value from 'FBab[0-9]+') as parent_fbab, p.is_obsolete
from feature b
join featureprop fp on fp.feature_id = b.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
left join feature p on p.uniquename = substring(fp.value from 'FBab[0-9]+')
where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - use balancer symbol and fullname%'
  and (p.feature_id is null or p.is_obsolete is true);

\echo === 4. the SM6 case: no rename flag on either balancer (expect 0 rows) ===
select f.uniquename, fp.value
from feature f
join featureprop fp on fp.feature_id = f.feature_id
join cvterm cvt on cvt.cvterm_id = fp.type_id
where f.uniquename in ('FBba0000039', 'FBba0000040')
  and cvt.name = 'internal_notes'
  and fp.value like 'FTA: Balancer - use balancer symbol and fullname%';

\echo === 5. the SM6 target strings: "SM6" only non-current, "Second Multiple 6" absent entirely ===
select s.name as value, f.uniquename, t.name as syn_type, fs.is_current
from synonym s
join cvterm t on t.cvterm_id = s.type_id
join feature_synonym fs on fs.synonym_id = s.synonym_id
join feature f on f.feature_id = fs.feature_id
where s.name in ('SM6', 'Second Multiple 6')
  and (f.uniquename like 'FBab%' or f.uniquename like 'FBba%')
order by 1, 2, 3;

\echo === 6. current names of the SM6 trio (parent has no current fullname) ===
select f.uniquename, t.name as syn_type, s.name as value
from feature f
join feature_synonym fs on fs.feature_id = f.feature_id
join synonym s on s.synonym_id = fs.synonym_id
join cvterm t on t.cvterm_id = s.type_id
where f.uniquename in ('FBab0004818', 'FBba0000039', 'FBba0000040')
  and t.name in ('symbol', 'fullname') and fs.is_current is true
group by 1, 2, 3
order by 1, 2, 3;

\echo === 7. rename flags are a strict subset of the merge flags (expect 24 / 38 / 24 / 0) ===
with renames as (
  select distinct f.feature_id from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - use balancer symbol and fullname%'
),
merges as (
  select distinct f.feature_id from feature f
  join featureprop fp on fp.feature_id = f.feature_id
  join cvterm cvt on cvt.cvterm_id = fp.type_id
  where cvt.name = 'internal_notes' and fp.value like 'FTA: Balancer - merge with parent%'
)
select (select count(*) from renames) as rename_flagged,
       (select count(*) from merges) as merge_flagged,
       (select count(*) from renames r where exists (select 1 from merges m where m.feature_id = r.feature_id)) as rename_also_merge,
       (select count(*) from renames r where not exists (select 1 from merges m where m.feature_id = r.feature_id)) as rename_only;

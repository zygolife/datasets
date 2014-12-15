.header on
.mode csv
select genus, species, count(*) from isolate group by genus,species;

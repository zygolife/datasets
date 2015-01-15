.header on
.mode csv
.output "query_results/genus_species_count.txt"
select genus, species, count(*) from isolate group by genus,species;
.output "query_results/genus_count.txt"
select genus, count(*) as count from isolate group by genus;
.output "query_results/type_strain.txt"
select genus,species,subspecies,collection,strain_id from isolate where type_strain='TRUE';

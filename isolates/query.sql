.header on
.mode csv
.output "query_results/genus_species_count.txt"
select genus, species, count(*) from isolate group by genus,species;

.output "query_results/genus_count.txt"
select genus, count(*) as count from isolate group by genus;

.output "query_results/type_strain.txt"
select genus,species,subspecies,collection,strain_id from isolate where type_strain='TRUE';

.output "query_results/genus_genomes_sampled.txt"
select isolate.genus, count(genomes.status) > 0 as any_genomes_sampled
FROM isolate LEFT OUTER JOIN genomes ON isolate.genus=genomes.genus group by isolate.genus;

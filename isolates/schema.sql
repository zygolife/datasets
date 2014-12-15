
CREATE TABLE isolate(
 --      isolate_primary_id	INTEGER	PRIMARY KEY AUTO_INCREMENT,
       strain_id		text,
       collection		text,
       genus			text,
       species			text,
       subvar_type		text, -- f. or var. or f.sp.
       subspecies		text,
       type_strain		integer,
       date_accessioned		text,
       source_locale		text,
       received_as		text,
       substrate		text,
       location_detail		text,
       country			text,
       nrrl_catalog		integer,
       nrrl_dead		integer,
       nrrl_lyo			integer
);

CREATE UNIQUE INDEX isolate_strain_idx ON isolate (strain_id, collection);
CREATE INDEX isolate_genus ON isolate (genus);
CREATE INDEX isolate_species ON isolate (species);

 

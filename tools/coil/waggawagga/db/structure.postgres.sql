CREATE TABLE public.genome_evaluations (
  id INTEGER PRIMARY KEY NOT NULL DEFAULT nextval('genome_evaluations_id_seq'::regclass),
  genome CHARACTER VARYING(255),
  entry_id CHARACTER VARYING(255),
  header TEXT,
  sequence TEXT,
  scores TEXT,
  "imagePath" CHARACTER VARYING(255),
  "windowSize" INTEGER,
  "segmentSize" INTEGER,
  created_at TIMESTAMP WITHOUT TIME ZONE NOT NULL,
  updated_at TIMESTAMP WITHOUT TIME ZONE NOT NULL,
  tool CHARACTER VARYING(255),
  taxonomy CHARACTER VARYING(1024),
  class_type CHARACTER VARYING(1024)
);

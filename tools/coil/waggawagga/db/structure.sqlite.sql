DROP TABLE IF EXISTS "genome_evaluations";
-- SQLite specific command to free deleted data
VACUUM;
CREATE TABLE "genome_evaluations" (
-- Rely on implicit SQLite col 'rowid' for AUTOINCREMENT
-- id INTEGER PRIMARY KEY AUTOINCREMENT,
  id INTEGER PRIMARY KEY,
  genome VARCHAR,
  entry_id VARCHAR,
  header TEXT,
  sequence TEXT,
  scores TEXT,
  "imagePath" VARCHAR,
  "windowSize" INT4,
  "segmentSize" INT4,
  tool VARCHAR,
  taxonomy VARCHAR,
  class_type VARCHAR,
  created_at DATETIME DEFAULT (DATETIME(CURRENT_TIMESTAMP, 'LOCALTIME')),
  updated_at DATETIME DEFAULT (DATETIME(CURRENT_TIMESTAMP, 'LOCALTIME'))
);
-- ### Notes
-- INSERT INTO X.TABLE(fieldname1, fieldname2) SELECT fieldname1, fieldname2 FROM Y.TABLE;
-- ### Merge SQLite3 databases
-- id,genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at
-- sqlite> attach 'development-1.sqlite3' as toMerge;
-- sqlite> INSERT INTO genome_evaluations SELECT * FROM toMerge.genome_evaluations;
-- sqlite> detach database toMerge;
-- sqlite> attach 'development-2.sqlite3' as toMerge;
-- sqlite> INSERT INTO genome_evaluations (genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at) SELECT genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at FROM toMerge.genome_evaluations;
-- sqlite> detach database toMerge;
-- sqlite> attach 'development-3.sqlite3' as toMerge;
-- sqlite> INSERT INTO genome_evaluations (genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at) SELECT genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at FROM toMerge.genome_evaluations;
-- sqlite> detach database toMerge;
-- sqlite> attach 'development-4.sqlite3' as toMerge;
-- sqlite> INSERT INTO genome_evaluations (genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at) SELECT genome,entry_id,header,sequence,scores,"imagePath","windowSize","segmentSize",tool,taxonomy,class_type,created_at,updated_at FROM toMerge.genome_evaluations;
-- sqlite> detach database toMerge;
-- sqlite> VACUUM;
-- Export of Data
-- .mode csv
-- .output seqlist_Arabidopsis_thaliana.TAIR10.pep.abinitio.csv
-- SELECT entry_id, sequence FROM genome_evaluations WHERE id IN (...)
-- .output
-- ### Export full database to MySQL
-- sqlite3 development-1.sqlite3 .dump > development-1.sqlite3.dump.sql
-- sqlite3 development-1.sqlite3 .dump | python sqlite_to_mysql_convert.py > dumped_data.sql
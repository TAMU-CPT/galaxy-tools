DROP TABLE IF EXISTS genome_evaluations;
CREATE TABLE IF NOT EXISTS `genome_evaluations` (
  id INTEGER PRIMARY KEY,
  genome VARCHAR(255),
  entry_id VARCHAR(255),
  gene VARCHAR(255)
  header TEXT,
  sequence TEXT,
  scores MEDIUMTEXT,
  `imagePath` VARCHAR(255),
  `windowSize` TINYINT,
  `segmentSize` TINYINT,
  tool VARCHAR(255),
  taxonomy VARCHAR(255),
  class_type VARCHAR(255),
  created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
  updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
);

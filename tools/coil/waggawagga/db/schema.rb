# encoding: UTF-8
# This file is auto-generated from the current state of the database. Instead
# of editing this file, please use the migrations feature of Active Record to
# incrementally modify your database, and then regenerate this schema definition.
#
# Note that this schema.rb definition is the authoritative source for your
# database schema. If you need to create the application database on another
# system, you should be using db:schema:load, not running all the migrations
# from scratch. The latter is a flawed and unsustainable approach (the more migrations
# you'll amass, the slower it'll run and the greater likelihood for issues).
#
# It's strongly recommended that you check this file into your version control system.

ActiveRecord::Schema.define(version: 20160523132506) do

  # These are extensions that must be enabled in order to support this database
  enable_extension "plpgsql"

  create_table "coiled_coil_components", force: true do |t|
    t.integer  "number"
    t.string   "chain"
    t.integer  "start"
    t.integer  "end"
    t.text     "register"
    t.string   "coilType"
    t.integer  "coiledCoil_id"
    t.integer  "pdbSequence_id"
    t.integer  "proteinSecondaryStructure_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "coiled_coil_components", ["coiledCoil_id"], name: "index_coiled_coil_components_on_coiledCoil_id", using: :btree
  add_index "coiled_coil_components", ["pdbSequence_id"], name: "index_coiled_coil_components_on_pdbSequence_id", using: :btree
  add_index "coiled_coil_components", ["proteinSecondaryStructure_id"], name: "index_coiled_coil_components_on_proteinSecondaryStructure_id", using: :btree

  create_table "coiled_coil_prediction_evaluations", force: true do |t|
    t.string   "tool"
    t.string   "parameters"
    t.integer  "tp"
    t.integer  "fp"
    t.integer  "tn"
    t.integer  "fn"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.integer  "pdbSequence_id"
  end

  add_index "coiled_coil_prediction_evaluations", ["pdbSequence_id"], name: "index_coiled_coil_prediction_evaluations_on_pdbSequence_id", using: :btree

  create_table "coiled_coil_predictions", force: true do |t|
    t.string   "tool"
    t.string   "parameters"
    t.text     "register"
    t.integer  "start"
    t.integer  "end"
    t.integer  "pdbSequence_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.float    "threshold"
    t.string   "version",        limit: 1024
  end

  add_index "coiled_coil_predictions", ["pdbSequence_id"], name: "index_coiled_coil_predictions_on_pdbSequence_id", using: :btree

  create_table "coiled_coils", force: true do |t|
    t.string   "coilType"
    t.integer  "pdb_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "coiled_coils", ["pdb_id"], name: "index_coiled_coils_on_pdb_id", using: :btree

  create_table "genome_evaluations", force: true do |t|
    t.string   "genome"
    t.string   "entry_id"
    t.text     "header"
    t.text     "sequence"
    t.text     "scores"
    t.string   "imagePath"
    t.integer  "windowSize"
    t.integer  "segmentSize"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.string   "tool"
    t.string   "taxonomy",    limit: 1024
  end

  create_table "matches", force: true do |t|
    t.integer  "coiled_coil_prediction_id"
    t.integer  "coiled_coil_component_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.float    "helixOverlap"
    t.float    "ccOverlap"
  end

  add_index "matches", ["coiled_coil_component_id"], name: "index_matches_on_coiled_coil_component_id", using: :btree
  add_index "matches", ["coiled_coil_prediction_id"], name: "index_matches_on_coiled_coil_prediction_id", using: :btree

  create_table "oligomerisation_predictions", force: true do |t|
    t.integer  "type"
    t.float    "score"
    t.string   "tool"
    t.string   "parameters"
    t.integer  "coiledCoilPrediction_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "oligomerisation_predictions", ["coiledCoilPrediction_id"], name: "index_oligomerisation_predictions_on_coiledCoilPrediction_id", using: :btree

  create_table "pdb_chain_components", force: true do |t|
    t.integer  "seqBegin"
    t.integer  "seqEnd"
    t.integer  "dbSeqBegin"
    t.integer  "dbSeqEnd"
    t.integer  "Protein_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.integer  "dbRef_id"
    t.integer  "pdbChain_id"
  end

  add_index "pdb_chain_components", ["Protein_id"], name: "index_pdb_chain_components_on_Protein_id", using: :btree
  add_index "pdb_chain_components", ["dbRef_id"], name: "index_pdb_chain_components_on_dbRef_id", using: :btree
  add_index "pdb_chain_components", ["pdbChain_id"], name: "index_pdb_chain_components_on_pdbChain_id", using: :btree

  create_table "pdb_chains", force: true do |t|
    t.string   "chain"
    t.integer  "pdbSequence_id"
    t.integer  "pdb_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.text     "dbrefSequence"
    t.integer  "incomplete"
  end

  add_index "pdb_chains", ["pdbSequence_id"], name: "index_pdb_chains_on_pdbSequence_id", using: :btree
  add_index "pdb_chains", ["pdb_id"], name: "index_pdb_chains_on_pdb_id", using: :btree

  create_table "pdb_import_errors", force: true do |t|
    t.string   "pdb"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "pdb_secondary_structures", force: true do |t|
    t.integer  "pdbSequence_id"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.string   "basicType"
    t.integer  "serNum"
    t.string   "structureId"
    t.string   "initResName"
    t.integer  "initSeqNum"
    t.string   "endResName"
    t.integer  "endSeqNum"
    t.integer  "structureType"
    t.integer  "sense"
    t.integer  "length"
    t.integer  "pdbChain_id"
  end

  add_index "pdb_secondary_structures", ["pdbChain_id"], name: "index_pdb_secondary_structures_on_pdbChain_id", using: :btree
  add_index "pdb_secondary_structures", ["pdbSequence_id"], name: "index_pdb_secondary_structures_on_pdbSequence_id", using: :btree

  create_table "pdb_sequences", force: true do |t|
    t.text     "sequence"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  create_table "pdbs", force: true do |t|
    t.string   "abbr"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.text     "definition"
  end

  add_index "pdbs", ["abbr"], name: "index_pdbs_on_abbr", unique: true, using: :btree

  create_table "protein_secondary_structures", force: true do |t|
    t.integer  "start"
    t.integer  "end"
    t.integer  "structureType"
    t.integer  "pdbSequence_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "protein_secondary_structures", ["pdbSequence_id"], name: "index_protein_secondary_structures_on_pdbSequence_id", using: :btree

  create_table "proteins", force: true do |t|
    t.string   "abbr"
    t.text     "sequence"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.string   "database"
    t.text     "definition"
  end

  add_index "proteins", ["abbr"], name: "index_proteins_on_abbr", unique: true, using: :btree

  create_table "sah_prediction_evaluations", force: true do |t|
    t.integer  "threshold"
    t.integer  "tp"
    t.integer  "fp"
    t.integer  "tn"
    t.integer  "fn"
    t.integer  "coiledCoilPrediction_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "sah_prediction_evaluations", ["coiledCoilPrediction_id"], name: "index_sah_prediction_evaluations_on_coiledCoilPrediction_id", using: :btree

  create_table "sah_score_evaluations", force: true do |t|
    t.integer  "minLenSahAa"
    t.integer  "SahAa"
    t.integer  "SahAaHelix"
    t.integer  "SahAaCc"
    t.integer  "SahAaMisc"
    t.datetime "created_at"
    t.datetime "updated_at"
    t.integer  "paramWindowSize"
    t.float    "paramThreshold"
    t.integer  "paramLenAa"
    t.string   "paramCcPredTool"
    t.integer  "SahAaPredCcHelix"
    t.integer  "SahAaPredCcHelixNotCcSocket"
    t.integer  "pdbSequence_id"
  end

  add_index "sah_score_evaluations", ["pdbSequence_id"], name: "index_sah_score_evaluations_on_pdbSequence_id", using: :btree

  create_table "sah_scores", force: true do |t|
    t.text     "scores"
    t.integer  "windowSize"
    t.string   "status"
    t.integer  "pdbSequence_id"
    t.datetime "created_at"
    t.datetime "updated_at"
  end

  add_index "sah_scores", ["pdbSequence_id"], name: "index_sah_scores_on_pdbSequence_id", using: :btree

end

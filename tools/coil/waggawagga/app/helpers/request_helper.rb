require 'csv'

module RequestHelper

	# SVG: Helper-method to generate n points on a circle (evenly distributed)
	#      Assumes that CSV-rows/cols are alphabetically-ordered (A->Y)
	# @param [String] filepath      Path to scoring-matrix (CSV)
	# @return [Array]
	def self.initialize_scoring_matrix(filepath)
		matrix = {'A' => {},	'C' => {},	'D' => {},	'E' => {},	'F' => {},
		          'G' => {},	'H' => {},	'I' => {},	'K' => {},	'L' => {},
		          'M' => {},	'N' => {},	'P' => {},	'Q' => {},	'R' => {},
		          'S' => {},	'T' => {},	'V' => {},	'W' => {},	'Y' => {}}
		aminoacids = matrix.keys

		# Read and parse CSV
		csv_arr = CSV.read(filepath, {:col_sep => ";"})
	  aminoacids.each_with_index { |aa, row_num|
		  csv_arr[row_num].each_with_index { |score, i| matrix[aa][aminoacids[i]] = score.to_f }
	  }

		return matrix
	end

	def initialize_scoring_matrix(filepath)
		return RequestHelper.initialize_scoring_matrix(filepath)
	end

end

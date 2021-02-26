# **********************************************************************
#  Copyright notice
#
#  (c) 2011-2017 Dominic Simm <dominic.simm@mpibpc.mpg.de>
#  All rights reserved
#
#  This file is part of Waggawagga.
#
#  Waggawagga is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  Waggawagga is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Waggawagga.  If not, see <http://www.gnu.org/licenses/>.
# **********************************************************************

# Changes 2011-2017 by Dominic Simm <dominic.simm@mpibpc.mpg.de>
# See the ChangeLog or git repository for details.

# require 'coiled_coil_sequence'

# Global object for each CLI result
class ShellCommand

	# Global variables $
	$lib_path = "#{Rails.root}/lib"
	$tmp_path = ""; $timestamp = ''
	$sequence = ''
	$oligomerisation_states = { 0 => 'antiparallel dimer', 1 => 'parallel dimer', 2 => 'dimer', 3 => 'trimer', 4 => 'tetramer' }

	# Class variables with automatic getter / setter
	attr_accessor :conf				# array: shell commands and patterns of coiledcoil tools
	attr_accessor :globals			# array: global variabls

	# Initializer for a new instance
	# @param tmp_path [string] 		absolute path to the tmp-directory 
	# @param tmp_filename [string] 	filename of the uploaded data (temporary file)
	# @param timestamp [int] 		  timestamp / jobId of the current request
	# @param window_size [int] 		chosen window-size for the coiledcoil tools
	# @param sequence [string]		submitted FASTA protein sequence
	# @param title [string]		    title of submitted FASTA protein sequence
	# @param id [int]				      identifier of the sequence [0,1,2...], for multiple FASTA-sequences submissions
	# @return void
	def initialize(tmp_path = '', tmp_filename = '', timestamp = '', window_size = '', sequence = '', title = '', id = '')
		
		# Setting globals
		$tmp_path = tmp_path; $timestamp = timestamp; $sequence = sequence; $title = title
		
		# Configuration hash-map
		@conf = {
			# window_size setting over cli-parameter
			:paircoil2 => {
					:path => "#{$lib_path}/paircoil2",
					:command => "export PAIRCOIL2='#{$lib_path}/paircoil2' && export PAIRCOIL_CONFIG=$PAIRCOIL2/.paircoil2_#{window_size} && cd $PAIRCOIL2 && ./paircoil2 -win #{window_size} #{$tmp_path}/#{tmp_filename} #{$tmp_path}/#{timestamp}_run#{id}_paircoil2.txt",
											# 1 M g 0.000 0.00
					:pattern => /([\d]+)[ ]+([A-Z])[ ]+([a-z])[ ]+([\d.]+)[ (]+([-]?[\d.inf]+)/,
          :rows => 5,       # number of rows in text-file (coiledcoil output)
          :max => 1.0,      # maximum value for coiledcoil score
          :tics => 0.1,     # point distance for y-axis in gnuplot
					:limit => 0.025, 	# limit for coiledcoil score = 0.025 | # 18.0
					:gnuId => 4,    	# specifies the coiledcoil-field in the text-file = 4 | 5
					:probId => 3,   	# specifies the coiledcoil-field in the regex-array (caution: some fields may be left out) = 3 | 4
					:comp => 0,     	# value must be compared by less (0) | more (1) to limit = 0 | 1
					:extra => 0,      # text-values must be extrapolated / complemented
					:trimere => 0,    # tool can calculate di-/trimere probabilities
				}, 
			:paircoil => {
			# window_size setting over config-file
					:path => "#{$lib_path}/paircoil",
					:command => "export PAIRCOIL='#{$lib_path}/paircoil' && export PAIRCOIL_CONFIG=$PAIRCOIL/paircoil_config_#{window_size} && cd $PAIRCOIL && ./paircoil #{$tmp_path}/#{tmp_filename} && cp #{$lib_path}/paircoil/result/test.log #{$tmp_path}/#{timestamp}_run#{id}_paircoil.txt",
											# 1 M g 0.000 0.00
					:pattern => /([\d.]+)@[ ]*([\d]+)-[ ]*([\d]+):([a-g])/,
          :rows => 4,       # number of rows in text-file (coiledcoil output)
          :max => 1.0,      # maximum value for coiledcoil score
          :tics => 0.1,     # point distance for y-axis in gnuplot
					:limit => 0.84,		# limit for coiledcoil score
					:gnuId => 4,		  # specifies the coiledcoil-field in the text-file
					:probId => 3,		  # specifies the coiledcoil-field in the regex-array (caution: some fields may be left out)
					:comp => 1,			  # value must be compared by less (0) | more (1) to limit
					:extra => 1,		  # text-values must be extrapolated / complemented
					:trimere => 0,		# tool can calculate di-/trimere probabilities
				}, 
			# window_size setting over cli-parameter
			:ncoils => {
					:path => "#{$lib_path}/ncoils",
					# -m abs_path_to_matrix  # Matrices = {'MTK' => 'old.mat', 'MTIDK' => 'new.mat'}
					# -w 					           # Weighted (a&d = b,c,e,d,g) or not
					# -win [14,21,28]		     # Window size, default 21
					# ./ncoils-osf -m $coilsroot/old.mat -win int -w  < input.file > output.file
					:command => "export COILSDIR='#{$lib_path}/ncoils' && cd $COILSDIR && ./ncoils-osf -m #{$lib_path}/ncoils/old.mat -win #{window_size} < #{$tmp_path}/#{tmp_filename} > #{$tmp_path}/#{timestamp}_run#{id}_ncoils.txt",
					# :pattern => /([\d]+)[ ]+([A-Z])[ ]+([a-z])[ ]+([\d.]+)[ ]+([\d.]+)[ (]+([\d.]+)[ ]+([\d.]+)/,
					:pattern => /([\d]+)[ ]+([A-Z])[ ]+([a-z])[ \d.]+([\d.]+)[ (]+([\d.]+)[ ]+([\d.]+)/,
          :rows => 6,
          :max => 5.0,
          :tics => 0.5,
					:limit => 0.5,
					:gnuId => 5,
					:probId => 4,
					:comp => 1,
					:extra => 0,
					:trimere => 0,
				}, 
			:marcoil => {
			# No window_size setting
					:path => "#{$lib_path}/marcoil",
					# :command => "cd #{$lib_path}/marcoil && ./marcoil +cs #{$tmp_path}/#{tmp_filename} && cp #{$lib_path}/marcoil/Outputs/ProbPerState #{$tmp_path}/marcoil_run_#{timestamp}.txt",
					:command => "cd #{$lib_path}/marcoil && ./marcoil +dls #{$tmp_path}/#{tmp_filename} && cp #{$lib_path}/marcoil/Outputs/ProbList #{$tmp_path}/#{timestamp}_run#{id}_marcoil.txt",
					# Long pattern with heptad position probabilities
					# :pattern => ([\d]+)([A-Z])[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[ ]+CC-ph: ([a-z])
					# Short pattern without heptad position probabilities
					:pattern => /([\d]+)[ ]+([A-Z])[ ]+([\d.]+)[\t ]+([a-g])/,
          :rows => 4,
          :max => 100.0,
          :tics => 10.0,
					:limit => 90.0,
					:gnuId => 4,
					:probId => 3,
					:comp => 1,
					:extra => 1,
					:trimere => 0,
				}, 
			:multicoil => {
			# window_size setting over config-file
					:path => "#{$lib_path}/multicoil",
					# Multicoil needs to be executable a symbolic link 'MULTICOIL' in the user-directory, which points to its absolute path (location)
					:command => "export MULTICOIL_CONFIG='#{$lib_path}/multicoil' && cd $MULTICOIL_CONFIG && ./multicoil #{$tmp_path}/#{tmp_filename} -config multicoil_config_#{window_size} && mv #{$lib_path}/multicoil/TEST_RUNS/#{tmp_filename}.out #{$tmp_path}/#{timestamp}_run#{id}_multicoil.txt && rm #{$lib_path}/multicoil/TEST_RUNS/#{tmp_filename}.*",
                    # 1 M g g g 0.000 0.000 0.000 0.000
					# :pattern => /([\d]+)[ ]+([A-Z])[ ]+([a-z])[ (]+([a-g]),([a-g])[) ]+([\d.]+)[ ]+([\d.]+)[\t ]+([\d.]+)/,
                    # 1 M g 0.000 0.000 0.000
					:pattern => /([\d]+)[ ]+([A-Z])[ ]+([a-z])[ (]+[a-g,]+[) ]+([\d.]+)[ ]+([\d.]+)[\t ]+([\d.]+)/,
          :rows => 6,
          :max => 1.0,
          :tics => 0.1,
					:limit => 0.25,
					:gnuId => 6,
					:probId => 3,
					:comp => 1,
					:extra => 0,
					:trimere => 1,
				},
			# window_size setting over cli-parameter
			:multicoil2 => {
					:path => "#{$lib_path}/multicoil2",
					:command => "export CLASSPATH='#{$lib_path}/multicoil2/bin:#{$lib_path}/multicoil2/bin/commons-cli-1.2.jar' && cd #{$lib_path}/multicoil2/bin && java -cp $CLASSPATH coilPred/CoilPredictor -sfile #{$tmp_path}/#{tmp_filename} -out #{$tmp_path}/#{timestamp}_run#{id}_multicoil2.txt",
					# Short pattern without heptad position probabilities
                    # 1 M g 0.000 0.000 0.000
					:pattern => /([\d]+)[\t ]+([A-Z])[\t ]+([a-z])[\t (]+[a-g],[a-g][)\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)/,
					# Long pattern with heptad position probabilities
					# :pattern => /([\d]+)[\t ]+([A-Z])[\t ]+([a-z])[\t (]+([a-g]),([a-g])[)\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+([\d.]+)[\t ]+[a:]+([\d.]+)[ ]+[b:]+([\d.]+)[ ]+[c:]+([\d.]+)[ ]+[d:]+([\d.]+)[ ]+[e:]+([\d.]+)[ ]+[f:]+([\d.]+)[ ]+[g:]+([\d.]+)[ ]+/
          :rows => 6,
          :max => 1.0,
          :tics => 0.1,
					:limit => 0.25,
					:gnuId => 4,
					:probId => 3,
					:comp => 1,
					:extra => 0,
					:trimere => 1,
				},
			# No window_size setting
			:scorer2 => {
					:path => "#{$lib_path}/scorer2",
					# :command => "cd #{$lib_path}/scorer2 && echo 'scorer2('Inputs/sequences.csv')' | cat scorer2.R - | R --slave",
					:command => "cd #{$lib_path}/scorer2 && cat scorer2.R run.R | R --slave",
					:pattern => /([\d-]+\.[\d]+)/,
					:limit => 0,
					:gnuId => 0,
					:probId => 0,
					:comp => 0,
					:extra => 0,
					:trimere => 1,
				},
			# No window_size setting
			:procoil => {
					:path => "#{$lib_path}/procoil",
					:command => "cd #{$lib_path}/procoil && echo 'procoil(\"Inputs/sequences.csv\", \"###PATH###\", \"###TOOL###\", \"###TIMESTAMP###\")' | cat procoil.R - | R --slave",
					# :command => "cd #{$lib_path}/procoil && cat procoil.R run.R | R --slave",
					:pattern => /([a-z]+):([\d.-]+)/,
					:limit => 0,
					:gnuId => 0,
					:probId => 0,
					:comp => 0,
					:extra => 0,
					:trimere => 1,
				},								
			# No window_size setting
			:logicoil => {
					:path => "#{$lib_path}/logicoil",
					:command => "cd #{$lib_path}/logicoil && cat main.R | R --slave",
					# :command => "cd #{$lib_path}/procoil && cat procoil.R run.R | R --slave",
					:pattern => /"(.+)" "([0-9\.]+)" "([0-9\.]+)" "([0-9\.]+)" "([0-9\.]+)" "([0-9]+)"/,
					:limit => 0,
					:gnuId => 0,
					:probId => 0,
					:comp => 0,
					:extra => 0,
					:trimere => 1,
				},
			# No window_size setting
			:sah => {
          # :path => "#{$lib_path}/marcoil",
          # :command => #"cd #{$lib_path}/marcoil && ./marcoil +dls #{$tmp_path}/#{tmp_filename} && cp #{$lib_path}/marcoil/Outputs/ProbList #{$tmp_path}/#{timestamp}_run#{id}_marcoil.txt",
          # :pattern => /([\d]+)[ ]+([A-Z])[ ]+([\d.]+)[\t ]+([a-g])/,
          :path => '',
          :command => '',
          :pattern => '',
          :rows => 4,
          :max => 1,
          :tics => 0.1,
          :limit => 0,
          :gnuId => 4,
          :probId => 3,
          :comp => 1,
					:extra => 1,
					:trimere => 0,
				}
		}
	end

	# Executes Logicoil and adds the result to all predictions
	# @param coiled_coil_sequence coiledcoil_sequence	current sequence object
	# @param string title 								sequence title
	# return void
	def logicoil(coiledcoil_sequence, title)
		input = "id,register,sequence\n"; result = []
		# Create Logicoil input data/file
		coiledcoil_sequence.domains.each_with_index { |domain, i|
			register = domain.coiledcoil_data.map { |data| data[2] }
			input << title + "_#{i}" + ',' + register.join('') + ',' + domain.sequence.join('') + "\n"
		}
		File.open(conf[:logicoil][:path] + "/sequences.csv", "w") do |fh|
			fh.puts(input)
		end
		# Execute Logicoil and extract result data
		system(conf[:logicoil][:command])
		result_logicoil = IO.read(conf[:logicoil][:path] + "/logicoil.out")
		result_logicoil.each_line { |line|
			line.scan(conf[:logicoil][:pattern]) { |data| result.push(data[1..4]) }
		}
		# Add Logicoil result to coiledcoil_sequence.domains
		if result.size > 0
			coiledcoil_sequence.domains.each_with_index { |domain, j|
				# Calculate maximum
				sortArr = result[j].sort
				max = result[j].rindex(sortArr[3])
				domain.oligomerisation['logicoil'] = [max < 2 ? max : max+1, sortArr[3].to_f]
			}
		end
	end

	# Executes Scorer2.0 and adds the result to all predictions
	# @param coiled_coil_sequence coiledcoil_sequence	current sequence object
	# @param string title 								sequence title
	# return void
	def scorer2(coiledcoil_sequence, title)
		input = "id,register,sequence\n"; result = []
		# Create Scorer2 input data/file
		coiledcoil_sequence.domains.each_with_index { |domain, i|
			register = domain.coiledcoil_data.map { |data| data[2] }
			input << title + "_#{i}" + ',' + register.join('') + ',' + domain.sequence.join('') + "\n"
		}
		File.open("#{$lib_path}/scorer2/Inputs/sequences.csv", "w") do |fh|
			fh.puts(input)
		end
		# Execute Scorer2 and extract result data
		system(conf[:scorer2][:command])
		result_scorer = IO.read("#{$lib_path}/scorer2/Outputs/scorer2.out")
		result_scorer.each_line { |line|
			line.scan(conf[:scorer2][:pattern]) { |data| result.push(data[0].to_f) }
		}
		# Add Scorer2 result to coiledcoil_sequence.domains
		if result.size > 0
			coiledcoil_sequence.domains.each_with_index { |domain, j|
				domain.oligomerisation['scorer2'] = [(result[j] < 0 ? 3 : 2), result[j].to_f]
			}
		end
	end

	# Executes PrOCoil and adds the result to all predictions
	# @param coiled_coil_sequence coiledcoil_sequence	current sequence object
	# @param string tool 								current cc-tool
	# return void
	def procoil(coiledcoil_sequence, tool)
		input = "sequence,register\n"; result = []
		# Create PrOCoil input data/file
		coiledcoil_sequence.domains.each { |domain|
			register = domain.coiledcoil_data.map { |data| data[2] }
			input << domain.sequence.join('') + ',' + register.join('') + "\n"
		}
		File.open("#{$lib_path}/procoil/Inputs/sequences.csv", "w") do |fh|
			fh.puts(input)
		end
		# Execute PrOCoil and extract result data
		conf[:procoil][:command].gsub!("###PATH###","#{$lib_path}/procoil")
		conf[:procoil][:command].gsub!("###TOOL###","#{tool}")
		conf[:procoil][:command].gsub!("###TIMESTAMP###","#{$timestamp}")
		system(conf[:procoil][:command])
		result_scorer = IO.read("#{$lib_path}/procoil/Outputs/procoil.out")
		result_scorer.each_line { |line|
			line.scan(conf[:procoil][:pattern]) { |data| result.push(data) }
		}
		# Add PrOCoil result to coiledcoil_sequence.domains
		if result.size > 0
			coiledcoil_sequence.domains.each_with_index { |domain, j|
				domain.oligomerisation['procoil'] = [result[j][0] == 'trimer' ? 3 : 2, result[j][1].to_f]
			}
		end
  end

	# Calculates the overall scores for the uploaded sequence with the predefined
	# window-size
	# @param CoiledCoilSequence coiledcoilSequence  current sequence-part object
	# @return Integer[] scores                      SAH window scores for each interval: [0;seq.end-$window_size]
  def sah_total_score(coiledcoilSequence)
    totalScore = 0.0; score = 0.0; scores = []; seq = coiledcoilSequence.full_sequence;
    startPos = seq.start; endPos = seq.end; iterations = endPos-$window_size;
    puts "#{$window_size} #{iterations}"

    if iterations > 0
		  (0..iterations).each{ |i|
        windowStart = startPos+i.to_i; windowStop = ($window_size+1)+i.to_i;
        seqPart = seq.subseq(windowStart, windowStop, true)
        score, summedScore, cumulPathScore = sah_window_score(seq, seqPart)
        scores.push(score)
        totalScore += score
      }
      return scores
    else
      return []
    end
  end

  # Calculates the overall scores for the uploaded sequence with the predefined
	# window-size
  # @param CoiledCoilSequence coiledcoilSequence  current sequence-part object
  # @return Integer[] scores                      SAH window scores for each interval: [0;seq.end]
  def sah_total_endless_score(coiledcoilSequence)
    totalScore = 0.0; scores = []; seq = coiledcoilSequence.full_sequence;
    startPos = seq.start; endPos = seq.end; half = $window_size.divmod(2); iterations = endPos;
    # Determine max. params for SAH
    maxScore = 0; maxSummedScore = 0; maxCumulPathScore = 0
    # puts "#{$window_size} #{iterations} #{startPos}-#{endPos}"
    # puts half.inspect

    # Calculate foreach 'aa' a SAH score
    if iterations > 0
		  (0..iterations-1).each{ |i|
			  # Assign scores to first window-element: |aa-----49------>|
			  # windowStart = (startPos)+i.to_i; windowStop = ($window_size+1)+i.to_i;
			  # Assign scores to mid window-element:   |<--24--aa--24-->||<--14--aa--13-->||<--10--aa--10-->||<--7--aa--6-->|
			  # windowStart = (startPos-(half[0]+1))+i.to_i; windowStop = (half[0]+1+half[1]-1)+i.to_i;
        # Assign scores to mid window-element:   |<--24--aa--24--+4>|
        # shift = shift + 4                      # Add last line (half heptad), to see bottom interactions
			  windowStart = (startPos-(half[0]+1))+i.to_i; windowStop = (half[0]+half[1]-1+4)+i.to_i;
			  # puts ""
			  # puts windowStop.to_s + "-" + windowStart.to_s + "=" + (windowStop-windowStart+1).to_s

        seqPart = seq.subseq(windowStart, windowStop, true)
        # puts seqPart.inspect
        # puts seqPart.coiledcoil_data.length
			  # puts seqPart.end.to_s + "-" + seqPart.start.to_s + "=" + (seqPart.end-seqPart.start).to_s
        score, summedScore, cumulPathScore = sah_window_score(seq, seqPart)
			  # puts "%.04f %.04f %.04f" % [score, summedScore, cumulPathScore]
			  # Evaluation: Show maximum values
			  maxScore = maxScore > score ? maxScore : score
        maxSummedScore = maxSummedScore > summedScore ? maxSummedScore : summedScore
        maxCumulPathScore = maxCumulPathScore > cumulPathScore ? maxCumulPathScore : cumulPathScore
        scores.push(score)
        totalScore += score
      }
		  if $verbose
			  puts "maxSummedScore = " + maxSummedScore.to_s + "; maxCumulPathScore = " + maxCumulPathScore.to_s + "; maxScore = " + maxScore.to_s
		  end
		  return scores
    else
      return []
    end
  end

  # Calculates for a passed sequence and specified window-range the SAH-Score
  # @param CoiledCoilSequence fullsequence  full_sequence object
  # @param CoiledCoilDomain coiledcoil      current sequence-part object
  # @return Float score                     SAH window score
  def sah_window_score(fullsequence, coiledcoil)
    svgObj = SVG.new

    # Prepare interaction data
    gridMat = svgObj.calculate_simple_sah_grid(coiledcoil, false)
    catGrid = svgObj.transpose_simple_grid(gridMat) # ordered by categories
    grid = svgObj.convert_to_interaction_grid(coiledcoil, catGrid) # interaction grid matrix
    intGraph = svgObj.convert_to_interaction_graph(grid) # interaction graph
    cumulPathScore = svgObj.calculate_cumulative_path_score(intGraph)

    # Score calculation
    summedScore, strong, weak = svgObj.calculate_grid_interaction_score(grid)
    # puts "%.04f %d %d" % [summedScore, strong, weak]
    normScore = svgObj.calculate_normalized_interaction_score(coiledcoil, summedScore.to_f, cumulPathScore.to_f)
    # MaxSAH calculation
    # puts "start: " + coiledcoil.start.to_s + " end: " + coiledcoil.end.to_s
    # puts "summedScore = " + summedScore.to_s + "; cumulPathScore = " + cumulPathScore.to_s

    return ("%.4f" % normScore).to_f, summedScore, cumulPathScore
  end

  # Extract table-data out of a tmp-file based on a pattern using regexp
  # @param string object		identifier of coiledcoil tool
  # @param string array	  	window size for the coiledcoil tool
  # return array				array containing the relevant coiledcoil data from CLI-tool result
  def add_sah_score(coiledcoil_sequence, result)
    scores = sah_total_endless_score(coiledcoil_sequence)
    result.each_with_index{ |vec, i|
	    # VARIANT 1
	    # Assign scores to first window-element: |aa-----49------>|
	    if i < scores.size
	      result[i].push(scores[i])
	    else
	      result[i].push(0.0000)
	    end
	    # VARIANT 2
	    # Assign scores to mid window-element:   |<--24--aa--24-->|
	    # mid = ($window_size/2).to_i
	    # if i < scores.size
  	  #    result[(i+mid)%scores.size].push(scores[i])
	    # else
	    #    result[(i+mid)%scores.size].push(0.0000)
	    # end
	    # VARIANT 3
	    # Assign scores to last window-element:  |<------49-----aa|
	    # last = ($window_size).to_i
	    # if i < scores.size
	    #    result[(i+last)%scores.size].push(scores[i])
	    # else
	    #    result[(i+last)%scores.size].push(0.0000)
	    # end
    }
  end

	# Extract table-data out of a tmp-file based on a pattern using regexp
	# @param string id			identifier of coiledcoil tool
	# @param string int			window size for the coiledcoil tool
	# return array				array containing the relevant coiledcoil data from CLI-tool result
	def run_command(tool, id)
		if conf.has_key?(:"#{tool}")
			filename = "#{$timestamp}_run#{id}_#{tool}"
			# Run tool
			result_coil = system(conf[:"#{tool}"][:command]) if tool != 'sah'
			# Extract information 
			result = regexp_extract(filename, conf[:"#{tool}"][:pattern])

			# 'Extrapolate' shortened result data
			if(conf[:"#{tool}"][:extra] == 1)
				result = extrapolate(tool, result)
      end
      # Prepare data for SAH score and merge with existing tool output (last row)
      coiledcoil_sequence = prepare_coildata(result, tool)

      # Calculate SAH window scores - 49 (default: $window_size) and 21, 14 for pure SAH
      puts tool
      add_sah_score(coiledcoil_sequence, result)
      if tool == 'sah'
	      restore_window = $window_size
        $window_size = 28
        add_sah_score(coiledcoil_sequence, result)
        $window_size = 21
        add_sah_score(coiledcoil_sequence, result)
        $window_size = 14
        add_sah_score(coiledcoil_sequence, result)
	      $window_size = restore_window
      end
      # puts "windowSize = " + $window_size.to_s

      # Store result to disk
      array_to_file(result, filename)

      # Load GNUplot-template, change it to cli-tool-format and execute
      make_gnu(filename, conf[:"#{tool}"][:gnuId], (conf[:"#{tool}"][:rows].to_i+1).to_s, conf[:"#{tool}"][:limit], conf[:"#{tool}"][:max], coiledcoil_sequence, conf[:"#{tool}"][:tics], tool)
      # prepare_coildata(shell_data, conf[:"#{tool}"][1], limit[:"#{tool}"], conf[:"#{tool}"][2])
      return result
    else
      return nil
    end
	end

	# Write extrapolated result-data to file
	# @param array resultdata	extrapolated result data
	# @param string filename	name of the tmp_file
	# return void
	def array_to_file(arr, tmp_filename)	
		str = ""
		arr.each{ |entry|
			str << entry.join("\t") + "\n"
		}
		File.open("#{$tmp_path}/#{tmp_filename}.txt", "w") do |fh|
			fh.puts("#{str}")
		end
	end

	# Extract table-data out of a tmp-file based on a pattern using regexp
	# @param string filename	name of the file
	# @param string pattern		for the search with regexp
	# return array				containing the relevant data
	def regexp_extract(filename, pattern)
    result = []
    if File.exists?("#{$tmp_path}/#{filename}.txt")
      tmp_txt = IO.read("#{$tmp_path}/#{filename}.txt")
      tmp_txt.each_line { |line|
        # line.match(pattern) { |data| result.push(data.captures) }
        line.scan(pattern) { |data| result.push(data) }
      }
    end
		return result
	end

	# Extrapolate extracted compressed table-data to full datasets
	# @param string id 			name of the current tool
	# @param array result		extracted 'compressed' table-data
	# return array				containing the relevant data
	def extrapolate(id, result)
		arr = []
		if id == 'paircoil'
			(1..$sequence.length).each{ |i|
				found = false
				result.each { |entry| 
					if i >= entry[1].to_i && i <= entry[2].to_i
						bias = entry[3].ord-97
						arr.push(["#{i}", $sequence[i-1], ((i-entry[1].to_i+bias)%7+97).chr, entry[0]])
						found = true
					end
				} 
				if found == false
					arr.push(["#{i}", $sequence[i-1], ((i-1)%7+97).chr, '0.00'])
				end
			}
			return arr
		end
		if id == 'marcoil'
			counter = 1; bias = 0; offset = 0
			result.each_with_index { |entry, i|
				while counter < entry[0].to_i do
					arr.push([counter.to_s, $sequence[counter-1], ((counter-offset+bias)%7+97).chr, entry[2]])					
					counter = counter+1 	
				end
				arr.push([entry[0], entry[1], entry[3], entry[2]])
				counter = counter+1 	
				bias = entry[3].ord-97; offset = entry[0].to_i
			}
			return arr
		end
		if id == 'sah'
      i = 1
      $sequence.each_char { |aa|
        arr.push([i.to_s, aa, ((i-1)%7+97).chr, 0.0])
        i += 1
      }
			return arr
		end
	end

	# Load GNUplot-template, set correct source-file and params, then execute
	# @param string filename	name of the gnu-templatefile
	# @param string ccprob		position of coiledcoil probability (p-score)
	# @param string sahRow		SAH row position
	# @param string limit	    upper limit for the plot-cmd
	# @param string max	    	Maximum value
	# @param string coiledcoil_sequence   	CoiledCoilSequence object
	# @param string tics		  y-axis division
	# @param string tool		  name of used tool
	# return void				      output of the gnuplot-script
	def make_gnu(filename, ccprob, sahRow, limit, max, coiledcoil_sequence, tics, tool)
		if limit > 0
			gnu_file_tmpl = IO.read("#{Rails.root}/lib/assets/display.gnu_tmpl")
			gnu_file_tmpl.gsub!("###CCPROB###","#{ccprob}")
		else
			gnu_file_tmpl = IO.read("#{Rails.root}/lib/assets/display.sah.gnu_tmpl")
			gnu_file_tmpl.gsub!("###SAH28###",(sahRow.to_i+1).to_s)
			gnu_file_tmpl.gsub!("###SAH21###",(sahRow.to_i+2).to_s)
			gnu_file_tmpl.gsub!("###SAH14###",(sahRow.to_i+3).to_s)
		end
		gnu_file_tmpl.gsub!("###FILE###","#{$tmp_path}/#{filename}.txt")
		gnu_file_tmpl.gsub!("###OUTPUT###","#{$tmp_path_image}/#{filename}.svg")
		gnu_file_tmpl.gsub!("###TITLE###","#{$title[0..32]}...")
		gnu_file_tmpl.gsub!("###TOOL###","#{tool}")
		gnu_file_tmpl.gsub!("###SAH###","#{sahRow}")
		gnu_file_tmpl.gsub!("###LIMIT###","#{limit}")
		gnu_file_tmpl.gsub!("###MAX###","#{max}")
		gnu_file_tmpl.gsub!("###LENGTH###","#{coiledcoil_sequence.seq_length}")
		gnu_file_tmpl.gsub!("###TICS###","#{tics}")
		gnu_file_tmpl.gsub!("###CURVE###", (tool=='paircoil2'? 'below':'above') )
		File.open("#{$tmp_path}/display.gnu", "w") do |fh|
			fh.puts(gnu_file_tmpl.to_s)
		end
		result_gnuplot = `cd #{$tmp_path} && chmod 775 ./display.gnu && ./display.gnu`
	end


	# Extract coiled coil-data out of a data field (CLI computation result)
	# @param array data			complete coil data information from CL-tool
	# 							array-data pattern => "value", [global position, AA, heptad-position, score, p-score]
	# @param string id 			Cli-Tool name
	# @return object			CoiledCoilSequence; containing the relevant coiled coil data
	# TODO: Need to be refactored => double code through comparison
	def prepare_coildata(data, id)
		match = false; entry = []; dimere = 0.0; trimere = 0.0; min_length = 14
		register = {"a" => "b", "b" => "c", "c" => "d", "d" => "e", "e" => "f", "f" => "g", "g" => "a"}
		compare_field = conf[:"#{id}"][:probId]
		limit = conf[:"#{id}"][:limit]
		comparator = conf[:"#{id}"][:comp]
		
		coiledcoils = CoiledCoilSequence.new
		coiledcoil = CoiledCoilDomain.new

		# print data
		# print "\n"

		# Create CoiledCoilDomain for whole data
		coiledcoil.start = 1
		coiledcoil.end = data.count
		data.each_with_index { |entry, i|
			# Detect breaks in register
			if i+1 < coiledcoil.end && data[i+1][2] != register[entry[2]]
				coiledcoil.breaks.push(i)
			end						
			coiledcoil.coiledcoil_data.push(entry)
			coiledcoil.sequence.push(entry[1])
		}
		coiledcoils.full_sequence = coiledcoil

		# Create CoiledCoilDomain for regions divided by limits
		bias = 0; coiledcoil = CoiledCoilDomain.new
		if comparator == 0   # Paircoil2
			data.each_with_index { |entry, i|
				if (entry[compare_field].to_f < limit)
					if match == false
						match = true
						coiledcoil = CoiledCoilDomain.new
						coiledcoil.start = entry[0].to_i
					end
				# 	if i+1 == coiledcoil.start
				# 		bias = (entry[2].ord-97)%7
				# 	end
				# 	# Heptad-Zuweisungen relativ zur ersten Aminosaure
				# 	case (i+1-coiledcoil.start+bias)%7
				# 		when 0; coiledcoil.a.elements.push(entry)
				# 		when 1; coiledcoil.b.elements.push(entry)
				# 		when 2; coiledcoil.c.elements.push(entry)
				# 		when 3; coiledcoil.d.elements.push(entry)
				# 		when 4; coiledcoil.e.elements.push(entry)
				# 		when 5; coiledcoil.f.elements.push(entry)
				# 		when 6; coiledcoil.g.elements.push(entry)
				# 	end
					# Heptad-Zuweisung: Aminosaure entsprechend ihrer Vorhersage
				# 	if entry[2] == 'a' && found == false
				# 		coiledcoil.start = entry[0].to_i
				# 		found = true
				# 	end
					case entry[2]
						when 'a'; coiledcoil.a.elements.push(entry)
						when 'b'; coiledcoil.b.elements.push(entry)
						when 'c'; coiledcoil.c.elements.push(entry)
						when 'd'; coiledcoil.d.elements.push(entry)
						when 'e'; coiledcoil.e.elements.push(entry)
						when 'f'; coiledcoil.f.elements.push(entry)
						when 'g'; coiledcoil.g.elements.push(entry)
					end
					# Detect breaks in register
					if i+1 < coiledcoils.full_sequence.end && data[i+1][2] != register[entry[2]]
						coiledcoil.breaks.push(i)
					end
					coiledcoil.coiledcoil_data.push(entry)
					coiledcoil.sequence.push(entry[1])
				else
					if match == true
						coiledcoil.end = entry[0].to_i-1
						if min_length < coiledcoil.coiledcoil_data.count
							coiledcoils.domains.push(coiledcoil)
							coiledcoils.counter += 1
						end
						match = false
					end
				end
			}
		else	# All other tools without Paircoil2
			data.each_with_index { |entry, i|
				if (entry[compare_field].to_f > limit)
					if match == false
						match = true
						coiledcoil = CoiledCoilDomain.new
						coiledcoil.start = entry[0].to_i
						dimere = 0; trimere = 0
					end
				# 	if i+1 == coiledcoil.start
				# 		bias = (entry[2].ord-97)%7
				# 	end
				# 	# Heptad-Zuweisungen relativ zur ersten Aminosaure
				# 	case (i+1-coiledcoil.start+bias)%7
				# 		when 0; coiledcoil.a.elements.push(entry)
				# 		when 1; coiledcoil.b.elements.push(entry)
				# 		when 2; coiledcoil.c.elements.push(entry)
				# 		when 3; coiledcoil.d.elements.push(entry)
				# 		when 4; coiledcoil.e.elements.push(entry)
				# 		when 5; coiledcoil.f.elements.push(entry)
				# 		when 6; coiledcoil.g.elements.push(entry)
				# 	end
					# Heptad-Zuweisung: Aminosaure entsprechend ihrer Vorhersage
					case entry[2]
						when 'a'; coiledcoil.a.elements.push(entry)
						when 'b'; coiledcoil.b.elements.push(entry)
						when 'c'; coiledcoil.c.elements.push(entry)
						when 'd'; coiledcoil.d.elements.push(entry)
						when 'e'; coiledcoil.e.elements.push(entry)
						when 'f'; coiledcoil.f.elements.push(entry)
						when 'g'; coiledcoil.g.elements.push(entry)
					end
					if conf[:"#{id}"][:trimere]	== 1
						dimere = dimere + entry[4].to_f
						trimere = trimere + entry[5].to_f
					end
					# Detect breaks in register
					if i+1 < coiledcoils.full_sequence.end && data[i+1][2] != register[entry[2]]
						coiledcoil.breaks.push(i)
					end
					coiledcoil.coiledcoil_data.push(entry)
					coiledcoil.sequence.push(entry[1])
				else
					if match == true
						coiledcoil.end = entry[0].to_i-1
						if conf[:"#{id}"][:trimere]	== 1
							isTrimer = dimere < trimere ? true : false
							coiledcoil.oligomerisation["#{id}"] = [(isTrimer ? 3 : 2), isTrimer ? trimere : dimere]
						end					
						if min_length < coiledcoil.coiledcoil_data.count
							coiledcoils.domains.push(coiledcoil)
							coiledcoils.counter += 1
						end
						match = false
					end
				end
			}
		end
		# Sets info for regions which end on the last sequence-character
		if match == true
			coiledcoil.end = data.count
			if conf[:"#{id}"][:trimere]	== 1
				isTrimer = dimere < trimere ? true : false
				coiledcoil.oligomerisation["#{id}"] = [(isTrimer ? 3 : 2), isTrimer ? trimere : dimere]
			end					
			if min_length < coiledcoil.coiledcoil_data.count
				coiledcoils.domains.push(coiledcoil)
				coiledcoils.counter += 1
			end
		end
		coiledcoils.seq_length = data.count

		# print coiledcoils.full_sequence.coiledcoil_data
		# print "\n"
		# print coiledcoils.full_sequence.
		# coiledcoil.get_categories().each { |item|
		# 	item.elements.each { |pos| 
		# 		print pos
		# 		print "\n"
		# 	}
		# 	print "\n"
		# }
		# coiledcoils.domains.each { |domain|
		# 	print domain.breaks
		# 	print "\n"
		# }

		return coiledcoils
	end

end

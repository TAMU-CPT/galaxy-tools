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

# Global object for each CLI result
class CoiledCoilSequence

	# Attributes
	attr_accessor :title			  # string: FASTA title/identifier of the Sequence (for displaying purposes)
	attr_accessor :definition	  # string: FASTA title/identifier of the Sequence (for displaying purposes)
	attr_accessor :name				  # string: Unique identifier of the Sequence (for internal/identification purposes)
	attr_accessor :cli_name			# string: Name of CLI tool
	attr_accessor :full_sequence	# CoiledCoilDomain: From complete uploaded sequence
	attr_accessor :counter			# int: number of CoiledCoilDomains
	attr_accessor :domains			# CoiledCoilDomain[]: Array of detected CoiledCoilDomain's
	attr_accessor :domain_model	# string: SVG
	attr_accessor :seq_length		# int: Sequence length
	attr_accessor :seq_text			# string: Sequence information (heptad)
	attr_accessor :limit			  # float:
	attr_accessor :filename			# string: Filename of gnu-plot 

	def initialize(title = '', definition = '', name = '', sequence = '')
		@title = title || ''
		@definition = definition || ''
		@name = name || ''
		@cli_name = ''
		@full_sequence = CoiledCoilDomain.new(sequence)
		@counter = 0
		@domains = []
		@domain_model = ''
		@seq_length = 0
		@seq_text = []
		@limit = 0.0
		@filename = ''
	end

	def plain
		print cli_name + "\n"
		print "\n"
		coiledcoil_data().each { |item|
			item.elements.each { |pos| 
				print pos
				print "\n"
			}
			print "\n"
		}
	end
	
end

class CoiledCoilDomain < CoiledCoilSequence

	@@count = 0									# int: Global sequence counter (class variable)
	attr_accessor :coiledcoil_data				# string[]: Coiled-coil sequence information
	attr_accessor :coiledcoil_env_data		# string[]: Coiled-coil sequence information
	attr_accessor :a, :b, :c, :d, :e, :f, :g 	# string[]: Categorized coiled-coil information
	attr_accessor :sequence, :start, :end # string, int, int: Sequence information
	attr_accessor :breaks						      # int[]: positions where register has breaks
	attr_accessor :reverse						    # int: reverseFlag [0=false, 1=true]
	attr_accessor :oligomerisation				# hash: oligomerisation state; according probability [oligoID, score]
	attr_accessor :cc_svg						      # string: SVG
	attr_accessor :sah_svg						    # string: SVG

	def initialize(sequence = '')
		@@count += 1
		@coiledcoil_data = construct_cc_data(sequence) || []
		@a = Category.new("a")
		@b = Category.new("b")
		@c = Category.new("c")
		@d = Category.new("d")
		@e = Category.new("e")
		@f = Category.new("f")
		@g = Category.new("g")
		@sequence = []
		@start = 0; @end = sequence.size || 0; @reverse = 0
		@breaks = []
		@oligomerisation = Hash.new()
		@cc_svg = @sah_svg = ''
	end

  # Helper method for initialization
  # @param string sequence    amino acid sequence
  # return array coiledcoil_data
  def construct_cc_data(sequence)
    arr = []; i = 1
    sequence.each_char { |aa|
      arr.push([i.to_s, aa, ((i-1)%7+97).chr, 0.0])
      i += 1
    }
    return arr
  end

	# Getter for categories
	def get_categories
	  [@a, @b, @c, @d, @e, @f, @g]
  	end

	# Assign positions to categories
	# @param array points		circle coordinates / positions
	def assign_positions(points)
		heptads = get_categories; i = 0
		points.each { |point|
			heptads[i].pos = point
			i = (i+4) % heptads.count
		}
	end

	# Extract sub-sequence {startCoord-endCoord} out of full_sequence
	# Produce categorized coiled-coil information
	# @param int startCoord
	# @param int endCoord
	# @param bool strict      subsequence is cutted strict by position or on starting A-heptad
	# @param bool gapMode     set gapMode (gapless = false, break-based = true)
	# return CoiledCoilDomain
	def subseq(startCoord, endCoord, strict, gapMode = true)
		bias = 0; found = false
		coiledcoil = CoiledCoilDomain.new
		if startCoord.to_i < endCoord.to_i
			coiledcoil.start = self.start <= startCoord.to_i ? startCoord.to_i : self.start
			coiledcoil.end = self.end >= endCoord.to_i ? endCoord.to_i : self.end
			# puts "Coords: " + coiledcoil.start.to_s + " " + coiledcoil.end.to_s

			## Extract data from coil_data for new coiledcoil-domain
			self.coiledcoil_data.each_with_index { |entry, i|
				if (i >= coiledcoil.start-1 && i <= coiledcoil.end-1)
					## Heptad-Zuweisung: CoiledCoil-Domäne startet mit erstem 'a' ab 'start'
					if strict == false
						if entry[2] == 'a' && found == false
							shift = entry[0].to_i-coiledcoil.start
							coiledcoil.start = entry[0].to_i
							coiledcoil.end   = coiledcoil.end+shift
							found = true
						end
						if found == true
							if gapMode == true
								## Zuordnung genau nach Vorhersage
								case entry[2]
									when 'a'; coiledcoil.a.elements.push(entry)
									when 'b'; coiledcoil.b.elements.push(entry)
									when 'c'; coiledcoil.c.elements.push(entry)
									when 'd'; coiledcoil.d.elements.push(entry)
									when 'e'; coiledcoil.e.elements.push(entry)
									when 'f'; coiledcoil.f.elements.push(entry)
									when 'g'; coiledcoil.g.elements.push(entry)
								end
							else
								## Heptad-Zuweisung im Heptad-Raster relativ zum Start
								case (i-coiledcoil.start+1)%7
									when 0; coiledcoil.a.elements.push(entry)
									when 1; coiledcoil.b.elements.push(entry)
									when 2; coiledcoil.c.elements.push(entry)
									when 3; coiledcoil.d.elements.push(entry)
									when 4; coiledcoil.e.elements.push(entry)
									when 5; coiledcoil.f.elements.push(entry)
									when 6; coiledcoil.g.elements.push(entry)
								end
							end
						end

					## Strict sequence cutting: Start on specified start position
					else
            if gapMode == true
		            ## Zuordnung genau nach Vorhersage
								case entry[2]
									when 'a'; coiledcoil.a.elements.push(entry)
									when 'b'; coiledcoil.b.elements.push(entry)
									when 'c'; coiledcoil.c.elements.push(entry)
									when 'd'; coiledcoil.d.elements.push(entry)
									when 'e'; coiledcoil.e.elements.push(entry)
									when 'f'; coiledcoil.f.elements.push(entry)
									when 'g'; coiledcoil.g.elements.push(entry)
								end
						else
								## Heptad-Zuweisung im Heptad-Raster relativ zum Start
		            if i == startCoord.to_i
									bias = (entry[2].ord-97)%7 # Biased from a
									print "Erste Aminosaure: #{entry} bias=#{bias.to_s} \n"
			          end
								case (i-startCoord.to_i+bias)%7
									when 0; coiledcoil.a.elements.push(entry)
									when 1; coiledcoil.b.elements.push(entry)
									when 2; coiledcoil.c.elements.push(entry)
									when 3; coiledcoil.d.elements.push(entry)
									when 4; coiledcoil.e.elements.push(entry)
									when 5; coiledcoil.f.elements.push(entry)
									when 6; coiledcoil.g.elements.push(entry)
								end
						end
					end
					coiledcoil.coiledcoil_data.push(entry)
					coiledcoil.sequence.push(entry[1])
				end
			}

			# print coiledcoil.sequence
			# print "\n"
			# coiledcoil.get_categories().each { |item|
			# 	item.elements.each { |pos|
			# 		print pos
			# 		print "\n"
			# 	}
			# 	print "\n"
			# }

			return coiledcoil
		end
	end
end

class Category
	attr_accessor :name, :pos, :elements

	def initialize(name)
		@elements = []
		@name = name
		@pos = [0, 0]
	end

end

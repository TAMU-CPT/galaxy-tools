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


class SVG

	# Class constants
	STRONG = 1.0; MEDIUM = 0.75; WEAK = 0.5; LOW = 0.25; PENALTY_ONE = -0.5; PENALTY_TWO = -0.75;
	UP_RIGHT = 1; UP_LEFT = 2; BOTTOM_RIGHT = 4; BOTTOM_LEFT = 8;

	# attr_accessor :skip_headers
	attr_accessor :globals			# array: global variabls

	# SVG: Helper-method to generate n points on a circle (evenly distributed)
	# @param [Array] points   number of points
	# @param [Integer] radius     radius of the circle
	# @param [Integer] rotate     in degree, change position of the points on the circle
	# @param [Integer] offset_x   offset of the circle in x-direction
	# @param [Integer] offset_y   offset of the circle in y-direction
	# @return [Array]         point, center of small circle
	def point_on_circle(number, radius, rotate = 0, offset_x = 0, offset_y = 0)
		point = ''
		deg = 2*Math::PI/number			# radian
		rotate = rotate*Math::PI/180	# convert to radian
		(1..1).each { |i|
			x_coord = Math::cos( (i - 0.2) * deg + rotate) * radius + offset_x
			y_coord = Math::sin( (i - 0.2) * deg + rotate) * radius + offset_y
			point = [x_coord, y_coord]
		}
		return point
	end

	# SVG: Helper-method to generate n points on a circle (evenly distributed)
	# @param [Array] points   number of points
	# @param [Integer] radius     radius of the circle
	# @param [Integer] rotate     in degree, change position of the points on the circle
	# @param [Integer] offset_x   offset of the circle in x-direction
	# @param [Integer] offset_y   offset of the circle in y-direction
	# @return [Array]         points, centers of small circles
	def points_on_circle(number, radius, rotate = 0, offset_x = 0, offset_y = 0)
		points = []
		deg = 2*Math::PI/number			# radian
		rotate = rotate*Math::PI/180	# convert to radian
		(1..number).each { |i|
			x_coord = Math::cos( (i - 0.2) * deg + rotate) * radius + offset_x
			y_coord = Math::sin( (i - 0.2) * deg + rotate) * radius + offset_y
			points.push([x_coord, y_coord])
		}
		return points
	end

	# SVG: Helper-method to generate n points on a circle (evenly distributed)
	# @param [Array] points   number of points
	# @param [Integer] radius     radius of the circle
	# @param [Integer] rotate     in degree, change position of the points on the circle
	# @param [Integer] offset_x   offset of the circle in x-direction
	# @param [Integer] offset_y   offset of the circle in y-direction
	# @return [Array] points  centers of small circles
	def points_on_spiral(number, radius, rotate = 0, offset_x = 0, offset_y = 0)
		points = []; spiral_offset = (16.to_f/number);
		deg = 2*Math::PI/number			# radian
		rotate = rotate*Math::PI/180	# convert to radian
		(1..number).each { |i|
			sho = i*spiral_offset
			x_coord = Math::cos( (i - 0.2) * deg + rotate) * (radius+sho) + offset_x
			y_coord = Math::sin( (i - 0.2) * deg + rotate) * (radius+sho) + offset_y
			points.push([x_coord, y_coord])
		}
		return points
	end

	# SVG: Helper-method to calculate next point (on a line)
	# @param [Array] point    position of source point
	# @param [Integer] center     center of main shape circle
	# @param [Integer] dist       distance to next point (pos./neg. value)
	# @return [Array] point   position of the next point according to the center
	def calc_nextpoint(point, center, distance)
		delta = [point[0]-center[0], point[1]-center[1]]
		slope = delta[1]/delta[0]
		ratio = distance.abs/Math.sqrt(delta[0]**2+delta[1]**2)
		interval = ratio*delta[0]
		distance < 0 ? sign = -1 : sign = 1
		return nextpoint = [point[0] + sign*interval, point[1] + sign*slope*interval]
	end

	# SVG: Helper-method to calculate the corresponding color of an AA
	# @param [String] aa     abbreviation of amino acid
	# @param [String] type   type of color-set
	# @return [String] color
	def calc_aa_color(aa, type)
			color_cc = ["black", "red", "blue", "green"]
			color_sah = ["black", "red", "blue", "black"]
			aminoacid = Hash["A" => 0, "C" => 3, "D" => 1, "E" => 1, "F" => 0, "G" => 3,
						 "H" => 2, "I" => 0, "K" => 2, "L" => 0, "M" => 0, "N" => 3, "P" => 3,
						 "Q" => 3, "R" => 2, "S" => 3, "T" => 3, "V" => 0, "W" => 0, "Y" => 0]
		if aa != nil && aa != "" && aa != '' && aminoacid[aa] != nil
			if type == "cc"
				return color_cc[aminoacid[aa]]
			else
				return color_sah[aminoacid[aa]]
			end
		else
			return 'black'
		end
	end

	# DEPRECATED
	# Helper-method: Calculates ALL strong and weak interactions for the current amino-acid
	# @param [CoiledCoilDomain] coiledcoil
	# @param [Integer] pos            Current sequence position
	# @param [String] element         Current amino-acid
	# @return [Integer, Integer]      Count of strong, weak interactions
	def check_all_sah_neighbors(coiledcoil, pos, element)
		strong = weak = 0; last = false
		sequence = coiledcoil.sequence
		len = coiledcoil.end-coiledcoil.start

		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		 bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		if pos-4 >= 0 && len-pos >= 4 	# for all interactions
			if element == 'E'
				strong += if ['K','R'].rindex(sequence[pos-4]) != nil then 1 else 0 end
				strong += if ['K','R'].rindex(sequence[pos+3]) != nil then 4 else 0 end
				if [1,5].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos-4], "sah") == "blue" then 1 else 0 end
				elsif [4,5].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos+3], "sah") == "blue" then 4 else 0 end
				end
				weak += if calc_aa_color(sequence[pos-3], "sah") == "blue" then 2 else 0 end
				weak += if calc_aa_color(sequence[pos+4], "sah") == "blue" then 8 else 0 end
			elsif ['K','R'].rindex(element) != nil
				strong += if sequence[pos-3] == 'E' then 2 else 0 end
				strong += if sequence[pos+4] == 'E' then 8 else 0 end
				if [2,10].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos-3], "sah") == "red" then 2 else 0 end
				elsif [8,10].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos+4], "sah") == "red" then 8 else 0 end
				end
				weak += if calc_aa_color(sequence[pos-4], "sah") == "red" then 1 else 0 end
				weak += if calc_aa_color(sequence[pos+3], "sah") == "red" then 4 else 0 end
			elsif element == 'Q'
				strong += if sequence[pos+4] == 'E' then 8 else 0 end
				if strong != 8
					weak += if calc_aa_color(sequence[pos+4], "sah") == "red" || sequence[pos+4] == "Q" then 8 else 0 end
				end
				weak += if calc_aa_color(sequence[pos-4], "sah") == "red" || sequence[pos-4] == "Q" then 1 else 0 end
				weak += if calc_aa_color(sequence[pos-3], "sah") == "red" || sequence[pos-3] == "Q" then 2 else 0 end
				weak += if calc_aa_color(sequence[pos+3], "sah") == "red" || sequence[pos+3] == "Q" then 4 else 0 end
			end
		end

		return strong, weak
	end

	# DEPRECATED
	# Helper-method: Calculates bottom (strong and weak) interactions for the current amino-acid
	# @param [CoiledCoilDomain] coiledcoil
	# @param [Integer] pos            Current sequence position
	# @param [String] element         Current amino-acid
	# @return [Integer, Integer]      Count of strong, weak interactions
	def old_rules_check_sah_neighbors(coiledcoil, pos, element)
		strong = 0; weak = 0; last = false
		sequence = coiledcoil.sequence
		len = coiledcoil.end-coiledcoil.start

		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		 bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		if len-pos >= 4  		# only bottom interactions are calculated
				# E -> K,R,Q
				if element == 'E'
				strong += if ['K','R'].rindex(sequence[pos+3]) != nil then 4 else 0 end
				if [4,5].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos+3], "sah") == "blue" || sequence[pos+3] == "Q" then 4 else 0 end
				end
				weak += if calc_aa_color(sequence[pos+4], "sah") == "blue" || sequence[pos+4] == "Q"  then 8 else 0 end
			elsif ['K','R'].rindex(element) != nil
				strong += if sequence[pos+4] == 'E' then 8 else 0 end
				if [8,10].rindex(strong) == nil
					weak += if calc_aa_color(sequence[pos+4], "sah") == "red" then 8 else 0 end
				end
				weak += if calc_aa_color(sequence[pos+3], "sah") == "red" then 4 else 0 end
			elsif element == 'Q'
				strong += if sequence[pos+4] == 'E' then 8 else 0 end
				if strong != 8
					weak += if calc_aa_color(sequence[pos+4], "sah") == "red" || sequence[pos+4] == "Q" then 8 else 0 end
				end
				weak += if calc_aa_color(sequence[pos+3], "sah") == "red" || sequence[pos+3] == "Q" then 4 else 0 end
			end
		end

		return strong, weak
	end

	# Helper-method: Calculates bottom (strong and weak) interactions for the current amino-acid
	# @param [CoiledCoilDomain] coiledcoil
	# @param [Integer] pos            Current sequence position
	# @param [String] element         Current amino-acid
	# @return [Integer, Integer]      Count of strong, weak interactions
	def check_sah_neighbors(coiledcoil, pos, element)
		strong = 0; weak = 0; last = false
		sequence = coiledcoil.sequence
		len = coiledcoil.end-coiledcoil.start

		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		 bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		if len-pos >= 4  		# only bottom interactions are calculated
			# E -> K,R (blue) or Q
			if element == 'E'
				weak += if calc_aa_color(sequence[pos+3], "sah") == "blue" || sequence[pos+3] == "Q" then 4 else 0 end
				weak += if calc_aa_color(sequence[pos+4], "sah") == "blue" || sequence[pos+4] == "Q" then 8 else 0 end
			# K,R -> E,D
			elsif ['K','R'].rindex(element) != nil
				strong += if ['E','D'].rindex(sequence[pos+3]) != nil then 4 else 0 end
				strong += if ['E','D'].rindex(sequence[pos+4]) != nil then 8 else 0 end
			# D -> K,R
			elsif element == 'D'
				weak += if ['K','R'].rindex(sequence[pos+3]) != nil then 4 else 0 end
				weak += if ['K','R'].rindex(sequence[pos+4]) != nil then 8 else 0 end
			# Q -> E,(Q)
			elsif element == 'Q'
				weak += if sequence[pos+3] == 'E' then 4 else 0 end
				weak += if sequence[pos+4] == 'E' then 8 else 0 end
				# weak += if calc_aa_color(sequence[pos+3], "sah") == "red" || sequence[pos+3] == "Q" then 4 else 0 end
				# weak += if calc_aa_color(sequence[pos+4], "sah") == "red" || sequence[pos+4] == "Q" then 8 else 0 end
			end
		end

		return strong, weak
	end

	# DEPRECATED
	# Helper-method: Check interactions of triangle neighbors in SAH-grid (parent, child_left, child_right)
	# @param [String] aa		the parent aminoacid
	# @param [String] nbl		the bottom left aminoacid
	# @param [String] nbr		the bottom right aminoacid
	# @return [Int, Int] 		strong, weak interaction indices
	def old_rules_classify_sah_grid_neighbors(aa, nbl, nbr)
		strong = 0; weak = 0;
		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		if aa == 'E'
			strong += if ['K','R'].rindex(nbr) != nil then 4 else 0 end
			if [4,5].rindex(strong) == nil
				weak += if calc_aa_color(nbr, "sah") == "blue" || nbr == "Q" then 4 else 0 end
			end
			weak += if calc_aa_color(nbl, "sah") == "blue" || nbl == "Q"  then 8 else 0 end
		elsif ['K','R'].rindex(aa) != nil
			strong += if nbl == 'E' then 8 else 0 end
			if [8,10].rindex(strong) == nil
				weak += if calc_aa_color(nbl, "sah") == "red" then 8 else 0 end
			end
			weak += if calc_aa_color(nbr, "sah") == "red" then 4 else 0 end
		elsif aa == 'Q'
			strong += if nbl == 'E' then 8 else 0 end
			if strong != 8
				weak += if calc_aa_color(nbl, "sah") == "red" || nbl == "Q" then 8 else 0 end
			end
			weak += if calc_aa_color(nbr, "sah") == "red" || nbr == "Q" then 4 else 0 end
		end

		return strong, weak
	end

	# Helper-method: Check interactions of triangle neighbors in SAH-grid (parent, child_left, child_right)
	# @param [String] aa		the parent aminoacid
	# @param [String] nbl		the bottom left aminoacid
	# @param [String] nbr		the bottom right aminoacid
	# @return [Int, Int] 		strong, weak interaction indices
	def classify_sah_grid_neighbors(aa, nbl, nbr)
		strong = 0; weak = 0;
		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		# E -> K,R (blue) or Q
		if aa == 'E'
				weak += if ['K','R','Q'].rindex(nbr) != nil then 4 else 0 end
				weak += if ['K','R','Q'].rindex(nbl) != nil then 8 else 0 end
		# K,R -> E,D
		elsif ['K','R'].rindex(aa) != nil
				strong += if ['E','D'].rindex(nbr) != nil then 4 else 0 end
				strong += if ['E','D'].rindex(nbl) != nil then 8 else 0 end
		# D -> K,R
		elsif aa == 'D'
				weak += if ['K','R'].rindex(nbr) != nil then 4 else 0 end
				weak += if ['K','R'].rindex(nbl) != nil then 8 else 0 end
		# Q -> E,(Q)
		elsif aa == 'Q'
				weak += if nbr == 'E' then 4 else 0 end
				weak += if nbl == 'E' then 8 else 0 end
				# weak += if calc_aa_color(nbr, "sah") == "red" || nbr == "Q" then 4 else 0 end
				# weak += if calc_aa_color(nbl, "sah") == "red" || nbl == "Q" then 8 else 0 end
		end

		return strong, weak
	end

	# Helper-method: Check interactions of triangle neighbors in SAH-grid (parent, child_left, child_right)
	# Version 3: 04.2014
	# @param [String] aa		the parent aminoacid
	# @param [String] nbl		the bottom left aminoacid
	# @param [String] nbr		the bottom right aminoacid
	# @return [Int, Int, Int] 		strong, middle, weak interaction indices
	def classify_greyscale_sah_grid_neighbors_v3(aa, nbl, nbr)
		strong = 0; middle = 0; weak = 0;
		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		# E -> K,R (blue) or Q
		if aa == 'E'
			weak += if ['K','R','Q'].rindex(nbr) != nil then 4 else 0 end
			weak += if ['K','R','Q'].rindex(nbl) != nil then 8 else 0 end
		# K,R -> E,D
		elsif ['K','R'].rindex(aa) != nil
			strong += if ['E','D'].rindex(nbr) != nil then 4 else 0 end
			strong += if ['E','D'].rindex(nbl) != nil then 8 else 0 end
		# D -> K,R
		elsif aa == 'D'
			weak += if ['K','R'].rindex(nbr) != nil then 4 else 0 end
			middle += if ['K','R'].rindex(nbl) != nil then 8 else 0 end
		# Q -> E,(Q)
		elsif aa == 'Q'
			weak += if nbr == 'E' then 4 else 0 end
			weak += if nbl == 'E' then 8 else 0 end
			# weak += if calc_aa_color(nbr, "sah") == "red" || nbr == "Q" then 4 else 0 end
			# weak += if calc_aa_color(nbl, "sah") == "red" || nbl == "Q" then 8 else 0 end
		end

		return strong, middle, weak
	end

	# Helper-method: Check interactions of triangle neighbors in SAH-grid (parent, child_left, child_right)
	# Version 4: 04.2014
	# @param [String] aa		the parent aminoacid
	# @param [String] nbl		the bottom left aminoacid
	# @param [String] nbr		the bottom right aminoacid
	# @return [Int, Int, Int] 		strong, middle, weak interaction indices
	def classify_greyscale_sah_grid_neighbors_v4(aa, nbl, nbr)
		strong = 0; medium = 0; weak = 0;

		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8

		# SCORING
		# Strong (score = 1.0)
		# E  =>  K/R : i+4
		# H/K/R => E : i+4
		# Medium (score = 0.75)
		# E  =>  K/R : i+3
		# D  =>  K/R : i+4
		# H/K/R => E : i+3
		# Weak   (score = 0.5)
		# E   =>   H : i+3
		# D  =>  K/R : i+3
		# K/R  =>  D : i+3
		# Q   =>   E : i+4, i+3
		# E   =>   Q : i+4, i+3

		# E -> K/R, H, Q
		if aa == 'E'
			strong += if ['K','R'].rindex(nbl) != nil then BOTTOM_LEFT else 0 end
			medium += if ['K','R'].rindex(nbr) != nil then BOTTOM_RIGHT else 0 end
			weak += if ['H','Q'].rindex(nbr) != nil then BOTTOM_RIGHT else 0 end
			weak += if                   nbl == 'Q' then BOTTOM_LEFT else 0 end
			# H/K/R -> E
		elsif ['H','K','R'].rindex(aa) != nil
			strong += if nbl == 'E' then BOTTOM_LEFT else 0 end
			medium += if nbr == 'E' then BOTTOM_RIGHT else 0 end
			# K/R -> D
		elsif ['K','R'].rindex(aa) != nil
			weak += if nbr == 'D' then BOTTOM_RIGHT else 0 end
			# D -> K/R
		elsif aa == 'D'
			medium += if ['K','R'].rindex(nbl) != nil then BOTTOM_LEFT else 0 end
			weak += if ['K','R'].rindex(nbr) != nil then BOTTOM_RIGHT else 0 end
			# Q -> E
		elsif aa == 'Q'
			weak += if nbl == 'E' then BOTTOM_LEFT else 0 end
			weak += if nbr == 'E' then BOTTOM_RIGHT else 0 end
		end

		return strong, medium, weak
	end

	# Helper-method: Check interactions of triangle neighbors in SAH-grid (parent, child_left, child_right)
	# Version 5: 27.09.2016
	# @param [String] aa		the parent aminoacid
	# @param [String] nbl		the bottom left aminoacid
	# @param [String] nbr		the bottom right aminoacid
	# Scoring-Matrices available as Gloabls
	# @param [Float] bottomRightMatrix  Bottom right interaction scoring-matrix
	# @param [Float] bottomLeftMatrix   Bottom left interaction scoring-matrix
	# @return [Int, Int, Int, Int] 		strong, middle, weak, helix_stability interaction indices
	def classify_greyscale_sah_grid_neighbors(aa, nbl, nbr) #, bottomRightMatrix, bottomLeftMatrix)
		strong = 0; medium = 0; weak = 0; helix_stability = 0;
		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8

		# Bottom right neighbor
		# puts aa + " " + nbr
		# puts $bottomRightMatrix.inspect # Does this value fit to the real meant value???
		# if !$bottomRightMatrix[aa][nbr].nil?
		case $bottomRightMatrix[aa][nbr]
			when STRONG
				   strong += BOTTOM_RIGHT
			when MEDIUM
					medium += BOTTOM_RIGHT
			when WEAK
				   weak += BOTTOM_RIGHT
			when LOW
					 helix_stability += BOTTOM_RIGHT
			when PENALTY_ONE
					 helix_stability += BOTTOM_RIGHT
			when PENALTY_TWO
					 helix_stability += BOTTOM_RIGHT

		end
		# end
		# Bottom left neighbor
		# puts aa + " " + nbl
		# puts $bottomLeftMatrix[aa][nbl] # Does this value fit to the real meant value???
		# if !$bottomLeftMatrix[aa][nbl].nil?
		case $bottomLeftMatrix[aa][nbl]
			when STRONG
				strong += BOTTOM_LEFT
			when MEDIUM
				medium += BOTTOM_LEFT
			when WEAK
				weak += BOTTOM_LEFT
			when LOW
				helix_stability += BOTTOM_LEFT
			when PENALTY_ONE
				helix_stability += BOTTOM_LEFT
			when PENALTY_TWO
				helix_stability += BOTTOM_LEFT
		end
		# end

		return strong, medium, weak, helix_stability
	end

	# Helper-method: Score triangle neighbors in SAH-grid
	# Version 3: 04.2014
	# @param [String] aa      the parent aminoacid
	# @param [String] nbl     the bottom left aminoacid
	# @param [String] nbr     the bottom right aminoacid
	# @return [Float] score   Summed neighbor score (nbl+nbr)
	def score_sah_grid_neighbors_v3(aa, nbl, nbr)
		score = 0.0
		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8
		# Strong (score=1.0)
		# K/R  => E/D : i+3 / i+4
		# Medium (score=0.75)
		# D  => K/R   : i+4
		# Weak   (score=0.5)
		# E  => K/R/Q : i+3 / i+4
		# D  => K/R   : i+3
		# Q  => E     : i+3 / i+4

		# E -> K/R/Q
		if aa == 'E'
			score += if ['K','R','Q'].rindex(nbr) != nil then 0.5 else 0 end
			score += if ['K','R','Q'].rindex(nbl) != nil then 0.5 else 0 end
		# K/R -> E/D
		elsif ['K','R'].rindex(aa) != nil
			score += if ['E','D'].rindex(nbr) then 1.0 else 0 end
			score += if ['E','D'].rindex(nbl) then 1.0 else 0 end
		# D -> K/R
		elsif aa == 'D'
			score += if ['K','R'].rindex(nbr) != nil then 0.5 else 0 end
			score += if ['K','R'].rindex(nbl) != nil then 0.75 else 0 end
		# Q -> E
		elsif aa == 'Q'
			score += if nbr == 'E' then 0.5 else 0 end
			score += if nbl == 'E' then 0.5 else 0 end
		end

		return score
	end

	# Helper-method: Score triangle neighbors in SAH-grid
	# Version 4: 28.09.2016
	# @param [String] aa      the parent aminoacid
	# @param [String] nbl     the bottom left aminoacid
	# @param [String] nbr     the bottom right aminoacid
	# @return [Float] score   Summed neighbor score (nbl+nbr)
	def score_sah_grid_neighbors_v4(aa, nbl, nbr)
		score = 0.0

		# Calculate strong/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8

		# SCORING
		# Strong (score = 1.0)
		# E  =>  K/R : i+4
		# H/K/R => E : i+4
		# Medium (score = 0.75)
		# E  =>  K/R : i+3
		# D  =>  K/R : i+4
		# H/K/R => E : i+3
		# Weak   (score = 0.5)
		# E   =>   H : i+3
		# D  =>  K/R : i+3
		# K/R  =>  D : i+3
		# Q   =>   E : i+4, i+3
		# E   =>   Q : i+4, i+3
		# Negative (score = -0.5)
		# L/I/V/F/Y/M => L/I/V/F/Y/M : i+4, i+3

		# E -> K/R, H, Q
		if aa == 'E'
			score += if ['K','R'].rindex(nbl) != nil then STRONG else 0 end
			score += if ['K','R'].rindex(nbr) != nil then MEDIUM else 0 end
			score += if ['H','Q'].rindex(nbr) != nil then WEAK else 0 end
			score += if                   nbl == 'Q' then WEAK else 0 end
		# H/K/R -> E
		elsif ['H','K','R'].rindex(aa) != nil
			score += if nbl == 'E' then STRONG else 0 end
			score += if nbr == 'E' then MEDIUM else 0 end
		# K/R -> D
		elsif ['K','R'].rindex(aa) != nil
			score += if nbr == 'D' then WEAK else 0 end
		# D -> K/R
		elsif aa == 'D'
			score += if ['K','R'].rindex(nbl) != nil then MEDIUM else 0 end
			score += if ['K','R'].rindex(nbr) != nil then WEAK else 0 end
		# Q -> E
		elsif aa == 'Q'
			score += if nbl == 'E' then WEAK else 0 end
			score += if nbr == 'E' then WEAK else 0 end
		# L/I/V/F/Y/M => L/I/V/F/Y/M
		elsif ['L','I','V','F','Y','M'].rindex(aa) != nil
			score += if ['L','I','V','F','Y','M'].rindex(nbl) != nil then NEGATIVE else 0 end
			score += if ['L','I','V','F','Y','M'].rindex(nbr) != nil then NEGATIVE else 0 end
		# Hydrophobe Interactions (additional)
		elsif ['A','G','P'].rindex(aa) != nil
			case aa
				when 'A'
					score += 1
				when 'G'
					score -= 2
				when 'P'
					score -= 3
			end
		end

		return score
	end

	# Helper-method: Score triangle neighbors in SAH-grid
	# Version 5: 27.09.2016
	# @param [String] aa      the parent aminoacid
	# @param [String] nbl     the bottom left aminoacid
	# @param [String] nbr     the bottom right aminoacid
	# Scoring-Matrices available as Gloabls
	# @param [Float] bottomRightMatrix  Bottom right interaction scoring-matrix
	# @param [Float] bottomLeftMatrix   Bottom left interaction scoring-matrix
	# @return [Float] score   Summed neighbor score (nbl+nbr)
	def score_sah_grid_neighbors(aa, nbl, nbr) #, bottomRightMatrix, bottomLeftMatrix)
		score = 0.0

		# Calculate strong/medium/weak interactions
		# Rules: top-right    = (i-4) = 1, top-left    = (i-3) = 2
		# 		   bottom-right = (i+3) = 4, bottom-left = (i+4) = 8

		# Scores from Scoring-Matrices (left/right)
		# puts $bottomLeftMatrix[aa].to_s
		# puts $bottomRightMatrix[aa].to_s
		# 	if !$bottomLeftMatrix[aa][nbl].nil? && !$bottomLeftMatrix[aa][nbr].nil?
		# 	puts "%.4f %.4f" % [$bottomLeftMatrix[aa][nbl], $bottomLeftMatrix[aa][nbr]]
		# end
		score += !$bottomLeftMatrix[aa][nbl].nil? ? $bottomLeftMatrix[aa][nbl] : 0.0
		score += !$bottomRightMatrix[aa][nbr].nil? ? $bottomRightMatrix[aa][nbr] : 0.0

		return score
	end

	# Helper-method to check if edge already exists
	# @param [Array] edges    String array of edges
	# @param [String] src     Source edge
	# @param [String] tar     Target edge
	# @return [Boolean]       Edge exists (true) or not (false)
	def edge_exists(edges, src, tar)
		if result = edges.assoc("#{src},#{tar}") != nil
			return true
		elsif result = edges.assoc("#{tar},#{src}") != nil
			return true
		else
			return false
		end
	end

	# Helper-method to convert interactionScore in number of edges
	# Makes use of 'Binary' number decomposition
	# @param [Integer] interactionScore   Computed edges score
	# @return [Integer] edges             Number of edges {0-4}
	def number_of_edges(interactionScore)
		edges = 0; div = 8
		while div > 0 do
			if (quot = interactionScore/div) >= 1
				interactionScore -= div
				edges += 1
			end
			div /= 2
		end
		return edges
	end

	# Helper-method to check drawing of specific SAH plot edge (specified by bin = {1,2,4,8})
	# @return [False] | string, int, int 		cssClass, interact, weak
	def draw_edge(bin, pos, interact, weak, category)
		cssClass = ''; result = false
		if interact-bin >= 0
			if weak-bin >= 0
				cssClass = ' potential'
				weak -= bin
			end
			if bin == 8 && category.name == 'f' then ; else result = cssClass end
			interact -= bin;
		end

		if result != false
			return [result, interact, weak]
		else
			return result
		end
	end

	# 'Abstract' helper method to check drawing of specific SAH plot edge (specified by bin = {1,2,4,8})
	# @param [Object] bin         identifier of testing edge
	# @param [Object] pos         position in sequence
	# @param [Object] interact    interaction code
	# @param [Object] weak        interaction code
	# @param [Object] category
	# @return [False] | string, int, int 		cssClass, interact, weak
	def draw_grid_edge(bin, pos, interact, weak, category)
		cssClass = ''; result = false
		if interact-bin >= 0
			if weak-bin >= 0
				cssClass = ' potential'
				weak -= bin
			end
			result = cssClass
			#if bin == 8 && category == 'f' then ; else result = cssClass end
			interact -= bin;
		end

		if result != false
			return [result, interact, weak]
		else
			return result
		end
	end

	# 'Abstract' helper method to check drawing of specific SAH plot edge (specified by bin = {1,2,4,8})
	# @param [Object] bin         identifier of testing edge
	# @param [Object] pos         position in sequence
	# @param [Object] interact    interaction code
	# @param [Object] weak        interaction code
	# @param [Object] middle      interaction code
	# @return [False] | string, int, int 		cssClass, interact, weak, middle
	def draw_greyscale_grid_edge(bin, pos, interact, weak, middle)
		cssClass = ''; result = false
		if interact-bin >= 0
			if weak-bin >= 0
				cssClass = ' potential'
				weak -= bin
			elsif middle-bin >= 0
				cssClass = ' more-potential'
				middle -= bin
			end
			result = cssClass
			#if bin == 8 && category == 'f' then ; else result = cssClass end
			interact -= bin;
		end

		if result != false
			return [result, interact, weak, middle]
		else
			return result
		end
	end

	# # DEPRECATED
	# # SVG: Helper-method to draw ALL interactions between AA in heptad-net plot
	# # @return [String]         shape - SVG markup (code)
	# def draw_all_interactions(coiledcoil, edges, category, entry, i, y, color, abs_start)
	# 	svg_data = ''; svg_data_f = ''; interact = 0
	# 	pos = entry[0].to_i; aa = entry[1]
	#
	# 	if ['D','E','K','R','Q'].rindex(aa) != nil
	# 		internal_pos = pos-coiledcoil.start
	# 		strong, weak = check_sah_neighbors(coiledcoil, internal_pos, aa)
	# 		interact = strong + weak
	# 		score = number_of_edges(strong) + number_of_edges(weak)*0.25
	#
	# 		if interact > 0
	# 			yEdge = -2+i*40
	# 			svg_data << '<g transform="rotate(45 36, '+ "#{y}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
	# 			# residue interaction bottom left (bottom)
	# 			# svg_data << draw_edge(8, pos, pos+4, strong, weak, edges, category)
	# 			bin = 8; tar = pos+4
	# 			if interact-bin >= 0 && !edge_exists(edges, pos, tar)
	# 				edges.push(["#{pos},#{tar}", pos, tar])
	# 				potential = if weak-bin >= 0 then " potential" else "" end
	# 				# print "Bottom left: " + potential + "\n"
	# 				if category.name != 'f'
	# 					svg_data << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
	# 				end
	# 				if potential != "" then weak -= bin end
	# 				# print "Weak: " + weak.to_s + "\n"
	# 				interact -= bin;
	# 			else
	# 				interact -= if interact-bin >= 0 then bin else 0 end;
	# 				weak -= if weak-bin >= 0 then bin else 0 end;
	# 			end
	# 			# residue interaction bottom right (right)
	# 			bin = 4; tar = pos+3
	# 			if interact-bin >= 0 && !edge_exists(edges, pos, tar)
	# 				edges.push(["#{pos},#{tar}", pos, tar])
	# 				potential = if weak-bin >= 0 then " potential" else "" end
	# 				# print "Bottom right: " + potential + "\n"
	# 				svg_data << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
	# 				if potential != "" then weak -= bin end
	# 				# print "Weak: " + weak.to_s + "\n"
	# 				interact -= bin;
	# 			else
	# 				interact -= if interact-bin >= 0 then bin else 0 end;
	# 				weak -= if weak-bin >= 0 then bin else 0 end;
	# 			end
	# 			# residue interaction top left (left)
	# 			bin = 2; tar = pos-3
	# 			if interact-bin >= 0 && !edge_exists(edges, pos, tar)
	# 				edges.push(["#{pos},#{tar}", pos, tar])
	# 				potential = if weak-bin >= 0 then " potential" else "" end
	# 				# print "Top left: " + potential + "\n"
	# 				if category.name != 'f'
	# 					svg_data << '<line x1="20" y1="58" x2="33" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
	# 				end
	# 				if potential != "" then weak -= bin end
	# 				# print "Weak: " + weak.to_s + "\n"
	# 				interact -= bin;
	# 			else
	# 				interact -= if interact-bin >= 0 then bin else 0 end;
	# 				weak -= if weak-bin >= 0 then bin else 0 end;
	# 			end
	# 			# residue interaction top right (top)
	# 			bin = 1; tar = pos-4
	# 			if interact-bin >= 0 && !edge_exists(edges, pos, tar)
	# 				edges.push(["#{pos},#{tar}", pos, tar])
	# 				potential = if weak-bin >= 0 then " potential" else "" end
	# 				# print "Top right: " + potential + "\n"
	# 				svg_data << '<line x1="41" y1="38" x2="41" y2="49" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
	# 				if potential != "" then weak -= bin end
	# 				# print "Weak: " + weak.to_s + "\n"
	# 				interact -= bin;
	# 			else
	# 				interact -= if interact-bin >= 0 then bin else 0 end;
	# 				weak -= if weak-bin >= 0 then bin else 0 end;
	# 			end
	# 			svg_data << '</g>'
	# 		end
	# 	end
	#
	# 	if category.name == 'f'
	# 		if interact-8 >= 0 && !edge_exists(edges, pos, pos+4)
	# 			# residue interaction bottom left (bottom)
	# 			edges.push(["#{pos},#{pos+4}", pos, pos+4])
	# 			potential = if weak-8 >= 0 then " potential" else "" end
	# 			svg_data_f << '<g transform="rotate(45 36, '+ "#{y}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
	# 			svg_data_f << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
	# 			svg_data_f << '</g>'
	# 			if potential != "" then weak -= 8 end
	# 			interact -= 8;
	# 		end
	# 		if interact-2 >= 0 && !edge_exists(edges, pos, pos-3)
	# 			# residue interaction top left (left)
	# 			edges.push(["#{pos},#{pos-3}", pos, pos-3])
	# 			potential = if weak-2 >= 0 then " potential" else "" end
	# 			svg_data_f << '<g transform="rotate(45 36, '+ "#{y}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
	# 			svg_data_f << '<line x1="25" y1="55" x2="36" y2="55" class="interaction' + potential + '"' + "#{color[:edge]}" + '/>'
	# 			# svg_data_f << '<line x1="20" y1="58" x2="33" y2="58" class="interaction' + potential + '"' + "#{color[:edge]}" + '/>'
	# 			svg_data_f << '</g>'
	# 			if potential != "" then weak -= 2 end
	# 			interact -= 2;
	# 		end
	# 	end
	#
	# 	return svg_data, svg_data_f, score
	# end

	# SVG: Helper-method to draw interactions between AA in heptad-net plot (only bottom ones)
	# @return [String]         shape - SVG markup (code)
	def draw_interactions(coiledcoil, category, entry, i, y, color)
		svg_data = ''; svg_data_f = ''; interact = weak = 0
		pos = entry[0].to_i; aa = entry[1]

		if ['D','E','K','R','Q'].rindex(aa) != nil
			domain_pos = pos-coiledcoil.start #
			strong, weak = check_sah_neighbors(coiledcoil, domain_pos, aa)
			interact = strong + weak
			score = number_of_edges(strong) + number_of_edges(weak)*0.25

			if interact > 0
				yEdge = -2+i*40
				svg_data << '<g transform="rotate(45 36, '+ "#{y}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
				# residue interaction bottom left (bottom)
				if (result = draw_edge(8, pos, interact, weak, category)) != false
					potential = result[0]; interact = result[1]; weak = result[2]
					svg_data << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
				end
				# residue interaction bottom right (right)
				if (result = draw_edge(4, pos, interact, weak, category)) != false
					potential = result[0]; interact = result[1]; weak = result[2]
					svg_data << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
				end

				svg_data << '</g>'
			end
		end

		if category.name == 'f'
			# residue interaction bottom right (right)
			if (result = draw_edge(4, pos, interact, weak, category)) != false
				potential = result[0]; interact = result[1]; weak = result[2]
				svg_data_f << '<g transform="rotate(45 36, '+ "#{y}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
				svg_data_f << '<line x1="46" y1="63" x2="46" y2="76" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
				svg_data_f << '</g>'
			end
		end

		return svg_data, svg_data_f, score
	end


	# SVG: Helper-method to calculate interaction-score for the heptad-net plot (only bottom connections)
	# param grid [array]    pre-ordered grid-array with aminoacid objects
	# return score [float]  interaction score for the heptad-net
	def calculate_grid_interaction_score(grid)
		interact = 0; strong = 0; weak = 0; score = 0; numRows = grid.size;
		heptadPos = { :a => 0, :b => 0, :c => 0, :d => 0, :e => 0, :f => 0, :g => -1 }
		total_strong = 0; total_middle = 0; total_weak = 0; newWeightScore = 0.0; hydrophob_score = 0;

		# Determine last 4 AA
		grid_helper = Array.new(grid.flatten(1))
		grid_helper.delete_if{ |el| el.nil? }.sort! { |a,b| a[0].to_i <=> b[0].to_i }

		# Travel by row
		grid.each_with_index { |row, y|

			# Debug
			#print row
			#print "\n"

			# Travel by column
			row.each_with_index { |entry, x|

				# TODO: Additional Last 4 AA !!! Has negative influence on short windows
				if entry != nil && y < numRows-1 && entry[0].to_i < grid_helper[-4][0].to_i
					aa = entry[1]; heptad = entry[2]; pos = entry[0].to_i

					# Version <= 3
					# if ['D','E','K','R','Q'].rindex(aa) != nil
					# Version 4
					if ['D','E','H','K','Q','R','L','I','V','F','Y','M','A','G','P'].rindex(aa) != nil
						if heptad != 'f'

							# Debug
							#print grid[y+1][x+1]
							#print " " + (y+1).to_s + "," + (x+1).to_s + " "
							#print grid[y+1][x-1]
							#print " " + (y-1).to_s + "," + (x+1).to_s + " "

							nbl = grid[y+1][x-1] == nil ? "" : grid[y+1][x-1][1]
							if heptad == 'c'
								nbr = grid[y][0] == nil ? "" : grid[y][0][1]
							else
								nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							end
							strong, middle, weak, helix_stability = classify_greyscale_sah_grid_neighbors(aa, nbl, nbr)
							newWeightScore += score_sah_grid_neighbors(aa, nbl, nbr)

							# Debug
							#print aa + "#{pos}: " + nbl  + " " + nbr # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"

						else

							# Debug
							#print grid[y+1][6]
							#print " " + (y+1).to_s + "," + (6).to_s + " "
							#print grid[y+2][x+1]
							#print " " + (y+2).to_s + "," + (x+1).to_s + " "

							nbl = (y+2) < numRows-1 && grid[y+2].size == 7 && grid[y+2][6] != nil ? grid[y+2][6][1] : ""
							nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							strong, middle, weak, helix_stability = classify_greyscale_sah_grid_neighbors(aa, nbl, nbr)
							newWeightScore += score_sah_grid_neighbors(aa, nbl, nbr)

							# Debug
							#print aa + "#{pos}: " + nbl  + " " + nbr + " (Heptad f)" # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"

						end
						total_strong += number_of_edges(strong); total_middle += number_of_edges(middle); total_weak += number_of_edges(weak);
						# interact = strong + weak
						# Old score
						# score += number_of_edges(strong) + number_of_edges(weak)*0.25

					# Hydrophobe Interactions (additional)
					# elsif ['V','I','M','L','Y','P','F','G','A'].rindex(aa) != nil
					# Version 3
					# elsif ['V','I','M','L','Y'].rindex(aa) != nil
					# Version 4
					# elsif ['A','G','P'].rindex(aa) != nil
					# 	case aa
					# 		when 'A'
					# 			hydrophob_score += 1
					# 		when 'G'
					# 			hydrophob_score -= 2
					# 		when 'P'
					# 			hydrophob_score -= 3
					# 	end
					end
				end
		} }
		# puts "New weighted score: "
		# puts newWeightScore
		# puts "Old score: "
		# puts score
		# puts "Hydrophob penalty: " + hydrophob_score.to_s

		# if (newWeightScore - hydrophob_score) >= 0
		# newWeightScore -= hydrophob_score
		# else
		# 		newWeightScore = 0
		# end

		return newWeightScore, total_strong, total_middle, total_weak
	end

	# SVG: Helper-method to calculate surrounding interaction-score for the heptad-net plot (only bottom connections)
	# @param [CoiledCoilSequence]   coiledcoilseq
	# @param [CoiledCoilDomain]     coiledcoil
	# @return [Float, Float, Float] Score of surrounding interactions
	def calculate_surrounding_interactions(coiledcoilseq, coiledcoil)
			expand = 4 # Expand the heptad-net for outlying interactions
		# 	puts "start, end = " + (coiledcoil.start-expand).to_s + ", " + (coiledcoil.end+expand).to_s
			env_coiledcoil = coiledcoilseq.subseq(coiledcoil.start-expand, coiledcoil.end+expand, true)
			gridMat = calculate_simple_sah_grid(env_coiledcoil, true)
			catGrid = transpose_simple_grid(gridMat) # ordered by categories
			grid = convert_to_interaction_grid(env_coiledcoil, catGrid) # start matrix
			envGridScore, envStrong, envWeak = calculate_grid_interaction_score(grid)
		# 	puts "strong, weak = " + envStrong.to_s + ", " + envWeak.to_s
			return envGridScore, envStrong, envWeak
	end

	# SVG: Helper-method to draw interactions between AA in heptad-net plot (only bottom ones)
	# @return [String, String]         shape - SVG markup (code)
	def draw_grid_interactions(coiledcoil, grid)
		svg_data_f = ''; color = { :edge => '' }
		svg_data = {:a => '', :b => '', :c => '', :d => '', :e => '', :f => '', :g => '' }
		interact = 0; strong = 0; weak = 0; score = 0; numRows = grid.size;
		heptadPos = { :a => 0, :b => 0, :c => 0, :d => 0, :e => 0, :f => 0, :g => -1 }
		total_strong = 0; total_weak = 0

		grid.each_with_index { |row, y|
			# print row
			# print "\n"
			row.each_with_index { |entry, x|

				if entry != nil && y < numRows-1
					aa = entry[1]; heptad = entry[2]; pos = entry[0].to_i
					i = heptadPos[:"#{heptad}"]+y/2
					yEdge  = -2+i*40; yCoord = 64+i*40 	# y-coordinate of text-node 'amino-acid'

					# Only take edges from parent nodes 'D','E','K','R','Q' into concern
					if ['E','D','K','R','Q'].rindex(aa) != nil
						if heptad != 'f'
							#print grid[y+1][x+1]
							#print " " + (y+1).to_s + "," + (x+1).to_s + " "
							#print grid[y+1][x-1]
							#print " " + (y-1).to_s + "," + (x+1).to_s + " "
							nbl = grid[y+1][x-1] == nil ? "" : grid[y+1][x-1][1]
							if heptad == 'c'
								nbr = grid[y][0] == nil ? "" : grid[y][0][1]
							else
								nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							end
							strong, weak = classify_sah_grid_neighbors(aa, nbl, nbr)

							#print aa + "#{pos}: " + nbl  + " " + nbr # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"
						else
							#print grid[y+1][6]
							#print " " + (y+1).to_s + "," + (6).to_s + " "
							#print grid[y+2][x+1]
							#print " " + (y+2).to_s + "," + (x+1).to_s + " "
							nbl = (y+2) < numRows-1 && grid[y+2].size == 7 && grid[y+2][6] != nil ? grid[y+2][6][1] : ""
							nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							strong, weak = classify_sah_grid_neighbors(aa, nbl, nbr)

							#print aa + "#{pos}: " + nbl  + " " + nbr + " (Heptad f)" # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"
						end
						total_strong += number_of_edges(strong); total_weak += number_of_edges(weak);
						interact = strong + weak
						# score += number_of_edges(strong) + number_of_edges(weak)*0.25
					end

					if interact > 0
						if heptad != 'f'
								svg_data[:"#{heptad}"] << "\n" + '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
								# residue interaction bottom left (bottom)
								if (result = draw_grid_edge(8, pos, interact, weak, heptad)) != false
									potential = result[0]; interact = result[1]; weak = result[2]
									# svg_data << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
									svg_data[:"#{heptad}"] << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								end
								# residue interaction bottom right (right)
								if (result = draw_grid_edge(4, pos, interact, weak, heptad)) != false
									potential = result[0]; interact = result[1]; weak = result[2]
									# svg_data << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
									svg_data[:"#{heptad}"] << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								end
								svg_data[:"#{heptad}"] << '</g>' + "\n"
						else
							# residue interaction bottom left (bottom) - for repeated f-group on the right
							if (result = draw_grid_edge(8, pos, interact, weak, heptad)) != false
								potential = result[0]; interact = result[1]; weak = result[2]
								svg_data_f << '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
								svg_data_f << '<line x1="46" y1="63" x2="46" y2="76" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data_f << '</g>'
							end
							# residue interaction bottom right
							if (result = draw_grid_edge(4, pos, interact, weak, heptad)) != false
								potential = result[0]; interact = result[1]; weak = result[2]
								svg_data[:"#{heptad}"] << "\n" + '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
								svg_data[:"#{heptad}"] << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data[:"#{heptad}"] << '</g>' + "\n"
							end
						end
					end
				end
			}
		}

		return svg_data, svg_data_f #, score, total_strong, total_weak
	end

	# SVG: Helper-method to draw interactions between AA in heptad-net plot (only bottom ones)
	# @return [String, String]         shape - SVG markup (code)
	def draw_greyscale_grid_interactions(coiledcoil, grid)
		svg_data_f = ''; color = { :edge => '' }
		svg_data = {:a => '', :b => '', :c => '', :d => '', :e => '', :f => '', :g => '' }
		interact = 0; strong = 0; middle = 0; weak = 0; score = 0; numRows = grid.size;
		heptadPos = { :a => 0, :b => 0, :c => 0, :d => 0, :e => 0, :f => 0, :g => -1 }
		reverseHeptadPos = { 0 => 'f', 1 => 'b', 2 => 'e', 3 => 'a', 4 => 'd', 5 => 'g', 6 => 'c'}
		total_strong = 0; total_middle = 0; total_weak = 0;

		grid.each_with_index { |row, y|
			# print row
			# print "\n"
			row.each_with_index { |entry, x|

				if entry != nil && y < numRows-1
					# aa = entry[1]; heptad = entry[2]; pos = entry[0].to_i
					aa = entry[1]; heptad = reverseHeptadPos[x]; pos = entry[0].to_i
					i = heptadPos[:"#{heptad}"]+y/2
					yEdge  = -2+i*40; yCoord = 64+i*40 	# y-coordinate of text-node 'amino-acid'

					# Only take edges from parent nodes 'D','E','K','R','Q' into concern
					if ['E','D','K','R','Q'].rindex(aa) != nil
						if heptad != 'f'
							#print grid[y+1][x+1]
							#print " " + (y+1).to_s + "," + (x+1).to_s + " "
							#print grid[y+1][x-1]
							#print " " + (y-1).to_s + "," + (x+1).to_s + " "
							nbl = grid[y+1][x-1] == nil ? "" : grid[y+1][x-1][1]
							if heptad == 'c'
								nbr = grid[y][0] == nil ? "" : grid[y][0][1]
							else
								nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							end
							strong, middle, weak, helix_stability = classify_greyscale_sah_grid_neighbors(aa, nbl, nbr)

							#print aa + "#{pos}: " + nbl  + " " + nbr # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"
						else
							#print grid[y+1][6]
							#print " " + (y+1).to_s + "," + (6).to_s + " "
							#print grid[y+2][x+1]
							#print " " + (y+2).to_s + "," + (x+1).to_s + " "
							nbl = (y+2) < numRows-1 && grid[y+2].size == 7 && grid[y+2][6] != nil ? grid[y+2][6][1] : ""
							nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1][1]
							strong, middle, weak, helix_stability = classify_greyscale_sah_grid_neighbors(aa, nbl, nbr)

							#print aa + "#{pos}: " + nbl  + " " + nbr + " (Heptad f)" # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"
						end
						total_strong += number_of_edges(strong); total_middle += number_of_edges(middle); total_weak += number_of_edges(weak);
						interact = strong + middle + weak
						# score += number_of_edges(strong) + number_of_edges(weak)*0.25
					end

					if interact > 0
						if heptad != 'f'
							svg_data[:"#{heptad}"] << "\n" + '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
							# residue interaction bottom left (bottom)
							if (result = draw_greyscale_grid_edge(8, pos, interact, weak, middle)) != false
								potential = result[0]; interact = result[1]; weak = result[2]; middle = result[3]
								# svg_data << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data[:"#{heptad}"] << '<line x1="41" y1="67" x2="41" y2="78" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
							end
							# residue interaction bottom right (right)
							if (result = draw_greyscale_grid_edge(4, pos, interact, weak, middle)) != false
								potential = result[0]; interact = result[1]; weak = result[2]; middle = result[3]
								# svg_data << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data[:"#{heptad}"] << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
							end
							svg_data[:"#{heptad}"] << '</g>' + "\n"
						else
							# residue interaction bottom left (bottom) - for repeated f-group on the right
							if (result = draw_greyscale_grid_edge(8, pos, interact, weak, middle)) != false
								potential = result[0]; interact = result[1]; weak = result[2]; middle = result[3]
								svg_data_f << '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
								svg_data_f << '<line x1="46" y1="63" x2="46" y2="76" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data_f << '</g>'
							end
							potential = ''
							# residue interaction bottom right
							if (result = draw_greyscale_grid_edge(4, pos, interact, weak, middle)) != false
								potential = result[0]; interact = result[1]; weak = result[2]; middle = result[3]
								svg_data[:"#{heptad}"] << "\n" + '<g transform="rotate(45 36, '+ "#{yCoord}" +') translate(-4,' + "#{yEdge}" + ')" id="pathX_interaction">'
								svg_data[:"#{heptad}"] << '<line x1="48" y1="58" x2="62" y2="58" class="interaction' + potential + '" ' + "#{color[:edge]}" + '/>'
								svg_data[:"#{heptad}"] << '</g>' + "\n"
							end
						end
					end
				end
			}
		}

		return svg_data, svg_data_f #, score, total_strong, total_weak
	end

	# DEPRECATED
	# SVG: Helper-method to create positioning shape
	# @param [Array] pos		top-left corner of main shape rec-angle
	# @param [Integer]   precise	flag, start-position should be chosen absolute or relative in coiledcoil pattern
	# @return [String]         shape - SVG markup (code)
	def draw_sah_aminoacids(coiledcoil)
		svg_tmpl = svg_tmpl_conn = IO.read("#{Rails.root}/lib/assets/sah_domain_elements.svg")
		categories = coiledcoil.get_categories
		sequence = coiledcoil.sequence
		totalScore = 0

		# DEBUG
		# print sequence.join + "\n"
		# categories.each { |item|
		# 	item.elements.each { |pos|
		# 		print pos
		# 		print "\n"
		# 	}
		# 	print "\n"
		# }

		# get (absolute) start position from category a
		abs_start = categories[0].elements[0][0].to_i
		categories.each { |category|
			svg_data = ''; svg_data_f = ''; k = 0
			# while category.elements[k][0].to_i < abs_start do
			# 	k += 1
			# end
			category.elements[(0+k)..(6+k)].each_with_index { |entry, i|
				pos = entry[0]; aa = entry[1]; cls = entry[2]
				y = 64+i*40 	# y-coordinate of text-node 'amino-acid'

				# regulate colors of amino-acids
				# color = { :node => calc_aa_color(aa, "sah"), :edge => '' }
				color = { :node => 'lightgrey', :edge => 'style="fill:#777; stroke:#777;"' }
				if entry[2] == category.name
					color[:node] = calc_aa_color(aa, "sah")
					color[:edge] = ''
				end

				# create amino-acid text-nodes
				if ['a','d'].rindex(category.name) != nil  
					svg_data <<	'<rect x="30" y="28" width="16" height="16" class="even" stroke-width="0px"/>'
					svg_data << '<a xlink:title="' + aa + pos + " [#{cls}]" '"><text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">' + aa + '</text></a>'
				else
					svg_data << '<a xlink:title="' + aa + pos + " [#{cls}]" '"><text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">' + aa + '</text></a>'
				end					
				if category.name == 'f'
					svg_data_f << '<text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">(' + aa + ')</text>'
				end
				# regulate sub-index of aa (only first)
				if abs_start == pos.to_i
					# svg_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey ' + "#{color[:node]}" + ' damntiny" x="45" y="'+ "#{y+4}" +'">' + pos + '</text>'
					svg_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey damntiny" x="45" y="'+ "#{y+4}" +'">' + pos + '</text>'
				end

				# draw connections
				svg_conn, svg_conn_f, score = draw_interactions(coiledcoil, category, entry, i, y, color)
				svg_data << svg_conn
				svg_data_f << svg_conn_f
				totalScore += score.to_f
			}

			name = "[CATEGORY_#{category.name.capitalize}]"
			svg_tmpl.gsub!(name, svg_data)
			if category.name == 'f'
				name = "[CATEGORY_(#{category.name.capitalize})]"
				svg_tmpl.gsub!(name, svg_data_f)
			end
		}

		# calculate score for subsequence in heptad net
		print "SAH length: " + ((coiledcoil.end+1-abs_start)*2-7).to_s + "\n"
		maxScore = (coiledcoil.end+1-abs_start)*2-7
		# print "Max. Score: " + maxScore.to_s + "\n"
		# print "Tot. Score: " + totalScore.to_s + "\n"
		totalScore /= maxScore

		return svg_tmpl, totalScore
	end

	# SVG: Helper-method to create positioning shape
	# @param [Array] pos		top-left corner of main shape rec-angle
	# @param [Integer]   precise	flag, start-position should be chosen absolute or relative in coiledcoil pattern
	# @return [String]         shape - SVG markup (code)
# 	def draw_category_grid_sah_aminoacids(coiledcoil)
# 		svg_tmpl = svg_tmpl_conn = IO.read("#{Rails.root}/lib/assets/sah_domain_elements.svg")
# 		position = Hash["a" => [3,0], "b" => [1,0], "c" => [6,1], "d" => [4,1], "e" => [2,1], "f" => [0,1], "g" => [5,2]]
# 		grid = Array.new(15) { Array.new(7) }; edges = []; totalScore = 0.0; skip = 0;
# 		svg_data = {:a => '', :b => '', :c => '', :d => '', :e => '', :f => '', :f => '', :g => ''}; svg_tmp_data_f = ''
# 		sequence = coiledcoil.sequence
# 		categories = coiledcoil.get_categories
#
# 		# DEBUG
# 		# print sequence.join + "\n"
# 		# categories.each { |item|
# 		# 	item.elements.each { |pos|
# 		# 		print pos
# 		# 		print "\n"
# 		# 	}
# 		# 	print "\n"
# 		# }
#
# 		# get (absolute) start position from category a
# 		abs_start = categories[0].elements[0][0].to_i
# 		categories.each { |category|
#
# 			puts JSON.dump(category)
#
# 			svg_tmp_data = ''; k = 0; skip = 0
# 			while category.elements[k][0].to_i < abs_start do
# 				k += 1
# 			end
# 			category.elements[(0+k)..(6+k)].each_with_index { |entry, i|
# 				pos = entry[0]; aa = entry[1]; cls = entry[2]
#
# 				# Calculate skip of node
# 				#puts (entry[0].to_i-7).to_s + " > " + category.elements[i+skip-1][0]
# 				if i>0 && entry[0].to_i-category.elements[i-1][0].to_i > 7
# 					skip += 1
# 				end
# 				#puts "skip: " + skip.to_s
# 				y = 64+(i+skip)*40 	# y-coordinate of text-node 'amino-acid'
#
# 				# regulate colors of amino-acids
# 				# color = { :node => calc_aa_color(aa, "sah"), :edge => '' }
# 				color = { :node => 'lightgrey', :edge => 'style="fill:#777; stroke:#777;"' }
# 				if entry[2] == category.name
# 					color[:node] = calc_aa_color(aa, "sah")
# 					color[:edge] = ''
# 				end
#
# 				# create amino-acid text-nodes
# 				svg_tmp_data << '<a xlink:title="' + aa + pos + " [#{cls}]" '"><text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">' + aa + '</text></a>'
# 				if category.name == 'f'
# 					svg_tmp_data_f << '<text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">(' + aa + ')</text>'
# 				end
# 				# regulate sub-index of aa (only first)
# 				if abs_start == pos.to_i
# 					#svg_tmp_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey ' + "#{color[:node]}" + ' damntiny" x="45" y="'+ "#{y+4}" +'">' + pos + '</text>'
# 					svg_tmp_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey damntiny" x="45" y="'+ "#{y+4}" +'">' + pos + '</text>'
# 				end
#
# 				# GRID representation: grid[row][col]
# 				if position[category.name][1]+(skip+k+i)*2 < 15
# 					grid[position[category.name][1]+(i+k+skip)*2][position[category.name][0]] = entry
# 				end
# 			}
# 			# print grid
# 			svg_data[:"#{category.name}"] = svg_tmp_data
# 		}
#
# 		# draw connections
# 		svg_conn, svg_conn_f, score = draw_grid_interactions(coiledcoil, grid)
# 		svg_conn.each_pair { |key, value|
# 			svg_data[:"#{key}"] << value
# 			name = "[CATEGORY_#{key.capitalize}]"
# 			svg_tmpl.gsub!(name, svg_data[:"#{key}"])
# 			if key == :f
# 				name = "[CATEGORY_(#{key.capitalize})]"
# 				svg_tmpl.gsub!(name, svg_tmp_data_f + svg_conn_f)
# 			end
# 		}
# 		totalScore = score.to_f
#
# 		# calculate score for subsequence in heptad net
# 		# print "SAH length: " + (coiledcoil.end+1-abs_start).to_s + "\n"
# 		maxScore = (coiledcoil.end+1-abs_start)*2-7
# 		# print "Max. Score: " + maxScore.to_s + "\n"
# 		# print "Tot. Score: " + totalScore.to_s + "\n"
# 		totalScore /= maxScore
#
# 		return svg_tmpl, totalScore
# 	end

	# SVG: Helper-method to create positioning shape
	# @param [Array] pos		top-left corner of main shape rec-angle
	# @param [Integer]   precise	flag, start-position should be chosen absolute or relative in coiledcoil pattern
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [String]         shape - SVG markup (code)
	def draw_grid_sah_aminoacids(coiledcoil, fullsequence, abs_start, gapMode)
		svg_tmpl = svg_tmpl_conn = IO.read("#{Rails.root}/lib/assets/sah_domain_elements.svg"); max_end_coord = 0
		position = Hash["a" => [3,0], "b" => [1,0], "c" => [6,1], "d" => [4,1], "e" => [2,1], "f" => [0,1], "g" => [5,2]]

		# grid = Array.new(15) { Array.new(7) }; edges = []; totalScore = 0.0;
		svg_data = {:a => '', :b => '', :c => '', :d => '', :e => '', :f => '', :f => '', :g => ''}; svg_tmp_data_f = ''

		categories = coiledcoil.get_categories

		# get (absolute) start position from category a
		#abs_start = categories[0].elements[0][0].to_i
		#abs_start = coiledcoil.start

		# Prepare interaction data
		shift = abs_start - coiledcoil.start - 1 # Add shift-difference to endCoord
		shift = shift + 4                        # Add last line / half heptad, to see bottom interactions
		seqPart = fullsequence.subseq(abs_start, coiledcoil.end+shift, true)

		if gapMode == false
				gridMat = calculate_simple_gapless_sah_grid(seqPart, true)
		else
				gridMat = calculate_simple_sah_grid(seqPart, true)
		end
		catGrid = transpose_simple_grid(gridMat) # ordered by categories
		grid = convert_to_interaction_grid(coiledcoil, catGrid) # start matrix
		intGraph = convert_to_interaction_graph(grid) # interaction graph
		cumulPathScore = calculate_cumulative_path_score(intGraph)

		# Render SVG-output
		catGrid.each_with_index{ |categoryVec, line_idx|
			category = categories[line_idx]; svg_tmp_data = '';
			categoryVec.each_with_index{ |entry, entry_idx|
				break if entry_idx >= 7
				if entry != nil
					max_end_coord = max_end_coord < entry[0].to_i ? entry[0].to_i : max_end_coord
					pos = entry[0]; aa = entry[1]; cls = entry[2]
					y = 64+(entry_idx)*40 	# y-coordinate of text-node 'amino-acid'

					# regulate colors of amino-acids
					color = { :node => 'lightgrey', :edge => 'style="fill:#777; stroke:#777;"' }
					## Color aminoacids in grey, when out of pattern
				# 	if cls == category.name
						color[:node] = calc_aa_color(aa, "sah") # use right aa-colors when in correct motif
						color[:edge] = ''
				# 	end

					# create amino-acid text-nodes
					if ['a','d'].rindex(category.name) != nil  
						svg_tmp_data <<	'<rect transform="rotate(10 36,'+ "#{y}" +')" x="34" y="' + "#{y-14}" + '" width="16" height="16" rx="4" class="core" />'
						svg_tmp_data << '<a xlink:title="' + aa + pos + " [#{cls}]" '"><text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">' + aa + '</text></a>'
					else
						svg_tmp_data << '<a xlink:title="' + aa + pos + " [#{cls}]" '"><text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">' + aa + '</text></a>'
					end					
					if category.name == 'f'
						svg_tmp_data_f << '<text transform="rotate(10 36,'+ "#{y}" +')" class="light_grey ' + "#{color[:node]}" + '" x="36" y="'+ "#{y}" +'">(' + aa + ')</text>'
					end
					# regulate sub-index of aa (only first)
					# if abs_start == pos.to_i
					  # Index beyond
					  # svg_tmp_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey ' + "#{color[:node]}" + ' damntiny" x="45" y="'+ "#{y+4}" +'">' + pos + '</text>'
					  # Index above
						svg_tmp_data << '<text transform="rotate(10 36,'+ "#{y+4}" +')" class="light_grey damntiny" x="48" y="'+
								"#{y-8}" +'">' + pos + '</text>'
					# end
				end
			}
			svg_data[:"#{category.name}"] = svg_tmp_data
		}

		# Validation output
		# puts "\n\n\nInteraction Grid\n"
		# grid.each{ |line|
		# 	print line
		# 	puts "\n"
		# }

		# Draw connections
		svg_conn, svg_conn_f = draw_greyscale_grid_interactions(coiledcoil, grid)
		svg_conn.each_pair { |key, value|
			svg_data[:"#{key}"] << value
			name = "[CATEGORY_#{key.capitalize}]"
			svg_tmpl.gsub!(name, svg_data[:"#{key}"])
			if key == :f
				name = "[CATEGORY_(#{key.capitalize})]"
				svg_tmpl.gsub!(name, svg_tmp_data_f + svg_conn_f)
			end
		}

		# Calculate new interaction score -> normalized Score (Interactions+Paths)
		summedScore, strong, middle, weak = calculate_grid_interaction_score(grid)
		# envSummedScore, envStrong, envWeak = calculate_surrounding_interactions(fullsequence, seqPart)
		puts "summedScore = " + summedScore.to_s + "; cumulPathScore = " + cumulPathScore.to_s
		normScore = calculate_normalized_interaction_score(seqPart, summedScore.to_f, cumulPathScore.to_f)

		# Show number of visible interactions
		return svg_tmpl, normScore, strong, middle, weak, max_end_coord
	end


	# SVG: Helper-method to calculate normalized Heptad-Net score
	# @param [CoiledCoilDomain] coiledcoil   CoiledCoilDomain object with prediction
	# @param [Float] summedScore             Summed score for the current heptad net window
	# return [Float]                         Normalized score over the window
	def calculate_normalized_interaction_score(coiledcoil, summedScore, cumulPathScore)
		# calculate score for visible subsequence in heptad net
		# print "SAH length: " + ((coiledcoil.end+1-coiledcoil.start)*2-7).to_s + "\n"
		# puts "#{coiledcoil.start} - #{coiledcoil.end} = #{coiledcoil.end-coiledcoil.start}"

		# Dependent on window-size
		heptadRepeats = $window_size.divmod(7)
		# Dependent on coiled-coil-size
		# heptadRepeats = (coiledcoil.end-coiledcoil.start).divmod(7)

		# Maximum interconnections of n-times full heptads (only down-pointing interactions)
		# strong = 1.0; middle = 0.75; weak = 0.5;
		#  10.5 = 3 x strong, 3 x medium, 3 x ll-network, 3 x rr-network
		# maxScore = (heptadRepeats[0]*10.5) + 2.5

		maxSumScoreValues = [73.5, 42.0, 31.5, 21.0]
		# maxSumScoreValues = [76.0, 44.5, 34.0, 23.5]
		# maxSumScoreValues = [75.25, 43.75, 33.25, 22.75]
		# Measured values
		# maxEnvScoreValues = [88.25, 57.75, 46.25, 36.75]

		# Maximum cumulPathScore derived from full-interconnected heptad-net: 49 | 28 | 21 | 14
		maxPathScoreValues = [73.5, 42.0, 31.5, 21.0]

		# Determine maximum values for the windows
		case heptadRepeats[0]
				when 6..8
						# maxScore     = maxEnvScoreValues[0]
						maxScore     = maxSumScoreValues[0]
						maxPathScore = maxPathScoreValues[0]
				when 4..5
						# maxScore     = maxEnvScoreValues[1]
						maxScore     = maxSumScoreValues[1]
						maxPathScore = maxPathScoreValues[1]
				when 3
						# maxScore     = maxEnvScoreValues[2]
						maxScore     = maxSumScoreValues[2]
						maxPathScore = maxPathScoreValues[2]
				when 0..2
						# maxScore     = maxEnvScoreValues[3]
						maxScore     = maxSumScoreValues[3]
						maxPathScore = maxPathScoreValues[3]
				else
						# maxScore     = maxEnvScoreValues[0]
						maxScore     = maxSumScoreValues[0]
						maxPathScore = maxPathScoreValues[0]
		end

		# Normalize both scores
		normCumulPathScore = (cumulPathScore / maxPathScore).to_f
		normNetScore       = (summedScore / maxScore).to_f
		normedScore        = ((normNetScore + normCumulPathScore)/2).to_f

		# puts "%d-%d: %.4f %.4f " % [coiledcoil.start, coiledcoil.end, summedScore, cumulPathScore]
		# puts "%.4f %.4f " % [normNetScore, normCumulPathScore]
		# puts "%.4f" % cumulPathScore

		return normedScore >= 0 ? normedScore : 0.0
	end


	# SVG: Helper-method to create positioning shape
	# @param CoiledCoilDomain coiledcoil      CoiledCoilDomain object with prediction
	# @param CoiledCoilDomain fullsequence	  Full CoiledCoilDomain object for whole sequence
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# return string svg_data                  SVG markup (code) to integrate in SAH-SVG-template
	def draw_sah_overview(coiledcoil, fullsequence, startCoord, gapMode)
		position = Hash["a" => 0, "b" => 1, "c" => 2, "d" => 3, "e" => 4, "f" => 5, "g" => 6]
		endCoord = coiledcoil.end; len = endCoord-startCoord; line_index = 0;
		# Set coords for overview
		heptadPos = position[fullsequence.coiledcoil_data[startCoord][2]]; preSeqSize = heptadPos-1+14
		windowStart = startCoord-preSeqSize > 0 ? startCoord-preSeqSize : 1; # start two rows before (14aa) with heptad 'a'
		windowEnd = startCoord+(7*16) >= fullsequence.end ? fullsequence.end : startCoord+(7*16)-1; # fill until end of graphic

		# Get overview cutaway from fullsequence with gapless/break-based aa distribution
		seqPart = fullsequence.subseq(windowStart, windowEnd, true)
		if gapMode == false
			gridMat = calculate_simple_gapless_sah_grid(seqPart, true)
		else
			gridMat = calculate_simple_sah_grid(seqPart, true)
		end
		heptads = gridMat.size

		# Render SVG-output
		svg_data = '<g transform="translate(0,16)">'
		# Render aminoacids depending on their categories
		gridMat.each_with_index{ |line, line_idx|
			line.each_with_index{ |node, node_idx|
				if node != nil
					line_index = line_idx if node[0].to_i == startCoord
					svg_data <<	'<a xlink:title="' + node[1] + node[0] + " [#{node[2]}]" '"><text class="tiny ' + calc_aa_color(node[1], 'sah') + '" x="' + (20+12*node_idx).to_s + '" y="24">' + node[1] + '</text></a>'
				end
			}
			svg_data <<	'</g><g transform="translate(0,' + (32+16*line_idx).to_s + ')">'
		}
		svg_data << '</g>'

		# Flags for window-range
		# Right-positioned start-flag
		# svg_data <<	'<rect x="102" y="28" rx="4" ry="4" width="26" height="16" class="even" stroke-width="0px"/>'
		# svg_data <<	'<rect x="102" y="28" width="4" height="16" class="even" stroke-width="0px"/>'
		# svg_data <<	'<text x="106" y="39" class="lightgrey damntiny">' + (windowStart).to_s + '</text>'
		# Left-positioned start-flag
		svg_data <<	'<rect x="-9"  y="28" rx="4" ry="4" width="26" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<rect x="13"  y="28" width="4" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<text x="-5" y="39" class="lightgrey damntiny">' + (windowStart).to_s + '</text>'
		# End flag
		y = 12+16*heptads
		svg_data <<	'<rect x="102" y="' + (y).to_s + '" rx="4" ry="4" width="26" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<rect x="102" y="' + (y).to_s + '" width="4" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<text x="106" y="' + (y+11).to_s + '" class="lightgrey damntiny">' + (windowEnd).to_s + '</text>'

		window_y = 26+16*line_index
		# Draw subset
		svg_data <<	'<line x1="104" y1="' + (window_y).to_s + '" x2="176" y2="38" class="lines" stroke-width="2px"/>'
		svg_data <<	'<line x1="104" y1="' + (window_y+114+1).to_s + '" x2="228" y2="326" class="lines" stroke-width="2px"/>'''
		svg_data <<	'<polygon points="104,' + (window_y+3).to_s + ' 170,44 216,308 104,' + (window_y+112).to_s + '" style="fill: url(#thin_stripes); fill-opacity: 0.2;"/>'
		# Frame and position for selected SAH-window labeling flag
		svg_data <<	'<rect x="-9"  y="' + (window_y+1).to_s + '" rx="4" ry="4" width="26" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<rect x="13"  y="' + (window_y+1).to_s + '" width="4" height="16" class="even" stroke-width="0px"/>'
		svg_data <<	'<text class="lightgrey damntiny" x="-5" y="' + (window_y+12).to_s + '">' + (startCoord).to_s + '</text>'
		svg_data << '<rect x="16"  y="' + (window_y).to_s + '" rx="4" ry="4" width="88" height="114" class="sahFrame"/>'

		return svg_data
	end

	# Helper-method to convert categorized grid-matrix (catMat) to interaction-graph
	# @param array  interactionGrid       data matrix to convert
	# return [Graph]                      Interaction graph
	def convert_to_interaction_graph(interactionGrid)
		grid = interactionGrid; numRows = interactionGrid.size
		weights = {:strong => 1, :weak => 0.25}; nbr = nil; nbl = nil
		graph = Graph.new

		# # Determine last 4 AA in seqPart (ordered)
		grid_helper = Array.new(grid.flatten(1))
		grid_helper.delete_if{ |el| el.nil? }.sort! { |a,b| a[0].to_i <=> b[0].to_i }

		# Iterate over grid-matrix
		grid.each_with_index { |row, y|
			row.each_with_index { |parent, x|
				# TODO: Additional Last 4 AA !!! Has negative influence on short windows
				if parent != nil && y < numRows-1 && parent[0].to_i < grid_helper[-4][0].to_i
					pos = parent[0].to_i; aa = parent[1]; heptad = parent[2]; prntLbl = aa + parent[0];

					# # Version <= 3 Only take edges from parent nodes 'D','E','K','R','Q' into concern
					# if ['D','E','K','R','Q'].rindex(aa) != nil
					# Version 4
					if ['D','E','H','K','Q','R','L','I','V','F','Y','M','A','G','P'].rindex(aa) != nil
						# Normal case
						if heptad != 'f'

							# Debug
							#print grid[y+1][x+1]
							#print " " + (y+1).to_s + "," + (x+1).to_s + " "
							#print grid[y+1][x-1]
							#print " " + (y-1).to_s + "," + (x+1).to_s + " "

							# Get neighbor positions
							nbl = grid[y+1][x-1] == nil ? "" : grid[y+1][x-1]
							if heptad == 'c'
								nbr = grid[y][0] == nil ? "" : grid[y][0]
							else
								nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1]
							end

							# Debug
							#print aa + "#{pos}: " + nbl  + " " + nbr # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"

						# Left edge-case: Doubled "f"-heptad-position
						else
							# Debug
							#print grid[y+1][6]
							#print " " + (y+1).to_s + "," + (6).to_s + " "
							#print grid[y+2][x+1]
							#print " " + (y+2).to_s + "," + (x+1).to_s + " "

							# Get neighbor positions
							nbl = (y+2) < numRows-1 && grid[y+2].size == 7 && grid[y+2][6] != nil ? grid[y+2][6] : ""
							nbr = grid[y+1][x+1] == nil ? "" : grid[y+1][x+1]

							# Debug
							#print aa + "#{pos}: " + nbl  + " " + nbr + " (Heptad f)" # + "\n"
							#print " s" + strong.to_s + " w" + weak.to_s + "\n"
						end
						# strong, weak = classify_sah_grid_neighbors(aa, nbl[1], nbr[1])
						strong, middle, weak, helix_stability = classify_greyscale_sah_grid_neighbors(aa, nbl[1], nbr[1])

						# Declaration: edge_br/l == 0 => no edge
						edge_bl = 0; edge_br = 0
						# edge_br = (strong == 12 || strong == 4) ? STRONG : 0;
						# edge_bl = (strong == 12 || strong == 8) ? STRONG : 0;
						# edge_br = (edge_br == 0 && (weak == 12 || weak == 4)) ? WEAK : edge_br;
						# edge_bl = (edge_bl == 0 && (weak == 12 || weak == 8)) ? WEAK : edge_bl;
						case strong
							when BOTTOM_RIGHT
								edge_br = STRONG
							when BOTTOM_LEFT
								edge_bl = STRONG
						end
						case middle
							when BOTTOM_RIGHT
								edge_br = MEDIUM
							when BOTTOM_LEFT
								edge_bl = MEDIUM
						end
						case weak
							when BOTTOM_RIGHT
								edge_br = WEAK
							when BOTTOM_LEFT
								edge_bl = WEAK
						end

						# Graph: Add parent edge
						graph.add_node(pos, aa, heptad, calc_aa_color(aa, 'sah'))
						# Connect nodes: parent and nbl
						if edge_bl > 0
							# adjacency_matrix[pos - 1][nbl[0].to_i - 1] = edge_bl
							# adjacency_matrix[nbl[0].to_i - 1][pos - 1] = edge_bl
							graph.add_node(nbl[0].to_i, nbl[1], nbl[2], calc_aa_color(nbl[1], 'sah'))
							graph.add_edge(prntLbl, nbl[1]+nbl[0], 'left', edge_bl)
						end
						# Connect nodes: parent and nbr
						if edge_br > 0
							# adjacency_matrix[pos - 1][nbr[0].to_i - 1] = edge_br
							# adjacency_matrix[nbr[0].to_i - 1][pos - 1] = edge_br
							graph.add_node(nbr[0].to_i, nbr[1], nbr[2], calc_aa_color(nbr[1], 'sah'))
							graph.add_edge(prntLbl, nbr[1]+nbr[0], 'right', edge_br)
						end

					end
				end
		} }

		return graph
	end

	# Helper-method: Calculate cumulative path score for passed interaction graph
	# @param [Graph]  graph     Interaction graph
	# return [Float]            Cumulative path score
	def calculate_cumulative_path_score(graph)
		if graph.is_a?(Graph)
				# graph.find_cliques
				return graph.get_paths
		else
				return 0.0
		end
	end

	# Helper-method to convert categorized grid-matrix (catMat) to interaction-grid
	# @param CoiledCoilDomain coiledcoil  CoiledCoilDomain object
	# @param array  catMat                data matrix to convert
	# return array  interMat              converted matrix - interaction matrix to calculate connections between AAs
	def convert_to_interaction_grid(coiledcoil, catMat)
		position = Hash["a" => [3,0], "b" => [1,0], "c" => [6,1], "d" => [4,1], "e" => [2,1], "f" => [0,1], "g" => [5,2]]
		counter = 0;
		interMat = Array.new(20) { Array.new(7) }
		categories = coiledcoil.get_categories

		catMat.each_with_index{ |categoryVec, line_idx|
			category = categories[line_idx]
			categoryVec.each_with_index{ |entry, entry_idx|
				# break if entry_idx >= 7
				if entry != nil
					# GRID representation: grid[row][col]
					if position[category.name][1]+(entry_idx)*2 < 20
						interMat[position[category.name][1]+(entry_idx)*2][position[category.name][0]] = entry
						counter += 1
					end
				end
			}
		}
		return interMat
	end

	# Helper-method to transpose a matrix
	# @param array  gridMat   data for the svg-file
	# return array  newMat    transposed matrix
	def transpose_simple_grid(gridMat)
		newMat = Array.new(7) { Array.new(20) };
		gridMat.each_with_index{ |line, line_idx|
			line.each_with_index{ |node, node_idx|
				newMat[node_idx][line_idx] = node if line_idx < 20
			}
		}
		return newMat
	end

	# Helper-method to calculate simple SAH overview grid-matrix
	# @param CoiledCoilDomain  coiledcoil     data for the svg-file
	# @param Boolean draw_mode                method called for displaying net-view [true: false]
	# return array gridMat                    shape - SVG markup (code)
	def calculate_simple_sah_grid(coiledcoil, draw_mode)
		position = Hash["a" => 0, "b" => 1, "c" => 2, "d" => 3, "e" => 4, "f" => 5, "g" => 6]
		size = 18; round = 0; lastPos = 0; start = true

		# Min size of 18 for the net view, else scalable
		if !draw_mode
			size = (coiledcoil.coiledcoil_data.size/7).ceil+10 if (coiledcoil.coiledcoil_data.size/7).ceil >= 18
		end
		grid = Array.new(size) { Array.new(7) };
		coiledcoil.coiledcoil_data.each_with_index{ |el, i|
			pos = position[el[2]]
			if lastPos >= pos && !start
				round = round + 1
			end
			break if round >= size # Break condition to guarantue no "index-out-of-bounds" error
			grid[round][pos.to_i] = el
			lastPos = pos; start = false
		}

		# Trim (empty array lines) and return
		return trim_array(grid)
	end

	# Helper-method to calculate simple SAH overview grid-matrix
	# @param CoiledCoilDomain  coiledcoil     data for the svg-file
	# @param Boolean  draw_mode                method called for displaying net-view [true: ]
	# return array gridMat                    shape - SVG markup (code)
	def calculate_simple_gapless_sah_grid(coiledcoil, draw_mode)
		position = Hash["a" => 0, "b" => 1, "c" => 2, "d" => 3, "e" => 4, "f" => 5, "g" => 6]
		size = 18; round = 0; start = true;

		# Min size of 18 for the net view, else scalable
		if !draw_mode
			size = (coiledcoil.coiledcoil_data.size/7).ceil+10 if (coiledcoil.coiledcoil_data.size/7).ceil >= 18
		end
		grid = Array.new(size) { Array.new(7) };
		coiledcoil.coiledcoil_data.each_with_index{ |el, i|
			if i%7 == 0 && !start
				round = round + 1
			end
			break if round >= 18 # Break condition to guarantuee no "index-out-of-bounds" error
			grid[round][i%7] = el
			start = false
		}

		# Trim (empty array lines) and return
		return trim_array(grid)
	end

	# Array: Helper-method to Trim empty array lines
	# @param [Array] array		array to trim
	# @return [Array]         (trimmed) array
	def trim_array(array)
		remove_indeces = Array.new();
		array.each_with_index { |line, index|
			if line != nil
				if line.all? {|x| x.nil?} # => true
					remove_indeces.unshift(index)
				end
			end
		}
		remove_indeces.each { |index| array.delete_at(index) }
		return array
	end

	# SVG: Helper-method to create positioning shape
	# @param [Array] pos		top-left corner of main shape rec-angle
	# @param [Integer]   width
	# @param [Integer]   height
	# @return [String]         shape - SVG markup (code)
	def draw_rect_shape(position, width, height)
		data = '<rect x="' + position[0].to_s + '" y="' + position[1].to_s + '" width="' + width.to_s + '" height="' + height.to_s + '" rx="4" ry="4" fill="white" stroke="black" stroke-width="1px"/>'
		return data
	end

	# SVG: Helper-method to create positioning shape
	# @param [CoiledCoilSequence]
	# @param [Array] pos             top-left corner of main shape rec-angle
	# @param [Integer]   width
	# @param [Integer]   height
	# @return [String]               shape - SVG markup (code)
	def draw_sequences(coiledcoils, position, width, height)
		data = ''; oligomerisation = ['&#8693;', '&#8648;', '2', '3', '4']
		global_width = width-6; last_end = 0
		domain_position = position[0]+6
#		data << '<path transform="translate(0, ' + (height+20).to_s + ')" d="M100,10 Q100,30 398,20 Q400,40 402,20 Q700,30 700,10 " width="' + global_width.to_s + '" stroke="black" fill="none" />'
#		data <<	'<text class="tiny" x="' + (width/2+70).to_s + '" y="' + (height+60).to_s + '" id="text" >' + coiledcoils.seq_length.to_s + ' aa</text>'
		# data <<	'<text class="tiny" x="' + (width+120).to_s + '" y="' + (height/2+24).to_s + '" id="text" style="text-anchor:start" >' + coiledcoils.seq_length.to_s + ' aa</text>'
		data <<	'<text class="tiny" x="' + (width+120).to_s + '" y="' + (height/2+24).to_s + '" style="text-anchor:start" >' + coiledcoils.seq_length.to_s + ' aa</text>'
		coiledcoils.domains.each_with_index { |coiledcoil, i|
			# print [coiledcoil.start, coiledcoil.end, coiledcoils.seq_length]
			length = coiledcoil.end.to_i - coiledcoil.start.to_i
			perc = length.abs.to_f/coiledcoils.seq_length
			if ((last_end+1).to_i < coiledcoil.start.to_i)
				# calc intermed. space
				intermediate = coiledcoil.start.to_i - (last_end+1).to_i
				domain_position += (intermediate.to_f/coiledcoils.seq_length*global_width).to_i
			end
			domain_width = (perc*global_width).to_i
			typeArr = coiledcoil.oligomerisation.map { |state| state[0] }
			type = typeArr[0].to_i <= 2 ? "2" : "3";
			# Use user selection for coiledcoildomain-image
			# data << '<a xlink:href="javascript:showdomain(\'' + coiledcoils.name + '\', \'' + coiledcoils.cli_name + '\', $(\'#coiledcoil_type\').attr(\'value\'), ' + i.to_s + ', ' + coiledcoil.start.to_s + ', ' + coiledcoil.end.to_s + ');" xlink:title="Show region: ' + coiledcoil.start.to_s + '-' + coiledcoil.end.to_s + '">'
			# Use prediction for coiledcoildomain-image
			data << '<a xlink:href="javascript:showdomain(\'' + coiledcoils.name + '\', \'' + coiledcoils.cli_name + '\', \'' + type + '\', ' + i.to_s + ', ' + coiledcoil.start.to_s + ', ' + coiledcoil.end.to_s + ', $(\'#coiledcoil_gap_mode\').val());" xlink:title="Show region: ' + coiledcoil.start.to_s + '-' + coiledcoil.end.to_s + '">'
			data << '<rect x="' + domain_position.to_s + '" y="' + (position[1]-5).to_s + '" width="' + (domain_width-3).to_s + '" height="' + (height+10).to_s + '" rx="4" ry="4" stroke-width="1px" id="' + coiledcoils.name + '_' + coiledcoil.start.to_s + '-' + coiledcoil.end.to_s + '" class="ccdomain" style="fill: url(#gradient);"/>'
			data << '</a>'

			# irregular breaks in register
			coiledcoil.breaks.each { |pos|
				xPos = ((pos-coiledcoil.start).to_f/coiledcoils.seq_length*global_width).to_i-4
				data << '<a xlink:title="Break-Position: ' + pos.to_s + '">'
				data << '<line x1="' + (domain_position+xPos).to_s + '" y1="' + (position[1]-10).to_s + '" x2="' + (domain_position+xPos).to_s + '" y2="' + (position[1]+height+18).to_s + '" style="stroke:rgb(255,0,0)" />'
				data << '</a>'
			}

			# Textures for domains
			if domain_width > 66
				data <<	'<text class="tiny" x="' + domain_position.to_s + '" y="' + (position[1]-10).to_s + '" id="text" >' + coiledcoil.start.to_s + '-' + coiledcoil.end.to_s + '</text>'
			end
			if (domain_width > 56 && coiledcoil.oligomerisation.size > 3 || domain_width > 46 && coiledcoil.oligomerisation.size == 3 || domain_width > 36 && coiledcoil.oligomerisation.size == 2 || domain_width > 26 && coiledcoil.oligomerisation.size == 1)
				shift = 0
				coiledcoil.oligomerisation.each_pair { |tool, value|
					data <<	'<circle class="budge" cx="' + (domain_position+12+shift).to_s + '" cy="' + (position[1]+8).to_s + '" r="8" />'
					if [$tools[1],$tools[2]].rindex(tool) != nil
						data << '<a xlink:title="' + tool.capitalize + ' (' + ($oligomerisation_states[value[0]].capitalize) + '): ' + (value[1]/length).round(2).to_s + ' (averaged cumulative Probability)">'
					else
						data << '<a xlink:title="' + tool.capitalize + ' (' + ($oligomerisation_states[value[0]].capitalize) + '): ' + value[1].round(2).to_s + ' (Score)">'
					end
					data <<	'<text class="tiny" x="' + (domain_position+8+shift).to_s + '" y="' + (position[1]+12).to_s + '" id="text" >' + oligomerisation[value[0]] + '</text>'
					data << '</a>'
					if coiledcoil.oligomerisation.size > 1
						shift += 12
					end
				}
			elsif coiledcoil.oligomerisation.size > 0
				title = ''
				coiledcoil.oligomerisation.each_pair { |key, value|
					if [$tools[1],$tools[2]].rindex(key) != nil
						title << key.capitalize + ' (' + ($oligomerisation_states[value[0]].capitalize) + '): ' + (value[1]/length).round(2).to_s + ' (averaged cumulative Probability)' + "\n"
					else
						title << key.capitalize + ' (' + ($oligomerisation_states[value[0]].capitalize) + '): ' + value[1].round(2).to_s + ' (Score)' + "\n"
					end
				}
				data << '<a xlink:title="' + title + '">'
				data <<	'<text class="tiny" x="' + (domain_position).to_s + '" y="' + (position[1]+12).to_s + '" id="text" style="font-style:italic; fill:#CCC; stroke:#060;" >&#149;</text>'
				data << '</a>'
			end
			domain_position += domain_width
			last_end = coiledcoil.end
		}
		return data
	end

	# SVG: Helper-method to create positioning shapes (outer circles)
	# @param [Array] center   center of main shape circle
	# @param [Integer] radius     initial radius of the main shape circle
	# @param [Integer] rotate     rotation angle of the spiral
	# @param [Integer] number     number of further outer shapes
	# @param [String] text    center text
	# @return [String]        circles - SVG markup (code)
	def draw_shapes(center, radius, rotate = 0, number, text)
		data = ''; data_text = ''; d = '';
		spiral_offset = 20; # spiral_offset
		so = 0; sqo = 5; # spiral_quarter_offset
		wheel_shift = (rotate>180) ? 2*radius : 0;
		text_shift = (rotate>180) ? 0 : 11;
		nextpoint = [center[0]+2, center[1]+radius]
		data << '<path class="shape" d="M 0,0 '
		(0..number-1).each { |i|
			nextpoint = calc_nextpoint(nextpoint, center, 20)
			so = i*spiral_offset
			d << 'A ' + (radius+so+sqo).to_s + ',' + (radius+so+sqo).to_s + ' 0 0 1 ' + (2*radius+so+2*sqo).to_s + ',0 '
			d << 'A ' + (radius+so+3*sqo).to_s + ',' + (-(radius+so+3*sqo)).to_s + ' 0 0 1 ' + (-(so+4*sqo)).to_s + ',0 ' if (i<=number-1)
			data_text << '<text transform="translate(-4, ' + ((text_shift+radius+8+40*i)).to_s + ')" class="lightgrey" x="' + center[0].to_s + '" y="' + center[1].to_s + '">' + (i+1).to_s + '</text>' + "\n\t" if(i < 5)
			data_text << '<text transform="" class="shapeText lightgrey" x="' + (center[0]+2).to_s + '" y="' + (center[1]+40).to_s + '">' + text + '</text>' + "\n\t"
		}
		data << d + '"/>'
		data = '<g transform="translate(' + (center[0]-radius+wheel_shift).to_s + ',' + (center[1]).to_s + ') rotate(' + rotate.to_s + ') ">' + data + '</g>' + data_text
		return data
	end

	# SVG: Helper-method to create positioning shapes (outer circles)
	# @param [Array] center	center of main shape circle
	# @param [Integer] radius     initial radius of the main shape circle
	# @param [Integer] number     number of further outer shapes
	# @param [String] text    center text
	# @return [String]         circles - SVG markup (code)
	def draw_XXX_shapes(center, radius, number, text)
		data = ''
		#nextpoint = [-(16+radius-40), 0]
		nextpoint = [center[0]+2, center[1]+radius]
		(0..number-1).each { |i|
			nextpoint = calc_nextpoint(nextpoint, center, 40)
			data << '<circle class="shape" cx="' + center[0].to_s + '" cy="' + center[1].to_s + '" r="' + (radius+40*i).to_s + '" />' + "\n\t"
			# Ring numbers straight horizontally to the left
			#data << '<text transform="translate(' + (-(16+radius+40*i)).to_s + ',0)" class="lightgrey" x="' + center[0].to_s + '" y="' + center[1].to_s + '">' + (i+1).to_s + '</text>' + "\n\t"
			#data << '<text transform="translate(' + nextpoint[0].to_s + ',' + nextpoint[1].to_s + ')" class="lightgrey" x="' + center[0].to_s + '" y="' + center[1].to_s + '">' + (i+1).to_s + '</text>' + "\n\t"
			# Ring numbers straight vertically to the bottom
			data << '<text transform="translate(0, ' + ((6+radius+40*i)).to_s + ')" class="lightgrey" x="' + center[0].to_s + '" y="' + center[1].to_s + '">' + (i+1).to_s + '</text>' + "\n\t"
			data << '<text transform="" class="shapeText lightgrey" x="' + (center[0]+2).to_s + '" y="' + (center[1]+40).to_s + '">' + text + '</text>' + "\n\t"
		}
		return data
	end

	# SVG: Helper-method to create svg-data with drawn circles and labels according to given points
	# @param [Object]			CoiledCoilDomain-object: positions and data (grouped by heptad-position)
	# @param [Array] center   center of main shape
	# @param [Integer] radius     radius of the circles
	# @return [String]         labeled circles - SVG markup
	def draw_circles(coiledcoil, center, radius)
		data = ''
		categories = coiledcoil.get_categories
		categories.each { |category|
			x = category.pos[0]; y = category.pos[1];
			if ['a','d'].rindex(category.name) != nil
				# data << '<g id="label" filter="url(#filter)"><circle class="heptads" style="fill: url(#gradient);" cx="' + x.to_s + '" cy="' + y.to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				data << '<g filter="url(#filter)"><circle class="heptads" style="fill: url(#gradient);" cx="' + x.to_s + '" cy="' + y.to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				label_pos = calc_nextpoint(category.pos, center, 2*(-radius))
				data << '<text transform="translate(-4,4)" class="label small" x="' + label_pos[0].to_s + '" y="' + label_pos[1].to_s + '" id="' + category.name + '">' + category.name + '</text>' + "\n\t"
			elsif ['e','f','g'].rindex(category.name) != nil
				circle_pos = calc_nextpoint(category.pos, center, 2*(-radius+26))
				label_pos = calc_nextpoint(category.pos, center, 2*(-radius+12))
				#data << '<g id="label" filter="url(#filter)"><circle class="outer-heptads" cx="' + x.to_s + '" cy="' + y.to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				# data << '<g id="label" filter="url(#filter)"><circle class="outer-heptads" style="fill: url(#gradient-inactive);" cx="' + circle_pos[0].to_s + '" cy="' + circle_pos[1].to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				data << '<g filter="url(#filter)"><circle class="outer-heptads" style="fill: url(#gradient-inactive);" cx="' + circle_pos[0].to_s + '" cy="' + circle_pos[1].to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				data << '<text transform="translate(-4,4)" class="label small" x="' + label_pos[0].to_s + '" y="' + label_pos[1].to_s + '" id="' + category.name + '">' + category.name + '</text>' + "\n\t"
			else
				# data << '<g id="label" filter="url(#filter)"><circle class="outer-heptads" style="fill: url(#gradient-inactive);" cx="' + x.to_s + '" cy="' + y.to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				data << '<g filter="url(#filter)"><circle class="outer-heptads" style="fill: url(#gradient-inactive);" cx="' + x.to_s + '" cy="' + y.to_s + '" r="' + radius.to_s + '" /></g>' + "\n\t"
				label_pos = calc_nextpoint(category.pos, center, 2*(-radius))
				# data << '<text transform="translate(-4,4)" class="label category" x="' + label_pos[0].to_s + '" y="' + label_pos[1].to_s + '" id="' + category.name + '">' + category.name + '</text>' + "\n\t"
				data << '<text transform="translate(-4,4)" class="label category" x="' + label_pos[0].to_s + '" y="' + label_pos[1].to_s + '">' + category.name + '</text>' + "\n\t"
			end
		}
		return data
	end

	# SVG: Helper-method to create svg-data with lines (connections between circles)
	# @param [Object]			CoiledCoilDomain-object: positions and data (grouped by heptad-position)
	# @param [Integer] radius 	radius of the circles
	# @return [String]       	Arrows - SVG markup
	def draw_connections(coiledcoil, center, radius)
		data = ''; width = 4.5
		categories = coiledcoil.get_categories
		categories.each { |category|
			startpoint = calc_nextpoint(category.pos, center, -radius)
			endindex = (categories.rindex(category)+1) % categories.count
			endpoint = calc_nextpoint(categories[endindex].pos, center, -radius)
			data <<
				'<g id="g' + category.pos[0].to_s + '">
					<path
						 d="M ' + startpoint[0].to_s + ',' + startpoint[1].to_s + ' L ' + endpoint[0].to_s + ',' + endpoint[1].to_s + '"
						 class="arrow dotted_arrow" style="marker-end: url(#arrow); stroke-width: ' + width.to_s + 'px;" />
				</g>' + "\n\t"
			width -= 0.5
		}
		return data
	end

	# SVG: Helper-method to create svg-data with coiled coil assignments
	# @param [Object] 		CoiledCoilDomain-object: positions and data (grouped by heptad-position)
	# @param [Array] center 	center-coords of the main circle (shape)
	# @param [Integer] radius 	radius of the small circles
	# @param [String] mode      identifier of oligomerization type
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [String] 		Arrows - SVG markup
	def draw_categories(coiledcoil, center, radius, mode, gapMode)
		data = ''; data_end = ''; data_path_end = ''
		## Make deep copy
		coiledcoil_copy = Marshal.load(Marshal.dump(coiledcoil))
		categories = coiledcoil_copy.get_categories
		## All categories should start in a row: a is the lowest position
		amino_start = categories[0].elements[0][0].to_i
		categories.each { |category|
			if category.elements[0][0].to_i < amino_start
					category.elements.delete_at(0)
			end
		}
		## Insert gaps, when shift/break in pattern
		# Prepare interaction data
		# shift = amino_start-coiledcoil.start # Add shift-difference to endCoord
		# seqPart = fullsequence.subseq(abs_start, coiledcoil.end+shift, true)
		seqPart = coiledcoil_copy
		if gapMode == true
				gridMat = calculate_simple_sah_grid(seqPart, false)
				# puts gridMat.inspect
				gridMat.delete_at(0) if gridMat[0][0] == nil
				# puts gridMat.inspect
				categories.each_with_index { |category, cat_pos|
						# puts "Category: Before"
						# puts category.inspect
						col = 0
						category.elements.each { |element|
								# puts "- " + gridMat[col].inspect
								# puts "- Grid: " + gridMat[col][cat_pos].inspect
								# puts "- Elem: " + element.inspect
								# Break, when all entries in grid column are empty, indicator that no further data follows
								break if gridMat[col] && gridMat[col].all?(&:nil?) # => Go to next category if all col-entries = nil
								if gridMat[col][cat_pos] != nil && element != nil && element[0].to_i == gridMat[col][cat_pos][0].to_i
								else
										if element != nil
												category.elements.insert(col, nil)
												col += 1 if gridMat[col][cat_pos] != nil
												# puts "Category: Insert"
												# puts category.elements.inspect
										end
								end
								col += 1
						}
						# puts "Category: After"
						# puts category.inspect
				}
				# puts ""
		end

		categories.each { |category|
		# 	puts category.inspect
			delta = [category.pos[0]-center[0], category.pos[1]-center[1]]
			delta_norm = [delta[0]/delta[0].abs, delta[1]/delta[0].abs]
			slope = delta[1]/delta[0];
			alternate = 0; box = 0; count = 0; distance = 0; increment = 40; skip = 0 	# distance dependent on font-size

			collide_shift = [0,0];
			if mode == 'double'
				if category.name == 'a'
					category.pos = categories[4].pos
					collide_shift = delta[1] > 0 ? [radius+66,radius-84] : [radius-100,radius+54] # ul : or
					increment = 48
				elsif category.name == 'd'
					category.pos = categories[6].pos
					collide_shift = delta[1] > 0 ? [radius-112,radius-86] : [radius+76,radius+54] # ur : ol
					increment = 50
				end
			end
			if ['e','f','g'].rindex(category.name) != nil
					#shift_factor = 10
					#collide_shift = [delta_norm[0]*shift_factor,delta_norm[1]*shift_factor]
					shift_one = calc_nextpoint(category.pos, center, radius+increment/2)
					shift_two = calc_nextpoint(category.pos, center, radius)
					collide_shift = [shift_one[0]-shift_two[0], shift_one[1]-shift_two[1]]
			end

			# Render fields on circular shapes
			category.elements.each_with_index { |entry, i|
				nextpoint = calc_nextpoint(category.pos, center, distance)
				distance += increment
				#if(slope.abs < 0.25)
				#	((distance/increment)%2 == 0) ? alternate = 10 : alternate = -10
				#end
				if i > 0 && i < 5
					data << '<g><rect class="transparent" transform="translate(-12,' + (-16 + alternate).to_s + ')" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + collide_shift[1]).to_s + '" width="24" height="28" rx="4" />'
			  	data_path_end << '</g>'
				end
			}

			# Render aminoacids: normal / reverse
			distance = 0
			elements = (coiledcoil.reverse == 1) ? category.elements.reverse : category.elements
			elements.each_with_index { |entry, i|
				count += 1; shift = 0
				nextpoint = calc_nextpoint(category.pos, center, distance)
				distance += increment
				#if(slope.abs < 0.25)
				#	((distance/increment)%2 == 0) ? alternate = 10 : alternate = -10
				#end
				next if entry == nil
				id = entry[1].to_s + entry[0].to_s

				# Draw the 5 inner positions
				if (count+skip) < 6
					if i < 1
						alternate = 0
					else
						# if entry[0].to_i-7 > elements[i-1][0].to_i  		# Skip node
						# 	skip += 1 if gapMode == true                  # Turn skipping off, when in gapless mode
						# # elsif entry[0].to_i-7 < elements[i-1][0].to_i 	# Repeat node
						# # 	skip -= 1; shift = delta[1] > 0 ? 16 : -16;
						# end
					end
					if (count+skip) < 6
						nextpoint = calc_nextpoint(category.pos, center, distance+(skip-1)*increment)
						alternate *= (-1) ** skip
						# data << '<text transform="translate(0,' + (8 + alternate.to_i).to_s + ')" class="' + (i==0 ? 'bold' : '') + ' badgelabel ' + calc_aa_color(entry[1].to_s,"cc").to_s + '" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + shift + collide_shift[1]).to_s + '" id="' + id + '" title="' + id + '">' + entry[1].to_s + '</text>'
						# data << '<text transform="translate(0,' + (8 + alternate.to_i).to_s + ')" class="' + (i==0 ? 'bold' : '') + ' badgelabel ' + calc_aa_color(entry[1].to_s,"cc").to_s + '" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + shift + collide_shift[1]).to_s + '" id="' + id + '">' + entry[1].to_s + '</text>'
						data << '<text transform="translate(0,' + (8 + alternate.to_i).to_s + ')" class="' + (i==0 ? 'bold' : '') + ' badgelabel ' + calc_aa_color(entry[1].to_s,"cc").to_s + '" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + shift + collide_shift[1]).to_s + '">' + entry[1].to_s + '</text>'
					end
				# 	(distance > increment) ? data << '</g>' + "\n\t" : data << "\n\t"

				# Draw the other categories on the shapes
				elsif count == 6
					data << data_path_end
					box = nextpoint
					data << '<a xlink:href="javascript:splitter(\'' + id + center[0].to_s + '\');" xlink:title="Show more aminoacids ..." style="text-decoration:none;">'
					data << '<rect id="r' + id + center[0].to_s + '" transform="translate(-12,' + (-16 + alternate.to_i).to_s + ')" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + collide_shift[1]).to_s + '" width="48" height="28" rx="4" style="stroke-dasharray:6, 4;" />'
					data << '<text id="t' + id + center[0].to_s + '" transform="translate(-8,' + (4 + alternate.to_i).to_s + ')" class="collect" x="' + (nextpoint[0] + collide_shift[0]).to_s + '" y="' + (nextpoint[1] + collide_shift[1]).to_s + '" style="fill-opacity:0; display:none">'
					data << '<tspan class="' + calc_aa_color(entry[1].to_s,"cc").to_s + '">' + id + '</tspan>'
					data_end << '</text></a>' + "\n\t"
				elsif count < 25
					if((count-6) % 4 == 0)
							data << ', <tspan x="' + (box[0] + collide_shift[0]).to_s + '" dy="20" class="' + calc_aa_color(entry[1].to_s,"cc").to_s + '">'
					else
							data << ', <tspan class="' + calc_aa_color(entry[1].to_s,"cc").to_s + '">'
					end
					data << id + '</tspan>'
				end
			}
			# Close text-tag only when outlying collect-textbox is added
			data << data_end
			data_end = ''; data_path_end = ''
		}
		return data
	end

	# Load SVG-template and fill it with data
	# @param [Object] 			CoiledCoilDomain-object: positions and data (grouped by heptad-position)
	# @param [Bool] noDomain		flag, are there any domains
	# @return [String]				output as SVG file
	def make_domain_model_svg(coiledcoils, noDomains)
		pos = [100, 24]; width = 600; height = 16;
		svg_data = draw_rect_shape(pos, width, height)
		if noDomains == true
			svg_data << draw_sequences(coiledcoils, pos, width, height)
		else
			svg_data << '<text transform="translate(0,100)" class="badgelabel " x="400" y="32">The analysed sequence has no coiled coil domains.</text>'
		end

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/domain_model.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		# svg_tmpl.gsub!('width="800" height="100" viewBox="0 0 800 100"', 'width="800" height="100" viewBox="0 0 800 64"')
		# svg_tmpl.gsub!('width="800" height="100" viewBox="0 0 800 100"', 'width="800" height="64" viewBox="0 0 800 64"')
		svg_tmpl.gsub!('width="800" height="100" viewBox="0 0 800 100"', 'width="800" height="48" viewBox="0 0 800 48"')

		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilDomain]  coiledcoil         data for the svg-file
	# @param [COiledCoilSequence] 	fullsequence		flag, start-position should be chosen absolute or relative in coiledcoil pattern
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]					                        output as SVG file
	def make_sah_domain(coiledcoil, fullsequence, gapMode)
		categories = coiledcoil.get_categories
		startCoord = categories[0].elements[0][0].to_i # Begin always on heptad-postion 'a'
		shift = startCoord - coiledcoil.start - 1      # Add shift-difference to endCoord
		shift = shift + 4                              # Add last line / half heptad, to see bottom interactions
		endCoord = coiledcoil.end+shift

		# Create SVG subparts content
		overview = draw_sah_overview(coiledcoil, fullsequence, startCoord, gapMode)
		svg_data, score, strong, middle, weak, end_coord = draw_grid_sah_aminoacids(coiledcoil, fullsequence, startCoord, gapMode)

		# Fill SVG-template
		svg_tmpl = IO.read("#{Rails.root}/lib/assets/sah_domain.svg")
		svg_tmpl.gsub!("###OVERVIEW###", overview)
		svg_tmpl.gsub!("###CONTENT###", svg_data)
		svg_tmpl.gsub!("###SCORE###", "%.4f" % score)
		svg_tmpl.gsub!('###STRONG###', strong.to_s)
		svg_tmpl.gsub!('###MIDDLE###', middle.to_s)
		svg_tmpl.gsub!('###WEAK###', weak.to_s)
		svg_tmpl.gsub!('###POSITION###', startCoord.to_s+"-"+end_coord.to_s)
		# viewbox: xCoord yCoord width height
		svg_tmpl.gsub!('height="467px" width="700px" viewBox="0 0 360 380">', 'height="467px" width="700px" viewBox="0 0 360 380">')

		return svg_tmpl, startCoord, end_coord
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilDomain] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_cc_svg(coiledcoils, gapMode)
		category_radius = 16; shape_radius = 120; center = [400, 400]
		rings = (size = (coiledcoils.end-coiledcoils.start)/7) < 5 ? size : 10
		points = points_on_spiral(7, shape_radius, 164, center[0], center[1])
		coiledcoils.assign_positions(points)
		svg_data =  draw_shapes(center, shape_radius, rings, '1')
		svg_data << draw_connections(coiledcoils, center, category_radius)
		svg_data << draw_circles(coiledcoils, center, category_radius)
		svg_data << draw_categories(coiledcoils, center, category_radius, '', gapMode)

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="800" height="500" viewBox="0 0 800 800"')
		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilDomain] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_dbl_cc_svg(coiledcoils, gapMode)
		category_radius = 16; shape_radius = 120
		rings = (size = (coiledcoils.end-coiledcoils.start)/7) < 5 ? size : 10
		center_pri = [360, 400]; center_sec = [820, 400]
		# draw second circle [right]
		# 		points_sec = points_on_circle(7, shape_radius, 164, center_sec[0], center_sec[1])
		points_sec = points_on_spiral(7, shape_radius, 164, center_sec[0], center_sec[1])
		coiledcoils.assign_positions(points_sec)
		svg_data =  draw_shapes(center_sec, shape_radius, 0, rings, '2')
		svg_data << draw_shapes(center_pri, shape_radius, 181, rings, '1')
		svg_data << draw_connections(coiledcoils, center_sec, category_radius)
		svg_data << draw_circles(coiledcoils, center_sec, category_radius)
		svg_data << draw_categories(coiledcoils, center_sec, category_radius, 'double', gapMode)
		# draw first circle [left]
# 		points_pri = points_on_circle(7, shape_radius, 344, center_pri[0], center_pri[1])
		points_pri = points_on_spiral(7, shape_radius, 344, center_pri[0], center_pri[1])
		coiledcoils.assign_positions(points_pri)
		svg_data << draw_connections(coiledcoils, center_pri, category_radius)
		svg_data << draw_circles(coiledcoils, center_pri, category_radius)
		svg_data << draw_categories(coiledcoils, center_pri, category_radius, 'double', gapMode)

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="850" height="500" viewBox="0 0 1200 800"')
		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilSequence] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_hetero_dbl_cc_svg(coiledcoils, gapMode)
		category_radius = 16; shape_radius = 120
		rings = (size = (coiledcoils.domains[0].end-coiledcoils.domains[0].start)/7) < 5 ? size : 10
		center_pri = [360, 400]; center_sec = [820, 400]
		# draw second circle [right]
		# 		points_sec = points_on_circle(7, shape_radius, 164, center_sec[0], center_sec[1])
		points_sec = points_on_spiral(7, shape_radius, 164, center_sec[0], center_sec[1])
		coiledcoils.domains[1].assign_positions(points_sec)
		svg_data =  draw_shapes(center_sec, shape_radius, 0, rings, '2')
		svg_data << draw_shapes(center_pri, shape_radius, 181, rings, '1')
		svg_data << draw_connections(coiledcoils.domains[1], center_sec, category_radius)
		svg_data << draw_circles(coiledcoils.domains[1], center_sec, category_radius)
		svg_data << draw_categories(coiledcoils.domains[1], center_sec, category_radius, 'double', gapMode)
		# draw first circle [left]
		# 		points_pri = points_on_circle(7, shape_radius, 344, center_pri[0], center_pri[1])
		points_pri = points_on_spiral(7, shape_radius, 344, center_pri[0], center_pri[1])
		coiledcoils.domains[0].assign_positions(points_pri)
		svg_data << draw_connections(coiledcoils.domains[0], center_pri, category_radius)
		svg_data << draw_circles(coiledcoils.domains[0], center_pri, category_radius)
		svg_data << draw_categories(coiledcoils.domains[0], center_pri, category_radius, 'double', gapMode)

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="850" height="500" viewBox="0 0 1200 800"')
		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# Requires a hendecad-prediction (11 categories)
	# @param [CoiledCoilDomain] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_dbl_11cc_svg(coiledcoils, gapMode)
		category_radius = 24; shape_radius = 150
		rings = (size = (coiledcoils.end-coiledcoils.start)/7) < 5 ? size : 5
		# draw first circle
		center_pri = [350, 400]
		points_pri = points_on_circle(11, shape_radius, 344, center_pri[0], center_pri[1])
		coiledcoils.assign_positions(points_pri)
		svg_data =  draw_shapes(center_pri, shape_radius, rings, '1')
		svg_data << draw_connections(coiledcoils, center_pri, category_radius)
		svg_data << draw_circles(coiledcoils, center_pri, category_radius)
		svg_data << draw_categories(coiledcoils, center_pri, category_radius, '', gapMode)
		# draw second circle
		# center_sec = [850, 400]
		# points_sec = points_on_circle(11, shape_radius, 164, center_sec[0], center_sec[1])
		# coiledcoils.assign_positions(points_sec)
		# svg_data << draw_shapes(center_sec, shape_radius, rings, '2')
		# svg_data << draw_connections(coiledcoils, center_sec, category_radius)
		# svg_data << draw_circles(coiledcoils, center_sec, category_radius)
		# svg_data << draw_categories(coiledcoils, center_sec, category_radius, 'double')

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="850" height="500" viewBox="0 0 1200 800"')
		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilDomain] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_trpl_cc_svg(coiledcoils, gapMode)
		category_radius = 16; shape_radius = 120
		rings = (size = (coiledcoils.end-coiledcoils.start)/7) < 5 ? size : 10
		# draw first circle
		center_pri = [600, 140]
		points_pri = points_on_circle(7, shape_radius, 75, center_pri[0], center_pri[1])
		coiledcoils.assign_positions(points_pri)
		svg_data =  draw_shapes(center_pri, shape_radius, rings, '1')
		svg_data << draw_connections(coiledcoils, center_pri, category_radius)
		svg_data << draw_circles(coiledcoils, center_pri, category_radius)
		svg_data << draw_categories(coiledcoils, center_pri, category_radius, 'triple', gapMode)
		# draw second circle
		center_sec = [380, 600]
		points_sec = points_on_circle(7, shape_radius, 300, center_sec[0], center_sec[1])
		coiledcoils.assign_positions(points_sec)
		svg_data << draw_shapes(center_sec, shape_radius, rings, '2')
		svg_data << draw_connections(coiledcoils, center_sec, category_radius)
		svg_data << draw_circles(coiledcoils, center_sec, category_radius)
		svg_data << draw_categories(coiledcoils, center_sec, category_radius, 'triple', gapMode)
		# draw third circle
		center_trd = [830, 600]
		points_trd = points_on_circle(7, shape_radius, 210, center_trd[0], center_trd[1])
		coiledcoils.assign_positions(points_trd)
		svg_data << draw_shapes(center_trd, shape_radius, rings, '3')
		svg_data << draw_connections(coiledcoils, center_trd, category_radius)
		svg_data << draw_circles(coiledcoils, center_trd, category_radius)
		svg_data << draw_categories(coiledcoils, center_trd, category_radius, 'triple', gapMode)

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="800" height="800" viewBox="60 -250 1000 1320"')

# => write svg-data to file
#		File.open("#{$tmp_path_image}/coiled_coil_run_#{timestamp}.svg", "w") do |fh|
#			fh.puts(svg_tmpl.to_s)
#		end
		return svg_tmpl
	end

	# Load SVG-template and fill it with data
	# @param [CoiledCoilDomain] coiledcoils		data of one coiled coil domain
	# @param [Boolean] gapMode 	indicator for [gapless=false/break_based=true] mode
	# @return [Void]								output as SVG file
	def make_hetero_trpl_cc_svg(coiledcoils, gapMode)
		category_radius = 16; shape_radius = 120
		rings = (size = (coiledcoils.domains[0].end-coiledcoils.domains[0].start)/7) < 5 ? size : 10
		# draw first circle
		center_pri = [600, 140]
		points_pri = points_on_spiral(7, shape_radius, 75, center_pri[0], center_pri[1])
		coiledcoils.domains[0].assign_positions(points_pri)
		svg_data =  draw_shapes(center_pri, shape_radius, rings, '1')
		svg_data << draw_connections(coiledcoils.domains[0], center_pri, category_radius)
		svg_data << draw_circles(coiledcoils.domains[0], center_pri, category_radius)
		svg_data << draw_categories(coiledcoils.domains[0], center_pri, category_radius, '', gapMode)
		# draw second circle
		center_sec = [380, 600]
		points_sec = points_on_spiral(7, shape_radius, 300, center_sec[0], center_sec[1])
		coiledcoils.domains[1].assign_positions(points_sec)
		svg_data << draw_shapes(center_sec, shape_radius, rings, '2')
		svg_data << draw_connections(coiledcoils.domains[1], center_sec, category_radius)
		svg_data << draw_circles(coiledcoils.domains[1], center_sec, category_radius)
		svg_data << draw_categories(coiledcoils.domains[1], center_sec, category_radius, '', gapMode)
		# draw third circle
		center_trd = [830, 600]
		points_trd = points_on_spiral(7, shape_radius, 210, center_trd[0], center_trd[1])
		coiledcoils.domains[2].assign_positions(points_trd)
		svg_data << draw_shapes(center_trd, shape_radius, rings, '3')
		svg_data << draw_connections(coiledcoils.domains[2], center_trd, category_radius)
		svg_data << draw_circles(coiledcoils.domains[2], center_trd, category_radius)
		svg_data << draw_categories(coiledcoils.domains[2], center_trd, category_radius, '', gapMode)

		svg_tmpl = IO.read("#{Rails.root}/lib/assets/coiled_coil.svg")
		svg_tmpl.gsub!("###CONTENT###",svg_data)
		svg_tmpl.gsub!('width="800" height="800" viewBox="0 0 800 800"', 'width="800" height="800" viewBox="60 -250 1000 1320"')

# => write svg-data to file
#		File.open("#{$tmp_path_image}/coiled_coil_run_#{timestamp}.svg", "w") do |fh|
#			fh.puts(svg_tmpl.to_s)
#		end
		return svg_tmpl
	end

	# Set the headers to SVG/XML
	# helper_method :set_content_type
	def set_content_type
		headers["Content-Type"] = "image/svg+xml"
	end

	# Remove the xml definition from an SVG file so it can be embedded in XHTML.
	# helper_method :strip_svg
	def strip_svg(svg)
		# return svg.gsub(/<.xml.*?>/, '').gsub(/<defs>.*<\/defs>/m, '').gsub(/<.DOCTYPE.*>/, '').gsub(/style='.*'/, '')
		return svg.gsub(/<defs.*>.*<\/defs>/m, '')
	end


	# def to_s
	# 	text  = ""
	# 	if !@skip_headers then
	# 	  text << %|<?xml version="1.0" standalone="no"?>\n|
	# 	  text << %|<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n|
	# 	end
	# 	text << %|<svg xmlns="http://www.w3.org/2000/svg" width="#{@width}" height="#{@height}" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:xhtml="http://www.w3.org/1999/xhtml"|
	# 	text << %| viewBox="#{@view_box}"| if @view_box
	# 	text << %|>\n|
	#
	# 	@scripts.each { |script|
	# 	  text << script.to_s
	# 	}
	#
	# 	unless @styles.empty?
	# 	  text << %|<defs>\n|
	# 	  text << %|<style type="text/css"><![CDATA[\n|
	# 	  text << @styles.collect { |define| define.to_s + "\n" }.join
	# 	  text << %|]]></style>\n|
	# 	  text << %|</defs>\n|
	# 	end
	#
	# 	text << %|<title>#{@title}</title>\n| if @title
	# 	text << %|<desc>#{@desc}</desc>\n|    if @desc
	# 	text << @elements.collect { |element| element.to_s + "\n" }.join
	# 	text << %|</svg>\n|
	# 	return text
	# end

end

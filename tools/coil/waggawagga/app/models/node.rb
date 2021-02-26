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

class Node

  # Class variable to count instances
  @@count = 0
  # Define public instance variables with getter and setter
  attr_accessor :number, :pos, :aa, :heptad, :type, :parent, :adj_list

  # Initializer; called on creation
  # @param [Integer] pos    sequence position of the current node
  # @param [String] aa      aminoacid of the current node
  # @param [String] heptad  heptad of the current node
  # @param [String] type    type of the current node (class color)
  # @param [String] parent  parent node of the current node (label)
  # @param [Array] adjList  adjacency list of the node; list of directly connected child nodes
  # @return [void]
  def initialize(pos=0, aa='', heptad='', type='', parent='', adjList=[])
    @@count += 1
    @number = @@count
    @pos = pos
    @aa = aa
    @heptad = heptad
    @type = type
    @parent = parent
    @adj_list = adjList
  end

  # @param [Integer] node-index
  # @return [Boolean]
  def has_child(index)
		  return @adj_list.select{ |edge| edge.dst_node_id == index } != nil
  end

  def Node.count
    return @@count
  end

end
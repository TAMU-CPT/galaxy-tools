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

class Edge

  # Define starting and end-point
  attr_accessor :label, :src_node_id, :dst_node_id, :direction, :weight

  # Initializer; called on creation
  # @param [Integer] src        start-point of the edge
  # @param [Integer] dst        end-point of the edge
  # @param [String] direction   direction of the edge
  # @param [Float] cost         weight of the edge
  # @return [void]
  def initialize( label='', src_node_id=nil, dst_node_id=nil, direction='', weight=nil )
    @label = label
    @src_node_id = src_node_id
    @dst_node_id = dst_node_id
    @direction = direction
    @weight = weight
  end

  def edge( dst_node_id, weight)
    @label = label
    @src_node_id = src_node_id
    @dst_node_id = dst_node_id
    @direction = direction
    @weight = weight
  end

end
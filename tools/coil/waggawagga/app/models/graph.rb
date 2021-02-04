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

require_relative 'node'
require_relative 'edge'

require 'json'
require 'pp'

class Graph

  # Initializer; called on creation
  def initialize
    # Private instance variables
    @labels = Hash.new
    @labels.default = nil
    @nodes = Array.new
		@visited = Array.new

    @edge_labels = Hash.new
    @edge_labels.default = nil
    @edges = Array.new
		@visited_edges = Array.new

    @pattern = Array.new
		@cum_path_score = 0
  end

  # @param [Integer] pos    sequence position of the current node
  # @param [String] aa      aminoacid of the current node
  # @param [String] heptad  heptad of the current node
  # @param [String] type    type of the current node (class color)
  # @return [void]
  def add_node(pos, aa, heptad, type)
    label = "#{aa}#{pos}"
    if @labels.has_key?(label)
        # raise NodeAlreadyDefinedException.new(true), "Node already defined"
        # puts "Exception: Node already defined"
				return
    end
    @nodes.push(Node.new(pos, aa, heptad, type))
    idx = @nodes.size - 1  # idx = Node.count
    @labels.store(label, idx)
  end

  # @param [String] label
  # @return [Node]
  def get_node(label)
    idx = @labels[label]
    if idx == nil
      # raise NoSuchElementException.new(true), "No such element"
      puts "Exception: No such element #{label}"
      puts @labels
	  	return nil
    end
    @nodes[idx]
  end

  # @param [String] label
  # @return [Int]
  def get_node_id(label)
    idx = @labels[label]
    if idx == nil
      # raise NoSuchElementException.new(true), "No such element"
      puts "Exception: No such element #{label}"
      puts @labels
	  	return nil
    end
    idx.to_i
  end

  # Helper-method to add an edge to current graph
  # @param [String] src   label of the source node
  # @param [String] dst   label of the destination node
  # @param [String] dir   direction of the edge
  # @param [Float] cost   cost/weight of the edge
  # @return [boolean]
  def add_edge(src, dst, dir, cost)
    if get_node_id(src) != nil && get_node_id(dst) != nil
        src_node = get_node(src)
        dst_node = get_node(dst)
        edge_label = '%s%d%s%d' % [src_node.aa, src_node.pos, dst_node.aa, dst_node.pos]
        edge = Edge.new(edge_label, get_node_id(src), get_node_id(dst), dir, cost)
        # Register new edge
        src_node.adj_list.push( edge )
        @edges.push( edge )
        idx = @edges.size - 1
        @edge_labels.store(edge_label, idx)
				# Add Parent node
				dst_node.parent = src
    else
				puts "Exception: One of the nodes src (#{src}) / dst (#{dst}) does not exist"
		end
  end

  # @param [String] label
  # @return [Edge]
  def get_edge(label)
	  idx = @edge_labels[label]
	  if idx == nil
		  # raise NoSuchElementException.new(true), "No such element"
		  puts "Exception: No such element #{label}"
		  puts @edge_labels
		  return nil
	  end
	  @edges[idx]
  end

  # @param [String] label
  # @return [Int]
  def get_edge_id(label)
	  idx = @edge_labels[label]
	  if idx == nil
		  # raise NoSuchElementException.new(true), "No such element"
		  puts "Exception: No such element #{label}"
		  puts @edge_labels
		  return nil
	  end
	  idx.to_i
  end

  # @param [String] node  label of the node
  # @return [Array]
  def get_edges(node)
    if get_node_id(node) != nil
	    # @nodes[get_node_id(node)].adj_list
    else
	    []
    end
  end

  # @param [String] node  label of the node
  # @return [Array]
  def get_parent(node)
    if get_node_id(node) != nil
	    @nodes[get_node_id(node)].parent
    else
	    []
    end
  end

  # @param [String] label  label of the node
  # @return [Boolean]
  def has_edges(label)
    if get_node_id(label) != nil
	    !@nodes[get_node_id(label)].adj_list.empty?
    else
	    false
		end
  end

  # @param [String] subpath   Path directions of connected edges
  # @return [Object]
  def calc_path_pattern(subpath)
		  path_score = 0
		  # puts "Transferred subpath: " + subpath
		  case subpath
				  when 'leftleft'
						  path_score = 2
				  when 'rightright'
						  path_score = 1.5
				  when 'leftright'
						  path_score = 0.2
				  when 'rightleft'
						  path_score = 0.2
				  else
						  path_score = 0
		  end
		  path_score
  end

  # @param [Array] edge_pair   Set of two edges
  # @return [boolean]
  def neighbored_dist(edge_pair)
		if edge_pair.size == 2
			  diff = (edge_pair.first.dst_node_id-edge_pair.last.dst_node_id).abs
				if [3,4].rindex(diff) != nil # diff is in set
						return true
				else
						return false
				end
		else
			false
		end
  end

  # @param [Array] edge_pair   Set of two edges
  # @return [boolean]
  def neighbored_edges(edge_pair)
		if edge_pair.size == 2
				node01 = @nodes[edge_pair.first.src_node_id]
				node02 = @nodes[edge_pair.first.dst_node_id]
				node03 = @nodes[edge_pair.last.src_node_id]
				node04 = @nodes[edge_pair.last.dst_node_id]
				# if node01.has_child(edge_pair.last.dst_node_id) || node02.has_child(edge_pair.first.dst_node_id)
				if node01 == node04 || node02 == node03
						# puts @labels.key(edge_pair.first.dst_node_id) + " and " + @labels.key(edge_pair.last.dst_node_id) + " are neighbored."
						return true
				else
						# puts @labels.key(edge_pair.first.dst_node_id) + " and " + @labels.key(edge_pair.last.dst_node_id) + " are not neighbored."
						# puts node01.adjList.inspect
						# puts node02.adjList.inspect
						return false
				end
		else
			false
		end
  end

  # @param [Array] edge_pair   Set of two edges
  # @return [boolean]
  def alternate_check(edge_pair)
		negative_aa = ['E', 'D']; positive_aa = ['K', 'R', 'H']
		network_nodes = []; state = ''
	  if edge_pair.size == 2
				network_nodes.push(@nodes[edge_pair.first.src_node_id])
				network_nodes.push(@nodes[edge_pair.last.src_node_id])
				network_nodes.push(@nodes[edge_pair.last.dst_node_id])
				network_nodes.each { |node|
					act_state = negative_aa.rindex(node.aa).nil? ? nil : 'negative'
					act_state = positive_aa.rindex(node.aa).nil? ? nil : 'positive'
					if act_state == state
						return false
					else
						state = act_state
					end
				}
			return true
		else
			return false
		end
  end

	# # @param [String] parent_label   unique label of the node
	# # @return [Object]
	# def check_node(parent_label)
	# 	  # puts "@visited:" + @visited.inspect
	# 	  if !@visited.include?(parent_label)
	# 			  # Variante 1 (findet keine ausreichend gute Lösung)
	# 			  # 1. Starte im Wurzelknoten und durchlaufe alle erreichbaren Knoten einzeln
	# 			  #   2. Verfolge teuersten Pfad bis zum letztmöglichen Knoten (Blatt)
	# 			  #   2. Steige auf und werte 2-Tupel-Pfad-Abschnitte aus
	# 			  #   2. Bei Gabelungen gehe erneut in die Tiefe bis zum nächsten Blatt
	# 			  #
  #    	    # Remember nodes visited
	# 	      @visited.push(parent_label)
	# 			  if has_edges(parent_label)
	# 					  parent_node = @nodes[@labels[parent_label]]
	# 					  puts "\n" + parent_label + " **"
	#
	# 					  # Detect higher score-pattern => change order of both children
	# 					  if !@pattern.empty?
	# 						  puts "adjList: " + JSON.dump(parent_node.adj_list)
	# 						  puts "subpathScore: " + calc_path_pattern(@pattern.last.direction+parent_node.adj_list.first.direction).to_s
	# 					  end
	# 						if !@pattern.empty? && calc_path_pattern(@pattern.last.direction+parent_node.adj_list.first.direction) < 1
	# 						    parent_node.adj_list.reverse!
	# 						    puts "adjList: " + JSON.dump(parent_node.adj_list)
	# 						    puts "subpathScore: " + calc_path_pattern(@pattern.last.direction+parent_node.adj_list.first.direction).to_s
	# 						end
	#
	# 					  # Visit children
	# 					  parent_node.adj_list.each { |child_edge|
	# 							  child_label = @labels.key(child_edge.dst_node_id)
	# 							  @pattern.push(child_edge)
	# 							  puts child_edge.direction.to_s + " neighbor: " + @labels.key(child_edge.dst_node_id)
	# 							  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_before_call: " + @pattern.inspect
	# 							  check_node(child_label)
	#
	# 							  # Check path pattern
	# 							  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_after_call: " + @pattern.inspect
	#
	# 							  # Last both pattern edges must be not visited, neighbored and alternating charged
	# 							  last_edges = @pattern.last(2)
	# 							  if last_edges.size == 2 and !@visited_edges.include?(last_edges[-1].label) and !@visited_edges.include?(last_edges[-2].label) and neighbored_edges(last_edges) and alternate_check(last_edges)
	# 								    edges = @pattern.pop(2)
	# 								    labels = edges.map {|edge| edge.label }
	# 								    subpath = edges.map {|edge| edge.direction}
	# 							      @cum_path_score += calc_path_pattern(subpath.join)
	#
	# 								    puts "subpath: " + JSON.dump(subpath) + labels.inspect
	# 										puts "cumPathScore = " + @cum_path_score.to_s
	# 										# Add edges to visited
	# 								    edges.each { |edge| @visited_edges.push(edge.label) }
	# 										# puts @visited_edges.inspect
	# 							  # Drop it
	# 							  else
	# 									  @pattern.pop if !@pattern.empty? && @pattern.last.dst_node_id > child_edge.dst_node_id
	# 							  end
	# 							  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_after_pop: " + @pattern.inspect
	# 							  puts "UP to " + parent_label
	# 					  }
	# 			  else # Current node has no more edges
	# 					  # Check path pattern
	# 					  # puts "Last element: " + label
	# 					  # puts "last_pattern: " + @pattern.inspect
	# 					  # if ( neighbored(@pattern.last(2)) )
	# 							#   subpath = @pattern.pop(2).map! {|edge| edge.direction}
	# 							#   puts "subpath: " + JSON.dump(subpath)
	# 							#   @cumPathScore += calcPathPattern(subpath.join())
	# 							#   puts "cumPathScore = " + @cumPathScore.to_s
	# 					  # else
	# 							#   @pattern.pop
	# 					  # end
	# 					  # puts "last_pattern_after_pop: " + @pattern.inspect
	# 						return nil
	# 				end
	# 	  else
	# 				# puts parent_label + " already visited"
	# 	  end
	# end

  # @param [String] parent_label   unique label of the node
  # @return [Object]
  def check_node(parent_label, direction)
	  # puts "@visited:" + @visited.inspect
	  if !@visited.include?(parent_label)
		  # Variante 2
		  # Annahme: Gerade Netzwerke (ll & rr) sind immer zu bevorzugen
		  # 2 getrennte Durchläufe: (kein Konfliktpotential)
		  #   1. Suche nur nach ll-Netzwerken (Merke alle Kanten)
		  #   2. Suche nur nach rr-Netzwerken (Merke alle Kanten)
		  #   3. Bestimme alle Kanten, subtrahiere gesehene (1./2.) und bestimme durch Sortierung
		  #      benachbarte lr oder rl-Netzwerke

		  # Remember nodes visited
		  @visited.push(parent_label)
		  if has_edges(parent_label)
			  parent_node = @nodes[@labels[parent_label]]
			  # puts "\n" + parent_label + " **"

			  # Visit children
			  parent_node.adj_list.each { |child_edge|
				  # Skip all edges with other oriented direction
				  next if child_edge.direction != direction

				  child_label = @labels.key(child_edge.dst_node_id)
				  @pattern.push(child_edge)
				  # puts child_edge.direction.to_s + " neighbor: " + @labels.key(child_edge.dst_node_id)
				  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_before_call: " + @pattern.inspect
				  check_node(child_label, direction)

				  # Check path pattern
				  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_after_call: " + @pattern.inspect

				  # Last both pattern edges must be not visited, neighbored and alternating charged
				  last_edges = @pattern.last(2)
				  if last_edges.size == 2 and !@visited_edges.include?(last_edges[-1].label) and !@visited_edges.include?(last_edges[-2].label) and neighbored_edges(last_edges) and alternate_check(last_edges)
					  edges = @pattern.pop(2)
					  # labels = edges.map {|edge| edge.label }
					  subpath = edges.map {|edge| edge.direction}
					  @cum_path_score += calc_path_pattern(subpath.join)

					  # puts "subpath: " + JSON.dump(subpath) + labels.inspect
					  # puts "cumPathScore = " + @cum_path_score.to_s
					  # Add edges to visited
					  edges.each { |edge| @visited_edges.push(edge.label) }
					  # puts @visited_edges.inspect
				  else
					  # Drop it
					  @pattern.pop if !@pattern.empty? && @pattern.last.dst_node_id > child_edge.dst_node_id
				  end
				  # puts @labels.key(child_edge.dst_node_id) + " nth_pattern_after_pop: " + @pattern.inspect
				  # puts "UP to " + parent_label
			  }
		  # Current node has no (more) edges
		  else
			  return nil
		  end
	  else
		  # puts parent_label + " already visited"
	  end
  end

  # # Find cliques helper
  # def clique_helper(parent_label)
		# # Remember nodes visited
		# @pattern.push(parent_label)
		# @visited.push(parent_label)
	 #  # Check nodes and straight edges
	 #  if has_edges(parent_label)
		#   parent_node = @nodes[@labels[parent_label]]
		#   parent_node.adj_list.each { |child_edge|
		# 	  child_label = @labels.key(child_edge.dst_node_id)
		# 		clique_helper(child_label)
		#   }
  #   end
  # end
  #
  # # Find cliques
  # def find_cliques
	 #  puts @labels
	 #  @visited = Array.new;
	 #  @labels.each_key { |root_label|
		#   @pattern = Array.new
		#   if !@visited.include?(root_label)
		# 	  clique_helper(root_label)
		#   end
		#   if @pattern.size > 1
		#     puts @pattern.to_s
	 #    end
	 #  }
  # end

  # Get paths in Heptad net by using Depth First Search (DFS)
  # @return [Float]
  def get_paths
		@visited_edges = Array.new
		# puts @labels
		# Go through all single labels
		@labels.each_key { |root_label|
				# Check nodes and straight edges
				@pattern = Array.new
				@visited = Array.new;
				check_node(root_label, 'left')
				# if @pattern.size > 1
				# 	puts @pattern.to_s
				# end
		}
		@labels.each_key { |root_label|
				@pattern = Array.new
				@visited = Array.new;
				check_node(root_label, 'right')
				# if @pattern.size > 1
				# 	puts @pattern.to_s
				# end
		}
		# puts @cum_path_score
		# Check around the corner edges
		set_diff = Set.new(@edge_labels.keys).-(Set.new(@visited_edges)).to_a
		set_diff_comp = Array.new(set_diff)
		# puts set_diff_comp.to_s
		set_diff_comp.each { |el|
			comp = set_diff.pop()
			set_diff.each { |cur|
				if neighbored_edges([get_edge(comp), get_edge(cur)]) and alternate_check([get_edge(comp), get_edge(cur)])
					@cum_path_score += calc_path_pattern([get_edge(comp).direction, get_edge(cur).direction].join)
				end
			}
		}
		# puts @cum_path_score
		@cum_path_score
  end

end

class NodeAlreadyDefinedException < RuntimeError
  attr :ok_to_retry
  def initialize(ok_to_retry)
    @ok_to_retry = ok_to_retry
  end
end


# # Testing
# g = Graph.new
# g.addNode(841, "R", "a", "0.999")
# g.addNode(842, "E", "b", "0.999")
# g.addNode(843, "R", "c", "0.999")
# g.addNode(844, "R", "d", "0.999")
# g.addNode(845, "E", "e", "0.999")
#
# g.addNode(846, "A", "e", "0.999")
# g.addNode(847, "E", "e", "0.999")
# g.addNode(848, "L", "e", "0.999")
# g.addNode(849, "R", "e", "0.999")
# g.addNode(850, "A", "e", "0.999")
#
# g.addNode(851, "Q", "e", "0.999")
# g.addNode(852, "Q", "e", "0.999")
# g.addNode(853, "E", "e", "0.999")
# g.addNode(854, "E", "e", "0.999")
# g.addNode(855, "E", "e", "0.999")
#
# g.addNode(856, "T", "e", "0.999")
# g.addNode(857, "R", "e", "0.999")
# g.addNode(858, "K", "e", "0.999")
# g.addNode(859, "Q", "e", "0.999")
# g.addNode(860, "Q", "e", "0.999")
#
# g.addNode(861, "E", "e", "0.999")
# g.addNode(862, "L", "e", "0.999")
# g.addNode(863, "E", "e", "0.999")
# g.addNode(864, "A", "e", "0.999")
# g.addNode(865, "L", "e", "0.999")
#
# g.addNode(866, "Q", "e", "0.999")
# g.addNode(867, "K", "e", "0.999")
# g.addNode(868, "S", "e", "0.999")
# g.addNode(869, "Q", "e", "0.999")
# g.addNode(870, "K", "e", "0.999")

# g.addNode(871, "E", "e", "0.999")
# g.addNode(872, "A", "e", "0.999")
# g.addNode(873, "E", "e", "0.999")
# g.addNode(874, "L", "e", "0.999")
# g.addNode(875, "T", "e", "0.999")

# g.addNode(876, "R", "e", "0.999")
# g.addNode(877, "E", "e", "0.999")
# g.addNode(878, "L", "e", "0.999")
# g.addNode(879, "E", "e", "0.999")
# g.addNode(880, "K", "e", "0.999")

# g.addNode(881, "Q", "e", "0.999")
# g.addNode(882, "K", "e", "0.999")
# g.addNode(883, "E", "e", "0.999")
# g.addNode(884, "N", "e", "0.999")
# g.addNode(885, "K", "e", "0.999")

# g.addNode(886, "Q", "e", "0.999")
# g.addNode(887, "V", "e", "0.999")
# g.addNode(888, "E", "e", "0.999")
# g.addNode(889, "E", "e", "0.999")
#
# # puts g.getNodeID('Q818')
#
# g.addEdge('R841', 'R844', 'right', 1.0)
# g.addEdge('R841', 'E845', 'left', 1.0)
# g.addEdge('E842', 'E845', 'right', 1.0)
# g.addEdge('E842', 'A846', 'left', 1.0)
#
# g.addEdge('R843', 'A846', 'right', 1.0)
# g.addEdge('R843', 'E847', 'left', 1.0)
# g.addEdge('R844', 'E847', 'right', 1.0)
# g.addEdge('R844', 'L848', 'left', 1.0)
#
# g.addEdge('E845', 'L848', 'right', 1.0)
# g.addEdge('E845', 'R849', 'left', 1.0)
# g.addEdge('A846', 'R849', 'right', 1.0)
# g.addEdge('A846', 'A850', 'left', 1.0)
#
# g.addEdge('E847', 'A850', 'right', 1.0)
# g.addEdge('E847', 'Q851', 'left', 1.0)
# g.addEdge('L848', 'Q851', 'right', 1.0)
# g.addEdge('L848', 'Q852', 'left', 1.0)
#
# g.addEdge('R849', 'Q852', 'right', 1.0)
# g.addEdge('R849', 'E853', 'left', 1.0)
# g.addEdge('A850', 'E853', 'right', 1.0)
# g.addEdge('A850', 'E854', 'left', 1.0)
#
# g.addEdge('Q851', 'E854', 'right', 1.0)
# g.addEdge('Q851', 'E855', 'left', 1.0)
# g.addEdge('Q852', 'E855', 'right', 1.0)
# g.addEdge('Q852', 'T856', 'left', 1.0)
#
# g.addEdge('E853', 'T856', 'right', 1.0)
# g.addEdge('E853', 'R857', 'left', 1.0)
# g.addEdge('E854', 'R857', 'right', 1.0)
# g.addEdge('E854', 'K858', 'left', 1.0)
#
# g.addEdge('E855', 'K858', 'right', 1.0)
# g.addEdge('E855', 'Q859', 'left', 1.0)
# g.addEdge('T856', 'Q859', 'right', 1.0)
# g.addEdge('T856', 'Q860', 'left', 1.0)
#
# g.addEdge('R857', 'Q860', 'right', 1.0)
# g.addEdge('R857', 'E861', 'left', 1.0)
# g.addEdge('K858', 'E861', 'right', 1.0)
# g.addEdge('K858', 'L862', 'left', 1.0)
#
# g.addEdge('Q859', 'L862', 'right', 1.0)
# g.addEdge('Q859', 'E863', 'left', 1.0)
# g.addEdge('Q860', 'E863', 'right', 1.0)
# g.addEdge('Q860', 'A864', 'left', 1.0)
#
# g.addEdge('E861', 'A864', 'right', 1.0)
# g.addEdge('E861', 'L865', 'left', 1.0)
# g.addEdge('L862', 'L865', 'right', 1.0)
# g.addEdge('L862', 'Q866', 'left', 1.0)
#
# g.addEdge('E863', 'Q866', 'right', 1.0)
# g.addEdge('E863', 'K867', 'left', 1.0)
# g.addEdge('A864', 'K867', 'right', 1.0)
# g.addEdge('A864', 'S868', 'left', 1.0)
#
# g.addEdge('L865', 'S868', 'right', 1.0)
# g.addEdge('L865', 'Q869', 'left', 1.0)
# g.addEdge('Q866', 'Q869', 'right', 1.0)
# g.addEdge('Q866', 'K870', 'left', 1.0)
#
# g.addEdge('K867', 'K870', 'right', 1.0)
# g.addEdge('K867', 'E871', 'left', 1.0)
# g.addEdge('S868', 'E871', 'right', 1.0)
# g.addEdge('S868', 'A872', 'left', 1.0)
#
# g.addEdge('Q869', 'A872', 'right', 1.0)
# g.addEdge('Q869', 'E873', 'left', 1.0)
# g.addEdge('K870', 'E873', 'right', 1.0)
# g.addEdge('K870', 'L874', 'left', 1.0)
#
# g.addEdge('E871', 'L874', 'right', 1.0)
# g.addEdge('E871', 'T875', 'left', 1.0)
# g.addEdge('A872', 'T875', 'right', 1.0)
# g.addEdge('A872', 'R876', 'left', 1.0)
#
# g.addEdge('E873', 'R876', 'right', 1.0)
# g.addEdge('E873', 'E877', 'left', 1.0)
# g.addEdge('L874', 'E877', 'right', 1.0)
# g.addEdge('L874', 'L878', 'left', 1.0)
#
# g.addEdge('T875', 'L878', 'right', 1.0)
# g.addEdge('T875', 'E879', 'left', 1.0)
# g.addEdge('R876', 'E879', 'right', 1.0)
# g.addEdge('R876', 'K880', 'left', 1.0)
#
# g.addEdge('E877', 'K880', 'right', 1.0)
# g.addEdge('E877', 'Q881', 'left', 1.0)
# g.addEdge('L878', 'Q881', 'right', 1.0)
# g.addEdge('L878', 'K882', 'left', 1.0)
#
# g.addEdge('E879', 'K882', 'right', 1.0)
# g.addEdge('E879', 'E883', 'left', 1.0)
# g.addEdge('K880', 'E883', 'right', 1.0)
# g.addEdge('K880', 'N884', 'left', 1.0)
#
# g.addEdge('Q881', 'N884', 'right', 1.0)
# g.addEdge('Q881', 'K885', 'left', 1.0)
# g.addEdge('K882', 'K885', 'right', 1.0)
# g.addEdge('K882', 'Q886', 'left', 1.0)
#
# g.addEdge('E883', 'Q886', 'right', 1.0)
# g.addEdge('E883', 'V887', 'left', 1.0)
# g.addEdge('N884', 'V887', 'right', 1.0)
# g.addEdge('N884', 'E888', 'left', 1.0)
#
# g.addEdge('K885', 'E888', 'right', 1.0)
# g.addEdge('K885', 'E889', 'left', 1.0)
# g.addEdge('Q886', 'E889', 'right', 1.0)
#
# puts g.getPaths()
# # print edges = g.getEdges('E819')
# # print edges[0].dst

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Pedro Cerqueira
    github: @pedrorvc
    
DESCRIPTION

    This script serves to create xml files contaning the information necessary
    for the execution of BRIG (Blast Ring Image Generator), reducing the time
    performing the tedious task of setting up all the information on the GUI 
    and provides a quick way to produce an image.
    
    The arguments for this script provide some (but not all) 
    of the available options in BRIG, which were the ones I used to change the most.
    
    USAGE:

    brigaid.py -q reference_sequence.fna -rfd path/to/reference/dir -od path/to/output/dir -of path/to/output/dir/output_file
                                -oi path/to/output/BRIG/output_image -t Image_title -a annotation_file.gbk --genes genes_of_interest.txt
                                --contig-order contig_order.tsv 
    
"""

import argparse
import csv
import os
import xml.etree.ElementTree as ET
from collections import OrderedDict
from xml.dom import minidom

from Bio import SeqIO
from matplotlib import cm


def listdir_fullpath(path):
    """ Gets the full path of the files from a directory

        Args:
            path (str): full path to a directory

        Returns:
            list containing the full path of every file contained in the input directory
    
    """
    
    return [os.path.join(path, f) for f in os.listdir(path)]


def ring_attributes(colour, name, position):
    """ Creates ring attributes.
    
        Args:
            colour (str): color of the ring.
            name (str): name of the ring.
            position (str): position of the ring.
            
        Returns:
            ring_attrs (dict): attributes of any regular ring of the BRIG xml.
    """
    
    ring_attrs = {"colour" : colour,
                 "name": name,
                 "position" : position,
                 "upperInt" : "90",
                 "lowerInt" : "70",
                 "legend" : "yes",
                 "size" : "30",
                 "labels" : "no",
                 "blastType" : "blastn"}
    
    return ring_attrs


def annotation_ring_attributes(position):
    """ Creates annotation ring attributes.
    
        Args:
            position (str): position of the ring.
            
        Returns:
            annotation_ring_attrs (dict): attributes of the annotation ring of the BRIG xml.
    """
    
    annotation_ring_attrs = {"colour" : '172,14,225',
                 "name": 'null',
                 "position" : position,
                 "upperInt" : "70",
                 "lowerInt" : "50",
                 "legend" : "yes",
                 "size" : "30",
                 "labels" : "no",
                 "blastType" : "blastn"}
    
    return annotation_ring_attrs


def create_feature_attrs(label, colour, decoration, start, stop):
    """ Create attributes for the Feature SubElements of the annotation ring.
    
        Args:
            label (str): name of the gene/CDS to annotate
            colour (str): colour of the decoration for the annotation
            decoration (str): shape of the gene/CDS to annotate, for example, 'clockwise-arrow'
            start (str): start of the gene/CDS to annotate
            stop (str): stop of the gene/CDS to annotate
            
        Results:
            feature_element_attrs (dict): attributes of the feature element.
            feature_range_element_attrs (dict): attributes of the feature range element
    """
    
    feature_element_attrs = {'label' : label,
                             'colour' : colour,
                             'decoration' : decoration}
    
    feature_range_element_attrs = {'start' : start,
                                   'stop' : stop}
    
    return feature_element_attrs, feature_range_element_attrs


def create_annotation_ring_tsv(annotation_ring, annotation_file):
    """ Uses a tsv file to annotate the reference genome.
    
        Args:
            annotation_ring: ElementTree SubElement object containing the 'ring' tag and its attributes.
            annotation_file (str): Full path to the file containing annotations for the reference genome.
            
            
    """
    
    with open(annotation_file) as tsvfile:
      reader = csv.DictReader(tsvfile, dialect='excel-tab')
      
      # Obtain the annotations from the file contents
      for row in reader:
          start = row['#START']
          stop = row['STOP']
          label = row['Label']
          colour = row['Colour']
          decoration = row['Decoration']
            
          # Create xml attributes
          feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, colour, decoration, start, stop)
          
          # Create xml elements
          feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
          feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
          
                
                
def annotation_ring_feature_elements_gbk_concat(annotation_ring, record, genome_size=False):
    """ Creates the annotation ring feature elements, using a concatenated Genbank annotation file.
    
        Args:
            annotation_ring: ElementTree SubElement object containing the 'ring' tag and its attributes.
            record (SeqRecord): Object of BioPython containing the information of the input Genbank.
            genome_size (bool): Size of genome. Integer when a Genbank divided by contigs is provided.
                                                Boolean (False) when a concatenated Genbank is provided.
                                                       
    """
    
    #if type(genome_size) == int:

    # Obtain the features of the Genbank file records
    for fea in record.features:
        
        # Get the start and end position of the genome
        # Also get the strand
        if fea.type == 'CDS':
            start = str(fea.location.start.position)
            end = str(fea.location.end.position)
            strand = fea.location.strand
            
            # Get the label of the gene or product
            if 'gene' in fea.qualifiers:
                label = str(fea.qualifiers['gene'][0])
            elif 'product' in fea.qualifiers:
                product = fea.qualifiers['product'][0]
                label = str(product)
            else:
                continue
            
            # Define the decoration of the annotation based on the strand
            if strand == -1:
                decoration = 'counterclockwise-arrow'
            
            elif strand == 1:
                decoration = 'clockwise-arrow'
                
            # Create xml attributes
            feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
                
            # Create xml elements
            feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
            feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)
        
        # If a genome size is provided, get the size of the records
        if type(genome_size) == int:
            
            if fea.type == 'source':
                size = fea.location.end.position

    try:
        size
        genome_size += size
        
        return genome_size
    
    except NameError:
        pass
    

    
            

def annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, record, genes, genome_size=False):
    """ Creates the annotation ring feature elements, using a concatenated Genbank annotation file 
        and specific gene annotations.
    
        Args:
            annotation_ring: ElementTree SubElement object containing the 'ring' tag and its attributes.
            record (SeqRecord): Object of BioPython containing the information of the input Genbank.
            genome_size (bool): Size of genome. Integer when a Genbank divided by contigs is provided.
                                                Boolean (False) when a concatenated Genbank is provided.
                                        
    """
            
    for f in record.features:
    
        if f.type == 'CDS':
            
            # Find the 'gene' tag and determine if the gene belongs to the specified genes to be annotated
            if 'gene' in f.qualifiers and f.qualifiers['gene'][0] in genes:
                label = f.qualifiers['gene'][0]
            elif 'product' in f.qualifiers and f.qualifiers['product'][0] in genes:
                product = f.qualifiers['product'][0]
                label = product
            else:
                continue
                
            # Determine the start, stop and strand of the gene
            start = str(f.location.start.position + genome_size)
            end = str(f.location.end.position + genome_size)
            strand = f.location.strand
                
            
            # Define the decoration of the annotation based on the strand
            if strand == -1:
                decoration = 'counterclockwise-arrow'
            
            elif strand == 1:
                decoration = 'clockwise-arrow'
                            
            # Create xml attributes
            feature_element_attrs, feature_range_element_attrs = create_feature_attrs(label, "black", decoration, start, end)
                
            # Create xml elements
            feature_element = ET.SubElement(annotation_ring, 'feature', attrib=feature_element_attrs)        
            feature_range_element = ET.SubElement(feature_element, 'featureRange', attrib=feature_range_element_attrs)

        # If a genome size is provided, get the size of the records
        if type(genome_size) == int:
            
            if f.type == "source":
                size = f.location.end.position
        
    try:
        size
        genome_size += size
        
        return genome_size
    
    except NameError:
        pass
      

def create_annotation_ring_gbk_concat(annotation_ring, annotation_file, genes_of_interest, records):
    """ Create annotation ring using a concatenated Genbank annotation file.
    
        Args:
            annotation_ring: ElementTree SubElement object containing the 'ring' tag and its attributes.
            annotation_file (str): Full path to the file containing annotations for the reference genome.
            genes_of_interest (str): Full path to the file containing the genes to search for in the Genbank file.
            records (SeqRecord): Object of BioPython containing the information of the input Genbank.
            
    """
    
    if genes_of_interest != []:
        
        # Get the genes to serach in the Genbank file
        with open(genes_of_interest, "r") as f:
            genes = f.readlines()
            genes = [gene.rstrip() for gene in genes]
        
        # Create feature elements of the annotation ring
        for seq_record in records:
            annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, seq_record, genes)

    else:
        
        for seq_record in records:
            annotation_ring_feature_elements_gbk_concat(annotation_ring, seq_record)
        

def create_annotation_ring_gbk_contigs(annotation_ring, annotation_file, records, genes_of_interest, contig_order):
    """ Create annotation ring using a Genbank annotation file divided by contigs.
    
        Args:
            annotation_ring: ElementTree SubElement object containing the 'ring' tag and its attributes.
            annotation_file (str): Full path to the file containing annotations for the reference genome.
            genes_of_interest (str): Full path to the file containing the genes to search for in the Genbank file.
            records (SeqRecord): Object of BioPython containing the information of the input Genbank.
            contig_order (str): Full path to the file containing the order of the contigs.

    """
    
    if contig_order != []:
        
        
        with open(contig_order) as tsvfile:
            reader = csv.DictReader(tsvfile, dialect='excel-tab')
            
            # Create an OrderedDict with the contents of the file
            # The keys are the order are a number representing the order of the contig
            # The values are the names of the contigs
            content_dict = OrderedDict()
            for r in reader:
                content_dict[r["order"]] = r["contig"]
    
        # Create an OrderedDict with the content of each contig
        # The keys are the names of the contigs
        # The values are SeqRecord objects from BipPython
        seq_records_dict = OrderedDict()
        for record in records:            
            seq_records_dict[record.id] = record
            
        if genes_of_interest != []:
            
            with open(genes_of_interest, "r") as f:
                genes = f.readlines()
                genes = [gene.rstrip() for gene in genes]
                
            genome_size = 0
            for i in range(1, len(records)+1):                                   
                ord_record = seq_records_dict[content_dict[str(i)]]
                                
                gsize = annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, ord_record, genes, genome_size)
                
                genome_size = gsize
        else:
            
            genome_size = 0
            for i in range(1, len(records)+1):                                   
                ord_record = seq_records_dict[content_dict[str(i)]]
                
                gsize = annotation_ring_feature_elements_gbk_concat(annotation_ring, ord_record, genome_size)
                
                genome_size = gsize

    else:
    
        if genes_of_interest != []:
            
            with open(genes_of_interest, "r") as f:
                genes = f.readlines()
                genes = [gene.rstrip() for gene in genes]

            
            for seq_record in records:
                annotation_ring_feature_elements_genes_of_interest_gbk_concat(annotation_ring, seq_record, genes)
    
        else:
            for seq_record in records:
                annotation_ring_feature_elements_gbk_concat(annotation_ring, seq_record)


def write_xml(root_elem, output_file):
    """ Writes a xml file.
    
        Args:
            root_elem is a ElementTree Element object containing all the information 
            required for the output file.
            output_file (str): full path to the output file
        
            
    """

    xml_file = ET.tostring(root_elem, encoding='utf8').decode('utf8')
    
    pretty_xml_file = minidom.parseString(xml_file).toprettyxml(indent='    ')
    
    output_file = output_file + ".xml"
    
    with open(output_file, "w") as f:
        f.write(pretty_xml_file)


####### Create xml elemnts
        
# Create root element
        
def create_root_element(blast_options, legend_position, query_file, 
                               output_folder, image_output_file, title, image_format):
    """
    Creates the root element of the xml file and its attributes.
    
    Args: 
        blast_options (str): additional options for blast, for example, -evalue or num_threads
        legend_position (str): position of the legend on the image
        query_file (str): full path to the query file
        output_folder (str): full path to the output folder
        image_output_file (str): full path to the image output file
        title (str): title of the output image
        image_format (str): format of the image output file
            
    Returns:
        root: ElementTree Element object containing the BRIG tag and its attributes
        
    """

    root_attrs = {"blastOptions" : blast_options,
                  "legendPosition" : legend_position,
                  "queryFile" : query_file,
                  "outputFolder" : output_folder,
                  "blastPlus" : "yes",
                  "outputFile" : os.path.join(output_folder, image_output_file),
                  "title" : title,
                  "imageFormat" : image_format,
                  "queryFastaFile" : query_file,
                  "cgXML" : os.path.join(output_folder + "/scratch", os.path.basename(query_file) + ".xml")}
    
    
    root = ET.Element('BRIG', attrib=root_attrs)
    
    return root

#### Create root children

# Create cgview_settings element
    
def create_cgview_settings_element(root, height, width):
    """ Creates the cgview_settings element of the xml file and its attributes.
    
        Args:
            root: ElementTree Element object containing the BRIG tag and its attributes.
            height (str): height of the output image in pixels
            width (str): width of the output image in pixels
        
        Returns:
            cgview_settings: ElementTree SubElement object containing the cgview settings tag and its attributes
    """

    cgview_settings_attrs = {"arrowheadLength" : "medium",
                             "backboneColor" : "black",
                             "backboneRadius" : "600",
                             "backboneThickness" : "medium",
                             "backgroundColor" : "white",
                             "borderColor" : "black",
                             "featureSlotSpacing" : "medium",
                             "featureThickness" : "30",
                             "giveFeaturePositions" : "false",
                             "globalLabel" : "true",
                             "height" : height,
                             "isLinear" : "false",
                             "labelFont" : "SansSerif,plain,25",
                             "labelLineLength" : "medium",
                             "labelLineThickness" : "medium",
                             "labelPlacementQuality" : "best",
                             "labelsToKeep" : "1000",
                             "longTickColor" : "black",
                             "minimumFeatureLength" : "medium",
                             "moveInnerLabelsToOuter" :"true",
                             "origin" : "12",
                             "rulerFont" : "SansSerif,plain,35",
                             "rulerFontColor" : "black",
                             "rulerPadding" : "40",
                             "rulerUnits" : "bases",
                             "shortTickColor" : "black",
                             "shortTickThickness" : "medium",
                             "showBorder" : "false",
                             "showShading" : "true",
                             "showWarning" : "false",
                             "tickDensity" : "0.2333",
                             "tickThickness" : "medium",
                             "titleFont" : "SansSerif,plain,45",
                             "titleFontColor" : "black",
                             "useColoredLabelBackgrounds" : "false",
                             "useInnerLabels" : "true",
                             "warningFont" : "Default,plain,35",
                             "warningFontColor" : "black",
                             "width" : width,
                             "zeroTickColor" : "black",
                             "tickLength" : "medium"}


    cgview_settings = ET.SubElement(root, 'cgview_settings', attrib=cgview_settings_attrs)
        
    return cgview_settings


# Create brig_settings element 
    
def create_brig_settings_element(root, java_memory):
    """ Creates the brig_settings element of the xml file and its attributes.
    
        Args:
            root: ElementTree Element object containing the BRIG tag and its attributes.
            java_memory (str): amount of memory (in bytes) java is allowed to use for BRIG
    
        Returns:
            brig_settings: ElementTree SubElement object containing the brig settings tag and its attributes

    """
    
    

    brig_settings_attrs = {"Ring1" : "172,14,225",
                           "Ring2" : "222,149,220",
                           "Ring3" : "161,221,231",
                           "Ring4" : "49,34,221",
                           "Ring5" : "116,152,226",
                           "Ring6" : "224,206,38",
                           "Ring7" : "40,191,140",
                           "Ring8" : "158,223,139",
                           "Ring9" : "226,38,122",
                           "Ring10" :"211,41,77",
                           "defaultUpper" : "70",
                           "defaultLower" : "50",
                           "defaultMinimum" : "50",
                           "genbankFiles" : "gbk,gb,genbank",
                           "fastaFiles" : "fna,faa,fas,fasta,fa",
                           "emblFiles" : "embl",
                           "blastLocation" : "",
                           "divider" : "3",
                           "multiplier" : "3",
                           "memory" : java_memory,
                           "defaultSpacer" : "0"}
    
    brig_settings = ET.SubElement(root, 
                                  "brig_settings", 
                                  attrib=brig_settings_attrs)
    
    return brig_settings


## Create special element

def create_special_element(root):
    """Creates the 'special' element of the xml file and its attributes
    
        Args:
            root: ElementTree Element object containing the BRIG tag and its attributes.
            
        Returns:
            gc_content_special: ElementTree SubElement object containing the 'special' tag and its attributes
            gc_skew_special: ElementTree SubElement object containing the 'special' tag and its attributes
                

    """
    
    gc_content_special = ET.SubElement(root, 'special', attrib={'value' : 'GC Content'})
    gc_skew_special = ET.SubElement(root, 'special', attrib={'value' : 'GC Skew'})
    
    return gc_content_special, gc_skew_special


# Create reference dir element

def create_reference_directory_element(root, reference_directory):
    """ Creates the 'reference directory' element of the xml file and its attributes.

        Args:
            root: ElementTree Element object containing the 'BRIG' tag and its attributes.
            reference_directory (str): full path to the reference directory that contains 
                                        the fasta files used to build the rings.  
    
        Returns:
            ref_file: ElementTree SubElement object containing the 'refFile' tag and its attributes

    
    """

    ref_dir = ET.SubElement(root, 
                            "refDir", 
                            attrib={"location" : reference_directory})
    
    # Obtain the full path for all the files in the directory
    ref_dir_list = listdir_fullpath(reference_directory)
    
    for f in ref_dir_list:
        ref_file = ET.SubElement(ref_dir, 
                                 "refFile", 
                                 attrib={"location" : f})
        
    return ref_file


# Create the ring where the annotations are defined

def create_annotation_ring(root, reference_directory, annotation_file, genes_of_interest, contig_order):
    """ Creates the ring that will contain the annotations for the reference genome.
    
        Args:
            root: ElementTree Element object containing the 'BRIG' tag and its attributes.
            reference_directory (str): full path to the reference directory that contains 
                                        the fasta files used to build the rings.  
            annotation_file (str): Full path to the file containing annotations for the reference genome.
            genes_of_interest (str): Full path to the file containing a list of specific genes.
            contig_order (str): Full path to the tab-delimited file containing the order of the contigs.
            
            
    """
    
    # Determine the position of the annotation ring, which will be the position after the last reference genome
    ring_position = len(os.listdir(reference_directory)) + 2
    
    # Create the annotation ring element
    annotation_ring = ET.SubElement(root, 'ring', attrib=annotation_ring_attributes(str(ring_position)))
    
    # Check for tab-delimited annotation file input 
    if list(SeqIO.parse(annotation_file, "genbank")) == []:
        create_annotation_ring_tsv(annotation_ring, annotation_file)
    
    else:
        
        # Get the records of the Genbank file
        records = [r for r in SeqIO.parse(annotation_file, "genbank")]
        
        ### Check if a contig order file has been provided
        
        if len(records) > 1:   # If more than 1 record exists, then the Genbank file is divided by contigs
            create_annotation_ring_gbk_contigs(annotation_ring, annotation_file, records, genes_of_interest, contig_order)
        else:
            create_annotation_ring_gbk_concat(annotation_ring, annotation_file, genes_of_interest, records)
            


## Create remaining rings
        
def create_ring_element(root, reference_directory, colormap):
    """ Creates the ring elements of the xml file, containing the position and color of the rings.
    
        Args: 
            root: ElementTree Element object containing the 'BRIG' tag and its attributes.
            
            reference_directory (str): full path to the reference directory that contains 
                                        the fasta files used to build the rings.  
            colormap (str): name of the colormap (available in matplotlib) to use for the color of the rings
        
        Returns:
            ring_number_element: ElementTree SubElement object containing the 'ring' tag and its attributes
            ring_sequence_element: ElementTree SubElement object containing the 'sequence' tag and its attributes
            
    """
    
    ref_dir_list = listdir_fullpath(reference_directory)    
    
    # Gets the colormap from matplotlib with as many colors as the number of files
    cmap = cm.get_cmap(colormap, len(ref_dir_list))

    list_colormap = cmap.colors.tolist()
    
    # Remove the fourth element (transparency) because it is not necessary
    colors_to_use = []
    for l in list_colormap:
        convert = [round(x * 255) for x in l]
        convert.pop()
        colors_to_use.append(convert)
    
    #reversed_colors_to_use = colors_to_use[::-1]
    
    # Check if the user provided an order for the rings
    has_digit = [os.path.basename(x).split("_")[0].isdigit() for x in ref_dir_list] 
    
    if True in has_digit:
        # Obtain the ring positions
        ring_positions = [os.path.basename(x).split("_")[0] for x in ref_dir_list]
        
        # Reverse sort the positions of the rings, because they will be created
        # in a descending order of their positions
        ring_positions.sort(reverse=True)
        ref_dir_list.sort(reverse=True)
        
        for ring in range(len(ref_dir_list)):
            
            # The ring positions start at 2 due to the special rings (GC Content and GC Skew)
            ring_position = int(ring_positions[ring]) + 1
            
            # Select a color for the ring
            ring_color = ",".join([str(e) for e in colors_to_use[ring]])
            
            # Define the name of the ring
            ring_name = os.path.basename(ref_dir_list[ring]).split("_")[1]
            
            # Create the xml elements
            ring_number_element = ET.SubElement(root, 
                                    'ring',
                                    ring_attributes(ring_color, ring_name, str(ring_position)))
            
            ring_sequence_element = ET.SubElement(ring_number_element, 
                                      "sequence", 
                                      attrib={"location" : ref_dir_list[ring]})
        
        
    else:
        # Sort files by lowercase
        ref_dir_list.sort(key=lambda y: y.lower())
 
        # The number of rings starts at 2 due to the GC Content and GC Skew
        ring_number = len(ref_dir_list) + 1
        for ring in range(len(ref_dir_list)):
            
            # Select a color for the ring
            ring_color = ",".join([str(e) for e in colors_to_use[ring]])
            
            # Define the name of the ring
            ring_name = os.path.basename(ref_dir_list[ring]).split("_")[0]
            
            # Create the xml elements
            ring_number_element = ET.SubElement(root, 
                                                'ring',
                                                ring_attributes(ring_color, ring_name, str(ring_number)))
            
            ring_sequence_element = ET.SubElement(ring_number_element, 
                                                  "sequence", 
                                                  attrib={"location" : ref_dir_list[ring]})
            
            ring_number -= 1
    
    return ring_number_element, ring_sequence_element
  
## Create special rings

def create_special_ring_element(root):
    """ Create the 'special' ring element and its attributes.
    
        Args:
            root: ElementTree Element object containing the 'BRIG' tag and its attributes.
            
        Returns:
            gc_content_location: ElementTree SubElement object containing the 'sequence' tag and its attributes
            gc_skew_location: ElementTree SubElement object containing the 'sequence' tag and its attributes
    """
    
    # Create ring attributes
    gc_content_ring_attrs = ring_attributes('225,0,0', "GC Content", "0")
    gc_skew_ring_attrs = ring_attributes('225,0,0', "GC Skew", "1")
    
    # Add ring element to root
    gc_skew_ring = ET.SubElement(root, 'ring', attrib=gc_skew_ring_attrs)
    gc_content_ring = ET.SubElement(root, 'ring', attrib=gc_content_ring_attrs)
    
    # Add sequence element to ring
    gc_content_location = ET.SubElement(gc_content_ring, 'sequence', attrib={'location' : 'GC Content'})
    gc_skew_location = ET.SubElement(gc_skew_ring, 'sequence', attrib={'location' : 'GC Skew'})

    return gc_content_location, gc_skew_location



def main(query_file, reference_directory, output_folder, output_xml, image_output_file, title, annotation_file,
         genes_of_interest, contig_order, blast_options, legend_position, image_format, height, width, java_memory, colormap):

    
    root = create_root_element(blast_options, legend_position, query_file, 
                               output_folder, image_output_file, title, image_format)
    
    cgview_settings = create_cgview_settings_element(root, height, width)
    
    brig_settings = create_brig_settings_element(root, java_memory)
    
    special = create_special_element(root)

    refdir = create_reference_directory_element(root, reference_directory)
        
    if annotation_file:
        create_annotation_ring(root, reference_directory, annotation_file, genes_of_interest, contig_order)
    
    rings = create_ring_element(root, reference_directory, colormap)
        
    special_ring = create_special_ring_element(root)
    
    write_xml(root, output_xml)

    print("\n File written to {}".format(output_xml))

def parse_arguments():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-q', '--query', type=str, required=True, dest='query_file',
                        help='Path to the query/reference FASTA file.')
    
    parser.add_argument('-rfd', '--ref_dir', type=str, required=True, dest='reference_directory',
                        help='Path to the directory where the FASTA files to compare against the reference are located.')
    
    parser.add_argument('-od', '--out_dir', type=str, required=True, dest='output_folder',
                        help='Path to the output directory for the results of BRIG.')
    
    parser.add_argument('-of', '--out_xml', type=str, required=True, dest='output_file',
                        help='Path to the output of this script.')

    parser.add_argument('-oi', '--out_img', type=str, required=True, dest='image_output_file',
                        help='Path to the output file of the resulting image of BRIG.')

    parser.add_argument('-t', '--title', type=str, required=True, dest='title',
                        help='Title of the resulting image from BRIG.')
    
    parser.add_argument('-a', '--annotation', type=str, required=False, dest='annotation_file', default=False,
                        help='File containing annotations for the reference genome. '
                        'The annoation file can be a tab-delimited file (.tsv) or a Genbank format file (.gbk, .gb)')
    
    parser.add_argument('--genes', type=str, required=False, dest='genes_of_interest', default=[],
                        help='File containing a list of specific genes (one gene per line) to search when a Genbank annotation file is provided. ')
   
    parser.add_argument('--contig_order', type=str, required=False, dest='contig_order', default=[],
                        help='Tab-delimited file containing the order of the contigs when a Genbank (divided by contigs) annotation file is provided. '
                             'Example:  order    contig '
                                          '1     Contig8')

    parser.add_argument('-b', '--blast_options', type=str, required=False, dest="blast_options", default="-evalue 0.001 -num_threads 6",
                        help='Options for running BLAST.')
    
    parser.add_argument('-l', '--legend_pos', type=str, required=False, dest="legend_position", default="middle-right",
                        help='Positon of the legend on the resulting image.'
                        'The options available are upper, center or lower, '
                        'paired with left, center or right')
    
    parser.add_argument('-if', '--image_format', type=str, required=False, dest="image_format", default="jpg",
                        help='Format of the resulting image file.'
                        'The available options are: jpg, png, svg or svgz.')

    parser.add_argument('-ht', '--height', type=str, required=False, dest="height", default="3000",
                        help='Height (in pixels) of the resulting image.')
    
    parser.add_argument('-wd', '--width', type=str, required=False, dest="width", default="3000",
                        help='Width (in pixels) of the resulting image.')

    parser.add_argument('-jm', '--java_memory', type=str, required=False, dest="java_memory", default="1500",
                        help='Amount of memory (in bytes) that Java is allowed to use for BRIG.')

    parser.add_argument('-cm', '--colormap', type=str, required=False, dest="colormap", default="viridis",
                        help='Colormap from matplotlib to use for the color of the rings. '
                        'The available options are: viridis, plasma, inferno, magma and cividis.'
                        'More options for colormaps at: '
                        'https://matplotlib.org/users/colormaps.html')

    args = parser.parse_args()
    
    return [args.query_file, args.reference_directory, args.output_folder, args.output_file,
            args.image_output_file, args.title, args.annotation_file, args.genes_of_interest, args.contig_order,
            args.blast_options, args.legend_position, args.image_format, args.height, args.width, args.java_memory, args.colormap]
    
if __name__ == '__main__':
    
    args = parse_arguments()
    main(args[0], args[1], args[2], args[3], args[4], args[5], args[6],
         args[7], args[8], args[9], args[10], args[11], args[12], args[13],
         args[14], args[15])

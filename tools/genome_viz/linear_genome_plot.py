from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicRecord
import matplotlib.pyplt as plt
from itertools import cycle
import numpy as np
import argparse

class CPTTranslator(BiopythonTranslator):
    """ 
    This is a customized translator from the dna_features_viewer module to fit Galaxy
    """

    global ignored_features_types
    global chop_block
    global label_fields
    global custom_name_colors
    global custom_feature_colors

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                color_specific = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in custom_name_colors.keys())
                if color_specific:
                    print(custom_name_colors[feature.qualifiers["product"][0]])
                    return custom_name_colors[feature.qualifiers["product"][0]]
                else:
                    return custom_feature_colors[feature.type]

    def compute_feature_label(self, feature): # remove the chop_blocks
        self.label_fields = label_fields
        if feature.type == "CDS":
            if "product" in feature.qualifiers:
                verify_chops = any(re.search(("(\\b"+str(item)+"\\b)"),feature.qualifiers["product"][0]) for item in chop_block)
                if verify_chops:
                    return None
                else:
                    return BiopythonTranslator.compute_feature_label(self, feature)

    def compute_filtered_features(self, features):
        return [
            feature for feature in features if feature.type not in ignored_features_types
        ]
        
    
    def compute_feature_legend_text(self, feature):
        return feature.type
    
    def compute_feature_box_color(self, feature):
        if feature.type == "CDS":
            return "white"
    
    def compute_featurebox_linewidth(self, feature):
        return 0


if __name__ == "__main__":
    pass

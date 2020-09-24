
from dna_features_viewer import BiopythonTranslator

class MyCustomTranslator(BiopythonTranslator):
    """ Custom translator """

    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return "blue"
        elif feature.type == "terminator":
            return "green"
        else:
            return "gold"
    
    def compute_feature_label(self, feature):
        if feature.type == 'restriction_site':
            return None
        elif feature.type == "CDS":
            return "CDS location"
        else:
            return BiopythonTranslator.compute_feature_label(self, feature)
    
    def compute_filtered_features(self, features):
        return [ 
            feature for feature in features
            if (feature.type != "restriction_site")
            or ("BamHI" in str(feature.qualifiers.get("label",'')))
        ]

if __name__ == "__main__":
    graphic_record = MyCustomTranslator().translate_record("test-data/miro.gb")
    ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=10)
    ax.figure.tight_layout()
    ax.figure.savefig("test-data/custom_plot.svg")
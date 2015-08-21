import copy

def feature_lambda(feature_list, test, test_kwargs, subfeatures=True):
    # Either the top level set of [features] or the subfeature attribute
    for feature in feature_list:
        if test(feature, **test_kwargs):
            if not subfeatures:
                feature_copy = copy.deepcopy(feature)
                feature_copy.sub_features = []
                yield feature_copy
            else:
                yield feature

        if hasattr(feature, 'sub_features'):
            for x in feature_lambda(feature.sub_features, test, test_kwargs, subfeatures=subfeatures):
                yield x

def feature_test_type(feature, **kwargs):
    return feature.type == kwargs['type']

def feature_test_quals(feature, **kwargs):
    for attribute_value in feature.qualifiers.get(kwargs['qualifier'], []):
        if attribute_value in kwargs['attribute_list']:
            return True
    return False

def get_id(feature=None, parent_prefix=None):
    result = ""
    if parent_prefix is not None:
        result += parent_prefix + '|'
    if 'locus_tag' in feature.qualifiers:
        result += feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        result += feature.qualifiers['gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    else:
        result += '%s_%s_%s' % (feature.location.start, feature.location.end,
                                feature.location.strand)
    return result


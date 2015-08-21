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

def feature_test_qual_value(feature, **kwargs):
    for attribute_value in feature.qualifiers.get(kwargs['qualifier'], []):
        if attribute_value in kwargs['attribute_list']:
            return True
    return False

def feature_test_quals(feature, **kwargs):
    for key in kwargs:
        if key not in feature.qualifiers:
            return False

        for value in kwargs[key]:
            if value not in feature.qualifiers[key]:
                return False
    return True

def feature_test_contains(feature, **kwargs):
    if 'index' in kwargs:
        return feature.location.start < kwargs['index'] < feature.location.end
    elif 'range' in kwargs:
        return feature.location.start < kwargs['range']['start'] < feature.location.end and \
                feature.location.start < kwargs['range']['end'] < feature.location.end
    else:
        raise RuntimeError('Must use index or range keyword')

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


def ensure_location_in_bounds(start=0, end=0, parent_length=0):
    # This prevents frameshift errors
    while start < 0:
        start += 3
    while end < 0:
        end += 3
    while start > parent_length:
        start -= 3
    while end > parent_length:
        end -= 3
    return (start, end)

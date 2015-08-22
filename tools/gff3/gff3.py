import copy

def feature_lambda(feature_list, test, test_kwargs, subfeatures=True):
    """Recursively search through features, testing each with a test function, yielding matches.

    GFF3 is a hierachical data structure, so we need to be able to recursively
    search through features. E.g. if you're looking for a feature with
    ID='bob.42', you can't just do a simple list comprehension with a test
    case. You don't know how deeply burried bob.42 will be in the feature tree. This is where feature_lambda steps in.

    :type feature_list: list
    :param feature_list: an iterable of features

    :type test: function reference
    :param test: a closure with the method signature (feature, **kwargs) where
                 the kwargs are those passed in the next argument. This
                 function should return True or False, True if the feature is
                 to be yielded as part of the main feature_lambda function, or
                 False if it is to be ignored. This function CAN mutate the
                 features passed to it (think "apply").

    :type test_kwargs: dictionary
    :param test_kwargs: kwargs to pass to your closure when it is called.

    :type subfeatures: boolean
    :param subfeatures: when a feature is matched, should just that feature be
                        yielded to the caller, or should the entire sub_feature
                        tree for that feature be included? subfeatures=True is
                        useful in cases such as searching for a gene feature,
                        and wanting to know what RBS/Shine_Dalgarno_sequences
                        are in the sub_feature tree (which can be accomplished
                        with two feature_lambda calls). subfeatures=False is
                        useful in cases when you want to process (and possibly
                        return) the entire feature tree, such as applying a
                        qualifier to every single feature.

    :rtype: yielded list
    :return: Yields a list of matching features.
    """
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
    """Test qualifier values.

    For every feature, check that at least one value in
    feature.quailfiers(kwargs['qualifier']) is in kwargs['attribute_list']
    """
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

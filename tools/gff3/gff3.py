import copy
import logging
logging.basicConfig(level=logging.WARN)
log = logging.getLogger()


def feature_lambda(feature_list, test, test_kwargs, subfeatures=True, parent=None):
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
        feature._parent = parent
        if test(feature, **test_kwargs):
            if not subfeatures:
                feature_copy = copy.deepcopy(feature)
                feature_copy.sub_features = []
                yield feature_copy
            else:
                yield feature

        if hasattr(feature, 'sub_features'):
            for x in feature_lambda(feature.sub_features, test, test_kwargs, subfeatures=subfeatures, parent=feature):
                yield x


def fetchParent(feature):
    if not hasattr(feature, '_parent') or feature._parent is None:
        return feature
    else:
        return fetchParent(feature._parent)


def feature_test_true(feature, **kwargs):
    return True


def feature_test_type(feature, **kwargs):
    if 'type' in kwargs:
        return feature.type == kwargs['type']
    elif 'types' in kwargs:
        return feature.type in kwargs['types']
    raise Exception("Incorrect feature_test_type call, need type or types")


def feature_test_qual_value(feature, **kwargs):
    """Test qualifier values.

    For every feature, check that at least one value in
    feature.quailfiers(kwargs['qualifier']) is in kwargs['attribute_list']
    """
    for attribute_value in feature.qualifiers.get(kwargs['qualifier'], []):
        if attribute_value in kwargs['attribute_list']:
            return True
    return False


def feature_test_location(feature, **kwargs):
    if 'strand' in kwargs:
        if feature.location.strand != kwargs['strand']:
            return False

    return feature.location.start <= kwargs['loc'] <= feature.location.end


def feature_test_quals(feature, **kwargs):
    """
    Example::

        a = Feature(qualifiers={'Note': ['Some notes', 'Aasdf']})

        # Check if a contains a Note
        feature_test_quals(a, {'Note': None})  # Returns True
        feature_test_quals(a, {'Product': None})  # Returns False

        # Check if a contains a note with specific value
        feature_test_quals(a, {'Note': ['ome']})  # Returns True

        # Check if a contains a note with specific value
        feature_test_quals(a, {'Note': ['other']})  # Returns False
    """
    for key in kwargs:
        if key not in feature.qualifiers:
            return False

        # Key is present, no value specified
        if kwargs[key] is None:
            return True

        # Otherwise there is a key value we're looking for.
        # so we make a list of matches
        matches = []
        # And check all of the feature qualifier valuse
        for value in feature.qualifiers[key]:
            # For that kwargs[key] value
            matches.append(kwargs[key] in value)

        # If none matched, then we return false.
        if not any(matches):
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
    elif 'Gene' in feature.qualifiers:
        result += feature.qualifiers['Gene'][0]
    elif 'product' in feature.qualifiers:
        result += feature.qualifiers['product'][0]
    elif 'Product' in feature.qualifiers:
        result += feature.qualifiers['Product'][0]
    elif 'Name' in feature.qualifiers:
        result += feature.qualifiers['Name'][0]
    else:
        return feature.id
        # Leaving in case bad things happen.
        # result += '%s_%s_%s_%s' % (
            # feature.id,
            # feature.location.start,
            # feature.location.end,
            # feature.location.strand
        # )
    return result


def get_gff3_id(gene):
    return gene.qualifiers.get('Name', [gene.id])[0]


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


def coding_genes(feature_list):
    for x in feature_lambda(feature_list, feature_test_type, {'type': 'gene'}, subfeatures=True):
        if len(list(feature_lambda(x.sub_features, feature_test_type, {'type': 'CDS'}, subfeatures=False))) > 0:
            yield x


def genes(feature_list, feature_type='gene'):
    """
    Simple filter to extract gene features from the feature set.
    """

    for x in feature_lambda(feature_list, feature_test_type,
                            {'type': feature_type},
                            subfeatures=True):
        yield x

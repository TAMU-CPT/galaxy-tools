def bigwig_add_header(bw_handle, identifier, name='SignalP'):
    bw_handle.write("track type=wiggle_0 name=%s-%s visibility=full\n" % (name, identifier))


def bigwig_store(bw_handle, chrom, data):
    bw_handle.write("variableStep chrom=%s span=1\n" % chrom)
    for position, value in enumerate(data):
        bw_handle.write('%s %.3f\n' % (position + 1, value))

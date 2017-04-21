import sys

# Like 'import json' / 'import yaml', except.. tab data.
def loads(str_data):
    return NotImplementedError()

def load(handle):
    return NotImplementedError()

def dump(data, handle=sys.stdout):
    for row in data:
        handle.write('%s\n' % '\t'.join(
            map(str, row)
        ))

def dumps(data):
    output = ""
    for row in data:
        output += '%s\n' % '\t'.join(
            map(str, row)
        )
    return output

def dump_line(row, handle=sys.stdout):
    dump([row], handle=handle)

def dumps_line(row):
    return dumps([row])

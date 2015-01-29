import argparse
import yaml
import os
import subprocess
from Bio import SeqIO
import copy
import hashlib
from Bio.Blast.Applications import NcbiblastnCommandline

BASE_PATH = os.path.join(os.sep, 'tmp', 'cpt.ani.blast')
GENOMES = os.path.join(BASE_PATH, 'genome')
CHUNKS = os.path.join(BASE_PATH, 'chunks')
BLASTDB = os.path.join(BASE_PATH, 'blastdb')
SUFFIX = '_nucl'

def serialize_sequences(fasta_file):
    """Store a fasta file with multiple sequences to a set of files with single sequences

    Sometimes users will pass files with 10, 100 sequences. We need to treat these as individual genomes
    """
    records = SeqIO.parse(fasta_file, 'fasta')
    data = []
    for record in records:
        seqhash = hashlib.md5(str(record.seq)).hexdigest()
        file_path = os.path.join(GENOMES, seqhash + SUFFIX +'.fa')

        # Only write out if new file
        if not os.path.exists(file_path):
            with open(file_path, 'w') as handle:
                handle.write('>%s\n%s\n' % (seqhash, record.seq))

        data.append((
                record.id,
                seqhash,
                file_path,
                len(record.seq),
            ))
    return data

def chunky_range(length, window_size, step_size):
    """Convenience function to get a set of chunks for a range

    Here's an example

        >>> chunky_range(1000, 100, 100)
        (10, [0, 100, 200, 300, 400, 500, 600, 700, 800, 900])

    TODO: refactor into just returning the list, and grabbing num_chunks from len()
    """
    nearest_multiple = int(length/step_size)*step_size
    num_chunks = nearest_multiple / step_size
    return (num_chunks, range(0, nearest_multiple, step_size))

def chunky(fasta_path, parent_id, window_size=1000, step_size=500):
    """Given a fasta file, chunk it and store those chunks to a "_chunks.fa" file.

    Returns the path to the file, the number of chunks, and the list of chunks.
    """

    with open(fasta_path, 'r') as handle:
        record = list(SeqIO.parse(handle, 'fasta'))[0]
        length = len(record.seq)

        chunk_path = os.path.join(CHUNKS, '%s_%s_%s_%s.fa' % (parent_id, window_size, step_size, SUFFIX))

        (num_chunks, chunk_range) = chunky_range(length, window_size, step_size)
        if not os.path.exists(chunk_path):
            with open(chunk_path, 'w') as handle:
                for i in chunk_range:
                    handle.write('>%s_%s\n%s\n' % (
                        parent_id, i, str(record.seq[i:i+window_size])
                        ))
        return (chunk_path, num_chunks, chunk_range)

def makeblastdb(files):
    """Create a blastdb from a set of files.
    """
    # Sort
    files = sorted(files)
    # Not relocatable.
    dbname = hashlib.md5(''.join(files)).hexdigest()
    merged_file_location = os.path.join(BLASTDB, dbname + SUFFIX)

    # Concatenate all files
    if not os.path.exists(merged_file_location + '.fa'):
        with open(merged_file_location + '.fa', 'w') as handle:
            for chunk_file in files:
                handle.write(open(chunk_file).read())

    # Create the blastdb
    if not os.path.exists(merged_file_location + '.nin'):
        command = ['makeblastdb', '-dbtype', 'nucl', '-out',
                merged_file_location, '-in', merged_file_location + '.fa']
        output = subprocess.check_output(command)
    return merged_file_location

def blastn(query, query_id, db):
    """Blast a query against a database

        >>> blastn_output_path = blastn('/path/to/query.fa', 'query-id', '/path/to/db')
    """
    path = db + query_id + SUFFIX + '.tsv'
    if not os.path.exists(path):
        #-outfmt 6 -word_size 11 -evalue 10 -gapopen 5 -gapextend 2 -reward 2 -penalty -3
        blastn_cline = NcbiblastnCommandline(query=query, db=db, evalue="1",
                                             outfmt=6, out=path, word_size=11,
                                             gapopen=5, gapextend=2, reward=2,
                                             penalty=-3)
        (stdout, stderr) = blastn_cline()
    return path

def ani_analysis(genomes, results_location, window_size=1000, ani_yaml=None):
    """Run analysis of an ANI internal data structure.

    The `genomes` structure looks like this:

        {'0c55417ff7136f14e693de66105ff9c6': {
        'chunks': '/tmp/cpt.ani.blast/chunks/0c55417ff7136f14e693de66105ff9c6_1000_500.fa',
        'hash': '0c55417ff7136f14e693de66105ff9c6',
        'id': '01',
        'path': '/tmp/cpt.ani.blast/genome/0c55417ff7136f14e693de66105ff9c6.fa'},

    The results_location refers to the output of BLASTN. Window size must be passed for doing calculations.
    """

    # Set up dict
    blast_data = {}
    for genome in genomes:
        for genome2 in genomes:
            if genome not in blast_data:
                blast_data[genome] = {}
                blast_data[genome]['_meta'] = copy.deepcopy(genomes[genome])
                # Useful
                blast_data[genome]['_meta']['chunk_size'] = window_size
                # Useless
                if 'chunk_range' in blast_data[genome]['_meta']:
                    del blast_data[genome]['_meta']['chunk_range']
            if genome2 not in blast_data[genome]:
                blast_data[genome][genome2] = {}

    # Clone dict keys
    results_map = copy.deepcopy(blast_data)

    with open(results_location, 'r') as handle:
        for line in handle.readlines():
            split_line = line.split('\t')
            (qseqid, sseqid) = split_line[0:2]

            (pident, length, mismatch, gapopen, qstart, qend,
                    sstart, send, evalue, bitscore) = map(float, split_line[2:])

            (qseqid_hash, qseqid_idx) = qseqid.split('_')
            (sseqid_hash, sseqid_idx) = sseqid.split('_')

            if qseqid_idx not in blast_data[qseqid_hash][sseqid_hash]:
                blast_data[qseqid_hash][sseqid_hash][qseqid_idx] = {}

            dice = (pident * length) / window_size

            if sseqid_idx not in blast_data[qseqid_hash][sseqid_hash][qseqid_idx]:
                blast_data[qseqid_hash][sseqid_hash][qseqid_idx][sseqid_idx] = dice
            else:
                if dice > blast_data[qseqid_hash][sseqid_hash][qseqid_idx][sseqid_idx]:
                    blast_data[qseqid_hash][sseqid_hash][qseqid_idx][sseqid_idx] = dice

    with open(ani_yaml, 'w') as handle:
        yaml.dump(blast_data, handle)

    for genome_query in genomes:
        for genome_subject in genomes:
            qs_results = []
            # For each region of the genome, we want to ensure that in the target genome there is at least one high quality match.
            for idx_query in map(str, genomes[genome_query]['chunk_range']):
                if idx_query in blast_data[genome_query][genome_subject]:
                    best_hit = 0
                    for x in blast_data[genome_query][genome_subject][idx_query]:
                        if blast_data[genome_query][genome_subject][idx_query][x] > best_hit:
                            best_hit = blast_data[genome_query][genome_subject][idx_query][x]
                    qs_results.append(best_hit)
                else:
                    qs_results.append(0)
            results_map[genome_query][genome_subject] = sum(qs_results)/len(qs_results)

    return results_map

def format_results(genomes, results_data):
    """Convenience function to turn list of lists into a TSV
    """
    id_map = {}
    for genome in genomes:
        id_map[genomes[genome]['hash']] = genomes[genome]['id']

    reference_keys = sorted(genomes.keys())

    tabular = [ '\t'.join(['From\To'] + [id_map[k] for k in reference_keys]) ]

    for key in reference_keys:
        row = [id_map[key]]
        for key2 in reference_keys:
            row.append('%0.2f' % results_data[key][key2])
        tabular.append('\t'.join(row))
    return '\n'.join(tabular)

def ani(fasta_files, window_size, step_size, *args, **kwd):
    """Main ANI method.

    Given a list of fasta files, a window and step size (please ensure that
    window_size==step_size), calculate average nucleotide identity.
    """
    complete_genomes = {}
    for fasta in fasta_files:
        output_files = serialize_sequences(fasta)
        for sequence in output_files:
            complete_genomes[sequence[1]] = {
                    'id': sequence[0],
                    'path': sequence[2],
                    'hash': sequence[1],
                    'size': sequence[3],
                    }

    for genome in complete_genomes:
        (chunks, chunk_count, chunk_range) = \
            chunky(complete_genomes[genome]['path'],
                    complete_genomes[genome]['hash'],
                    window_size=window_size,
                    step_size=step_size)

        complete_genomes[genome]['chunks'] = chunks
        complete_genomes[genome]['chunk_count'] = chunk_count
        complete_genomes[genome]['chunk_range'] = chunk_range

    chunked_files = [complete_genomes[x]['chunks'] for x in complete_genomes]
    blastdb = makeblastdb(chunked_files)
    results_location = blastn(blastdb + '.fa', 'self', blastdb)
    ani_results = ani_analysis(complete_genomes, results_location, window_size=window_size, ani_yaml=kwd['ani_yaml_out'])

    results = format_results(complete_genomes, ani_results)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Average Nucleotide Identity')
    parser.add_argument('fasta_files', type=file, nargs='+', help='input fasta genomes')
    parser.add_argument('--ani_output', help='Output')
    parser.add_argument('--ani_yaml_out', help='ANI YAML Output')

    parser.add_argument('--window_size', type=int, default=1000, help='window size for splitting the genome')


    # Ensure paths exist
    for path in (GENOMES, CHUNKS, BLASTDB):
        if not os.path.exists(path):
            os.makedirs(path)

    parsed_args = parser.parse_args()
    args = vars(parsed_args)
    # Deprecate step_size option and ensure step_size==window_size
    args['step_size'] = args['window_size']
    # Execute
    results = ani(**args)
    # STDOUT

    with open(parsed_args.ani_output, 'w') as output_file:
        output_file.write(results)

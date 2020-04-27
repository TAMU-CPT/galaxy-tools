from cpt import OrfFinder
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio import Seq

data = "debug/cody-genome.fa"

#isps = OrfFinder(args.table, args.ftype, args.ends, args.isp_min_len, args.strand)
aadebug = open("aa-debug.fa","w")
ntdebug = open("nt-debug.fa","w")
bedbug = open("debug.bed","w")
gffbug = open("debug.gff3","w")
strand = OrfFinder(11,"CDS","closed",60,"both")
#table = CodonTable.ambiguous_generic_by_id[11]
table = CodonTable.unambiguous_dna_by_id[11]
print(strand)
print(table)
#strand.locate(data,ntdebug,aadebug,bedbug,gffbug)
#isps.locate(args.fasta_file,args.out_isp_nuc,args.out_isp_prot,args.out_isp_bed,args.out_isp_gff3,)

print(table.start_codons)
print(table.stop_codons)
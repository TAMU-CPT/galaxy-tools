#!/usr/bin/env python
import os
import re
import json
import glob
import shutil
import argparse
import subprocess
import numpy as np
from enum import Enum
from Bio import SeqIO
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from cpt_convert_mga_to_gff3 import mga_to_gff3
from gff3 import feature_lambda
from safe_reopen import extract_gff3_regions, gaps
from jinja2 import Environment, FileSystemLoader
import logging

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(name="autoreopen")

# TODO: This tool depends on PY2K ONLY TOOLS. THIS WILL(MAY?) NOT FUNCTINO UNDER PY3.
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "reopen-data")
env = Environment(
    loader=FileSystemLoader(SCRIPT_DIR), trim_blocks=True, lstrip_blocks=True
)


class AutoreopenException(Exception):
    """Class to base exceptions off of. This allows us to catch normal program
    exceptions without catching worse things as well."""

    pass


class Evidence(Enum):
    NONE = 0
    Assumption = 1
    BLAST = 2
    PhageTerm = 3


class PhageType(Enum):
    ShortTR = 1
    LongTR = 2
    Prime3Cos = 3
    Prime5Cos = 4
    PacHeadful = 6
    T4Headful = 7
    MuLike = 8
    Unknown = 9


class BlastBasedReopening(object):
    """Re-open genomes based on blast"""

    def __init__(
        self,
        db="test-data/canonical_nucl",
        genome=None,
        genome_real=None,
        protein=False,
        data_dir=None,
        reopen_blastN=None,
    ):
        self.reopen_blastN = reopen_blastN
        self.data_dir = data_dir
        self.genome = genome
        self.genome_nucleotide_real = genome_real
        self.protein = protein
        self.db = db  # 'test-data/canonical_nucl'
        self.canonical_map = {
            "NC_025442.1": "Chi",
            "NC_005880.2": "Phage K",
            "NC_008720.1": "N4",
            "NC_004775.1": "epsilon15",
            "NC_005282.1": "Felix01",
            "NC_001416.1": "lambda",
            "NC_000929.1": "Mu",
            "NC_001901.1": "N15",
            "NC_005856.1": "P1",
            "NC_002371.2": "P22",
            "NC_001895.1": "P2",
            "NC_001609.1": "P4",
            "NC_011048.1": "phi29",
            "NC_005045.1": "phiKMV",
            "NC_004831.2": "SP6",
            "NC_011421.1": "SPO1",
            "NC_004166.2": "SPP1",
            "NC_005833.1": "T1",
            "NC_003298.1": "T3",
            "NC_000866.4": "T4",
            "NC_005859.1": "T5",
            "NC_001604.1": "T7",
        }

    def blast2Xmfa(self, blast_results, name="nucl"):
        tpl = env.get_template("autoreopen_template_plot.html")
        with open(os.path.join(self.data_dir, "%s.html" % name), "w") as html:
            html.write(tpl.render(filename="%s.json" % name).encode("utf-8"))

        with open(os.path.join(self.data_dir, "%s.json" % name), "w") as handle:
            json.dump(
                {
                    "xmfa": "%s_regions.json" % name,
                    "fasta": [
                        {"path": "n.a", "name": "SubjectGenome", "length": 200000},
                        {"path": "n.a", "name": "QueryGenome", "length": 200000},
                    ],
                    "gff3": ["n.a", "n.a"],
                },
                handle,
                sort_keys=True,
                indent=2,
            )

        with open(os.path.join(self.data_dir, "%s_regions.json" % name), "w") as handle:
            region_data = []
            for hit in blast_results:
                region_data.append(
                    [
                        {
                            "rid": "1:fake",
                            "comment": "",
                            "start": hit["qstart"],
                            "seq": "",
                            "id": "1",
                            "strand": 1 if hit["sstart"] > hit["send"] else -1,
                            "end": hit["send"],
                        },
                        {
                            "rid": "2:fake",
                            "comment": "",
                            "start": hit["qstart"],
                            "seq": "",
                            "id": "2",
                            "strand": 1 if hit["qstart"] > hit["qend"] else -1,
                            "end": hit["send"],
                        },
                    ]
                )

            json.dump(region_data, handle, sort_keys=True, indent=2)

    def evaluate(self):
        if self.protein:
            blast_results = self.getBlastP(self.genome)
            self.blast2Xmfa(blast_results, "prot")
        else:
            blast_results = self.getBlastN(self.genome)
            self.blast2Xmfa(blast_results, "nucl")

        if len(blast_results) == 0:
            return {
                "location": None,
                "orientation": None,
                "results": [],
                "reopened_file": None,
            }

        orientation = "+"
        current_genome_file = None
        if not self.protein:
            # If not a protein file, we check the orientation. We probably
            # should do this for protein too, but ... yikes.
            genome_seq = SeqIO.read(self.genome, "fasta")
            current_genome_file = self.genome

            # If more than half the hits are 'backwards'
            if self.shouldReverse(blast_results) > 0.5:
                orientation = "-"
                new_genome_file = os.path.join(
                    self.data_dir,
                    os.path.basename(self.genome).replace(".fa", ".revcom.fa"),
                )
                rec = genome_seq.reverse_complement(id=True, description=True)
                rec.description += " autoreopen.revcom"
                SeqIO.write([rec], new_genome_file, "fasta")
                current_genome_file = new_genome_file
                blast_results = self.getBlastN(new_genome_file)

        # TODO: Filter out just ONE genome + best results
        # blast_results = blast_results

        # ref_genome_hits = [self.avg(hit['sstart'], hit['send']) for hit in blast_results]
        # our_genome_hits = [self.avg(hit['qstart'], hit['qend']) for hit in blast_results]
        # Now we sort relative to ref genome
        # hits = zip(our_genome_hits, ref_genome_hits)
        # hits = sorted(hits, key=lambda x: x[1])
        # Then we re-construct *_genome_hits from the sorted/zipped version
        # ref_genome_hits = [x[1] for x in hits]
        # our_genome_hits = [x[0] for x in hits]
        # These are now two independent lists, sorted relative to ref genome. We
        # can shift our_genome_hits without affecting ref genome ordering.
        # (score, (index, m, c)) = self.scoreResults(ref_genome_hits, our_genome_hits)

        # Now on to plotting the real data.
        totalData = [
            (hit["qstart"], hit["qend"], hit["sstart"], hit["send"])
            for hit in blast_results
        ]
        totalData = sorted(totalData, key=lambda x: self.avg(x[2], x[3]))

        if self.protein:
            # In order to use this we need to re-open according to our new end,
            # and then re-call genes / export to pfa / re-blast. This is SUPER
            # ugly and SUPER hacky.
            new_genome_file_reopened = self.genome_nucleotide_real
            if self.genome_nucleotide_real:
                customPhageReopener = PhageReopener(
                    self.genome_nucleotide_real, None, None, data_dir=self.data_dir
                )
                protein_fasta = customPhageReopener._orfCalls()
                # This gets a new protein_fasta file, which we can re-blast and re-analyse
                updated_blast_results = self.getBlastP(protein_fasta)
                self.blast2Xmfa(updated_blast_results, "prot2")
        else:
            genome_seq = SeqIO.read(current_genome_file, "fasta")
            new_genome_file_reopened = os.path.join(
                self.data_dir,
                os.path.basename(self.genome).replace(".fa", ".reopened.fa"),
            )
            with open(new_genome_file_reopened, "w") as handle:
                loc = totalData[0][0]
                rec = genome_seq[loc:] + genome_seq[0:loc]
                rec.description += " autoreopen.reopen(%s)" % loc
                SeqIO.write([rec], handle, "fasta")
                # copy the first time only.
                if not os.path.exists(self.reopen_blastN):
                    shutil.copy(
                        os.path.join(new_genome_file_reopened),
                        os.path.join(self.reopen_blastN),
                    )

            updated_blast_results = self.getBlastN(new_genome_file_reopened)
            self.blast2Xmfa(updated_blast_results, "nucl2")

        return {
            "location": totalData[0],
            "orientation": orientation,
            "results": blast_results,
            "reopened_file": new_genome_file_reopened,
        }

    def getBlastN(self, path):
        blast_results = (
            subprocess.check_output(
                [
                    "blastn",
                    "-db",
                    self.db,
                    "-num_threads",
                    "4",
                    "-query",
                    path,
                    "-outfmt",
                    "6 sseqid evalue pident qstart qend sstart send",
                ]
            )
            .strip()
            .split("\n")
        )
        blast_results = [x.split("\t") for x in blast_results if len(x.strip()) > 0]
        if len(blast_results) == 0:
            log.warning("Empty blastN")
            return []

        blast_results = [
            {
                "sseqid": self.canonical_map[sseqid],
                "evalue": float(evalue),
                "pident": float(pident),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
            }
            for (sseqid, evalue, pident, qstart, qend, sstart, send) in blast_results
            if float(evalue) < 0.001 and float(pident) > 80
        ]

        sseqids = [x["sseqid"] for x in blast_results]
        if len(set(sseqids)) > 1:
            log.warning("DO not know how to handle. This WILL MISBEHAVE")
        return blast_results

    def getBlastP(self, path):
        cmd = [
            "blastp",
            "-db",
            self.db,
            "-num_threads",
            "4",
            "-query",
            path,
            "-outfmt",
            "6 qseqid sseqid evalue pident qstart qend sstart send",
        ]
        blast_results = subprocess.check_output(cmd).strip().split("\n")
        blast_results = [x for x in blast_results if len(x.strip()) > 0]

        if len(blast_results) == 0:
            log.warning("Empty blastP")
            return []

        blastp_locmap = (
            open(os.path.join(SCRIPT_DIR, "test-data", "merged.pfa.idmap"), "r")
            .read()
            .strip()
            .split("\n")
        )

        def extremes(data):
            mi = None
            ma = None

            for d in data:
                d = int(d)
                if mi is None:
                    mi = d
                if ma is None:
                    ma = d

                ma = max(ma, d)
                mi = min(mi, d)
            return mi, ma

        blastp_locmap2 = {}
        for x in blastp_locmap:
            (key, value) = x.split("\t")
            (mi, ma) = extremes([y for y in re.split("[^0-9]+", value) if len(y) > 0])
            if "complement" in value:
                blastp_locmap2[key] = (ma, mi)
            else:
                blastp_locmap2[key] = (mi, ma)

        input_locmap = {}
        for tmprec in SeqIO.parse(path, "fasta"):
            desc = tmprec.description
            locstr = desc[desc.index("Location=") + 10 : -2]
            locdata = re.split("[^0-9+-]+", locstr)

            if locdata[2] == "+":
                input_locmap[tmprec.id] = (int(locdata[0]), int(locdata[1]))
            else:
                input_locmap[tmprec.id] = (int(locdata[1]), int(locdata[0]))

        blast_results = [x.split("\t") for x in blast_results]
        blast_results = [
            {
                "qseqid": qseqid[qseqid.rindex("_") + 1 :],
                "sseqid": self.canonical_map[sseqid[4 : sseqid.index("_prot_")]],
                "evalue": float(evalue),
                "pident": float(pident),
                "qstart": input_locmap[qseqid][0],
                "qend": input_locmap[qseqid][1],
                "sstart": blastp_locmap2[sseqid][0],
                "send": blastp_locmap2[sseqid][1],
            }
            for (
                qseqid,
                sseqid,
                evalue,
                pident,
                qstart,
                qend,
                sstart,
                send,
            ) in blast_results
            if float(evalue) < 0.001 and float(pident) > 70
        ]

        sseqid_set = set([x["sseqid"] for x in blast_results])
        if len(sseqid_set) > 1:
            log.warning(
                "DO not know how to handle. This WILL MISBEHAVE. "
                + ",".join(sseqid_set)
            )
        return blast_results

    def shouldReverse(self, blast_results):
        should_reverse = 0
        total = len(blast_results)

        for hit in blast_results:
            qdir = 1 if hit["qstart"] > hit["qend"] else -1
            sdir = 1 if hit["sstart"] > hit["send"] else -1
            if (qdir > 0 and sdir < 0) or (qdir < 0 and sdir > 0):
                should_reverse += 1

        return float(should_reverse) / (total + 1)

    @classmethod
    def rotate(cls, data, index):
        return data[index:] + data[0:index]

    @classmethod
    def avg(cls, a, b):
        return int(float(a + b) / 2)

    def scoreResults(self, ref_genome_hits, our_genome_hits):

        best_score = None
        best_params = None

        # This no longer needs to be re-sorted.
        for i in range(len(our_genome_hits)):
            hits = zip(self.rotate(our_genome_hits, i), ref_genome_hits)
            x_position = np.array([hit[0] for hit in hits])
            y_position = np.array([hit[1] for hit in hits])
            (m, cx), r2, c, d, e = np.polyfit(x_position, y_position, 1, full=True)
            r2 = r2[0]

            # Now we look to minimize R**2. Could also optimize for m == 1?
            if best_score is None:
                best_score = r2
                best_params = (i, m, cx)

            if r2 < best_score:
                best_score = r2
                best_params = (i, m, cx)

        return best_score, best_params


class PhageReopener:
    def __init__(
        self,
        fasta,
        fastq1,
        fastq2,
        closed=False,
        data_dir=None,
        html=None,
        reopen_phageTerm=None,
        reopen_TerL=None,
        reopen_blastN=None,
    ):
        # fasta, fastq1, fastq2 are all file handles
        self.base_name = os.path.join(
            data_dir, os.path.basename(fasta).replace(".fa", "")
        )
        self.html = html
        self.rec_file = open(fasta, "r")
        self.data_dir = data_dir
        self.fasta = SeqIO.read(fasta, "fasta")
        self.fq1 = fastq1
        self.fq2 = fastq2
        self.closed = closed

        self.reopen_phageTerm = reopen_phageTerm
        self.reopen_TerL = reopen_TerL
        self.reopen_blastN = reopen_blastN

    def _orfCalls(self):
        fnmga = self.base_name + ".mga"
        if not os.path.exists(fnmga):
            log.warn(
                "%s does not exist, calling genes in %s", fnmga, self.rec_file.name
            )
            # Run MGA
            subprocess.check_call(
                ["mga_linux_x64", "-s", self.rec_file.name], stdout=open(fnmga, "w")
            )

        # Convert to gff3
        fn = self.base_name + ".mga.gff3"
        self.mga_gff3 = fn
        with open(fnmga, "r") as handle, open(fn, "w") as output:
            self.rec_file.seek(0)
            for result in mga_to_gff3(handle, self.rec_file):
                # Store gFF3 data in self in order to access later.
                self.mga_rec = result
                GFF.write([result], output)

        # Process a feature id -> feature table in mem.
        self.featureDict = {}
        for f in feature_lambda(
            self.mga_rec.features, lambda x: True, {}, subfeatures=True
        ):
            self.featureDict[f.qualifiers["ID"][0]] = f

        # Extract
        fnfa = self.base_name + ".mga.fa"
        self.fnfa = fnfa
        subprocess.check_call(
            [
                "python2",
                os.path.join(SCRIPT_DIR, os.pardir, "gff3", "gff3_extract_sequence.py"),
                "--feature_filter",
                "CDS",
                self.rec_file.name,
                fn,
            ],
            stdout=open(fnfa, "w"),
        )

        # Translate
        fnpfa = self.base_name + ".mga.pfa"
        self.fnpfa = fnpfa
        subprocess.check_call(
            [
                "python2",
                os.path.join(SCRIPT_DIR, os.pardir, "fasta", "fasta_translate.py"),
                "--table",
                "11",
                "--strip_stops",
                "--target",
                "protein",
                fnfa,
            ],
            stdout=open(fnpfa, "w"),
        )
        return fnpfa

    def _runPhageTerm(self):
        fn = self.base_name + ".json"
        if not os.path.exists(fn):
            subprocess.check_call(
                [
                    "python2",
                    os.path.join(
                        SCRIPT_DIR, os.pardir, "external", "phageterm", "PhageTerm.py"
                    ),
                    "-c",
                    "6",
                    "-f",
                    self.fq1.name,
                    "-n",
                    self.fasta.id,
                    "-r",
                    self.rec_file.name,
                    "-p",
                    self.fq2.name,
                    "-s",
                    "30",
                ]
            )
            for f in glob.glob(self.fasta.id + "*"):
                shutil.move(f, self.data_dir)

            shutil.move(
                os.path.join(self.data_dir, self.fasta.id + ".json"),
                self.base_name + ".json",
            )

            shutil.copy(
                os.path.join(self.data_dir, self.fasta.id + "_PhageTerm_report.pdf"),
                os.path.join(self.data_dir, "report.pdf"),
            )

        # Take the sequence, need to SAFELY re-open it.
        # This requires generating gene calls
        input_seq = os.path.join(self.data_dir, self.fasta.id + "_sequence.fasta")
        tmpfile = os.path.join(self.data_dir, self.fasta.id + "_sequence.mga")

        if os.path.exists(input_seq):
            subprocess.check_call(
                ["mga_linux_x64", "-s", input_seq], stdout=open(tmpfile, "w")
            )
            # Ok, with gene calls done, we can now convert to gff3

            tmpfile2 = os.path.join(self.data_dir, self.fasta.id + "_sequence.mga.gff3")
            with open(tmpfile, "r") as handle, open(tmpfile2, "w") as output:
                for result in mga_to_gff3(handle, open(input_seq, "r")):
                    # Store gFF3 data in self in order to access later.
                    self.mga_rec = result
                    GFF.write([result], output)

            # Now with GFF3 + fasta, we can safe_reopen
            subprocess.check_call(
                [
                    "python2",
                    os.path.join(SCRIPT_DIR, "safe_reopen.py"),
                    input_seq,
                    tmpfile2,
                ],
                stdout=open(self.reopen_phageTerm, "w"),
            )

        with open(fn, "r") as handle:
            return json.load(handle)

    def _analysePhageTerm(self, results):
        phtype = None
        opening = results["reopening"]

        if results["P_class"] == "COS (3')":
            phtype = PhageType.Prime3Cos
        elif results["P_class"] == "COS (5')":
            phtype = PhageType.Prime5Cos
        elif results["P_class"] == "DTR (long)":
            phtype = PhageType.LongTR
        elif results["P_class"] == "DTR (short)":
            phtype = PhageType.ShortTR
        elif results["P_class"] == "Headful (pac)":
            phtype = PhageType.PacHeadful
        elif results["P_class"] == "Mu-like":
            phtype = PhageType.MuLike
        else:
            phtype = PhageType.Unknown

        if phtype == PhageType.Unknown:
            opening = ["Unknown"]

        return (results["phagename"], phtype, opening, Evidence.PhageTerm)

    def _safeOpeningLocationForFeature(self, feature):
        """Given a feature, find a 'safe' location to re-open the genome.

        Safe here means not in the middle of a gene. Thus, we search upstream
        in the genome until we find a place not covered by genes.
        """
        occupied_regions = extract_gff3_regions([self.mga_gff3])
        # Only acceptable places to open a genome.
        gaps_in_data = gaps(occupied_regions[self.fasta.id])

        best = None

        for gap in sorted(gaps_in_data, key=lambda x: x[0]):
            if feature.location.strand > 0:
                # If the end of the gap is upstream
                if gap[1] < feature.location.start:
                    if not best:
                        best = gap

                    if (
                        feature.location.start - gap[1]
                        < feature.location.start - best[1]
                    ):
                        best = gap
            else:
                # if beginning of gap is upstream
                if gap[0] > feature.location.end:
                    if not best:
                        best = gap

                    if gap[0] - feature.location.end < best[0] - feature.location.end:
                        best = gap
        # Ok, we /should/ have a 'best' gap (But I think I'm missing an edge
        # case since this isn't a circular data structure...)
        # Next, get midpoint
        mid = sum(best) / 2
        return mid

    def _upstreamFeatures(self, feature, distance=2000):
        features = sorted(self.mga_rec.features, key=lambda x: x.location.start)
        # 2 kb. I like 2.
        if feature.location.strand > 0:
            wantedFeatures = [
                x
                for x in features
                if feature.location.start - distance
                < x.location.start
                <= feature.location.start
            ]
        else:
            wantedFeatures = [
                x
                for x in features
                if feature.location.end
                < x.location.end
                <= feature.location.end + distance
            ]

        return wantedFeatures

    def _blast(self, pfa):
        blast_results = subprocess.check_output(
            [
                "blastp",
                "-query",
                pfa,
                "-db",
                os.path.join(SCRIPT_DIR, "test-data", "merged"),
                "-outfmt",
                "6 sseqid evalue pident qseqid",
                "-evalue",
                "0.001",
            ]
        ).split("\n")

        def blast_filter(d):
            if len(d.strip()) == 0:
                return False
            q = d.split("\t")
            q[2] = float(q[2])
            return q[2] > 50

        # Filter only good quality hits (>50%)
        blast_results = filter(blast_filter, blast_results)
        # Pull out the type
        blast_result = []
        for hit in blast_results:
            (sseqid, evalue, pident, qseqid) = hit.split("\t")
            nseqid = sseqid
            if "_" in nseqid:
                nseqid = nseqid[: nseqid.index("_")]
            if "|" in nseqid:
                nseqid = nseqid[: nseqid.index("|")]
            blast_result.append((nseqid, qseqid))

        # reduce to set
        blast_result = list(set(blast_result))
        ph_type = PhageType.Unknown
        # Handle results
        if len(blast_result) == 1:
            blast_type, protein_name = blast_result[0]

            if blast_type == "3primeCos":
                ph_type = PhageType.Prime3Cos
            elif blast_type == "5primeCos":
                ph_type = PhageType.Prime5Cos
            elif blast_type == "long-tr":
                ph_type = PhageType.LongTR
            elif blast_type == "short-tr":
                ph_type = PhageType.ShortTR
            elif blast_type == "pac-headful":
                ph_type = PhageType.PacHeadful
            elif blast_type == "t4-headful":
                ph_type = PhageType.T4Headful
            else:
                # log.warning("PhageType %s??", blast_type)
                ph_type = PhageType.Unknown

            f = self.featureDict[protein_name]
            safeLoc = self._safeOpeningLocationForFeature(f)

            # now to re-open there.
            if not os.path.exists(self.reopen_TerL):
                rec = self.fasta[safeLoc:] + self.fasta[0:safeLoc]
                rec.description += " blast.terL.reopen(%s)" % safeLoc
                SeqIO.write([rec], self.reopen_TerL, "fasta")

            return (
                blast_results,
                ph_type,
                ["forward" if f.strand > 0 else "reverse", safeLoc],
                f,
                self._upstreamFeatures(f),
            )
        else:
            log.warning("%s blast results for this protein", len(blast_results))
            return blast_results, ph_type, ["Unknown", "Unknown"], None, None

    def _guessTerS(self, upstreamFeatures):
        # Given a TerL feature and the features upstream of that, let's take a
        # stab at the TerS by interpro scan searching the three genes upstream
        # for a DNA Binding domain
        if not upstreamFeatures:
            return

        tmpfile = os.path.join(self.data_dir, "ters.pfa")
        # Write out our features
        with open(tmpfile, "w") as handle:
            recs = []
            for feature in upstreamFeatures:
                rec = SeqRecord(
                    feature.extract(self.fasta).seq, id=feature.qualifiers["ID"][0]
                )
                rec.seq = rec.seq.translate(table=11)
                recs.append(rec)

            SeqIO.write(recs, handle, "fasta")

    def giveEvidence(self):
        kwargs = {}
        # Naive annotation, we run these NO MATTER WHAT.
        self.protein_fasta = self._orfCalls()
        (
            blast_results,
            blast_type,
            blast_reopen_location,
            blast_feature,
            blast_upstraem,
        ) = self._blast(self.protein_fasta)
        # Cannot guess terS because bad bad
        # terS = self._guessTerS(blast_upstraem)

        # Try their tool
        try:
            results = self._runPhageTerm()
            (
                phage_name,
                phageterm_type,
                phageterm_location,
                phageterm_evidence,
            ) = self._analysePhageTerm(results)
            phageterm = {"type": phageterm_type.name, "location": phageterm_location}
            kwargs["PhageTerm"] = phageterm
            kwargs["Name"] = phage_name
        except Exception:
            phageterm = {"type": "PHAGETERM FAILED", "location": ["Unknown"]}
            kwargs["PhageTerm"] = phageterm
            kwargs["Name"] = "Unknown"

        # Next we failover to blast results of terminases.
        blast = {
            "type": blast_type.name,
            "location": blast_reopen_location,
            "TerL": blast_feature,
            "TerS": None,
            "results": [x.split() for x in blast_results],
        }
        kwargs["BLAST"] = blast

        # Then we do alignment relative to canonical
        bbr_nucl = BlastBasedReopening(
            genome=self.rec_file.name,
            db=os.path.join(SCRIPT_DIR, "test-data/canonical_nucl"),
            protein=False,
            data_dir=self.data_dir,
            reopen_blastN=self.reopen_blastN,
        )
        kwargs["canonical_nucl"] = bbr_nucl.evaluate()
        blastp_genome = kwargs["canonical_nucl"]["reopened_file"]

        # If no location because no blastN hits.
        if kwargs["canonical_nucl"]["location"]:
            kwargs["canonical_nucl"]["location"] = kwargs["canonical_nucl"]["location"][
                0
            ]
        else:
            blastp_genome = self.rec_file.name  # noqa

        bbr_prot = BlastBasedReopening(
            genome=self.fnpfa,
            db=os.path.join(SCRIPT_DIR, "test-data/canonical_prot"),
            protein=True,
            genome_real=kwargs["canonical_nucl"]["reopened_file"],
            data_dir=self.data_dir,
        )
        kwargs["canonical_prot"] = bbr_prot.evaluate()

        # Failing that, if the genome is closed
        if self.closed:
            # we assume it is headful
            # yield (PhageType.PacHeadful, None, Evidence.Assumption)
            pass
        else:
            # Otherwise we truly have no clue
            # yield (PhageType.Unknown, None, Evidence.NONE)
            pass

        tpl = env.get_template("autoreopen_template.html")
        self.html.write(tpl.render(**kwargs).encode("utf-8"))

        shutil.copy(
            os.path.join(SCRIPT_DIR, "autoreopen_mauved3.js"),
            os.path.join(self.data_dir, "mauved3.js"),
        )


def main():
    parser = argparse.ArgumentParser(description="Automatic re-opening tool")
    parser.add_argument("fasta", type=argparse.FileType("r"))
    parser.add_argument("fastq1", type=argparse.FileType("r"))
    parser.add_argument("fastq2", type=argparse.FileType("r"))
    parser.add_argument("--closed", action="store_true")
    parser.add_argument("--data_dir", help="Directory for HTML files to go.")
    parser.add_argument("--html", type=argparse.FileType("w"))

    parser.add_argument("--reopen_phageTerm", default="asdf")
    parser.add_argument("--reopen_TerL", default="bsdf")
    parser.add_argument("--reopen_blastN", default="csdf")
    args = parser.parse_args()

    # create new IntronFinder object based on user input
    kwargs = vars(args)
    kwargs["fasta"] = args.fasta.name
    pr = PhageReopener(**kwargs)
    pr.giveEvidence()


if __name__ == "__main__":
    main()

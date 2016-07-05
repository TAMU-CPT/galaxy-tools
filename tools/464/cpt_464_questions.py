#!/usr/bin/env python
import os
import json
import bleach
import argparse
from string import Template
import logging
logging.basicConfig(level=logging.INFO)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='464 Questionnaire', epilog="""Fill out questions, run as many times as you need!

These questions are intended to reinforce key concepts that should be part of your GenomeA papers.""")
    parser.add_argument('--name', help='Phage Name')
    parser.add_argument('--seq_cov', help='Sequence Coverage')
    parser.add_argument('--morph', help='Morphology', choices=['Myo', 'Sipho', 'Podo'])
    parser.add_argument('--host', help='Host')
    parser.add_argument('--host_range', help='Host Range')

    parser.add_argument('--sample_type', help='Isolation Sample type')
    parser.add_argument('--date', help='Isolation Date')
    parser.add_argument('--location', help='Isolation Location')
    parser.add_argument('--isolator', help='Isolator')

    parser.add_argument('--annotator', help='Annotator (You!)')
    parser.add_argument('--last_annotated', help='Date of Last Annotation')

    parser.add_argument('--genome_length', help='Genome Length', type=int)
    parser.add_argument('--cdss', help='# of CDSs', type=int)
    parser.add_argument('--cdss_assigned', help='# of CDSs with assigned functions', type=int)
    parser.add_argument('--cdss_hypo_novel', help='# of Hypothetical Novels', type=int)
    parser.add_argument('--cdss_hypo_consv', help='# of Hypothetical Conserveds', type=int)
    parser.add_argument('--coding_density', help='Coding Density', type=float)
    parser.add_argument('--pigs', help='PIGS Score', type=float)
    parser.add_argument('--gc', help='GC Content', type=float)
    parser.add_argument('--trna', help='# of tRNAs', type=int)
    parser.add_argument('--terminators', help='# of terminators', type=int)

    parser.add_argument('--missing_core', help='Are there any core genes missing or out of place? ' +
                        'Identify and explain.')
    parser.add_argument('--no_phage_homos', help='Do you have any gene with assigned function that' +
                        ' has no phage sequence homolog? We are interested in "bacteria only" hits' +
                        '--this does not include your hypothetical novel genes. Justify the inclus' +
                        'ion of bacteria-only annotations.')
    parser.add_argument('--unusual', nargs='+', help='What unusual or unexpected genes or genomic ' +
                        'features did you find? (e.g., morons, introns, inteins, invertible sequen' +
                        'ces, etc. and genes that have unique features')
    parser.add_argument('--no_sd', nargs='+', help='Do you have any genes without a legitimate (se' +
                        'quence and spacing) Shine-Dalgarno? If so, for each gene, justify retaini' +
                        'ng them in the genome (i.e. does it have a conserved gene product?).')
    parser.add_argument('--gaps', nargs='+', help='Justify any intergenic gap larger than 150 bp. ' +
                        '(explain why you did not annotate genes?)')
    parser.add_argument('--closed', help='Is your genome closed? How did you determine this? (link' +
                        ' to gel images & analysis text , including primer sequence).')
    parser.add_argument('--reopened', help='Has your phage genome been re-opened? If so, where and' +
                        ' why? If not, why not?')
    parser.add_argument('--sim', help='Is your phage similar to another phage? Which one? How simi' +
                        'lar? How was the similarity determined?')
    parser.add_argument('--gc_skew', help='Does your phage exhibit GC skew? Can you approximate an' +
                        ' origin of replication? (Open your genome in Artemis, go to "Graph"->"GC ' +
                        'Deviation". Use the right scroll bar to adjust the window to 5000. Scroll' +
                        ' through your genome to see if and where your deviation line crosses from' +
                        ' positive to negative.) Does this origin fit with the location of DNA rep' +
                        'lication genes?')
    parser.add_argument('--trna_genes', help='Do you have any tRNA genes?')

    parser.add_argument('--email', help="Email address")

    args = parser.parse_args()

    try:
        report_path = os.path.join('/', 'home', 'galaxy', '464-reports', args.email)
        with open(report_path, 'w') as handle:
            handle.write(json.dumps(vars(args)))
    except:
        pass

    s = Template("""
                 <html><head>
                 <link href="//maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css" rel="stylesheet">
                 </head>
                 <body>
                 <div class="container">
                    <h1>$name</h1>

                    <h3>Basic information</h3>
                    <table class="table table-striped">
                        <tr><td>Sequence Coverage</td><td>$seq_cov</td></tr>
                        <tr><td>Morphology</td><td>$morph</td></tr>
                        <tr><td>Host</td><td>$host</td></tr>
                        <tr><td>Host Range</td><td>$host_range</td></tr>
                    </table>

                    <h3>Isolation Information</h3>
                    <table class="table table-striped">
                        <tr><td>Sample Type</td><td>$sample_type</td></tr>
                        <tr><td>Date</td><td>$date</td></tr>
                        <tr><td>Location</td><td>$location</td></tr>
                        <tr><td>Isolator</td><td>$isolator</td></tr>
                    </table>

                    <h3>Annotation Properties</h3>
                    <table class="table table-striped">
                        <tr><td>Annotator</td><td><a href="mailto:$email">$annotator</a></td></tr>
                        <tr><td>Last Annotated</td><td>$last_annotated</td></tr>
                    </table>

                    <h3>Genome Properties</h3>
                    <table class="table table-striped">
                        <tr><td>Genome Length</td><td>$genome_length</td></tr>
                        <tr><td># of CDSs</td><td>$cdss</td></tr>
                        <tr><td># of Assigned CDSs</td><td>$cdss_assigned</td></tr>
                        <tr><td># of Hypothetical Novel CDSs</td><td>$cdss_hypo_novel</td></tr>
                        <tr><td># of Hypothetical Conserved CDSs</td><td>$cdss_hypo_consv</td></tr>
                        <tr><td>Coding Density</td><td>$coding_density</td></tr>
                        <tr><td>PIGS (Percent InterGenic Sequence)</td><td>$pigs</td></tr>
                        <tr><td>GC</td><td>$gc</td></tr>
                        <tr><td>tRNA</td><td>$trna</td></tr>
                        <tr><td>terminators</td><td>$terminators</td></tr>
                    </table>

                    <h3>Genome Questions</h3>
                    <table class="table table-striped">
                        <tr><td>Are there any core genes missing or out of place? Identify and explain.</td></tr>
                        <tr><td>$missing_core</td></tr>
                        <tr><td>Do you have any gene with assigned function that has no phage sequence homolog?
                        We are interested in "bacteria only" hits--this does not include your hypothetical novel
                         genes. Justify the inclusion of bacteria-only annotations.</td></tr>
                        <tr><td>$no_phage_homos</td></tr>
                        <tr><td>What unusual or unexpected genes or genomic features did you find? (e.g., morons,
                         introns, inteins, invertible sequences, etc. and genes that have unique features)</td></tr>
                        <tr><td>$unusual</td></tr>
                        <tr><td>Do you have any genes without a legitimate (sequence and spacing) Shine-Dalgarno?
                         If so, for each gene, justify retaining them in the genome (i.e. does it have a
                         conserved gene product?).</td></tr>
                        <tr><td>$no_sd</td></tr>
                        <tr><td>Justify any intergenic gap larger than 150 bp. (explain why you did not annotate
                         genes?)</td></tr>
                        <tr><td>$gaps</td></tr>
                        <tr><td>Is your genome closed? How did you determine this? (link to gel images &
                         analysis text , including primer sequence).</td></tr>
                        <tr><td>$closed</td></tr>
                        <tr><td>Has your phage genome been re-opened? If so, where and why? If not, why not?</td></tr>
                        <tr><td>$reopened</td></tr>
                        <tr><td>Is your phage similar to another phage? Which one? How similar? How was the
                         similarity determined?</td></tr>
                        <tr><td>$sim</td></tr>
                        <tr><td>Does your phage exhibit GC skew? Can you approximate an origin of replication?
                         (Open your genome in Artemis, go to "Graph"/"GC Deviation". Use the right scroll bar
                         to adjust the window to 5000. Scroll through your genome to see if and where your deviation
                         line crosses from positive to negative.) Does this origin fit with the location of DNA
                         replication genes?</td></tr>
                        <tr><td>$gc_skew</td></tr>
                        <tr><td>Do you have any tRNA genes?</td></tr>
                        <tr><td>$trna_genes</td></tr>
                    </table>
                </div>
                 </body></html>
                 """)

    fixed_vars = {}
    vargs = vars(args)
    for x in vargs.keys():
        if isinstance(vargs[x], list):
            fixed_vars[x] = [bleach.clean(y) for y in vargs[x]]
        else:
            fixed_vars[x] = bleach.clean(vargs[x])

    fixed_vars['unusual'] = '<ul><li>' + '</li><li>'.join(fixed_vars['unusual']) + '</li></ul>'
    fixed_vars['no_sd'] = '<ul><li>' + '</li><li>'.join(fixed_vars['no_sd']) + '</li></ul>'
    fixed_vars['gaps'] = '<ul><li>' + '</li><li>'.join(fixed_vars['gaps']) + '</li></ul>'

    print s.substitute(**fixed_vars)

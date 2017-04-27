from BCBio import GFF
import glob

# If they find this file in our git repo... (i.e. no need to obfuscate because
# no student will ever look here.)
ref = list(GFF.parse(open('assessment2-answers.gff3')))[0]
ref = ref[213:213 + 4113]
start_bounds = [f.location.start for f in ref.features]
end_bounds = [f.location.end for f in ref.features]
fMap = {f.location.start: f for f in ref.features}
total = len(ref.features)
scores = []

for f in glob.glob("*-464-*.gff"):
    with open(f, 'r') as handle:
        data = GFF.parse(f)
    for rec in data:
        correct = 0
        incorrect = 0
        owner = None
        for f in rec.features:
            if not owner:
                owner = f.qualifiers['owner'][0]

            if f.location.start in start_bounds:
                # at least partial credit they will receive
                if f.location.end == fMap[f.location.start].location.end:
                    correct += 1
                else:
                    correct += 0.5
            else:
                incorrect += 1
        scores.append((owner, correct, incorrect))

print "Given that we have jason.gill@tamu.edu's annotation set mixed in, we can curve based on his ability to annotate this portion"
max_score = max([x[1] for x in scores])
curve = total - max_score
print 'Curve: ', curve
print

for (owner, correct, incorrect) in sorted(scores, key=lambda x: x[0]):
    print "%30s c=%0.1f i=%d score=(%0.1f / %d) = %d" % (owner, correct, incorrect, correct, total, (100 * float(correct + curve) / total))

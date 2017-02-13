# Running 

## Prep

First copy r1.fq and r2.fq from the CPT server

```
scp <username>@cpt.tamu.edu:r*.fq reopen-data/
```

There are new pip dependencies,

```
pip install -r ../../requirements.txt
```

This tool does *awful* things, so you'll want to be running under python2.

## Ready to run

```
python autoreopen.py reopen-data/mis.fa reopen-data/r1.fq reopen-data/r2.fq
```


# Outline

Generally in three phases:

- identify type
- identify location
- reopen at location

## Identify Type Phase (done)

- if Positive result = (3', 5' COS, TR)
	- if TR (double coverage)
		- RETURN TR + location
	- else if 3'/5':
		- RETURN COS SITE + location
	/*: Open at the cos site.*/
- else naive annotation, find terminase, analyse.
	- if success:
		- RETURN type (+ assume correct) + location
- else if closed:
	- RETURN "ASSUME PAC"
- else (unclosed):
	- RETURN "UNKNOWN: AS_IS"


## Identify Location Phase (not done)

def alignToCanonical():
- blastn against the canonical phage database
- if there are good hits (> 50% percent identity? Pulling this metric out
  of thin air), then <SOMEHOW> re-open relative to that genome.


- if location from TerL Analysis or from PhageTerm
	- reopen there
- if we can alignToCanonical():
	- reopen there
	- if the alignment was to T1, then we rev-com the genome after re-opening relative to parent.
- else:
	- if closed (i.e. Pac): tell the user "Unknown Pac headful", no idea where to re-open
	- else (unclosed): tell the user "Completely Unknown"

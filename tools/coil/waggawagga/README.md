    # **********************************************************************
    #  Copyright notice
    #
    #  (c) 2011-2017 Dominic Simm <dominic.simm@mpibpc.mpg.de>
    #  All rights reserved
    #
    #  This file is part of Waggawagga.
    #
    #  Waggawagga is free software: you can redistribute it and/or modify
    #  it under the terms of the GNU General Public License as published by
    #  the Free Software Foundation, either version 2 of the License, or
    #  (at your option) any later version.
    #
    #  Waggawagga is distributed in the hope that it will be useful,
    #  but WITHOUT ANY WARRANTY; without even the implied warranty of
    #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #  GNU General Public License for more details.
    #
    #  You should have received a copy of the GNU General Public License
    #  along with Waggawagga.  If not, see <http://www.gnu.org/licenses/>.
    # **********************************************************************

WaggaWagga
==========

Please cite:
Simm D., Hatje K. and Kollmar M. (2014) Waggawagga: comparative visualization
of coiled-coil predictions and detection of stable single α-helices (SAH domains).
Bioinformatics. 31(5):767-769. doi://10.1093/bioinformatics/btu700

In a nutshell, Waggawagga-CLI is the command-line version of the web-application Waggawagga <http://waggawagga.motorprotein.de> for the SAH prediction in proteins. 

Installation
------------
The standalone version 'wagggawagga-cli' is designed to run instantly out-of-the-box. So after unpacking the downloaded tarball, you can directly start using 'Waggawagga'.

Besides the actual program, it comes with an own mobile SQLite database, where it stores your analysed sequence data for direct or later use. For a quickstart to see, if `waggawagga-cli` is working correctly on your system. Start with the following example, the produced output should look like this:

    ➜  ./waggawagga-cli -a example
    Import:
    - Create session-directory <$HOME>/waggawagga-cli/results/example_myosin_seqs
    - FASTA-file: test/example_myosin_seqs.fas
    - Found 4 sequences to import (4/4 sequences with len>=49aa)

    Progress: 100.00% (37.43 sec.)

    - Imported into SQLite-database under identifier 'genome' = "example_myosin_seqs"

    Evaluation:
    Calculation for window-size: 14
    All Found SAHs:  2
    Max SAH per seq: 2
    - Filter data
    Calculation for window-size: 21
    [...]

    Find your filtered results in
    - <$HOME>/waggawagga-cli/results/example_myosin_seqs

Quickstart
----------
Before you can start using 'waggawagga-cli', you need to import first some protein sequence data into the piggybacked SQLite database to work with. Essentially you simply need a multiple FASTA file carrying your sequences and a self-chosen identifier for your data to be found later in the database. For an example call have a look at the below listed example calls:

    # Analyse your sequences (FASTA)
    ➜  ./waggawagga-cli 'path/to/your/fasta.fas'
    
    # Analyse your sequences (FASTA) - More advanced with GNUplots and clean-up
    ➜  ./waggawagga-cli -a complete -f 'path/to/your/fasta.fas'
    
For more detailed instructions please have a look at the command line help or the documentation.

    ➜  ./waggawagga-cli --help

Troubleshooting
---------------
We are using the Traveling-Ruby framework [<https://github.com/phusion/traveling-ruby>] for packaging our Ruby-on-Rails-based *Waggawagga-CLI*-tool for different plaforms. Therefore we are currently not able to provide packages for others than the following platforms: Linux-x86 32-Bit, Linux-x86 64-Bit and Mac OS X / macOS (later 10.7.x).

If you get errors like the following on startup, there is an incompatibility between your system and the chosen package: 

    ➜  ./waggawagga-cli -a example
    [...]/lib/ruby/bin/ruby: line 6: [...]/ruby: No such file or directory
    [...]/lib/ruby/bin/ruby: line 6: [...]/ruby: Cannot execute binary file
    [...]/lib/ruby/bin/ruby: line 6: [...]/ruby: Exec format error: Binary file not executable.

Then try another version of 'waggawagga-cli':
 
- waggawagga-cli-0.5.4-linux-x86.tar.gz
- waggawagga-cli-0.5.4-linux-x86_64.tar.gz
- waggawagga-cli-0.5.4-osx.tar.gz

                                                       
Rights and Restrictions
-----------------------
Waggawagga-CLI may be used under a GNU General Public License. Using Waggawagga by non-academics requires permission. 

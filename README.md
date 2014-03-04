PileupModule
============
This module can be used to parse pileup file. 

Files
-----
* setup.py
* parsePileup.py
* README.md

Installing pileupmodule
---------------------

    python setup.py install


Using pileupmodule
------------------
After installation, type:

    from pileupmodule import ParsePileup
    p = ParsePileup()
    p.run(argument)

Output
------
Output is a column delimited file

    SCAFFOLD_ID	POSITION	BASE*	DEPTH	BASECOUNT(A C G T)***

#*Top 2 most common bases will be encapsulated in []
#**Repetition of A, C, G, T base count for N individuals from pileup file

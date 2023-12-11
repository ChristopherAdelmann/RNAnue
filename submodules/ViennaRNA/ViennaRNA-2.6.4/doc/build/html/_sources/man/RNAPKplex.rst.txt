#########
RNAPKplex
#########

:program:`RNAPKplex` - manual page for RNAPKplex 2.6.4

Synopsis
--------

.. code:: bash

    RNAPKplex [OPTION]...

DESCRIPTION
-----------

RNAPKplex 2.6.4

predicts RNA secondary structures including pseudoknots

Computes RNA secondary structures by first making two sequence intervals
accessible and unpaired using the algorithm of RNAplfold and then calculating
the energy of the interaction of those two intervals. The algorithm uses
O(n^2*w^4) CPU time and O(n*w^2) memory space.
The algorithm furthermore always considers dangle=2 model.


It  also  produces a PostScript file with a plot of the pseudoknot-free
secondary structure graph, in which the bases  forming  the  pseuodknot
are marked red.

Sequences are read in a simple text format where each sequence occupies
a single line. Each sequence may be preceded by a line of the form
.. code::

    > name
    

to assign a name to the sequence. If a name is given in the input, the
PostScript file "name.ps" is produced for the structure graph.  Other-
wise  the  file  name defaults to PKplex.ps. Existing files of the same
name will be overwritten.
The input format is similar to fasta except that  even  long  sequences
may  not  be  interrupted  by  line  breaks,  and  the header lines are
optional.  The program will continue to read new sequences until a line
consisting  of  the  single  character @ or an end of file condition is
encountered.

.. option:: -h, --help

    Print help and exit

.. option:: --detailed-help

    Print help, including all details and hidden options, and exit

.. option:: --full-help

    Print help, including hidden options, and exit

.. option:: -V, --version

    Print version and exit

.. option:: -v, --verbose

    Be verbose.


    *(default=off)*

I/O Options:
^^^^^^^^^^^^



    Command line options for input and output (pre-)processing

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.

.. option:: -c, --cutoff=FLOAT

    Report only base pairs with an average probability > cutoff in the dot plot.


    *(default="0.01")*

.. option:: -e, --energyCutoff=DOUBLE

    Energy cutoff or pseudoknot initiation cost. Minimum energy gain of a pseudoknot interaction for it to be returned. Pseudoknots with smaller energy gains are rejected.


    *(default="-8.10")*

.. option:: -s, --subopts=DOUBLE

    print suboptimal structures whose energy difference of the pseudoknot to the optimum pseudoknot is smaller than the given value.


    *(default="0.0")*


    NOTE: The final energy of a structure is calculated as the sum of the
    pseudoknot interaction energy, the penalty for initiating a  pseudoknot and
    the energy of the pseudoknot-free part of the structure. The :option:`-s` option only
    takes the pseudoknot interaction energy into account, so the final energy
    differences may be bigger than the specified value *(default=0.)*.

Energy Parameters:
^^^^^^^^^^^^^^^^^^



    Energy parameter sets can be adapted or loaded from user-provided input files

.. option:: -T, --temp=DOUBLE

    Rescale energy parameters to a temperature of temp C. Default is 37C.


    *(default="37.0")*

.. option:: -P, --paramFile=paramfile

    Read energy parameters from paramfile, instead of using the default parameter set.


    Different sets of energy parameters for RNA and DNA should accompany your
    distribution.
    See the RNAlib documentation for details on the file format. The placeholder
    file name ``DNA`` can be used to load DNA parameters without the need to
    actually specify any input file.

.. option:: -4, --noTetra

    Do not include special tabulated stabilizing energies for tri-, tetra- and hexaloop hairpins.


    *(default=off)*


    Mostly for testing.

.. option:: --salt=DOUBLE

    Set salt concentration in molar (M). Default is 1.021M.

Model Details:
^^^^^^^^^^^^^^



    Tweak the energy model and pairing rules additionally using the following
    parameters

.. option:: --noLP

    Produce structures without lonely pairs (helices of length 1).


    *(default=off)*


    For partition function folding this only disallows pairs that can only occur
    isolated. Other pairs may still occasionally occur as helices of length 1.

.. option:: --noGU

    Do not allow GU pairs.


    *(default=off)*

.. option:: --noClosingGU

    Do not allow GU pairs at the end of helices.


    *(default=off)*

.. option:: --nsp=STRING

    Allow other pairs in addition to the usual AU,GC,and GU pairs.


    Its argument is a comma separated list of additionally allowed pairs. If the
    first character is a "-" then AB will imply that AB and BA are allowed
    pairs, e.g. :option:`--nsp=`"-GA"  will allow GA and AG pairs. Nonstandard pairs are
    given 0 stacking energy.

.. option:: --helical-rise=FLOAT

    Set the helical rise of the helix in units of Angstrom.


    *(default="2.8")*


    Use with caution! This value will be re-set automatically to 3.4 in case DNA
    parameters are loaded via :option:`-P` DNA and no further value is provided.

.. option:: --backbone-length=FLOAT

    Set the average backbone length for looped regions in units of Angstrom.


    *(default="6.0")*


    Use with caution! This value will be re-set automatically to 6.76 in case DNA
    parameters are loaded via :option:`-P` DNA and no further value is provided.

REFERENCES
----------

*If you use this program in your work you might want to cite:*

R. Lorenz, S.H. Bernhart, C. Hoener zu Siederdissen, H. Tafer, C. Flamm, P.F. Stadler and I.L. Hofacker (2011),
"ViennaRNA Package 2.0",
Algorithms for Molecular Biology: 6:26

I.L. Hofacker, W. Fontana, P.F. Stadler, S. Bonhoeffer, M. Tacker, P. Schuster (1994),
"Fast Folding and Comparison of RNA Secondary Structures",
Monatshefte f. Chemie: 125, pp 167-188

R. Lorenz, I.L. Hofacker, P.F. Stadler (2016),
"RNA folding with hard and soft constraints",
Algorithms for Molecular Biology 11:1 pp 1-13

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Wolfgang Beyer

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
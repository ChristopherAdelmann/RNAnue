#########
RNA2Dfold
#########

:program:`RNA2Dfold` - manual page for RNA2Dfold 2.6.4

Synopsis
--------

.. code:: bash

    RNA2Dfold [OPTION]...

DESCRIPTION
-----------

RNA2Dfold 2.6.4

Compute MFE structure, partition function and representative sample structures
of k,l neighborhoods

The program partitions the secondary structure space into (basepair)distance
classes according to two fixed reference structures. It expects a sequence and
two secondary structures in dot-bracket notation as its inputs. For each
distance class, the MFE representative, Boltzmann probabilities and Gibbs free
energy is computed. Additionally, a stochastic backtracking routine allows one
to produce samples of representative suboptimal secondary structures from each
partition

.. option:: -h, --help

    Print help and exit

.. option:: --detailed-help

    Print help, including all details and hidden options, and exit

.. option:: --full-help

    Print help, including hidden options, and exit

.. option:: -V, --version

    Print version and exit

I/O Options:
^^^^^^^^^^^^



    Command line options for input and output (pre-)processing

.. option:: -j, --numThreads=INT

    Set the number of threads used for calculations (only available when compiled with OpenMP support)

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.
    The Minimum free energy (MFE) and a structure representative are calculated
    in any case.

.. option:: -p, --partfunc

    calculate partition function and thus, Boltzmann probabilities and Gibbs free energy


    *(default=off)*

.. option:: --stochBT=INT

    backtrack a certain number of Boltzmann samples from the appropriate k,l neighborhood(s)

.. option:: --neighborhood=<k>:<l>

    backtrack structures from certain k,l-neighborhood only, can be specified multiple times (<k>:<l>,<m>:<n>,...)

.. option:: -K, --maxDist1=INT

    maximum distance to first reference structure


    If this value is set all structures that exhibit a basepair distance greater
    than maxDist1 will be thrown into a distance class denoted by K=L=-1

.. option:: -L, --maxDist2=INT

    maximum distance to second reference structure


    If this value is set all structures that exhibit a basepair distance greater
    than maxDist1 will be thrown into a distance class denoted by K=L=-1

.. option:: -S, --pfScale=DOUBLE

    In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows).


    *(default="1.07")*


    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for
    long sequences.

.. option:: --noBT

    do not backtrack structures, calculate energy contributions only


    *(default=off)*

.. option:: -c, --circ

    Assume a circular (instead of linear) RNA molecule.


    *(default=off)*

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

.. option:: -d, --dangles=INT

    How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops


    (possible values="0", "2" default="2")


    With :option:`-d2` dangling energies will be added for the bases adjacent to a helix on
    both sides in any case. The option :option:`-d0` ignores dangling ends altogether
    (mostly for debugging).

.. option:: --noGU

    Do not allow GU pairs.


    *(default=off)*

.. option:: --noClosingGU

    Do not allow GU pairs at the end of helices.


    *(default=off)*

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

R. Lorenz, C. Flamm, I.L. Hofacker (2009),
"2D Projections of RNA folding Landscapes",
GI, Lecture Notes in Informatics, German Conference on Bioinformatics 2009: 157, pp 11-20

M. Zuker, P. Stiegler (1981),
"Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information",
Nucl Acid Res: 9, pp 133-148

J.S. McCaskill (1990),
"The equilibrium partition function and base pair binding probabilities for RNA secondary structures",
Biopolymers: 29, pp 1105-1119

I.L. Hofacker and P.F. Stadler (2006),
"Memory Efficient Folding Algorithms for Circular RNA Secondary Structures",
Bioinformatics

D. Adams (1979),
"The hitchhiker's guide to the galaxy",
Pan Books, London

The calculation of mfe structures is based on dynamic programming algorithm originally developed by M. Zuker and P. Stiegler. The partition function algorithm is based on work by J.S. McCaskill.

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
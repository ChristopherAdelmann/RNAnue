########
RNApdist
########

:program:`RNApdist` - manual page for RNApdist 2.6.4

Synopsis
--------

.. code:: bash

    RNApdist [OPTION]...

DESCRIPTION
-----------

RNApdist 2.6.4

Calculate distances between thermodynamic RNA secondary structures ensembles

This program reads RNA sequences from stdin and calculates structure distances
between the thermodynamic ensembles of their secondary structures.

To do this the partition function and matrix of base pairing probabilities is
computed for each sequence. The probability matrix is then condensed into a
vector holding for each base the probabilities of being unpaired, paired
upstream, or paired downstream, respectively. These profiles are compared
by a standard alignment algorithm.

The base pair probabilities are also saved as postscript "dot plots" (as in
RNAfold) in the files  "name_dp.ps", where name is the name of the sequence,
or a number if unnamed.

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

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.

.. option:: -X, --compare=p|m|f|c

    Specify the comparison directive. *(default="p")*


    Possible arguments for this option are: :option:`-Xp` compare the structures pairwise
    (p), i.e. first with 2nd, third with 4th etc.
    :option:`-Xm` calculate the distance matrix between all structures. The output is
    formatted as a lower triangle matrix.
    :option:`-Xf` compare each structure to the first one.
    :option:`-Xc` compare continuously, that is i-th with (i+1)th structure.

.. option:: -B, --backtrack[=<filename>]

    Print an "alignment" with gaps of the profiles. The aligned structures are written to <filename>, if specified.


    *(default="none")*


    Within the profile output, the following symbols will be used:

.. option:: ()

    essentially upstream (downstream) paired bases

.. option:: {}

    weakly upstream (downstream) paired bases

.. option:: |

    strongly paired bases without preference

.. option:: ,

    weakly paired bases without preference

.. option:: .

    essentially unpaired bases.


    If <filename> is not specified, the output is written to stdout, unless the


    "-Xm" option is set in which case "backtrack.file" is used.

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
    See the RNAlib documentation for details on the file format. When passing the
    placeholder file name "DNA", DNA parameters are loaded without the need to
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

    set energy model for treatment of dangling bases.


    (possible values="0", "2" default="2")

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
    pairs.
    e.g. RNAfold :option:`-nsp` :option:`-GA`  will allow GA and AG pairs. Nonstandard pairs are
    given 0 stacking energy.

.. option:: -e, --energyModel=INT

    Set energy model.


    Rarely used option to fold sequences from the artificial ABCD... alphabet,
    where A pairs B, C-D etc.  Use the energy parameters for GC (:option:`-e` 1) or AU (:option:`-e`
    2) pairs.

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

S. Bonhoeffer, J.S. McCaskill, P.F. Stadler, P. Schuster (1993),
"RNA multi-structure landscapes",
Euro Biophys J:22, pp 13-24

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Peter F Stadler, Ivo L Hofacker, Sebastian Bonhoeffer.

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
#######
RNApaln
#######

:program:`RNApaln` - manual page for RNApaln 2.6.4

Synopsis
--------

.. code:: bash

    RNApaln [OPTION]...

DESCRIPTION
-----------

RNApaln 2.6.4

RNA alignment based on sequence base pairing propensities

Uses string-alignment techniques to perform fast pairwise structural alignments
of RNAs. Similar to RNApdist secondary structure is incorporated in an
approximate manner by computing base pair probabilities, which are then reduced
to a vector holding the probability that a base is paired upstream, downstream,
or remains unpaired. Such pair propsensity vectors can then be compared using
standard alignment algorithms. In contrast to RNApdist, RNApaln performs
similarity (instead of distance) alignments, considers both sequence and
structure information, and uses affine (rather than linear) gap costs. RNApaln
can perform semi-local alignments by using free end gaps, a true local
alignment mode is planned.

The same approach has since been used in the StraL program from Gerhard
Steeger's group. Since StraL has optimized parameters and a multiple alignment
mode, it be be currently the better option.

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

.. option:: -B, --printAlignment[=filename]

    Print an "alignment" with gaps of the

.. option:: profiles

    The aligned structures are written to filename, if specified Otherwise output is written to stdout, unless the :option:`-Xm` option is set in which case "backtrack.file" is used.


    *(default="stdout")*


    The following symbols are used:

.. option:: (

    ) essentially upstream (downstream) paired bases

.. option:: {

    } weakly upstream (downstream) paired bases

.. option:: |

    strongly paired bases without preference

.. option:: ,

    weakly paired bases without preference

.. option:: .

    essentially unpaired bases.

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.

.. option:: -X, --mode=pmfc

    Set the alignment mode to be used.


    The alignment mode is passed as a single character value. The following
    options are available:
    ``p`` - Compare the structures pairwise, that is first with 2nd, third with 4th
    etc. This is the default.

.. option:: ``m``

    - Calculate the distance matrix between all structures. The output is


    formatted as a lower triangle matrix.


    ``f`` - Compare each structure to the first one.


    ``c`` - Compare continuously, that is i-th with (i+1)th structure.

.. option:: --gapo=open

    Set the gap open penalty

.. option:: --gape=ext

    Set the gap extension penalty

.. option:: --seqw=w

    Set the weight of sequence (compared to structure) in the scoring function.

.. option:: --endgaps

    Use free end-gaps


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

    How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.


    *(default="2")*


    With :option:`-d1` only unpaired bases can participate in at most one dangling end.
    With :option:`-d2` this check is ignored, dangling energies will be added for the bases
    adjacent to a helix on both sides in any case; this is the default for mfe
    and partition function folding (:option:`-p`).
    The option :option:`-d0` ignores dangling ends altogether (mostly for debugging).
    With :option:`-d3` mfe folding will allow coaxial stacking of adjacent helices in
    multi-loops. At the moment the implementation will not allow coaxial stacking
    of the two interior pairs in a loop of degree 3 and works only for mfe
    folding.


    Note that with :option:`-d1` and :option:`-d3` only the MFE computations will be using this
    setting while partition function uses :option:`-d2` setting, i.e. dangling ends will be
    treated differently.

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

Bonhoeffer S, McCaskill J S, Stadler P F, Schuster P (1993),
"RNA multi-structure landscapes",
Euro Biophys J: 22, pp 13-24

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Peter F Stadler, Ivo L Hofacker, Sebastian Bonhoeffer

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
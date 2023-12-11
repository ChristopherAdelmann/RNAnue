##########
RNAinverse
##########

:program:`RNAinverse` - manual page for RNAinverse 2.6.4

Synopsis
--------

.. code:: bash

    RNAinverse [OPTION]...

DESCRIPTION
-----------

RNAinverse 2.6.4

Find RNA sequences with given secondary structure

The program searches for sequences folding into a predefined structure, thereby
inverting the folding algorithm. Target structures (in bracket notation) and
starting sequences for the search are read alternately from stdin.
Characters in the start sequence other than "AUGC" (or the alphabet specified
with :option:`-a`) will be treated as wild cards and replaced by a random character. Any
lower case characters in the start sequence will be kept fixed during the
search. If necessary, the sequence will be elongated to the length of the
structure. Thus a string of "N"s as well as a blank line specify a random
start sequence.
For each search the best sequence found and its Hamming distance to the start
sequence are printed to stdout. If the the search was unsuccessful, a structure
distance to the target is appended.
The :option:`-Fp` and :option:`-R` options can modify the output format, see commandline options
below.
The program will continue to read new structures and sequences until a line
consisting of the single character "@" or an end of file condition is
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

    In conjunction with a negative value supplied to :option:`-R`, print the last subsequence and substructure for each unsuccessful search.


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.

.. option:: -F, --function=mp

    Use minimum energy (:option:`-Fm`), partition function folding (:option:`-Fp`) or both (:option:`-Fmp`).


    *(default="m")*


    In partition function mode, the probability of the target structure
    exp(:option:`-E`(S)/kT)/Q is maximized. This probability is written in brackets after
    the found sequence and Hamming distance. In most cases you'll want to use the
    :option:`-f` option in conjunction with :option:`-Fp`, see below.

.. option:: -f, --final=FLOAT

    In combination with :option:`-Fp` stop search when sequence is found with E(s)-F is smaller than final, where F=-kT*ln(Q).

.. option:: -R, --repeat[=INT]

    Search repeatedly for the same structure. If an argument is supplied to this option it must follow the option flag immediately. E.g.: :option:`-R5`


    *(default="1")*


    If repeats is negative search until :option:`--repeats` exact solutions are found, no
    output is done for unsuccessful searches. Be aware, that the program will not
    terminate if the target structure can not be found.
    If no value is supplied with this option, the default value is used.

.. option:: -a, --alphabet=ALPHABET

    Find sequences using only nucleotides from a given alphabet.

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

    How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops


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

D.H. Turner, N. Sugimoto, S.M. Freier (1988),
"RNA structure prediction",
Ann Rev Biophys Biophys Chem: 17, pp 167-192

M. Zuker, P. Stiegler (1981),
"Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information",
Nucl Acid Res: 9, pp 133-148

J.S. McCaskill (1990),
"The equilibrium partition function and base pair binding probabilities for RNA secondary structures",
Biopolymers: 29, pp 1105-1119

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

EXAMPLES
--------


To search 5 times for sequences forming a simple hairpin structure interrupted by one GA mismatch call

.. code::

    $ RNAinverse -R 5
    



and enter the lines

.. code::

    (((.(((....))).)))
    NNNgNNNNNNNNNNaNNN
    


AUTHOR
------


Ivo L Hofacker

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
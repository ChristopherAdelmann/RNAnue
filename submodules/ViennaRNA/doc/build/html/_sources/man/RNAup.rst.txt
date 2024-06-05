#####
RNAup
#####

:program:`RNAup` - manual page for RNAup 2.6.4

Synopsis
--------

.. code:: bash

    RNAup [OPTION]...

DESCRIPTION
-----------

RNAup 2.6.4

Calculate the thermodynamics of RNA-RNA interactions

RNAup calculates the thermodynamics of RNA-RNA interactions, by decomposing the
binding into two stages. (1) First the probability that a potential binding
sites remains unpaired (equivalent to the free energy needed to open the site)
is computed. (2) Then this accessibility is combined with the interaction
energy to obtain the total binding energy. All calculations are done by
computing partition functions over all possible conformations.


RNAup provides two different modes: By default RNAup computes accessibilities,
in terms of the free energies needed to open a region (default length 4). It
prints the region of highest accessibility and its opening energy to stdout,
opening energies for all other regions are written to a file.

.br
In interaction mode the interaction between two RNAs is calculated. It is
invoked if the input consists of two sequences concatenated with an ``&``,
or if the options -X[pf] or -b are given. Unless the -b option is specified
RNAup assumes that the longer RNA is a structured target sequence
while the shorter one is an unstructured small RNA.
.br
Additionally, for every position along the target sequence we write the best
free energy of binding for an interaction that includes this position to the
the output file.
Output to stdout consists of the location and free energy, dG,
for the optimal region of interaction. The binding energy dG is also split into
its components the interaction energy dGint and the opening energy dGu_l (and
possibly dGu_s for the shorter sequence).
.br
In addition we print the optimal interaction structure as computed by RNAduplex
for this region. Note that it can happen that the RNAduplex computed optimal
interaction does not coincide with the optimal RNAup region. If the two
predictions don't match the structure string is replaced by a run of "."
and a message is written to stderr.
.br

Each sequence should be in 5`` to 3`` direction. If the sequence is preceded
by a line of the form
.. code::

    > name
    


the output file "name_ux_up.out" is produced, where the "x" in "ux" is the
value set by the -u option. Otherwise the file name defaults to
RNA_ux_up.out. The output is concatenated if a file with the same name exists.
.br

RNA sequences are read from stdin as strings of characters. White space and
newline within a sequence cause an error! Newline is used to separate
sequences. The program will continue to read new sequences until a line
consisting of the single character @ or an end of file condition is
encountered.

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

.. option:: -o, --no_output_file

    Do not produce an output file.


    *(default=off)*

.. option:: --no_header

    Do not produce a header with the command line parameters used in the outputfile.


    *(default=off)*

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.

.. option:: -u, --ulength=length

    Specify the length of the unstructured region in the output.


    *(default="4")*


    The probability of being unpaired is plotted on the right border of the
    unpaired region. You can specify up to 20 different length values: use "-"
    to specify a range of continuous values (e.g. :option:`-u` 4-8) or specify a list of
    comma separated values (e.g. :option:`-u` 4,8,15).

.. option:: -c, --contributions=SHIME

    Specify the contributions listed in the output. *(default="S")*


    By default only the full probability of being unpaired is plotted. The :option:`-c`
    option allows one to get the different contributions (c) to the probability
    of being unpaired: The full probability of being unpaired ("S" is the sum
    of the probability of being unpaired in the exterior loop ("E"), within a
    hairpin loop ("H"), within an interior loop ("I") and within a multiloop
    ("M"). Any combination of these letters may be given.

Calculations of RNA-RNA interactions:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. option:: -w, --window=INT

    Set the maximal length of the region of interaction.


    *(default="25")*

.. option:: -b, --include_both

    Include the probability of unpaired regions in both (b) RNAs.


    *(default=off)*


    By default only the probability of being unpaired in the longer RNA (target)
    is used.

.. option:: -5, --extend5=INT

    Extend the region of interaction in the target to some residues on the 5' side.


    The underlying assumption is that it is favorable for an interaction if not
    only the direct region of contact is unpaired but also a few residues 5'

.. option:: -3, --extend3=INT

    Extend the region of interaction in the target to some residues on the 3' side.


    The underlying assumption is that it is favorable for an interaction if not
    only the direct region of contact is unpaired but also a few residues 3'

.. option:: --interaction_pairwise

    Activate pairwise interaction mode. *(default=off)*


    The first sequence interacts with the 2nd, the third with the 4th etc. If
    activated, two interacting sequences may be given in a single line separated
    by "&" or each sequence may be given on an extra line.

.. option:: --interaction_first

    Activate interaction mode using first sequence only.


    *(default=off)*


    The interaction of each sequence with the first one is calculated (e.g.
    interaction of one mRNA with many small RNAs). Each sequence has to be given
    on an extra line

.. option:: -S, --pfScale=DOUBLE

    In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows).


    *(default="1.07")*


    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for
    long sequences.

Structure Constraints:
^^^^^^^^^^^^^^^^^^^^^^



    Command line options to interact with the structure constraints feature of
    this program

.. option:: -C, --constraint

    Apply structural constraint(s) during prediction.


    *(default=off)*


    The program first reads the sequence(s), then a dot-bracket like string
    containing constraints on the structure. The following symbols are
    recognized:


    ``.`` ... no constraint for this base


    ``x`` ... the base is unpaired


    ``<`` ... the base pairs downstream, i.e. i is paired with j > i


    ``>`` ... the base pairs upstream, i.e. i is paired with j < i


    ``()`` ... base i pairs with base j


    ``|`` ... the corresponding base has to be paired intermolecularily (only for


    interaction mode)

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

.. option:: --saltInit=DOUBLE

    Provide salt correction for duplex initialization (in kcal/mol).

Model Details:
^^^^^^^^^^^^^^



    Tweak the energy model and pairing rules additionally using the following
    parameters

.. option:: -d, --dangles=INT

    Specify "dangling end" model for bases adjacent to helices in free ends and multi-loops.


    *(default="2")*


    With :option:`-d2` dangling energies will be added for the bases adjacent to a helix on
    both sides in any case.


    The option :option:`-d0` ignores dangling ends altogether (mostly for debugging).

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

U. Mueckstein, H. Tafer, J. Hackermueller, S.H. Bernhart, P.F. Stadler, and I.L. Hofacker (2006),
"Thermodynamics of RNA-RNA Binding",
Bioinformatics: 22(10), pp 1177-1182

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

EXAMPLES
--------


.B Output to stdout:

In Interaction mode RNAup prints the most favorable interaction energy
between the two sequences to stdout. The most favorable interaction energy
(dG) depends on the position in the longer sequence (region [i,j]) and the
position in the shorter sequence (region[k,l]): dG[i,j;k,l].  dG[i,j;k,l] is the
largest contribution to dG[i,j] = sum_kl dG[i,j;k,l] which is given in the
output file: therefore dG[i,j;k,l] <= dG[i,j].

.. code::

    ``....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8``
    > franz
    GGAGUAGGUUAUCCUCUGUU
    > sissi
    AGGACAACCU
    dG = dGint + dGu_l
    (((((.((((&)))).)))))   6,15  :   1,10  (-6.66 = -9.89 + 3.23)
    AGGUUAUCCU&AGGACAACCU
    RNAup output in file: franz_sissi_w25_u3_4_up.out
    


where the result line contains following information

.. code::

    RNAduplex results       [i,j]     [k,l]    dG = dGint + dGu_l
    (((((.((((&)))).)))))   6,15   :  1,10     (-6.66=-9.89+3.23)
    



.RD
.B Output to file:

Output to file contains a header including date, the command line of the
call to RNAup, length and names of the input sequence(s) followed
by the sequence(s). The first sequence is the target sequence.
Printing of the header can be turned off using the -nh option.

The line directly after the header gives the column names for the output:

.. code::

    position     dGu_l for -u 3      dGu_l for -u 4       dG
    #     pos      u3S       u3H       u4S       u4H        dG
    


where all information refers to the target sequence. The dGu_l column contains
information about the -u value (u=3 or u=4) and the contribution to the free
energy to open all structures "S" or only hairpin loops "H", see option -c.
NA means that no results is possible (e.g. column u3S row 2: no region of
length 3 ending at position 2 exists).

.. code::

    #  Thu Apr 10 09:15:11 2008
    #  RNAup -u 3,4 -c SH -b
    #  20 franz
    #  GGAGUAGGUUAUCCUCUGUU
    #  10 sissi
    #  AGGACAACCU
    #     pos      u3S       u3H       u4S       u4H        dG
    1        NA        NA        NA        NA    -1.540
    2        NA        NA        NA        NA    -1.540
    3     1.371        NA        NA        NA    -1.217
    4     1.754     5.777     1.761        NA    -1.393
    5     1.664     3.140     1.811     5.800    -1.393
    



If the -b option is selected position and dGu_s values for the shorter sequence
are written after the information for the target sequence.

AUTHOR
------


Ivo L Hofacker, Peter F Stadler, Ulrike Mueckstein, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
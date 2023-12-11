###########
RNALalifold
###########

:program:`RNALalifold` - manual page for RNALalifold 2.6.4

Synopsis
--------

.. code:: bash

    RNALalifold [options] <file1.aln>

DESCRIPTION
-----------

RNALalifold 2.6.4

calculate locally stable secondary structures for a set of aligned RNAs

reads aligned RNA sequences from stdin or file.aln and calculates locally
stable RNA secondary structure with a maximal base pair span. For a sequence of
length n and a base pair span of L the algorithm uses only O(n+L*L) memory and
O(n*L*L) CPU time. Thus it is practical to "scan" very large genomes for
short RNA


    structures.

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

.. option:: -q, --quiet

    Be quiet. *(default=off)*


    This option can be used to minimize the output of additional information and
    non-severe warnings which otherwise might spam stdout/stderr.

I/O Options:
^^^^^^^^^^^^



    Command line options for input and output (pre-)processing

.. option:: -f, --input-format=C|S|F|M

    File format of the input multiple sequence alignment (MSA).


    If this parameter is set, the input is considered to be in a particular file
    format. Otherwise, the program tries to determine the file format
    automatically, if an input file was provided in the set of parameters. In
    case the input MSA is provided in interactive mode, or from a terminal (TTY),
    the programs default is to assume CLUSTALW format.
    Currently, the following formats are available: ClustalW (``C``), Stockholm 1.0
    (``S``), FASTA/Pearson (``F``), and MAF (``M``).

.. option:: --csv

    Create comma separated output (csv)


    *(default=off)*

.. option:: --aln[=prefix]

    Produce output alignments and secondary structure plots for each hit found.


    This option tells the program to produce, for each hit, a colored and
    structure annotated (sub)alignment and secondary structure plot in PostScript
    format. It also adds the subalignment hit into a multi-Stockholm formatted
    file "RNALalifold_results.stk". The postscript output file names are
    "aln_start_end.eps" and "ss_start_end.eps". All files will be created in
    the current directory. The optional argument string can be used to set a
    specific prefix that is used to name the output files. The file names then
    become "prefix_aln_start_end.eps", "prefix_ss_start_end.eps", and
    "prefix.stk". Note: Any special characters in the prefix will be replaced
    by the filename delimiter, hence there is no way to pass an entire directory
    path through this option yet. (See also the "--filename-delim" parameter)

.. option:: --aln-stk[=prefix]

    Add hits to a multi-Stockholm formatted output file.


    *(default="RNALalifold_results")*


    The default file name used for the output is "RNALalifold_results.stk".
    Users may change the filename to "prefix.stk" by specifying the prefix as
    optional argument. The file will be create in the current directory if it
    does not already exist. In case the file already exists, output will be
    appended to it. Note: Any special characters in the prefix will be replaced
    by the filename delimiter, hence there is no way to pass an entire directory
    path through this option yet. (See also the "--filename-delim" parameter)

.. option:: --mis

    Output "most informative sequence" instead of simple consensus: For each column of the alignment output the set of nucleotides with frequency greater than average in IUPAC notation.


    *(default=off)*

.. option:: --split-contributions

    Split the free energy contributions into separate parts


    *(default=off)*


    By default, only the total energy contribution for each hit is returned.
    Using this option, this contribution is split into individual parts, i.e. the
    Nearest Neighbor model energy, the covariance pseudo energy, and if
    applicable, a remaining pseudo energy derived from special constraints, such
    as probing signals like SHAPE.

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

.. option:: --auto-id

    Automatically generate an ID for each alignment.


    *(default=off)*


    The default mode of RNALalifold is to automatically determine an ID from the
    input alignment if the input file format allows to do that. Alignment IDs
    are, for instance, usually given in Stockholm 1.0 formatted input. If this
    flag is active, RNALalifold ignores any IDs retrieved from the input and
    automatically generates an ID for each alignment.

.. option:: --id-prefix=STRING

    Prefix for automatically generated IDs (as used in output file names).


    *(default="alignment")*


    If this parameter is set, each alignment will be prefixed with the provided
    string. Hence, the output files will obey the following naming scheme:
    "prefix_xxxx_ss.ps" (secondary structure plot), "prefix_xxxx_dp.ps"
    (dot-plot), "prefix_xxxx_aln.ps" (annotated alignment), etc. where xxxx is
    the alignment number beginning with the second alignment in the input. Use
    this setting in conjunction with the :option:`--continuous-ids` flag to assign IDs
    beginning with the first input alignment.

.. option:: --id-delim=CHAR

    Change the delimiter between prefix and increasing number for automatically generated IDs (as used in output file names).


    *(default="_")*


    This parameter can be used to change the default delimiter "_" between the
    prefix string and the increasing number for automatically generated ID.

.. option:: --id-digits=INT

    Specify the number of digits of the counter in automatically generated alignment IDs.


    *(default="4")*


    When alignments IDs are automatically generated, they receive an increasing
    number, starting with 1. This number will always be left-padded by leading
    zeros, such that the number takes up a certain width. Using this parameter,
    the width can be specified to the users need. We allow numbers in the range
    [1:18].

.. option:: --id-start=LONG

    Specify the first number in automatically generated alignment IDs.


    *(default="1")*


    When alignment IDs are automatically generated, they receive an increasing
    number, usually starting with 1. Using this parameter, the first number can
    be specified to the users requirements. Note: negative numbers are not
    allowed.
    Note: Setting this parameter implies continuous alignment IDs, i.e. it
    activates the :option:`--continuous-ids` flag.

.. option:: --filename-delim=CHAR

    Change the delimiting character used in sanitized filenames.


    *(default="ID-delimiter")*


    This parameter can be used to change the delimiting character used while
    sanitizing filenames, i.e. replacing invalid characters. Note, that the
    default delimiter ALWAYS is the first character of the "ID delimiter" as
    supplied through the :option:`--id-delim` option. If the delimiter is a whitespace
    character or empty, invalid characters will be simply removed rather than
    substituted. Currently, we regard the following characters as illegal for use
    in filenames: backslash ``\``, slash ``/``, question mark ``?``, percent sign ``%``,
    asterisk ``*``, colon ``:``, pipe symbol ``|``, double quote ``"``, triangular
    brackets ``<`` and ``>``.

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.
    The Minimum free energy (MFE) and a structure representative are calculated
    in any case.

.. option:: -L, --maxBPspan=INT

    Set the maximum allowed separation of a base pair to span. I.e. no pairs (i,j) with j-i>span will be allowed.


    *(default="70")*

.. option:: --threshold=DOUBLE

    Energy threshold in kcal/mol per nucleotide above which secondary structure hits are omitted in the output.


    *(default="-0.1")*

.. option:: -g, --gquad

    Incoorporate G-Quadruplex formation into the structure prediction algorithm.


    *(default=off)*

Structure Constraints:
^^^^^^^^^^^^^^^^^^^^^^



    Command line options to interact with the structure constraints feature of
    this program

.. option:: --shape=file1,file2

    Use SHAPE reactivity data to guide structure predictions.


    Multiple shapefiles for the individual sequences in the alignment may be
    specified  as a comma separated list. An optional association of particular
    shape files to a specific  sequence in the alignment can be expressed by
    prepending the sequence number to the filename,  e.g.
    "5=seq5.shape,3=seq3.shape" will assign the reactivity values from file
    seq5.shape to  the fifth sequence in the alignment, and the values from file
    seq3.shape to sequence 3. If  no assignment is specified, the reactivity
    values are assigned to corresponding sequences in  the order they are given.

.. option:: --shapeMethod=D[mX][bY]

    Specify the method how to convert SHAPE reactivity data to pseudo energy contributions.


    *(default="D")*


    Currently, the only data conversion method available is that of to Deigan et
    al 2009.  This method is the default and is recognized by a capital ``D`` in
    the provided parameter, i.e.:  :option:`--shapeMethod=`"D" is the default setting.
    The slope ``m`` and the intercept ``b`` can be set to a  non-default value if
    necessary. Otherwise m=1.8 and b=-0.6 as stated in the paper mentionen
    before.  To alter these parameters, e.g. m=1.9 and b=-0.7, use a  parameter
    string like this: :option:`--shapeMethod=`"Dm1.9b-0.7". You may also provide only one
    of the two  parameters like: :option:`--shapeMethod=`"Dm1.9" or
    :option:`--shapeMethod=`"Db-0.7".

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

.. option:: --cfactor=DOUBLE

    Set the weight of the covariance term in the energy function


    *(default="1.0")*

.. option:: --nfactor=DOUBLE

    Set the penalty for non-compatible sequences in the covariance term of the energy function


    *(default="1.0")*

.. option:: -R, --ribosum_file=ribosumfile

    use specified Ribosum Matrix instead of normal


    energy model.


    Matrixes to use should be 6x6 matrices, the order of the terms is ``AU``, ``CG``,
    ``GC``, ``GU``, ``UA``, ``UG``.

.. option:: -r, --ribosum_scoring

    use ribosum scoring matrix. *(default=off)*


    The matrix is chosen according to the minimal and maximal pairwise identities
    of the sequences in the file.

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

Plotting:
^^^^^^^^^



    Command line options for changing the default behavior of structure layout
    and pairing probability plots

.. option:: --aln-EPS[=prefix]

    Produce colored and structure annotated subalignment for each hit.


    The default file name used for the output is "aln_start_end.eps" where
    "start" and "end" denote the first and last column of the subalignment
    relative to the input (1-based). Users may change the filename to
    "prefix_aln_start_end.eps" by specifying the prefix as optional argument.
    Files will be create in the current directory. Note: Any special characters
    in the prefix will be replaced by the filename delimiter, hence there is no
    way to pass an entire directory path through this option yet. (See also the
    "--filename-delim" parameter)

.. option:: --aln-EPS-cols=INT

    Number of columns in colored EPS alignment output.


    *(default="60")*


    A value less than 1 indicates that the output should not be wrapped at all.

.. option:: --aln-EPS-ss[=prefix]

    Produce colored consensus secondary structure plots in PostScript format.


    The default file name used for the output is "ss_start_end.eps" where
    "start" and "end" denote the first and last column of the subalignment
    relative to the input (1-based). Users may change the filename to
    "prefix_ss_start_end.eps" by specifying the prefix as optional argument.
    Files will be create in the current directory. Note: Any special characters
    in the prefix will be replaced by the filename delimiter, hence there is no
    way to pass an entire directory path through this option yet. (See also the
    "--filename-delim" parameter)

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

I.L. Hofacker, B. Priwitzer, and P.F. Stadler (2004),
"Prediction of Locally Stable RNA Secondary Structures for Genome-Wide Surveys",
Bioinformatics: 20, pp 186-190

Stephan H. Bernhart, Ivo L. Hofacker, Sebastian Will, Andreas R. Gruber, and Peter F. Stadler (2008),
"RNAalifold: Improved consensus structure prediction for RNA alignments",
BMC Bioinformatics: 9, pp 474


*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Ivo L Hofacker, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
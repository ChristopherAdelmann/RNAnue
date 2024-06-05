#######
RNAplot
#######

:program:`RNAplot` - manual page for RNAplot 2.6.4

Synopsis
--------

.. code:: bash

    RNAplot [OPTIONS] [<input0>] [<input1>]...

DESCRIPTION
-----------

RNAplot 2.6.4

Draw RNA Secondary Structures

The program reads (aligned) RNA sequences and structures in the format as
produced by RNAfold or Stockholm 1.0 and produces drawings of the secondary
structure graph.
Coordinates for the structure graphs are produced using either E. Bruccoleri's
naview routines, or a simple radial layout method.
For aligned sequences and consensus structures (:option:`--msa` option) the graph may be
annotated by covariance information. Additionally, a color-annotated EPS
alignment figure can be produced, similar to that obtained by RNAalifold and
RNALalifold.
If the sequence was preceded by a FASTA header, or if the multiple sequence
alignment contains an ID field, these IDs will be taken as names for the output
file(s): "name_ss.ps" and "name_aln.ps". Otherwise "rna.ps" and
"aln.ps" will be used. This behavior may be over-ruled by explicitly setting
a filename prefix using the :option:`--auto-id` option.
Existing files of the same name will be overwritten.

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

.. option:: -i, --infile=<filename>

    Read a file instead of reading from stdin.


    The default behavior of RNAplot is to read input from stdin or the file(s)
    that follow(s) the RNAplot command. Using this parameter the user can specify
    input file names where data is read from. Note, that any additional files
    supplied to RNAplot are still processed as well.

.. option:: -a, --msa

    Input is multiple sequence alignment in Stockholm 1.0 format. *(default=off)*


    Using this flag indicates that the input is a multiple sequence alignment
    (MSA) instead of (a) single sequence(s). Note, that only STOCKHOLM format
    allows one to specify a consensus structure. Therefore, this is the only
    supported MSA format for now!

.. option:: --mis

    Output "most informative sequence" instead of simple consensus *(default=off)*


    For each column of the alignment output this is the set of nucleotides with
    frequency greater than average in IUPAC notation.

.. option:: -j, --jobs[=number]

    Split batch input into jobs and start processing in parallel using multiple threads. *(default="0")*


    Default processing of input data is performed in a serial fashion, i.e. one
    sequence at a time. Using this switch, a user can instead start the
    computation for many sequences in the input in parallel. RNAplot will create
    as many parallel computation slots as specified and assigns input sequences
    of the input file(s) to the available slots. Note, that this increases memory
    consumption since input alignments have to be kept in memory until an empty
    compute slot is available and each running job requires its own dynamic
    programming matrices. A value of 0 indicates to use as many parallel threads
    as computation cores are available.

.. option:: -o, --output-format=ps|gml|xrna|svg

    Specify output format. *(default="ps")*


    Available formats are: PostScript (``ps``), Graph Meta Language (``gml``),
    Scalable Vector Graphics (``svg``), and XRNA save file (``xrna``). Output
    filenames will end in ".ps" ".gml" ".svg" ".ss", respectively.

.. option:: --pre=string

    Add annotation macros to postscript file, and add the postscript code in "string" just before the code to draw the structure. This is an easy way to add annotation.

.. option:: --post=string

    Same as :option:`--pre` but in contrast to adding the annotation macros. E.g to mark position 15 with circle use :option:`--post=`"15 cmark".

.. option:: --auto-id

    Automatically generate an ID for each sequence. *(default=off)*


    The default mode of RNAfold is to automatically determine an ID from the
    input sequence data if the input file format allows to do that. Sequence IDs
    are usually given in the FASTA header of input sequences. If this flag is
    active, RNAfold ignores any IDs retrieved from the input and automatically
    generates an ID for each sequence. This ID consists of a prefix and an
    increasing number. This flag can also be used to add a FASTA header to the
    output even if the input has none.

.. option:: --id-prefix=STRING

    Prefix for automatically generated IDs (as used in output file names).


    *(default="sequence")*


    If this parameter is set, each sequence will be prefixed with the provided
    string. Hence, the output files will obey the following naming scheme:
    "prefix_xxxx_ss.ps" (secondary structure plot), "prefix_xxxx_dp.ps"
    (dot-plot), "prefix_xxxx_dp2.ps" (stack probabilities), etc. where xxxx is
    the sequence number. Note: Setting this parameter implies :option:`--auto-id`.

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
    [1:18]. This option implies :option:`--auto-id`.

.. option:: --id-start=LONG

    Specify the first number in automatically generated IDs.


    *(default="1")*


    When sequence IDs are automatically generated, they receive an increasing
    number, usually starting with 1. Using this parameter, the first number can
    be specified to the users requirements. Note: negative numbers are not
    allowed.
    Note: Setting this parameter implies to ignore any IDs retrieved from the
    input data, i.e. it activates the :option:`--auto-id` flag.

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

.. option:: --filename-full

    Use full FASTA header to create filenames. *(default=off)*


    This parameter can be used to deactivate the default behavior of limiting
    output filenames to the first word of the sequence ID. Consider the following
    example: An input with FASTA header ``>NM_0001 Homo Sapiens some gene`` usually
    produces output files with the prefix "NM_0001" without the additional data
    available in the FASTA header, e.g. "NM_0001_ss.ps" for secondary structure
    plots. With this flag set, no truncation of the output filenames is done,
    i.e. output filenames receive the full FASTA header data as prefixes. Note,
    however, that invalid characters (such as whitespace) will be substituted by
    a delimiting character or simply removed, (see also the parameter option
    :option:`--filename-delim`).

Plotting:
^^^^^^^^^



    Command line options for changing the default behavior of structure layout
    and pairing probability plots

.. option:: --covar

    Annotate covariance of base pairs in consensus structure.


    *(default=off)*

.. option:: --aln

    Produce a colored and structure annotated alignment in PostScript format in the file "aln.ps" in the current directory.


    *(default=off)*

.. option:: --aln-EPS-cols=INT

    Number of columns in colored EPS alignment output.


    *(default="60")*


    A value less than 1 indicates that the output should not be wrapped at all.

.. option:: -t, --layout-type=INT

    Choose the plotting layout algorithm. *(default="1")*


    Select the layout algorithm that computes the nucleotide coordinates.
    Currently, the following algorithms are available:


    ``0``: simple radial layout


    ``1``: Naview layout (Bruccoleri et al. 1988)


    ``2``: circular layout


    ``3``: RNAturtle (Wiegreffe et al. 2018)


    ``4``: RNApuzzler (Wiegreffe et al. 2018)

.. option:: --noOptimization

    Disable the drawing space optimization of RNApuzzler.


    *(default=off)*

.. option:: --ignoreExteriorIntersections

    Ignore intersections with the exterior loop


    within the RNA-tree.


    *(default=off)*

.. option:: --ignoreAncestorIntersections

    Ignore ancestor intersections within the


    RNA-tree.


    *(default=off)*

.. option:: --ignoreSiblingIntersections

    Ignore sibling intersections within the


    RNA-tree.


    *(default=off)*

.. option:: --allowFlipping

    Allow flipping of exterior loop branches to resolve exterior branch intersections.


    *(default=off)*

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


Ivo L Hofacker, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
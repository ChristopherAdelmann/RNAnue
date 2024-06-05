########
RNALfold
########

:program:`RNALfold` - manual page for RNALfold 2.6.4

Synopsis
--------

.. code:: bash

    RNALfold [OPTION]...

DESCRIPTION
-----------

RNALfold 2.6.4

calculate locally stable secondary structures of RNAs

Compute locally stable RNA secondary structure with a maximal base pair span.
For a sequence of length n and a base pair span of L the algorithm uses only
O(n+L*L) memory and O(n*L*L) CPU time. Thus it is practical to "scan" very
large genomes for short RNA structures.
Output consists of a list of secondary structure components of size <= L, one
entry per line. Each output line contains the predicted local structure its
energy in kcal/mol and the starting position of the local structure.

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

.. option:: -i, --infile=filename

    Read a file instead of reading from stdin


    The default behavior of RNALfold is to read input from stdin. Using this
    parameter the user can specify an input file name where data is read from.

.. option:: -o, --outfile[=filename]

    Print output to file instead of stdout.


    This option may be used to write all output to output files rather than
    printing to stdout. The number of output files created for batch input
    (multiple sequences) depends on three conditions: (i) In case an optional
    filename is given as parameter argument, a single file with the specified
    filename will be written into. If the optional argument is omitted, (ii)
    FASTA input or an active :option:`--auto-id` switch will write to multiple files that
    follow the naming scheme "prefix.lfold". Here, "prefix" is taken from the
    sequence id as specified in the FASTA header. Lastly, (iii) single-line
    sequence input without FASTA header will be written to a single file
    "RNALfold_output.lfold". In case an output file already exists, any output
    of the program will be appended to it.
    Since the filename argument is optional, it must immediately follow the short
    option flag to not be mistaken as new parameter to the program. For instance
    \``-ornafold.out\`` will write to a file "rnafold.out".
    Note: Any special characters in the filename will be replaced by the filename
    delimiter, hence there is no way to pass an entire directory path through
    this option yet. (See also the "--filename-delim" parameter)

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

.. option:: --auto-id

    Automatically generate an ID for each sequence. *(default=off)*


    The default mode of RNALfold is to automatically determine an ID from the
    input sequence data if the input file format allows to do that. Sequence IDs
    are usually given in the FASTA header of input sequences. If this flag is
    active, RNALfold ignores any IDs retrieved from the input and automatically
    generates an ID for each sequence. This ID consists of a prefix and an
    increasing number. This flag can also be used to add a FASTA header to the
    output even if the input has none.

.. option:: --id-prefix=STRING

    Prefix for automatically generated IDs (as used in output file names).


    *(default="sequence")*


    If this parameter is set, each sequence will be prefixed with the provided
    string. Hence, the output files will obey the following naming scheme:
    "prefix_xxxx.lfold" where xxxx is the sequence number. Note: Setting this
    parameter implies :option:`--auto-id`.

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

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.
    The Minimum free energy (MFE) and a structure representative are calculated
    in any case.

.. option:: -L, --span=INT

    Set the maximum distance between any two pairing nucleotides.


    *(default="150")*


    This option specifies the window length L and therefore the upper limit for
    the distance between the bases i and j of any pair (i, j), i.e. (j - i + 1)
    <= L.

.. option:: -z, --zscore[=DOUBLE]

    Limit the output to predictions with a Z-score below a threshold.


    *(default="-2")*


    This option activates z-score regression using a trained SVM. Any predicted
    structure that exceeds the specified threshold will be ommited from the
    output.
    Since the Z-score threshold is given as a negative number, it must
    immediately preceed the short option to not be mistaken as a separate
    argument, e.g. :option:`-z-2`.9 sets the threshold to a value of :option:`-2`.9

.. option:: --zscore-pre-filter

    Apply the z-score filtering in the forward recursions.


    *(default=off)*


    The default mode of z-score filtering considers the entire structure space to
    decide whether or not a locally optimal structure at any position i is
    reported or not. When using this post-filtering step, however, alternative
    locally optimal structures


    starting at i with higher energy but lower z-score can be easily missed. The


    pre-filter


    option restricts the structure space already in the forward recursions, such


    that


    only optimal solution among those candidates that satisfy the z-score


    threshold are considered. Therefore, good results according to the z-score
    threshold criterion are less likely to be superseded by results with better
    energy but worse z-score. Note, that activating this switch results in higher
    computation time which scales linear in the window length.

.. option:: --zscore-report-subsumed

    Report subsumed structures if their z-score is less than that of the enclosing structure.


    *(default=off)*


    In default mode, RNALfold only reports locally optimal structures if they are
    no constituents of another, larger structure with less free energy. In
    z-score mode, however, such a larger structure may have a higher z-score,
    thus may be less informative than the smaller substructure. Using this switch
    activates reporting both, the smaller and the larger structure if the z-score
    of the smaller is lower than that of the larger.

.. option:: -b, --backtrack-global

    Backtrack a global MFE structure. *(default=off)*


    Instead of just reporting the locally stable secondary structure a global MFE
    structure can be constructed that only consists of locally optimal
    substructures. This switch activates a post-processing step that takes the
    locally optimal structures to generate the global MFE structure which
    constitutes the MFE value reported in the last line. The respective global
    MFE structure is printed just after the inut sequence part on the last line,
    preceding the global MFE score.
    Note, that this option implies :option:`-o`/--outfile since the locally optimal
    structures must be read after the regular prediction step! Also note, that
    using this option in combination with :option:`-z`/--zscore implies :option:`--zscore-pre-filter`
    to ensure proper construction of the global MFE structure!

.. option:: -g, --gquad

    Incoorporate G-Quadruplex formation into the structure prediction algorithm.


    *(default=off)*

Structure Constraints:
^^^^^^^^^^^^^^^^^^^^^^



    Command line options to interact with the structure constraints feature of
    this program

.. option:: --shape=filename

    Use SHAPE reactivity data to guide structure predictions.

.. option:: --shapeMethod=method

    Select SHAPE reactivity data incorporation strategy.


    *(default="D")*


    The following methods can be used to convert SHAPE reactivities into pseudo
    energy contributions.


    ``D``: Convert by using the linear equation according to Deigan et al 2009.


    Derived pseudo energy terms will be applied for every nucleotide involved in
    a stacked pair. This method is recognized by a capital ``D`` in the provided
    parameter, i.e.: :option:`--shapeMethod=`"D" is the default setting. The slope ``m``
    and the intercept ``b`` can be set to a non-default value if necessary,
    otherwise m=1.8 and b=-0.6. To alter these parameters, e.g. m=1.9 and b=-0.7,
    use a parameter string like this: :option:`--shapeMethod=`"Dm1.9b-0.7". You may also
    provide only one of the two parameters like: :option:`--shapeMethod=`"Dm1.9" or
    :option:`--shapeMethod=`"Db-0.7".


    ``Z``: Convert SHAPE reactivities to pseudo energies according to Zarringhalam


    et al 2012. SHAPE reactivities will be converted to pairing probabilities by
    using linear mapping. Aberration from the observed pairing probabilities will
    be penalized during the folding recursion. The magnitude of the penalties can
    affected by adjusting the factor beta (e.g. :option:`--shapeMethod=`"Zb0.8").


    ``W``: Apply a given vector of perturbation energies to unpaired nucleotides


    according to Washietl et al 2012. Perturbation vectors can be calculated by
    using RNApvmin.

.. option:: --shapeConversion=method

    Select method for SHAPE reactivity conversion.


    *(default="O")*


    This parameter is useful when dealing with the SHAPE incorporation according
    to Zarringhalam et al. The following methods can be used to convert SHAPE
    reactivities into the probability for a certain nucleotide to be unpaired.


    ``M``: Use linear mapping according to Zarringhalam et al.
    ``C``: Use a cutoff-approach to divide into paired and unpaired nucleotides
    (e.g. "C0.25")
    ``S``: Skip the normalizing step since the input data already represents
    probabilities for being unpaired rather than raw reactivity values
    ``L``: Use a linear model to convert the reactivity into a probability for
    being unpaired (e.g. "Ls0.68i0.2" to use a slope of 0.68 and an intercept
    of 0.2)
    ``O``: Use a linear model to convert the log of the reactivity into a
    probability for being unpaired (e.g. "Os1.6i-2.29" to use a slope of 1.6
    and an intercept of :option:`-2`.29)

.. option:: --commands=filename

    Read additional commands from file


    Commands include hard and soft constraints, but also structure motifs in
    hairpin and interior loops that need to be treeted differently. Furthermore,
    commands can be set for unstructured and structured domains.

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

.. option:: -m, --modifications[=STRING]

    Allow for modified bases within the RNA sequence string.


    *(default="7I6P9D")*


    Treat modified bases within the RNA sequence differently, i.e. use
    corresponding energy corrections and/or pairing partner rules if available.
    For that, the modified bases in the input sequence must be marked by their
    corresponding one-letter code. If no additional arguments are supplied, all
    available corrections are performed. Otherwise, the user may limit the
    modifications to a particular subset of modifications, resp. one-letter
    codes, e.g. :option:`-mP6` will only correct for pseudouridine and m6A bases.


    Currently supported one-letter codes and energy corrections are:


    ``7``: 7-deaza-adenonsine (7DA)


    ``I``: Inosine


    ``6``: N6-methyladenosine (m6A)


    ``P``: Pseudouridine


    ``9``: Purine (a.k.a. nebularine)


    ``D``: Dihydrouridine

.. option:: --mod-file=STRING

    Use additional modified base data from JSON file.

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

I.L. Hofacker, B. Priwitzer, and P.F. Stadler (2004),
"Prediction of Locally Stable RNA Secondary Structures for Genome-Wide Surveys",
Bioinformatics: 20, pp 186-190


*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Ivo L Hofacker, Peter F Stadler, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.

SEE ALSO
--------


RNAplfold(1) RNALalifold(1)
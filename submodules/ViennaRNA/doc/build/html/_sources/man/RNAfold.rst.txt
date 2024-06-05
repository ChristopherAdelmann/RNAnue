#######
RNAfold
#######

:program:`RNAfold` - manual page for RNAfold 2.6.4

Synopsis
--------

.. code:: bash

    RNAfold [OPTIONS] [<input0.fa>] [<input1.fa>]...

DESCRIPTION
-----------

RNAfold 2.6.4

Calculate minimum free energy secondary structures and partition function of
RNAs

The program reads RNA sequences, calculates their minimum free energy (mfe)
structure and prints the mfe structure in bracket notation and its free energy.
If not specified differently using commandline arguments, input is accepted
from stdin or read from an input file, and output printed to stdout. If the :option:`-p`
option was given it also computes the partition function (pf) and base pairing
probability matrix, and prints the free energy of the thermodynamic ensemble,
the frequency of the mfe structure in the ensemble, and the ensemble diversity
to stdout.

It also produces PostScript files with plots of the resulting secondary
structure graph and a "dot plot" of the base pairing matrix.
The dot plot shows a matrix of squares with area proportional to the pairing
probability in the upper right half, and one square for each pair in the
minimum free energy structure in the lower left half. For each pair i-j with
probability p>10E-6 there is a line of the form

i  j  sqrt(p)  ubox

in the PostScript file, so that the pair probabilities can be easily extracted.

Sequences may be provided in a simple text format where each sequence occupies
a single line. Output files are named "rna.ps" and "dot.ps". Existing files
of the same name will be overwritten.

It is also possible to provide sequence data in FASTA format. In this case, the
first word of the FASTA header will be used as prefix for output file names.
PostScript files "prefix_ss.ps" and "prefix_dp.ps" are produced for the
structure and dot plot, respectively. Note, however, that once FASTA input was
provided all following sequences must be in FASTA format too.

To avoid problems with certain operating systems and/or file systems the prefix
will ALWAYS be sanitized! This step substitutes any special character in the
prefix by a filename delimiter. See the :option:`--filename-delim` option to change the
delimiting character according to your requirements.

The program will continue to read new sequences until a line consisting of the
single character ``@`` or an end of file (EOF) condition is encountered.

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

    Read a file instead of reading from stdin.


    The default behavior of RNAfold is to read input from stdin or the file(s)
    that follow(s) the RNAfold command. Using this parameter the user can specify
    input file names where data is read from. Note, that any additional files
    supplied to RNAfold are still processed as well.

.. option:: -o, --outfile[=filename]

    Print output to file instead of stdout.


    This option may be used to write all output to output files rather than
    printing to stdout. The default filename is "RNAfold_output.fold" if no
    FASTA header precedes the input sequences and the :option:`--auto-id` feature is
    inactive. Otherwise, output files with the scheme "prefix.fold" are
    generated, where the "prefix" is taken from the sequence id, e.g. the FASTA
    header. The user may specify a single output file name for all data generated
    from the input by supplying a filename as argument following immediately
    after this parameter.
    In case a file with the same filename already exists, any output of the
    program will be appended to it. Note: Any special characters in the filename
    will be replaced by the filename delimiter, hence there is no way to pass an
    entire directory path through this option (yet). (See also the
    "--filename-delim" parameter)

.. option:: -j, --jobs[=number]

    Split batch input into jobs and start processing in parallel using multiple threads. A value of 0 indicates to use as many parallel threads as computation cores are available.


    *(default="0")*


    Default processing of input data is performed in a serial fashion, i.e. one
    sequence at a time. Using this switch, a user can instead start the
    computation for many sequences in the input in parallel. RNAfold will create
    as many parallel computation slots as specified and assigns input sequences
    of the input file(s) to the available slots. Note, that this increases memory
    consumption since input alignments have to be kept in memory until an empty
    compute slot is available and each running job requires its own dynamic
    programming matrices.

.. option:: --unordered

    Do not try to keep output in order with input while parallel processing is in place.


    *(default=off)*


    When parallel input processing (:option:`--jobs` flag) is enabled, the order in which
    input is processed depends on the host machines job scheduler. Therefore, any
    output to stdout or files generated by this program will most likely not
    follow the order of the corresponding input data set. The default of RNAfold
    is to use a specialized data structure to still keep the results output in
    order with the input data. However, this comes with a trade-off in terms of
    memory consumption, since all output must be kept in memory for as long as no
    chunks of consecutive, ordered output are available. By setting this flag,
    RNAfold will not buffer individual results but print them as soon as they
    have been computated.

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

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

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.
    The Minimum free energy (MFE) and a structure representative are calculated
    in any case.

.. option:: -p, --partfunc[=INT]

    Calculate the partition function and base pairing probability matrix.


    *(default="1")*


    In addition to the MFE structure we print a coarse representation of the pair
    probabilities in form of a pseudo bracket notation followed by the ensemble
    free energy. This notation makes use of the letters ``.``, ``,``, ``|``, ``{``, ``}``,
    ``(``, and ``)`` denoting bases that are essentially unpaired, weakly paired,
    strongly paired without preference, weakly upstream (downstream) paired, or
    strongly up- (down-)stream paired bases, respectively. On the next line the
    centroid structure as derived from the pair probabilities together with its
    free energy and distance to the ensemble is shown. Finally it prints the
    frequency of the mfe structure, and the structural diversity (mean distance
    between the structures in the ensemble).
    See the description of ``vrna_pf()`` and ``mean_bp_dist()`` and ``vrna_centroid()``
    in the RNAlib documentation for details.
    Note that unless you also specify :option:`-d2` or :option:`-d0`, the partition function and mfe
    calculations will use a slightly different energy model. See the discussion
    of dangling end options below.


    An additionally passed value to this option changes the behavior of partition
    function calculation:
    :option:`-p0` Calculate the partition function but not the pair probabilities, saving
    about 50% in runtime. This prints the ensemble free energy ``dG=-kT ln(Z)``.
    :option:`-p2` Compute stack probabilities, i.e. the probability that a pair ``(i,j)`` and
    the immediately interior pair ``(i+1,j-1)`` are formed simultaneously in
    addition to pair probabilities. A second postscript dot plot named
    "name_dp2.ps", or "dot2.ps" (if the sequence does not have a name), is
    produced that contains pair probabilities in the upper right half and stack
    probabilities in the lower left.

.. option:: --betaScale=DOUBLE

    Set the scaling of the Boltzmann factors. *(default="1.")*


    The argument provided with this option is used to scale the thermodynamic
    temperature in the Boltzmann factors independently from the temperature of
    the individual loop energy contributions. The Boltzmann factors then become
    ``exp(- dG/(kT*betaScale))`` where ``k`` is the Boltzmann constant, ``dG`` the free
    energy contribution of the state and ``T`` the absolute temperature.

.. option:: -S, --pfScale=DOUBLE

    In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows).


    *(default="1.07")*


    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for
    long sequences.

.. option:: --MEA[=gamma]

    Compute MEA (maximum expected accuracy) structure.


    *(default="1.")*


    The expected accuracy is computed from the pair probabilities: each base pair
    ``(i,j)`` receives a score ``2*gamma*p_ij`` and the score of an unpaired base is
    given by the probability of not forming a pair. The parameter gamma tunes the
    importance of correctly predicted pairs versus unpaired bases. Thus, for
    small values of gamma the MEA structure will contain only pairs with very
    high probability. Using :option:`--MEA` implies :option:`-p` for computing the pair
    probabilities.

.. option:: -c, --circ

    Assume a circular (instead of linear) RNA molecule.


    *(default=off)*

.. option:: --ImFeelingLucky

    Return exactly one stochastically backtracked structure.


    *(default=off)*


    This function computes the partition function and returns exactly one
    secondary structure stochastically sampled from the Boltzmann equilibrium
    according to its probability in the ensemble

.. option:: --bppmThreshold=cutoff

    Set the threshold/cutoff for base pair probabilities included in the postscript output.


    *(default="1e-5")*


    By setting the threshold the base pair probabilities that are included in the
    output can be varied. By default only those exceeding ``1e-5`` in probability
    will be shown as squares in the dot plot. Changing the threshold to any other
    value allows for increase or decrease of data.

.. option:: -g, --gquad

    Incoorporate G-Quadruplex formation into the structure prediction algorithm.


    *(default=off)*

Structure Constraints:
^^^^^^^^^^^^^^^^^^^^^^



    Command line options to interact with the structure constraints feature of
    this program

.. option:: --maxBPspan=INT

    Set the maximum base pair span.


    *(default="-1")*

.. option:: -C, --constraint[=filename]

    Calculate structures subject to constraints. *(default="")*


    The program reads first the sequence, then a string containing constraints on
    the structure encoded with the symbols:


    ``.`` (no constraint for this base)


    ``|`` (the corresponding base has to be paired


    ``x`` (the base is unpaired)


    ``<`` (base i is paired with a base j>i)


    ``>`` (base i is paired with a base j<i)


    and matching brackets ``(`` ``)`` (base i pairs base j)


    With the exception of ``|``, constraints will disallow all pairs conflicting
    with the constraint. This is usually sufficient to enforce the constraint,
    but occasionally a base may stay unpaired in spite of constraints. PF folding
    ignores constraints of type ``|``.

.. option:: --batch

    Use constraints for multiple sequences. *(default=off)*


    Usually, constraints provided from input file only apply to a single input
    sequence. Therefore, RNAfold will stop its computation and quit after the
    first input sequence was processed. Using this switch, RNAfold processes
    multiple input sequences and applies the same provided constraints to each of
    them.

.. option:: --canonicalBPonly

    Remove non-canonical base pairs from the structure constraint.


    *(default=off)*

.. option:: --enforceConstraint

    Enforce base pairs given by round brackets ``(`` ``)`` in structure constraint.


    *(default=off)*

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

.. option:: --motif=SEQUENCE,STRUCTURE,ENERGY

    Specify stabilizing energy of a ligand binding


    to a particular sequence/structure motif.


    Some ligands binding to RNAs require and/or induce particular sequence and
    structure motifs. For instance they bind to an interior loop, or small
    hairpin loop. If for such cases a binding free energy is known, the binding
    and therefore stabilizing effect of the ligand can be included in the folding
    recursions. Interior loop motifs are specified as concatenations of 5`` and 3``
    motif, separated by an ``&`` character.


    Energy contributions must be specified in kcal/mol.


    See the manpage for an example usage of this option.

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

Plotting:
^^^^^^^^^



    Command line options for changing the default behavior of structure layout
    and pairing probability plots

.. option:: --noPS

    Do not produce postscript drawing of the mfe structure.


    *(default=off)*

.. option:: --noDP

    Do not produce dot-plot postscript file containing base pair or stack probabilitities.


    *(default=off)*


    In combination with the :option:`-p` option, this flag turns-off creation of individual
    dot-plot files. Consequently, computed base pair probability output is
    omitted but centroid and MEA structure prediction is still performed.

.. option:: -t, --layout-type=INT

    Choose the layout algorithm. *(default="1")*


    Select the layout algorithm that computes the nucleotide coordinates.
    Currently, the following algorithms are available:


    ``0``: simple radial layout


    ``1``: Naview layout (Bruccoleri et al. 1988)


    ``2``: circular layout


    ``3``: RNAturtle (Wiegreffe et al. 2018)


    ``4``: RNApuzzler (Wiegreffe et al. 2018)

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

M. Zuker, P. Stiegler (1981),
"Optimal computer folding of large RNA sequences using thermodynamic and auxiliary information",
Nucl Acid Res: 9, pp 133-148

J.S. McCaskill (1990),
"The equilibrium partition function and base pair binding probabilities for RNA secondary structures",
Biopolymers: 29, pp 1105-1119

I.L. Hofacker & P.F. Stadler (2006),
"Memory Efficient Folding Algorithms for Circular RNA Secondary Structures",
Bioinformatics

A.F. Bompfuenewerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L. Hofacker, P.F. Stadler, S. Will (2007),
"Variations on {RNA} Folding and Alignment: Lessons from Benasque",
J. Math. Biol.

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

EXAMPLES
--------


Single line sequence input and calculation of partition function and MEA structure

.. code::

    $ RNAfold --MEA -d2 -p
    


The program will then prompt for sequence input. Using the example sequence
"CGACGTAGATGCTAGCTGACTCGATGC" and pressing ENTER the output of the program will be
similar to

.. code::

    CGACGUAGAUGCUAGCUGACUCGAUGC
    (((.((((.......)).)))))....
    minimum free energy =  -1.90 kcal/mol
    (((.((((.......))},})))....
    free energy of ensemble =  -2.86 kcal/mol
    (((.(.((.......))..)))).... {  0.80 d=2.81}
    (((.((((.......))).)))).... { -1.90 MEA=22.32}
    frequency of mfe structure in ensemble 0.20997; ensemble diversity 4.19
    



Here, the first line just repeats the sequence input. The second line contains a
MFE structure in dot bracket notation followed by the minimum free energy. After
this, the pairing probabilities for each nucleotide are shown in a pseudo dot-bracket
notation followed by the free energy of ensemble. The next two lines show the centroid
structure with its free energy and its distance to the ensemble as well as the MEA structure,
its free energy and the maximum expected accuracy, respectively. The last line finally
contains the frequency of the MFE representative in the complete ensemble of secondary
structures and the ensemble diversity. For further details about the calculation and
interpretation of the given output refer to the reference manual of RNAlib.

Since version 2.0 it is also possible to provide FASTA file sequence input. Assume
you have a file containing two sequences in FASTA format, e.g

.. code::

    $ cat sequences.fa
    >seq1
    CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
    GAACGAUCUAUAACACGACUUCACUCUU
    >seq2
    GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUA
    UAACACGACUUCACUCUU
    


In order to compute the MFE for the two sequences the user can use the following
command

.. code::

    $ RNAfold < sequences.fa
    


which would result in an output like this

.. code::

    >seq1
    CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
    .((.(((...((((..(((((........)))))))))...))).))................... ( -5.40)
    >seq2
    GAAUGACCCGAUAACCCCGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
    .......((((..............))))........................... ( -2.00)
    


CONSTRAINT EXAMPLES
-------------------


Secondary structure constraints may be given in addition to the sequence information, too.
Using the first sequence of the previous example and restricting the nucleotides of the
outermost helix to be unpaired, i.e. base pairs (2,47) and (3,46) the input file should
have the following form

.. code::

    $ cat sequence_unpaired.fa
    >seq1
    CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUG
    GAACGAUCUAUAACACGACUUCACUCUU
    .xx...................................
    .......xx...................
    


Calling RNAfold with the structure constraint option -C it shows the following result

.. code::

    $ RNAfold -C < sequence_unpaired.fa
    >seq1
    CGGCUCGCAACAGACCUAUUAGUUUUACGUAAUAUUUGGAACGAUCUAUAACACGACUUCACUCUU
    ....(((...((((..(((((........)))))))))...)))...................... ( -4.20)
    


This represents the minimum free energy and a structure representative of the RNA
sequence given that nucleotides 2,3,46 and 47 must not be involved in any base pair.
For further information about constrained folding refer to the details of the -C option
and the reference manual of RNAlib.

Since version 2.2 the ViennaRNA Package distinguishes hard and soft constraints.
As a consequence, structure predictions are easily amenable to a versatile set of constraints,
such as maximal base pair span, incorporation of SHAPE reactivity data, and RNA-ligand binding
to hairpin, or interior loop motifs.

*Restricting the maximal span of a base pair*

A convenience commandline option allows you to easily limit the distance (j - i + 1) between
two nucleotides i and j that form a basepair. For instance a limit of 600nt can be accomplished
using:

.. code::

    $ RNAfold --maxBPspan 600
    


*Guide structure prediction with SHAPE reactivity data*

Use SHAPE reactivity data to guide secondary structure prediction:

.. code::

    $ RNAfold --shape=reactivities.dat < sequence.fa
    


where the file reactivities.dat is a two column text file with sequence positions (1-based)
and normalized reactivity values (usually between 0 and 2. Missing values may be left out,
or assigned a negative score:

.. code::

    $ cat reactivities.dat
    9    -999       # No reactivity information
    10   -999
    11   0.042816   # normalized SHAPE reactivity
    12   0          # also a valid SHAPE reactivity
    15   0.15027    # Missing data for pos. 13-14
    ...
    42   0.16201
    


Note, that RNAfold will only process the first sequence in the input file, when provided
with SHAPE reactivity data!

*Complex structure constraints and grammar extensions*

Structure constraints beyond those that can be expressed with a pseudo-dot bracket notation
may be provided in a so-called command file:

.. code::

    $ RNAfold --commands=constraints.txt < sequence.fa
    


The command file syntax is a generalization of constraints as used in
UNAfold/mfold. Each line starts with a one or two letter command followed
by command parameters. For structure constraints, this amounts to a single
command character followed by three or four numbers. In addition, optional
auxiliary modifier characters may be used to limit the constraint to specific
loop types. For base pair specific constraints,
we currently distinguish pairs in exterior loops (E), closing pairs of hairpin
loops (H), closing (I) and enclosed (i) pairs of interior loops, and closing (M)
and enclosed (m) pairs of multibranch loops. Nucleotide-wise constraints may be
limited to their loop context using the corresponding uppercase characters. The
default is to apply a constraint to all (A) loop types. Furthermore, pairing
constraints for single nucleotides may be limited to upstream (U), or downstream (D)
orientation. The command file specification is as follows:

.. code::

    F i 0 k   [TYPE] [ORIENTATION] # Force nucleotides i...i+k-1 to be paired
    F i j k   [TYPE] # Force helix of size k starting with (i,j) to be formed
    P i 0 k   [TYPE] # Prohibit nucleotides i...i+k-1 to be paired
    P i j k   [TYPE] # Prohibit pairs (i,j),...,(i+k-1,j-k+1)
    P i-j k-l [TYPE] # Prohibit pairing between two ranges
    C i 0 k   [TYPE] # Nucleotides i,...,i+k-1 must appear in context TYPE
    C i j k          # Remove pairs conflicting with (i,j),...,(i+k-1,j-k+1)
    E i 0 k e        # Add pseudo-energy e to nucleotides i...i+k-1
    E i j k e        # Add pseudo-energy e to pairs (i,j),...,(i+k-1,j-k+1)
    UD m e    [LOOP] # Add ligand binding to unstructured domains with motif
    # m and binding free energy e
    
    # [LOOP]        = { E, H, I, M, A }
    # [TYPE]        = [LOOP] + { i, m }
    # [ORIENTATION] = { U, D }
    


Again, RNAfold by default only processes the first sequence in the input sequence
when provided with constraints in a command file. To apply the exact same constraints
to each of the input sequences in a multi FASTA file, use the batch mode commandline
option:

.. code::

    $ RNAfold --constraint=constraints.txt --batch < sequences.fa
    


*Ligand binding contributions to specific hairpin/interior loop motifs*

A convenience function allows one to specify a hairping/interior loop motif where a ligand
is binding with a particular binding free energy dG.
Here is an example that adds a theophylline binding motif. Free energy contribution of
this motif of dG=-9.22kcal/mol is derived from k_d=0.32umol/l, taken from Jenison et al.
1994. Although the structure motif consists of a symmetric interior loop of size 6,
followed by a small helix of 3 basepairs, and a bulge of 3 nucleotides, the entire
structure can still be represented by one interior loop.
See the below mofif description where the ``&`` character splits the motif into a 5' and
a 3' part. The first line gives the sequences motif, the second line shows the actual
structure motif of the aptamer pocket, and the third line is the interior loop motif
that fully encapsulates the theophylline aptamer:

.. code::

    GAUACCAG&CCCUUGGCAGC
    (...((((&)...)))...)
    (......(&).........)
    


To use the above information in the folding recursions of RNAfold, one only needs to
provide the motif itself, and binding free energy:

.. code::

    $ RNAfold --motif="GAUACCAG&CCCUUGGCAGC,(...((((&)...)))...),-9.22" < sequences.fa
    


Adding the --verbose option to the above call of RNAfold also prints the sequence
position of each motif found in the MFE structure. In case interior-loop like motifs
are provided, two intervals are printed denoting the 5`` and 3`` part, respectively.

*Ligand binding contributions to unpaired segments of the RNA structure*

The extension of the RNA folding grammar with unstructured domains allows for an easy
incorporation of ligands that bind to unpaired stretches of an RNA structure. To
model such interactions only two parameters are required: (i) a sequence motif in
IUPAC notation that specifies where the ligand binds to, and (ii) a binding free
energy that can be derived from the association/dissociation constant of the ligand.
With these two parameters in hand, the modification of RNAfold to include the competition
of regular intramolecular base pairing and ligand interaction is as easy as writing
a simple command file of the form:

.. code::

    UD m e    [LOOP]
    


where m is the motif string in upper-case IUPAC notation, and e the binding free energy
in kcal/mol and optional loop type restriction [LOOP]. See also the command file specification as defined above.

For instance, having a protein with a 4-nucleotide footprint binding ``AAAA``, a
binding free energy e = -5.0 kcal/mol, and a binding restriction to exterior- and
multibranch loops results in a command file:

.. code::

    $ cat commands.txt
    UD AAAA -5.0  ME
    


and the corresponding call to RNAfold to compute MFE and equilibrium probabilities becomes:

.. code::

    $ RNAfold --commands=commands.txt -p < sequence.fa
    


The resulting MFE plot will be annotated to display the binding site(s) of the ligand,
and the base pair probability dot-plot is extended to include the probability that
a particular nucleotide is bound by the ligand.



POST-TRANSCRIPTIONAL MODIFICATION EXAMPLES
------------------------------------------


Many RNA molecules harbor (post-transcriptional) modifications. These modified base
often change the pairing behavior or energy contribution for the loops they are part of.
To accommodate for that effect (to a certain degree) one may use additional correcting
energy parameters for loops with the respective modified bases. In literature, a few
stacking- and some terminal mismatch energies can be found. Some of them are already provided
within the ViennaRNA Package. The --modification and --mod-file command line parameters can be
used to apply these parameters in the predictions. While the former allows one to select a
subset of implemented modified base corrections, the latter enables the prediction programs
to read energy parameters for modified bases from a user-provided JSON file.

Consider, for instance, the following tRNA sequence with dihydrouridines and pseudouridines
annotated by their respective one-letter codes D and P:

.. code::

    $ cat tRNAphe.fa
    >tRNAphe
    GCCGAAAUAGCUCAGDDGGGAGAGCGPPAGACUGAAGAPCUAAAGGDCCCUGGUPCGAUCCCGGGUUUCGGCACCA
    


Now, a prediction that includes support for the destabilizing effect of D and the stabilizing
effects of P within base pair stacks can be done as follows:

.. code::

    $ RNAfold --modifications=DP < tRNAphe.fa
    >tRNAphe
    GCCGAAAUAGCUCAGDDGGGAGAGCGPPAGACUGAAGAPCUAAAGGDCCCUGGUPCGAUCCCGGGUUUCGGCACCA
    (((((((..((((........)))).(((((.......)))))....(((.((......)).)))))))))).... (-23.37)
    


AUTHOR
------


Ivo L Hofacker, Walter Fontana, Sebastian Bonhoeffer, Peter F Stadler, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
#########
RNAsubopt
#########

:program:`RNAsubopt` - manual page for RNAsubopt 2.6.4

Synopsis
--------

.. code:: bash

    RNAsubopt [OPTION]...

DESCRIPTION
-----------

RNAsubopt 2.6.4

calculate suboptimal secondary structures of RNAs

Reads RNA sequences from stdin and (in the default :option:`-e` mode) calculates all
suboptimal secondary structures within a user defined energy range above the
minimum free energy (mfe). It prints the suboptimal structures in dot-bracket
notation followed by the energy in kcal/mol to stdout. Be careful, the number
of structures returned grows exponentially with both sequence length and energy
range.

Alternatively, when used with the :option:`-p` option, RNAsubopt produces Boltzmann
weighted samples of secondary structures.

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


    The default behavior of RNAsubopt is to read input from stdin. Using this
    parameter the user can specify an input file name where data is read from.

.. option:: -o, --outfile[=filename]

    Print output to file instead of stdout.


    This option may be used to write all output to output files rather than
    printing to stdout. The default filename is "RNAsubopt_output.sub" if no
    FASTA header precedes the input sequences and the :option:`--auto-id` feature is
    inactive. Otherwise, output files with the scheme "prefix.sub" are
    generated, where the "prefix" is taken from the sequence id. The user may
    specify a single output file name for all data generated from the input by
    supplying an optional string as argument to this parameter. In case a file
    with the same filename already exists, any output of the program will be
    appended to it. Note: Any special characters in the filename will be replaced
    by the filename delimiter, hence there is no way to pass an entire directory
    path through this option yet. (See also the "--filename-delim" parameter)

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

.. option:: --auto-id

    Automatically generate an ID for each sequence. *(default=off)*


    The default mode of RNAsubopt is to automatically determine an ID from the
    input sequence data if the input file format allows to do that. Sequence IDs
    are usually given in the FASTA header of input sequences. If this flag is
    active, RNAsubopt ignores any IDs retrieved from the input and automatically
    generates an ID for each sequence. This ID consists of a prefix and an
    increasing number. This flag can also be used to add a FASTA header to the
    output even if the input has none.

.. option:: --id-prefix=STRING

    Prefix for automatically generated IDs (as used in output file names).


    *(default="sequence")*


    If this parameter is set, each sequences' FASTA id will be prefixed with the
    provided string. FASTA ids then take the form ">prefix_xxxx" where xxxx is
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
    available in the FASTA header, e.g. "NM_0001.sub". With this flag set, no
    truncation of the output filenames is performed, i.e. output filenames
    receive the full FASTA header data as prefixes. Note, however, that invalid
    characters (such as whitespace) will be substituted by a delimiting character
    or simply removed, (see also the parameter option :option:`--filename-delim`).

Algorithms:
^^^^^^^^^^^



    Select the algorithms which should be applied to the given RNA sequence(s).

.. option:: -e, --deltaEnergy=range

    Compute suboptimal structures with energy in a certain range of the optimum (kcal/mol).


    Default is calculation of mfe structure only.

.. option:: --deltaEnergyPost=range

    Only print structures with energy within range of the mfe after post reevaluation of energies.


    Useful in conjunction with :option:`-logML`, :option:`-d1` or :option:`-d3`: while the :option:`-e` option specifies
    the range before energies are re-evaluated, this option specifies the maximum
    energy after re-evaluation.

.. option:: -s, --sorted

    Sort the suboptimal structures by energy and lexicographical order.


    *(default=off)*


    Structures are first sorted by energy in ascending order. Within groups of
    the same energy, structures are then sorted in ascending in lexicographical
    order of their dot-bracket notation. See the :option:`--en-only` flag to deactivate
    this second step. Note that sorting is done in memory, thus it can easily
    lead to exhaution of RAM! This is especially true if the number of structures
    produced becomes large or the RNA sequence is rather long. In such cases
    better use an external sort method, such as UNIX "sort".

.. option:: --en-only

    Only sort structures by free energy. *(default=off)*


    In combination with :option:`--sorted`, this flag deactivates the second sorting
    criteria and sorts structures solely by their free energy instead of
    additionally sorting by lexicographic order in each energy band. This might
    save some time during the sorting process in situations where lexicographic
    order is not required.

.. option:: -p, --stochBT=number

    Randomly draw structures according to their probability in the Boltzmann ensemble.


    Instead of producing all suboptimals in an energy range, produce a random
    sample of suboptimal structures, drawn with probabilities equal to their
    Boltzmann weights via stochastic backtracking in the partition function. The
    :option:`-e` and :option:`-p` options are mutually exclusive.

.. option:: --stochBT_en=number

    Same as "--stochBT" but also print free energies and probabilities of the backtraced structures.

.. option:: --betaScale=DOUBLE

    Set the scaling of the Boltzmann factors. *(default="1.")*


    The argument provided with this option is used to scale the thermodynamic
    temperature in the Boltzmann factors independently from the temperature of
    the individual loop energy contributions. The Boltzmann factors then become
    ``exp(- dG/(kT*betaScale))`` where ``k`` is the Boltzmann constant, ``dG`` the free
    energy contribution of the state and ``T`` the absolute temperature.

.. option:: -N, --nonRedundant

    Enable non-redundant sampling strategy.


    *(default=off)*

.. option:: -S, --pfScale=DOUBLE

    In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows).


    *(default="1.07")*


    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for
    long sequences.

.. option:: -c, --circ

    Assume a circular (instead of linear) RNA molecule.


    *(default=off)*

.. option:: -D, --dos

    Compute density of states instead of secondary structures.


    *(default=off)*


    This option enables the evaluation of the number of secondary structures in
    certain energy bands around the MFE.

.. option:: -z, --zuker

    Compute Zuker suboptimals instead of all suboptimal structures within an energy band around the MFE.


    *(default=off)*

.. option:: -g, --gquad

    Incoorporate G-Quadruplex formation. *(default=off)*


    No support of G-quadruplex prediction for stochastic backtracking and
    Zuker-style suboptimals yet).

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
    sequence. Therefore, RNAsubopt will stop its computation and quit after the
    first input sequence was processed. Using this switch, RNAsubopt processes
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

.. option:: --logML

    Recompute energies of structures using a logarithmic energy function for multi-loops before output. *(default=off)*


    This option does not effect structure generation, only the energies that are
    printed out. Since logML lowers energies somewhat, some structures may be
    missing.

.. option:: --nsp=STRING

    Allow other pairs in addition to the usual AU,GC,and GU pairs.


    Its argument is a comma separated list of additionally allowed pairs. If the
    first character is a "-" then AB will imply that AB and BA are allowed
    pairs, e.g. :option:`--nsp=`"-GA"  will allow GA and AG pairs. Nonstandard pairs are
    given 0 stacking energy.

.. option:: --energyModel=INT

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

S. Wuchty, W. Fontana, I. L. Hofacker and P. Schuster (1999),
"Complete Suboptimal Folding of RNA and the Stability of Secondary Structures",
Biopolymers: 49, pp 145-165

M. Zuker (1989),
"On Finding All Suboptimal Foldings of an RNA Molecule",
Science 244.4900, pp 48-52

Y. Ding, and C.E. Lawrence (2003),
"A statistical sampling algorithm for RNA secondary structure prediction",
Nucleic Acids Research 31.24, pp 7280-7301

*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Ivo L Hofacker, Stefan Wuchty, Walter Fontana, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
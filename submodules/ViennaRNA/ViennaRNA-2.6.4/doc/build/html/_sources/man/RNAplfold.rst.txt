#########
RNAplfold
#########

:program:`RNAplfold` - manual page for RNAplfold 2.6.4

Synopsis
--------

.. code:: bash

    RNAplfold [OPTION]...

DESCRIPTION
-----------

RNAplfold 2.6.4

calculate locally stable secondary structure - pair probabilities

Computes local pair probabilities for base pairs with a maximal span of L. The
probabilities are averaged over all windows of size L that contain the base
pair. For a sequence of length n and a window size of L the algorithm uses only
O(n+L*L) memory and O(n*L*L) CPU time. Thus it is practical to "scan" very
large genomes for short stable RNA structures.


Output consists of a dot plot in postscript file, where the averaged pair probabilities
can easily be parsed and visually inspected.

The -u option makes i possible to compute the probability that a stretch of x consequtive
nucleotides is unpaired, which is useful for predicting possible binding sites. Again
this probability is averaged over all windows containing the region.

.B WARNING! Output format changed!!

The output is a plain text matrix containing on each line a position i followed by the
probability that i is unpaired, [i-1..i] is unpaired [i-2..i] is unpaired and so on to
the probability that [i-x+1..i] is unpaired.

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

.. option:: -c, --cutoff=FLOAT

    Report only base pairs with an average probability larger than ``cutoff`` in the dot plot.


    *(default="0.01")*

.. option:: -o, --print_onthefly

    Save memory by printing out everything during computation.


    *(default=off)*


    NOTE: activated per default for sequences over 1M bp.

.. option:: -O, --opening_energies

    Switch output from probabilities to their logarithms.


    *(default=off)*


    This is NOT exactly the mean energies needed to unfold the respective stretch
    of bases! (implies :option:`--ulength` option).

.. option:: --plex_output

    Create additional output files for RNAplex.


    *(default=off)*

.. option:: -b, --binaries

    Output accessibility profiles in binary format. *(default=off)*


    The binary files produced by RNAplfold do not need to be parsed by RNAplex,


    so that they are directly loaded into memory. This is useful when large
    sequences have to be searched for putative hybridization sites. Another
    advantage of the binary format is the 50% file size decrease.

.. option:: --noconv

    Do not automatically substitute nucleotide "T" with "U".


    *(default=off)*

.. option:: --auto-id

    Automatically generate an ID for each sequence. *(default=off)*


    The default mode of RNAplfold is to automatically determine an ID from the
    input sequence data if the input file format allows to do that. Sequence IDs
    are usually given in the FASTA header of input sequences. If this flag is
    active, RNAplfold ignores any IDs retrieved from the input and automatically
    generates an ID for each sequence. This ID consists of a prefix and an
    increasing number. This flag can also be used to add a FASTA header to the
    output even if the input has none.

.. option:: --id-prefix=STRING

    Prefix for automatically generated IDs (as used in output file names).


    *(default="sequence")*


    If this parameter is set, each sequences' FASTA id will be prefixed with the
    provided string. FASTA ids then take the form ">prefix_xxxx" where xxxx is
    the sequence number. Hence, the output files will obey the following naming
    scheme: "prefix_xxxx_dp.ps" (dot-plot), "prefix_xxxx_lunp" (unpaired
    probabilities), etc. Note: Setting this parameter implies :option:`--auto-id`.

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



    Select and change parameters of (additional) algorithms which should be
    included in the calculations.

.. option:: -W, --winsize=size

    Average the pair probabilities over windows of given size.


    *(default="70")*

.. option:: -L, --span=size

    Set the maximum allowed separation of a base pair to span.


    By setting the maximum base pair span no pairs (i,j) with j-i > span will be
    allowed. Defaults to winsize if parameter is omitted.

.. option:: -u, --ulength=length

    Compute the mean probability that regions of length 1 to a given length are unpaired.


    *(default="31")*


    Output is saved in a ``_lunp`` file.

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

    Specify "dangling end" model for bases adjacent to helices in free ends and multi-loops.


    *(default="2")*


    With :option:`-d2` dangling energies will be added for the bases adjacent to a helix on
    both sides in any case while :option:`-d0` ignores dangling ends altogether (mostly for
    debugging).

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

S. H. Bernhart, U. Mueckstein, and I.L. Hofacker (2011),
"RNA Accessibility in cubic time",
Algorithms Mol Biol. 6: 3.

S. H. Bernhart, I.L. Hofacker, and P.F. Stadler (2006),
"Local Base Pairing Probabilities in Large RNAs",
Bioinformatics: 22, pp 614-615

A.F. Bompfuenewerer, R. Backofen, S.H. Bernhart, J. Hertel, I.L. Hofacker, P.F. Stadler, S. Will (2007),
"Variations on RNA Folding and Alignment: Lessons from Benasque",
J. Math. Biol.


*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

AUTHOR
------


Stephan H Bernhart, Ivo L Hofacker, Peter F Stadler, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.

SEE ALSO
--------


RNALfold(1)
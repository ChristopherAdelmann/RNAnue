########
RNApvmin
########

:program:`RNApvmin` - manual page for RNApvmin 2.6.4

Synopsis
--------

.. code:: bash

    RNApvmin [options] <file.shape>

DESCRIPTION
-----------

RNApvmin 2.6.4

Calculate a perturbation vector that minimizes discripancies between predicted
and observed pairing probabilities

The program reads a RNA sequence from stdin and uses an iterative minimization
process to calculate a perturbation vector that minimizes the discripancies
between predicted pairing probabilites and observed pairing probabilities
(deduced from given shape reactivities). Experimental data is read from a given
SHAPE file and normalized to pairing probabilities. The experimental data has
to be provided in a multiline plain text file where each line has the format
``[position] [nucleotide] [absolute shape reactivity]`` (e.g. ``3 A 0.7``). The
objective function used for the minimization may be weighted by choosing
appropriate values for sigma and tau.

The minimization progress will be written to stderr. Once the minimization has
terminated, the obtained perturbation vector is written to stdout.

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

.. option:: -j, --numThreads=INT

    Set the number of threads used for calculations.

Algorithms:
^^^^^^^^^^^



    Select additional algorithms which should be included in the calculations.
    The Minimum free energy (MFE) and a structure representative are calculated
    in any case.

.. option:: --shapeConversion=STRING

    Specify the method used to convert SHAPE reactivities to pairing probabilities.


    *(default="O")*


    The following methods can be used to convert SHAPE reactivities into the
    probability for a certain nucleotide to be unpaired.


    ``M``: Use linear mapping according to Zarringhalam et al. 2012


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

.. option:: --tauSigmaRatio=DOUBLE

    Ratio of the weighting factors tau and sigma. *(default="1.0")*


    A high ratio will lead to a solution as close as possible to the experimental
    data, while a low ratio will lead to results close to the thermodynamic
    prediction without guiding pseudo energies.

.. option:: --objectiveFunction=INT

    The energies of the perturbation vector and the discripancies between predicted and observed pairing probabilities contribute to the objective function. This parameter defines, which function is used to process the contributions before summing them up. 0 square 1 absolute.


    *(default="0")*

.. option:: --sampleSize=INT

    The iterative minimization process requires to evaluate the gradient of the objective function.


    *(default="1000")*


    A sample size of 0 leads to an analytical evaluation which scales as O(N^4).
    Choosing a sample size >0 estimates the gradient by sampling the given number
    of sequences from the ensemble, which is much faster.

.. option:: -N, --nonRedundant

    Enable non-redundant sampling strategy.


    *(default=off)*

.. option:: --intermediatePath=STRING Write an output file for each iteration of the

    minimization process.


    Each file contains the used perturbation vector and the score of the
    objective function. The number of the iteration will be appended to the given
    path.

.. option:: --initialVector=DOUBLE

    Specify the vector of initial pertubations. *(default="0")*


    Defines the initial perturbation vector which will be used as starting vector
    for the minimization process. The value 0 results in a null vector. Every
    other value x will be used to populate the initial vector with random numbers
    from the interval [-x,x].

.. option:: --minimizer=ENUM

    Set the minimizing algorithm used for finding an appropriate perturbation vector.

.. option:: (possible values="conjugate_fr",

    "conjugate_pr", "vector_bfgs", "vector_bfgs2", "steepest_descent", "default" default="default")


    The default option uses a custom implementation of the gradient descent
    algorithms while all other options represent various algorithms implemented
    in the GNU Scientific Library. When the GNU Scientific Library can not be
    found, only the default minimizer is available.

.. option:: --initialStepSize=DOUBLE

    The initial stepsize for the minimizer methods.


    *(default="0.01")*

.. option:: --minStepSize=DOUBLE

    The minimal stepsize for the minizimer methods.


    *(default="1e-15")*

.. option:: --minImprovement=DOUBLE

    The minimal improvement in the default minizimer method that has to be surpassed to considered a new result a better one.


    *(default="1e-3")*

.. option:: --minimizerTolerance=DOUBLE

    The tolerance to be used in the GSL minimizer


    methods.


    *(default="1e-3")*

.. option:: -S, --pfScale=DOUBLE

    In the calculation of the pf use scale*mfe as an estimate for the ensemble free energy (used to avoid overflows).


    *(default="1.07")*


    The default is 1.07, useful values are 1.0 to 1.2. Occasionally needed for
    long sequences.

Structure Constraints:
^^^^^^^^^^^^^^^^^^^^^^



    Command line options to interact with the structure constraints feature of
    this program

.. option:: --maxBPspan=INT

    Set the maximum base pair span.


    *(default="-1")*

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

S. Washietl, I.L. Hofacker, P.F. Stadler, M. Kellis (2012)
"RNA folding with soft constraints: reconciliation of probing data and thermodynamics secondary structure prediction"
Nucl Acids Res: 40(10), pp 4261-4272


*The energy parameters are taken from:*

D.H. Mathews, M.D. Disney, D. Matthew, J.L. Childs, S.J. Schroeder, J. Susan, M. Zuker, D.H. Turner (2004),
"Incorporating chemical modification constraints into a dynamic programming algorithm for prediction of RNA secondary structure",
Proc. Natl. Acad. Sci. USA: 101, pp 7287-7292

D.H Turner, D.H. Mathews (2009),
"NNDB: The nearest neighbor parameter database for predicting stability of nucleic acid secondary structure",
Nucleic Acids Research: 38, pp 280-282

EXAMPLES
--------


RNApvmin acceptes a SHAPE file and a corresponding nucleotide sequence, which is read form stdin.

.. code::

    RNApvmin sequence.shape < sequence.fasta > sequence.pv
    


The normalized SHAPE reactivity data has to be stored in a text file, where each line contains the position
and the reactivity for a certain nucleotide ([position] [nucleotide] [SHAPE reactivity]).

.. code::

    1 A 1.286
    2 U 0.383
    3 C 0.033
    4 C 0.017
    ...
    ...
    98 U 0.234
    99 G 0.885
    


The nucleotide information in the SHAPE file is optional and will be used to cross check the given input sequence if present.
If SHAPE reactivities could not be determined for every nucleotide, missing values can simply be omited.

The progress of the minimization will be printed to stderr. Once a solution was found, the calculated perturbation vector
will be print to stdout and can then further be used to constrain RNAfold's MFE/partition function calculation by applying
the perturbation energies as soft constraints.

.. code::

    RNAfold --shape=sequence.pv --shapeMethod=W < sequence.fasta
    


AUTHOR
------


Dominik Luntzer, Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
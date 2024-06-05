##########
RNAparconv
##########

:program:`RNAparconv` - manual page for RNAparconv 2.6.4

Synopsis
--------

.. code:: bash

    RNAparconv [options] [<input file>] [<output file>]

DESCRIPTION
-----------

RNAparconv 2.6.4

Convert energy parameter files from ViennaRNA 1.8.4 to 2.0 format

Converts energy parameter files from "old" ViennaRNAPackage 1.8.4 format to
the new format used since ViennaRNAPackage 2.0.
The Program reads a valid energy parameter file or valid energy parameters from
stdin and prints the converted energy parameters to stdout or a specified
output file. Per default, the converted output file contains the whole set of
energy parameters used throughout ViennaRNAPackage 1.8.4. The user can specify
sets of energy parameters that should not be included in the output.

.. option:: -h, --help

    Print help and exit

.. option:: --full-help

    Print help, including hidden options, and exit

.. option:: -V, --version

    Print version and exit

I/O Options:
^^^^^^^^^^^^



    Command line options for input and output (pre-)processing

.. option:: -i, --input=filename

    Specify an input file name. If argument is missing the energy parameter input can be supplied via ``stdin``.

.. option:: -o, --output=filename

    Specify an output file name. If argument is missing the converted energy parameters are printed to ``stdout``.

.. option:: --vanilla

    Print just as much as needed to represent the given energy parameters data set. This option overrides all other output settings!


    *(default=off)*

.. option:: --dump

    Just dump Vienna 1.8.4 energy parameters in format used since 2.0. This option skips any energy parameter input!


    *(default=off)*

.. option:: --silent

    Print just energy parameters and appropriate comment lines but suppress all other output


    *(default=off)*

.. option:: --without-HairpinE

    Do not print converted hairpin energies and enthalpies


    *(default=off)*

.. option:: --without-StackE

    Do not print converted stacking energies and enthalpies


    *(default=off)*

.. option:: --without-IntE

    Do not print converted interior loop energies, enthalpies and asymetry factors


    *(default=off)*

.. option:: --without-BulgeE

    Do not print converted bulge loop energies and enthalpies


    *(default=off)*

.. option:: --without-MultiE

    Do not print converted multi loop energies and enthalpies


    *(default=off)*

.. option:: --without-MismatchE

    Do not print converted exterior loop mismatch energies and enthalpies


    *(default=off)*

.. option:: --without-MismatchH

    Do not print converted hairpin mismatch energies and enthalpies


    *(default=off)*

.. option:: --without-MismatchI

    Do not print converted interior loop mismatch energies and enthalpies


    *(default=off)*

.. option:: --without-MismatchM

    Do not print converted multi loop mismatch energies and enthalpies


    *(default=off)*

.. option:: --without-Dangle5

    Do not print converted 5' dangle energies and enthalpies


    *(default=off)*

.. option:: --without-Dangle3

    Do not print converted 3' dangle energies and enthalpies


    *(default=off)*

.. option:: --without-Misc

    Do not print converted Misc energies and enthalpies (TerminalAU, DuplexInit, lxc)


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


Ronny Lorenz

REPORTING BUGS
--------------


If in doubt our program is right, nature is at fault.
Comments should be sent to rna@tbi.univie.ac.at.
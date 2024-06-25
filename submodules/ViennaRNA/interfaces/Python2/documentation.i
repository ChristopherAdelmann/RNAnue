
// File: index.xml

// File: struct__struct__en.xml


%feature("docstring") _struct_en "

Data structure for energy_of_move()  

Attributes
----------
energy : int  

structure : short *  

C++ includes: ViennaRNA/move_set.h
";

// File: structconstrain.xml


%feature("docstring") constrain "

constraints for cofolding  

Attributes
----------
indx : int *  

ptype : char *  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structCOORDINATE.xml


%feature("docstring") COORDINATE "

this is a workarround for the SWIG Perl Wrapper RNA plot function that returns an array of type
COORDINATE  

Attributes
----------
X : float  

Y : float  

C++ includes: ViennaRNA/plotting/layouts.h
";

// File: structduplexT.xml


%feature("docstring") duplexT "

Data structure for RNAduplex.  

Attributes
----------
i : int  

j : int  

end : int  

structure : char *  

energy : double  

energy_backtrack : double  

opening_backtrack_x : double  

opening_backtrack_y : double  

offset : int  

dG1 : double  

dG2 : double  

ddG : double  

tb : int  

te : int  

qb : int  

qe : int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structdupVar.xml


%feature("docstring") dupVar "

Data structure used in RNApkplex.  

Attributes
----------
i : int  

j : int  

end : int  

pk_helix : char *  

structure : char *  

energy : double  

offset : int  

dG1 : double  

dG2 : double  

ddG : double  

tb : int  

te : int  

qb : int  

qe : int  

inactive : int  

processed : int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structenergy__corrections.xml


%feature("docstring") energy_corrections "

Attributes
----------
enc : short *  

strands : size_t  

ptypes : size_t  

stack_diff : int  

dangle5_diff : int  

dangle3_diff : int  

mismatch_diff : int  

terminal_diff : int  
";

%feature("docstring") energy_corrections::vrna_array "
energy_corrections::vrna_array";

// File: structinteract.xml


%feature("docstring") interact "

interaction data structure for RNAup  

Attributes
----------
Pi : double *  
    probabilities of interaction  

Gi : double *  
    free energies of interaction  

Gikjl : double  
    full free energy for interaction between [k,i] k<i in longer seq and [j,l] j<l in shorter seq  

Gikjl_wo : double  
    Gikjl without contributions for prob_unpaired.  

i : int  
    k<i in longer seq  

k : int  
    k<i in longer seq  

j : int  
    j<l in shorter seq  

l : int  
    j<l in shorter seq  

length : int  
    length of longer sequence  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structLIST.xml


%feature("docstring") LIST "

Attributes
----------
count : int  

head : LST_BUCKET *  

z : LST_BUCKET *  

hz : LST_BUCKET  
";

// File: structLST__BUCKET.xml


%feature("docstring") LST_BUCKET "

Attributes
----------
next : struct LST_BUCKET *  
";

// File: structnode.xml


%feature("docstring") node "

Data structure for RNAsnoop (fold energy list)  

Attributes
----------
k : int  

energy : int  

next : struct node *  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structPostorder__list.xml


%feature("docstring") Postorder_list "

Postorder data structure.  

Attributes
----------
type : int  

weight : int  

father : int  

sons : int  

leftmostleaf : int  

C++ includes: ViennaRNA/dist_vars.h
";

// File: structpu__contrib.xml


%feature("docstring") pu_contrib "

contributions to p_u  

Attributes
----------
H : double **  
    hairpin loops  

I : double **  
    interior loops  

M : double **  
    multi loops  

E : double **  
    exterior loop  

length : int  
    length of the input sequence  

w : int  
    longest unpaired region  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structpu__out.xml


%feature("docstring") pu_out "

Collection of all free_energy of beeing unpaired values for output.  

Attributes
----------
len : int  
    sequence length  

u_vals : int  
    number of different -u values  

contribs : int  
    [-c \"SHIME\"]  

header : char **  
    header line  

u_values : double **  
    (the -u values * [-c \"SHIME\"]) * seq len  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structsnoopT.xml


%feature("docstring") snoopT "

Data structure for RNAsnoop.  

Attributes
----------
i : int  

j : int  

u : int  

structure : char *  

energy : float  

Duplex_El : float  

Duplex_Er : float  

Loop_E : float  

Loop_D : float  

pscd : float  

psct : float  

pscg : float  

Duplex_Ol : float  

Duplex_Or : float  

Duplex_Ot : float  

fullStemEnergy : float  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structswString.xml


%feature("docstring") swString "

Some other data structure.  

Attributes
----------
type : int  

sign : int  

weight : float  

C++ includes: ViennaRNA/dist_vars.h
";

// File: structTree.xml


%feature("docstring") Tree "

Tree data structure.  

Attributes
----------
postorder_list : Postorder_list *  

keyroots : int *  

C++ includes: ViennaRNA/dist_vars.h
";

// File: structTwoDfold__vars.xml


%feature("docstring") TwoDfold_vars "

Variables compound for 2Dfold MFE folding.  

.. deprecated:: 2.6.4
    This data structure will be removed from the library soon! Use RNA.fold_compound() and the
    corresponding functions RNA.fold_compound_TwoD(), RNA.mfe_TwoD(), and
    RNA.fold_compound_free() instead!  

Attributes
----------
P : vrna_param_t *  
    Precomputed energy parameters and model details.  

do_backtrack : int  
    Flag whether to do backtracing of the structure(s) or not.  

ptype : char *  
    Precomputed array of pair types.  

sequence : char *  
    The input sequence  

S : short *  

S1 : short *  
    The input sequences in numeric form.  

maxD1 : unsigned int  
    Maximum allowed base pair distance to first reference.  

maxD2 : unsigned int  
    Maximum allowed base pair distance to second reference.  

mm1 : unsigned int *  
    Maximum matching matrix, reference struct 1 disallowed.  

mm2 : unsigned int *  
    Maximum matching matrix, reference struct 2 disallowed.  

my_iindx : int *  
    Index for moving in quadratic distancy dimensions.  

temperature : double  

referenceBPs1 : unsigned int *  
    Matrix containing number of basepairs of reference structure1 in interval [i,j].  

referenceBPs2 : unsigned int *  
    Matrix containing number of basepairs of reference structure2 in interval [i,j].  

bpdist : unsigned int *  
    Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j].  

reference_pt1 : short *  

reference_pt2 : short *  

circ : int  

dangles : int  

seq_length : unsigned int  

E_F5 : int ***  

E_F3 : int ***  

E_C : int ***  

E_M : int ***  

E_M1 : int ***  

E_M2 : int ***  

E_Fc : int **  

E_FcH : int **  

E_FcI : int **  

E_FcM : int **  

l_min_values : int **  

l_max_values : int **  

k_min_values : int *  

k_max_values : int *  

l_min_values_m : int **  

l_max_values_m : int **  

k_min_values_m : int *  

k_max_values_m : int *  

l_min_values_m1 : int **  

l_max_values_m1 : int **  

k_min_values_m1 : int *  

k_max_values_m1 : int *  

l_min_values_f : int **  

l_max_values_f : int **  

k_min_values_f : int *  

k_max_values_f : int *  

l_min_values_f3 : int **  

l_max_values_f3 : int **  

k_min_values_f3 : int *  

k_max_values_f3 : int *  

l_min_values_m2 : int **  

l_max_values_m2 : int **  

k_min_values_m2 : int *  

k_max_values_m2 : int *  

l_min_values_fc : int *  

l_max_values_fc : int *  

k_min_values_fc : int  

k_max_values_fc : int  

l_min_values_fcH : int *  

l_max_values_fcH : int *  

k_min_values_fcH : int  

k_max_values_fcH : int  

l_min_values_fcI : int *  

l_max_values_fcI : int *  

k_min_values_fcI : int  

k_max_values_fcI : int  

l_min_values_fcM : int *  

l_max_values_fcM : int *  

k_min_values_fcM : int  

k_max_values_fcM : int  

E_F5_rem : int *  

E_F3_rem : int *  

E_C_rem : int *  

E_M_rem : int *  

E_M1_rem : int *  

E_M2_rem : int *  

E_Fc_rem : int  

E_FcH_rem : int  

E_FcI_rem : int  

E_FcM_rem : int  

compatibility : vrna_fold_compound_t *  

C++ includes: ViennaRNA/2Dfold.h
";

// File: structTwoDpfold__vars.xml


%feature("docstring") TwoDpfold_vars "

Variables compound for 2Dfold partition function folding.  

.. deprecated:: 2.6.4
    This data structure will be removed from the library soon! Use RNA.fold_compound() and the
    corresponding functions RNA.fold_compound_TwoD(), RNA.pf_TwoD(), and RNA.fold_compound_free()
    instead!  

Attributes
----------
alloc : unsigned int  

ptype : char *  
    Precomputed array of pair types.  

sequence : char *  
    The input sequence  

S : short *  

S1 : short *  
    The input sequences in numeric form.  

maxD1 : unsigned int  
    Maximum allowed base pair distance to first reference.  

maxD2 : unsigned int  
    Maximum allowed base pair distance to second reference.  

temperature : double  

init_temp : double  

scale : FLT_OR_DBL *  

pf_scale : FLT_OR_DBL  

pf_params : vrna_exp_param_t *  

my_iindx : int *  
    Index for moving in quadratic distancy dimensions.  

jindx : int *  
    Index for moving in the triangular matrix qm1.  

reference_pt1 : short *  

reference_pt2 : short *  

referenceBPs1 : unsigned int *  
    Matrix containing number of basepairs of reference structure1 in interval [i,j].  

referenceBPs2 : unsigned int *  
    Matrix containing number of basepairs of reference structure2 in interval [i,j].  

bpdist : unsigned int *  
    Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j].  

mm1 : unsigned int *  
    Maximum matching matrix, reference struct 1 disallowed.  

mm2 : unsigned int *  
    Maximum matching matrix, reference struct 2 disallowed.  

circ : int  

dangles : int  

seq_length : unsigned int  

Q : FLT_OR_DBL ***  

Q_B : FLT_OR_DBL ***  

Q_M : FLT_OR_DBL ***  

Q_M1 : FLT_OR_DBL ***  

Q_M2 : FLT_OR_DBL ***  

Q_c : FLT_OR_DBL **  

Q_cH : FLT_OR_DBL **  

Q_cI : FLT_OR_DBL **  

Q_cM : FLT_OR_DBL **  

l_min_values : int **  

l_max_values : int **  

k_min_values : int *  

k_max_values : int *  

l_min_values_b : int **  

l_max_values_b : int **  

k_min_values_b : int *  

k_max_values_b : int *  

l_min_values_m : int **  

l_max_values_m : int **  

k_min_values_m : int *  

k_max_values_m : int *  

l_min_values_m1 : int **  

l_max_values_m1 : int **  

k_min_values_m1 : int *  

k_max_values_m1 : int *  

l_min_values_m2 : int **  

l_max_values_m2 : int **  

k_min_values_m2 : int *  

k_max_values_m2 : int *  

l_min_values_qc : int *  

l_max_values_qc : int *  

k_min_values_qc : int  

k_max_values_qc : int  

l_min_values_qcH : int *  

l_max_values_qcH : int *  

k_min_values_qcH : int  

k_max_values_qcH : int  

l_min_values_qcI : int *  

l_max_values_qcI : int *  

k_min_values_qcI : int  

k_max_values_qcI : int  

l_min_values_qcM : int *  

l_max_values_qcM : int *  

k_min_values_qcM : int  

k_max_values_qcM : int  

Q_rem : FLT_OR_DBL *  

Q_B_rem : FLT_OR_DBL *  

Q_M_rem : FLT_OR_DBL *  

Q_M1_rem : FLT_OR_DBL *  

Q_M2_rem : FLT_OR_DBL *  

Q_c_rem : FLT_OR_DBL  

Q_cH_rem : FLT_OR_DBL  

Q_cI_rem : FLT_OR_DBL  

Q_cM_rem : FLT_OR_DBL  

compatibility : vrna_fold_compound_t *  

C++ includes: ViennaRNA/2Dpfold.h
";

// File: structvrna__alignment__s.xml


%feature("docstring") vrna_alignment_s "

Attributes
----------
n_seq : unsigned int  

sequences : vrna_seq_t *  

gapfree_seq : char **  

gapfree_size : unsigned int *  

genome_size : unsigned long long *  

start : unsigned long long *  

orientation : unsigned char *  

a2s : unsigned int **  
";

// File: structvrna__array__header__s.xml


%feature("docstring") vrna_array_header_s "

The header of an array.  

Attributes
----------
num : size_t  
    The number of elements in an array.  

size : size_t  
    The actual capacity of an array.  

C++ includes: ViennaRNA/datastructures/array.h
";

// File: structvrna__basepair__s.xml


%feature("docstring") vrna_basepair_s "

Base pair data structure used in subopt.c.  

Attributes
----------
i : int  

j : int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structvrna__bp__stack__s.xml


%feature("docstring") vrna_bp_stack_s "

Base pair stack element.  

Attributes
----------
i : unsigned int  

j : unsigned int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structvrna__color__s.xml


%feature("docstring") vrna_color_s "

Attributes
----------
hue : float  

sat : float  

bri : float  
";

// File: structvrna__cpair__s.xml


%feature("docstring") vrna_cpair_s "

this datastructure is used as input parameter in functions of PS_dot.c  

Attributes
----------
i : int  

j : int  

mfe : int  

p : float  

hue : float  

sat : float  

type : int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structvrna__data__linear__s.xml


%feature("docstring") vrna_data_linear_s "

Attributes
----------
position : unsigned int  

value : float  

color : vrna_color_t  
";

// File: structvrna__dimer__conc__s.xml


%feature("docstring") vrna_dimer_conc_s "

Data structure for concentration dependency computations.  

Attributes
----------
Ac_start : double  
    start concentration A  

Bc_start : double  
    start concentration B  

ABc : double  
    End concentration AB.  

AAc : double  

BBc : double  

Ac : double  

Bc : double  

C++ includes: ViennaRNA/concentrations.h
";

// File: structvrna__dimer__pf__s.xml


%feature("docstring") vrna_dimer_pf_s "

Data structure returned by RNA.fold_compound.pf_dimer()  

Attributes
----------
F0AB : double  
    Null model without DuplexInit.  

FAB : double  
    all states with DuplexInit correction  

FcAB : double  
    true hybrid states only  

FA : double  
    monomer A  

FB : double  
    monomer B  

C++ includes: ViennaRNA/part_func.h
";

// File: structvrna__dotplot__auxdata__t.xml


%feature("docstring") vrna_dotplot_auxdata_t "

Attributes
----------
comment : char *  

title : char *  

top : vrna_data_lin_t **  

top_title : char **  

bottom : vrna_data_lin_t **  

bottom_title : char **  

left : vrna_data_lin_t **  

left_title : char **  

right : vrna_data_lin_t **  

right_title : char **  
";

// File: structvrna__elem__prob__s.xml


%feature("docstring") vrna_ep_t "

Data structure representing a single entry of an element probability list (e.g. list of pair
probabilities)  

See Also
--------
RNA.plist(), RNA.fold_compound.plist_from_probs(), RNA.db_from_plist(),  
RNA.PLIST_TYPE_BASEPAIR, RNA.PLIST_TYPE_GQUAD, RNA.PLIST_TYPE_H_MOTIF, RNA.PLIST_TYPE_I_MOTIF,
RNA.PLIST_TYPE_UD_MOTIF, RNA.PLIST_TYPE_STACK  

Attributes
----------
i : int  
    Start position (usually 5' nucleotide that starts the element, e.g. base pair)  

j : int  
    End position (usually 3' nucleotide that ends the element, e.g. base pair)  

p : float  
    Probability of the element.  

type : int  
    Type of the element.  

C++ includes: ViennaRNA/utils/structures.h
";

// File: structvrna__exp__param__s.xml


%feature("docstring") vrna_exp_param_t "

The data structure that contains temperature scaled Boltzmann weights of the energy parameters.  

Attributes
----------
id : int  
    An identifier for the data structure.  

    .. deprecated:: 2.6.4
        This attribute will be removed in version 3  

expstack : double  

exphairpin : double  

expbulge : double  

expinternal : double  

expmismatchExt : double  

expmismatchI : double  

expmismatch23I : double  

expmismatch1nI : double  

expmismatchH : double  

expmismatchM : double  

expdangle5 : double  

expdangle3 : double  

expint11 : double  

expint21 : double  

expint22 : double  

expninio : double  

lxc : double  

expMLbase : double  

expMLintern : double  

expMLclosing : double  

expTermAU : double  

expDuplexInit : double  

exptetra : double  

exptri : double  

exphex : double  

Tetraloops : char  

expTriloop : double  

Triloops : char  

Hexaloops : char  

expTripleC : double  

expMultipleCA : double  

expMultipleCB : double  

expgquad : double  

expgquadLayerMismatch : double  

gquadLayerMismatchMax : int  

kT : double  

pf_scale : double  
    Scaling factor to avoid over-/underflows.  

temperature : double  
    Temperature used for loop contribution scaling.  

alpha : double  
    Scaling factor for the thermodynamic temperature.  

    This allows for temperature scaling in Boltzmann factors independently from the energy
    contributions. The resulting Boltzmann factors are then computed by :math:`e^{-E/(\\alpha \\cdot
    K \\cdot T)}`  

model_details : vrna_md_t  
    Model details to be used in the recursions.  

param_file : char  
    The filename the parameters were derived from, or empty string if they represent the default.  

expSaltStack : double  

expSaltLoop : double  

SaltLoopDbl : double  

SaltMLbase : int  

SaltMLintern : int  

SaltMLclosing : int  

SaltDPXInit : int  

C++ includes: ViennaRNA/params/basic.h
";

// File: structvrna__fc__s.xml


%feature("docstring") vrna_fold_compound_t "

The most basic data structure required by many functions throughout the RNAlib.  

Note
----
Please read the documentation of this data structure carefully! Some attributes are only available
for specific types this data structure can adopt.  

Warnings
--------
Reading/Writing from/to attributes that are not within the scope of the current type usually result
in undefined behavior!  

See Also
--------
RNA.fold_compound().type, RNA.fold_compound(), RNA.fold_compound_comparative(),
RNA.fold_compound_free(), RNA.FC_TYPE_SINGLE, RNA.FC_TYPE_COMPARATIVE  

**SWIG Wrapper Notes**
  
    This data structure is wrapped as class `fold_compound` with several related functions attached
    as methods.  

    A new `fold_compound` can be obtained by calling one of its constructors:  

    *   `fold_compound(seq)` - Initialize with a single sequence, or two concatenated sequences
        separated by an ampersand character `&` (for cofolding)  
    *   `fold_compound(aln)` - Initialize with a sequence alignment *aln* stored as a list of
        sequences (with gap characters).  

    The resulting object has a list of attached methods which in most cases directly correspond to
    functions that mainly operate on the corresponding `C` data structure:  

    *   `type()` - Get the type of the *fold_compound* (See RNA.fc_type)  
    *   `length()` - Get the length of the sequence(s) or alignment stored within the
        `fold_compound`.  

    See, e.g.  :py:class:`RNA.fold_compound` in the :doc:`/api_python`.  

Attributes
----------
type : const vrna_fc_type_e  
    The type of the RNA.fold_compound().  

    Currently possible values are RNA.FC_TYPE_SINGLE, and RNA.FC_TYPE_COMPARATIVE  

    Warnings
    --------
    Do not edit this attribute, it will be automagically set by the corresponding get() methods for
    the RNA.fold_compound(). The value specified in this attribute dictates the set of other
    attributes to use within this data structure.  

length : unsigned int  
    The length of the sequence (or sequence alignment)  

cutpoint : int  
    The position of the (cofold) cutpoint within the provided sequence. If there is no cutpoint,
    this field will be set to -1.  

strand_number : unsigned int *  
    The strand number a particular nucleotide is associated with.  

strand_order : unsigned int *  
    The strand order, i.e. permutation of current concatenated sequence.  

strand_order_uniq : unsigned int *  
    The strand order array where identical sequences have the same ID.  

strand_start : unsigned int *  
    The start position of a particular strand within the current concatenated sequence.  

strand_end : unsigned int *  
    The end (last) position of a particular strand within the current concatenated sequence.  

strands : unsigned int  
    Number of interacting strands.  

nucleotides : vrna_seq_t *  
    Set of nucleotide sequences.  

alignment : vrna_msa_t *  
    Set of alignments.  

hc : vrna_hc_t *  
    The hard constraints data structure used for structure prediction.  

matrices : vrna_mx_mfe_t *  
    The MFE DP matrices.  

exp_matrices : vrna_mx_pf_t *  
    The PF DP matrices  

params : vrna_param_t *  
    The precomputed free energy contributions for each type of loop.  

exp_params : vrna_exp_param_t *  
    The precomputed free energy contributions as Boltzmann factors  

iindx : int *  
    DP matrix accessor  

jindx : int *  
    DP matrix accessor  

stat_cb : vrna_recursion_status_f  
    Recursion status callback (usually called just before, and after recursive computations in the
    library.  

    See Also
    --------
    RNA.recursion_status(), RNA.fold_compound.add_callback()  

auxdata : void *  
    A pointer to auxiliary, user-defined data.  

    See Also
    --------
    RNA.fold_compound.add_auxdata(), RNA.fold_compound().free_auxdata  

free_auxdata : vrna_auxdata_free_f  
    A callback to free auxiliary user data whenever the fold_compound itself is free'd.  

    See Also
    --------
    RNA.fold_compound().auxdata, RNA.auxdata_free()  

domains_struc : vrna_sd_t *  
    Additional structured domains.  

domains_up : vrna_ud_t *  
    Additional unstructured domains.  

aux_grammar : vrna_gr_aux_t *  
    Additional decomposition grammar rules.  

sequence : char *  
    The input sequence string.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_SINGLE  

sequence_encoding : short *  
    Numerical encoding of the input sequence.  

    See Also
    --------
    RNA.sequence_encode()  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_SINGLE  

encoding5 : short *  

encoding3 : short *  

sequence_encoding2 : short *  

ptype : char *  
    Pair type array.  

    Contains the numerical encoding of the pair type for each pair (i,j) used in MFE, Partition
    function and Evaluation computations.  

    Note
    ----
    This array is always indexed via jindx, in contrast to previously different indexing between mfe
    and pf variants!  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_SINGLE  

    See Also
    --------
    RNA.idx_col_wise(), RNA.ptypes()  

ptype_pf_compat : char *  
    ptype array indexed via iindx  

    .. deprecated:: 2.6.4
        This attribute will vanish in the future! It's meant for backward compatibility only!  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_SINGLE  

sc : vrna_sc_t *  
    The soft constraints for usage in structure prediction and evaluation.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_SINGLE  

sequences : char **  
    The aligned sequences.  

    Note
    ----
    The end of the alignment is indicated by a NULL pointer in the second dimension  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

n_seq : unsigned int  
    The number of sequences in the alignment.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

cons_seq : char *  
    The consensus sequence of the aligned sequences.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

S_cons : short *  
    Numerical encoding of the consensus sequence.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

S : short **  
    Numerical encoding of the sequences in the alignment.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

S5 : short **  
    S5[s][i] holds next base 5' of i in sequence s.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

S3 : short **  
    Sl[s][i] holds next base 3' of i in sequence s.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

Ss : char **  

a2s : unsigned int **  

pscore : int *  
    Precomputed array of pair types expressed as pairing scores.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

pscore_local : int **  
    Precomputed array of pair types expressed as pairing scores.  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

pscore_pf_compat : short *  
    Precomputed array of pair types expressed as pairing scores indexed via iindx.  

    .. deprecated:: 2.6.4
        This attribute will vanish in the future!  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

scs : vrna_sc_t **  
    A set of soft constraints (for each sequence in the alignment)  

    Warnings
    --------
    Only available if  

        type==RNA.FC_TYPE_COMPARATIVE  

oldAliEn : int  

maxD1 : unsigned int  
    Maximum allowed base pair distance to first reference.  

maxD2 : unsigned int  
    Maximum allowed base pair distance to second reference.  

reference_pt1 : short *  
    A pairtable of the first reference structure.  

reference_pt2 : short *  
    A pairtable of the second reference structure.  

referenceBPs1 : unsigned int *  
    Matrix containing number of basepairs of reference structure1 in interval [i,j].  

referenceBPs2 : unsigned int *  
    Matrix containing number of basepairs of reference structure2 in interval [i,j].  

bpdist : unsigned int *  
    Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j].  

mm1 : unsigned int *  
    Maximum matching matrix, reference struct 1 disallowed.  

mm2 : unsigned int *  
    Maximum matching matrix, reference struct 2 disallowed.  

window_size : int  
    window size for local folding sliding window approach  

ptype_local : char **  
    Pair type array (for local folding)  

zscore_data : vrna_zsc_dat_t  
    Data structure with settings for z-score computations.  

 : union vrna_fc_s  

C++ includes: ViennaRNA/fold_compound.h
";

/*
 Common data fields 
*/

/*
 User-defined data fields 
*/

/*
 Secondary Structure Decomposition (grammar) related data fields 
*/

/*
 Data fields available for single/hybrid structure prediction 
*/

/*
 Data fields for consensus structure prediction 
*/

/*
 Additional data fields for Distance Class Partitioning 
*/

/*
These data fields are typically populated with meaningful data only if used in the context of
Distance Class Partitioning  

*/

/*
 Additional data fields for local folding 
*/

/*
These data fields are typically populated with meaningful data only if used in the context of local
folding  

*/

// File: structvrna__gr__aux__s.xml


%feature("docstring") vrna_gr_aux_s "

Attributes
----------
cb_proc : vrna_grammar_cond_f  
    A callback for pre- and post-processing of auxiliary grammar rules.  

cb_aux_f : vrna_grammar_rule_f  

cb_aux_c : vrna_grammar_rule_f  

cb_aux_m : vrna_grammar_rule_f  

cb_aux_m1 : vrna_grammar_rule_f  

cb_aux : vrna_grammar_rule_f_aux  

cb_aux_exp_f : vrna_grammar_rule_f_exp  

cb_aux_exp_c : vrna_grammar_rule_f_exp  

cb_aux_exp_m : vrna_grammar_rule_f_exp  

cb_aux_exp_m1 : vrna_grammar_rule_f_exp  

cb_aux_exp : vrna_grammar_rule_f_aux_exp  

data : void *  

free_data : vrna_grammar_data_free_f  
";

// File: structvrna__hc__s.xml


%feature("docstring") vrna_hc_t "

The hard constraints data structure.  

The content of this data structure determines the decomposition pattern used in the folding
recursions. Attribute 'matrix' is used as source for the branching pattern of the decompositions
during all folding recursions. Any entry in matrix[i,j] consists of the 6 LSB that allows one to
distinguish the following types of base pairs:  

*   in the exterior loop (RNA.CONSTRAINT_CONTEXT_EXT_LOOP)  
*   enclosing a hairpin (RNA.CONSTRAINT_CONTEXT_HP_LOOP)  
*   enclosing an interior loop (RNA.CONSTRAINT_CONTEXT_INT_LOOP)  
*   enclosed by an exterior loop (RNA.CONSTRAINT_CONTEXT_INT_LOOP_ENC)  
*   enclosing a multi branch loop (RNA.CONSTRAINT_CONTEXT_MB_LOOP)  
*   enclosed by a multi branch loop (RNA.CONSTRAINT_CONTEXT_MB_LOOP_ENC)  

The four linear arrays 'up_xxx' provide the number of available unpaired nucleotides (including
position i) 3' of each position in the sequence.  

See Also
--------
RNA.fold_compound.hc_init(), RNA.hc_free(), RNA.CONSTRAINT_CONTEXT_EXT_LOOP,
RNA.CONSTRAINT_CONTEXT_HP_LOOP,
RNA.CONSTRAINT_CONTEXT_INT_LOOP, RNA.CONSTRAINT_CONTEXT_MB_LOOP,
RNA.CONSTRAINT_CONTEXT_MB_LOOP_ENC  

Attributes
----------
type : vrna_hc_type_e  

n : unsigned int  

state : unsigned char  

mx : unsigned char *  

matrix_local : unsigned char **  

 : union vrna_hc_s  

up_ext : int *  
    A linear array that holds the number of allowed unpaired nucleotides in an exterior loop.  

up_hp : int *  
    A linear array that holds the number of allowed unpaired nucleotides in a hairpin loop.  

up_int : int *  
    A linear array that holds the number of allowed unpaired nucleotides in an interior loop.  

up_ml : int *  
    A linear array that holds the number of allowed unpaired nucleotides in a multi branched loop.  

f : vrna_hc_eval_f  
    A function pointer that returns whether or not a certain decomposition may be evaluated.  

data : void *  
    A pointer to some structure where the user may store necessary data to evaluate its generic hard
    constraint function.  

free_data : vrna_auxdata_free_f  
    A pointer to a function to free memory occupied by auxiliary data.  

    The function this pointer is pointing to will be called upon destruction of the RNA.hc(), and
    provided with the RNA.hc().data pointer that may hold auxiliary data. Hence, to avoid leaking
    memory, the user may use this pointer to free memory occupied by auxiliary data.  

depot : vrna_hc_depot_t *  

C++ includes: ViennaRNA/constraints/hard.h
";

// File: structvrna__hc__up__s.xml


%feature("docstring") vrna_hc_up_s "

A single hard constraint for a single nucleotide.  

Attributes
----------
position : int  
    The sequence position (1-based)  

strand : int  

options : unsigned char  
    The hard constraint option  

C++ includes: ViennaRNA/constraints/hard.h
";

// File: structvrna__heat__capacity__s.xml


%feature("docstring") vrna_heat_capacity_s "

A single result from heat capacity computations.  

See Also
--------
RNA.fold_compound.heat_capacity()  

Attributes
----------
temperature : float  
    The temperature in C.  

heat_capacity : float  
    The specific heat at this temperature in Kcal/(Mol * K)  

C++ includes: ViennaRNA/heat_capacity.h
";

// File: structvrna__ht__entry__db__t.xml


%feature("docstring") vrna_ht_entry_db_t "

Default hash table entry.  

See Also
--------
RNA.ht_init(), RNA.ht_db_comp(), RNA.ht_db_hash_func(), RNA.ht_db_free_entry()  

Attributes
----------
structure : char *  
    A secondary structure in dot-bracket notation  

energy : float  
    The free energy of `structure`  

C++ includes: ViennaRNA/datastructures/hash_tables.h
";

// File: structvrna__hx__s.xml


%feature("docstring") vrna_hx_s "

Data structure representing an entry of a helix list.  

Attributes
----------
start : unsigned int  

end : unsigned int  

length : unsigned int  

up5 : unsigned int  

up3 : unsigned int  

C++ includes: ViennaRNA/utils/structures.h
";

// File: structvrna__md__s.xml


%feature("docstring") vrna_md_t "

The data structure that contains the complete model details used throughout the calculations.  

For convenience reasons, we provide the type name RNA.md() to address this data structure without
the use of the struct keyword  

See Also
--------
RNA.md.reset(), set_model_details(), RNA.md_update(), RNA.md()  

**SWIG Wrapper Notes**
    This data structure is wrapped as an object `md` with multiple related functions attached as
    methods.  

    A new set of default parameters can be obtained by calling the constructure of `md:`  

    *   `md()`-- Initialize with default settings  

    The resulting object has a list of attached methods which directly correspond to functions that
    mainly operate on the corresponding *C* data structure:  

    *   `reset()` - RNA.md.reset()  
    *   `set_from_globals()` - set_model_details()  
    *   `option_string()` - RNA.md.option_string()  

    Note, that default parameters can be modified by directly setting any of the following global
    variables. Internally, getting/setting default parameters using their global variable
    representative translates into calls of the following functions, therefore these wrappers for
    these functions do not exist in the scripting language interface(s):  
  

global variable  

`C getter`  

`C setter`  

temperature  

RNA.md_defaults_temperature_get()  

RNA.md_defaults_temperature()  

dangles  

RNA.md_defaults_dangles_get()  

RNA.md_defaults_dangles()  

betaScale  

RNA.md_defaults_betaScale_get()  

RNA.md_defaults_betaScale()  

tetra_loop  

this is an alias of *special_hp*  


special_hp  

RNA.md_defaults_special_hp_get()  

RNA.md_defaults_special_hp()  

noLonelyPairs  

this is an alias of *noLP*  


noLP  

RNA.md_defaults_noLP_get()  

RNA.md_defaults_noLP()  

noGU  

RNA.md_defaults_noGU_get()  

RNA.md_defaults_noGU()  

no_closingGU  

this is an alias of *noGUclosure*  


noGUclosure  

RNA.md_defaults_noGUclosure_get()  

RNA.md_defaults_noGUclosure()  

logML  

RNA.md_defaults_logML_get()  

RNA.md_defaults_logML()  

circ  

RNA.md_defaults_circ_get()  

RNA.md_defaults_circ()  

gquad  

RNA.md_defaults_gquad_get()  

RNA.md_defaults_gquad()  

uniq_ML  

RNA.md_defaults_uniq_ML_get()  

RNA.md_defaults_uniq_ML()  

energy_set  

RNA.md_defaults_energy_set_get()  

RNA.md_defaults_energy_set()  

backtrack  

RNA.md_defaults_backtrack_get()  

RNA.md_defaults_backtrack()  

backtrack_type  

RNA.md_defaults_backtrack_type_get()  

RNA.md_defaults_backtrack_type()  

do_backtrack  

this is an alias of *compute_bpp*  


compute_bpp  

RNA.md_defaults_compute_bpp_get()  

RNA.md_defaults_compute_bpp()  

max_bp_span  

RNA.md_defaults_max_bp_span_get()  

RNA.md_defaults_max_bp_span()  

min_loop_size  

RNA.md_defaults_min_loop_size_get()  

RNA.md_defaults_min_loop_size()  

window_size  

RNA.md_defaults_window_size_get()  

RNA.md_defaults_window_size()  

oldAliEn  

RNA.md_defaults_oldAliEn_get()  

RNA.md_defaults_oldAliEn()  

ribo  

RNA.md_defaults_ribo_get()  

RNA.md_defaults_ribo()  

cv_fact  

RNA.md_defaults_cv_fact_get()  

RNA.md_defaults_cv_fact()  

nc_fact  

RNA.md_defaults_nc_fact_get()  

RNA.md_defaults_nc_fact()  

sfact  

RNA.md_defaults_sfact_get()  

RNA.md_defaults_sfact()  

Attributes
----------
temperature : double  
    The temperature used to scale the thermodynamic parameters.  

betaScale : double  
    A scaling factor for the thermodynamic temperature of the Boltzmann factors.  

pf_smooth : int  
    A flat specifying whether energies in Boltzmann factors need to be smoothed.  

dangles : int  
    Specifies the dangle model used in any energy evaluation (0,1,2 or 3)  

    If set to 0 no stabilizing energies are assigned to bases adjacent to helices in free ends and
    multiloops (so called dangling ends). Normally (dangles = 1) dangling end energies are assigned
    only to unpaired bases and a base cannot participate simultaneously in two dangling ends. In the
    partition function algorithm RNA.fold_compound.pf() these checks are neglected. To provide comparability
    between free energy minimization and partition function algorithms, the default setting is 2.
    This treatment of dangling ends gives more favorable energies to helices directly adjacent to
    one another, which can be beneficial since such helices often do engage in stabilizing
    interactions through co-axial stacking.  
    If set to 3 co-axial stacking is explicitly included for adjacent helices in multiloops. The
    option affects only mfe folding and energy evaluation (RNA.mfe() and RNA.eval_structure()), as
    well as suboptimal folding (RNA.subopt()) via re-evaluation of energies. Co-axial stacking with
    one intervening mismatch is not considered so far. Note, that some function do not implement all
    dangle model but only a subset of (0,1,2,3). In particular, partition function algorithms can
    only handle 0 and 2. Read the documentation of the particular recurrences or energy evaluation
    function for information about the provided dangle model.  

special_hp : int  
    Include special hairpin contributions for tri, tetra and hexaloops.  

noLP : int  
    Only consider canonical structures, i.e. no 'lonely' base pairs.  

noGU : int  
    Do not allow GU pairs.  

noGUclosure : int  
    Do not allow loops to be closed by GU pair.  

logML : int  
    Use logarithmic scaling for multiloops.  

circ : int  
    Assume RNA to be circular instead of linear.  

gquad : int  
    Include G-quadruplexes in structure prediction.  

uniq_ML : int  
    Flag to ensure unique multi-branch loop decomposition during folding.  

energy_set : int  
    Specifies the energy set that defines set of compatible base pairs.  

backtrack : int  
    Specifies whether or not secondary structures should be backtraced.  

backtrack_type : char  
    Specifies in which matrix to backtrack.  

compute_bpp : int  
    Specifies whether or not backward recursions for base pair probability (bpp) computation will be
    performed.  

nonstandards : char  
    contains allowed non standard bases  

max_bp_span : int  
    maximum allowed base pair span  

min_loop_size : int  
    Minimum size of hairpin loops.  

    The default value for this field is TURN, however, it may be 0 in cofolding context.  

window_size : int  
    Size of the sliding window for locally optimal structure prediction.  

oldAliEn : int  
    Use old alifold energy model.  

ribo : int  
    Use ribosum scoring table in alifold energy model.  

cv_fact : double  
    Co-variance scaling factor for consensus structure prediction.  

nc_fact : double  
    Scaling factor to weight co-variance contributions of non-canonical pairs.  

sfact : double  
    Scaling factor for partition function scaling.  

rtype : int  
    Reverse base pair type array.  

alias : short  
    alias of an integer nucleotide representation  

pair : int  
    Integer representation of a base pair.  

pair_dist : float  
    Base pair dissimilarity, a.k.a. distance matrix.  

salt : double  
    Salt (monovalent) concentration (M) in buffer.  

saltMLLower : int  
    Lower bound of multiloop size to use in loop salt correction linear fitting.  

saltMLUpper : int  
    Upper bound of multiloop size to use in loop salt correction linear fitting.  

saltDPXInit : int  
    User-provided salt correction for duplex initialization (in dcal/mol). If set to 99999 the
    default salt correction is used. If set to 0 there is no salt correction for duplex
    initialization.  

saltDPXInitFact : float  
      

helical_rise : float  
      

backbone_length : float  
      

C++ includes: ViennaRNA/model.h
";

// File: structvrna__move__s.xml


%feature("docstring") vrna_move_t "

An atomic representation of the transition / move from one structure to its neighbor.  

An atomic transition / move may be one of the following:  

*   a **base pair insertion**,  
*   a **base pair removal**, or  
*   a **base pair shift** where an existing base pair changes one of its pairing partner.  

These moves are encoded by two integer values that represent the affected 5' and 3' nucleotide
positions. Furthermore, we use the following convention on the signedness of these encodings:  

*   both values are positive for *insertion moves*  
*   both values are negative for *base pair removals*  
*   both values have different signedness for *shift moves*, where the positive value indicates the
    nucleotide that stays constant, and the others absolute value is the new pairing partner  

Note
----
A value of 0 in either field is used as list-end indicator and doesn't represent any valid move.  

Attributes
----------
pos_5 : int  
    The (absolute value of the) 5' position of a base pair, or any position of a shifted pair.  

pos_3 : int  
    The (absolute value of the) 3' position of a base pair, or any position of a shifted pair.  

next : vrna_move_t *  
    The next base pair (if an elementary move changes more than one base pair), or `NULL` Has to be
    terminated with move 0,0.  

C++ includes: ViennaRNA/landscape/move.h
";

// File: structvrna__multimer__pf__s.xml


%feature("docstring") vrna_multimer_pf_s "

Attributes
----------
F_connected : double  
    Fully connected ensemble (incl. DuplexInititiation and rotational symmetry correction.  

F_monomers : double *  
    monomers  

num_monomers : size_t  
    Number of monomers.  
";

// File: structvrna__mx__mfe__s.xml


%feature("docstring") vrna_mx_mfe_s "

Minimum Free Energy (MFE) Dynamic Programming (DP) matrices data structure required within the
RNA.fold_compound().  

Attributes
----------
type : const vrna_mx_type_e  
    Type of the DP matrices  

length : unsigned int  
    Length of the sequence, therefore an indicator of the size of the DP matrices.  

strands : unsigned int  
    Number of strands  

c : int *  
    Energy array, given that i-j pair.  

f5 : int *  
    Energy of 5' end.  

f3 : int *  
    Energy of 3' end.  

fms5 : int **  
    Energy for connected interstrand configurations.  

fms3 : int **  
    nergy for connected interstrand configurations  

fML : int *  
    Multi-loop auxiliary energy array.  

fM1 : int *  
    Second ML array, only for unique multibrnach loop decomposition.  

fM2 : int *  
    Energy for a multibranch loop region with exactly two stems, extending to 3' end.  

ggg : int *  
    Energies of g-quadruplexes.  

Fc : int  
    Minimum Free Energy of entire circular RNA.  

FcH : int  
    Minimum Free Energy of hairpin loop cases in circular RNA.  

FcI : int  
    Minimum Free Energy of internal loop cases in circular RNA.  

FcM : int  
    Minimum Free Energy of multibranch loop cases in circular RNA.  

c_local : int **  
    Energy array, given that i-j pair.  

f3_local : int *  
    Energy of 5' end.  

fML_local : int **  
    Multi-loop auxiliary energy array.  

ggg_local : int **  
    Energies of g-quadruplexes.  

E_F5 : int ***  

l_min_F5 : int **  

l_max_F5 : int **  

k_min_F5 : int *  

k_max_F5 : int *  

E_F3 : int ***  

l_min_F3 : int **  

l_max_F3 : int **  

k_min_F3 : int *  

k_max_F3 : int *  

E_C : int ***  

l_min_C : int **  

l_max_C : int **  

k_min_C : int *  

k_max_C : int *  

E_M : int ***  

l_min_M : int **  

l_max_M : int **  

k_min_M : int *  

k_max_M : int *  

E_M1 : int ***  

l_min_M1 : int **  

l_max_M1 : int **  

k_min_M1 : int *  

k_max_M1 : int *  

E_M2 : int ***  

l_min_M2 : int **  

l_max_M2 : int **  

k_min_M2 : int *  

k_max_M2 : int *  

E_Fc : int **  

l_min_Fc : int *  

l_max_Fc : int *  

k_min_Fc : int  

k_max_Fc : int  

E_FcH : int **  

l_min_FcH : int *  

l_max_FcH : int *  

k_min_FcH : int  

k_max_FcH : int  

E_FcI : int **  

l_min_FcI : int *  

l_max_FcI : int *  

k_min_FcI : int  

k_max_FcI : int  

E_FcM : int **  

l_min_FcM : int *  

l_max_FcM : int *  

k_min_FcM : int  

k_max_FcM : int  

E_F5_rem : int *  

E_F3_rem : int *  

E_C_rem : int *  

E_M_rem : int *  

E_M1_rem : int *  

E_M2_rem : int *  

E_Fc_rem : int  

E_FcH_rem : int  

E_FcI_rem : int  

E_FcM_rem : int  

 : union vrna_mx_mfe_s  

C++ includes: ViennaRNA/dp_matrices.h
";

/*
 Common fields for MFE matrices 
*/

/*
 Default DP matrices 
*/

/*
Note
----
These data fields are available if  

*/

/*
 Local Folding DP matrices using window approach 
*/

/*
Note
----
These data fields are available if  

*/

/*
 Distance Class DP matrices 
*/

/*
Note
----
These data fields are available if  

*/

// File: structvrna__mx__pf__s.xml


%feature("docstring") vrna_mx_pf_s "

Partition function (PF) Dynamic Programming (DP) matrices data structure required within the
RNA.fold_compound().  

Attributes
----------
type : const vrna_mx_type_e  
    Type of the DP matrices  

length : unsigned int  
    Size of the DP matrices (i.e. sequence length)  

scale : FLT_OR_DBL *  
    Boltzmann factor scaling  

expMLbase : FLT_OR_DBL *  
    Boltzmann factors for unpaired bases in multibranch loop  

q : FLT_OR_DBL *  

qb : FLT_OR_DBL *  

qm : FLT_OR_DBL *  

qm1 : FLT_OR_DBL *  

probs : FLT_OR_DBL *  

q1k : FLT_OR_DBL *  

qln : FLT_OR_DBL *  

G : FLT_OR_DBL *  

qo : FLT_OR_DBL  

qm2 : FLT_OR_DBL *  

qho : FLT_OR_DBL  

qio : FLT_OR_DBL  

qmo : FLT_OR_DBL  

q_local : FLT_OR_DBL **  

qb_local : FLT_OR_DBL **  

qm_local : FLT_OR_DBL **  

pR : FLT_OR_DBL **  

qm2_local : FLT_OR_DBL **  

QI5 : FLT_OR_DBL **  

q2l : FLT_OR_DBL **  

qmb : FLT_OR_DBL **  

G_local : FLT_OR_DBL **  

Q : FLT_OR_DBL ***  

l_min_Q : int **  

l_max_Q : int **  

k_min_Q : int *  

k_max_Q : int *  

Q_B : FLT_OR_DBL ***  

l_min_Q_B : int **  

l_max_Q_B : int **  

k_min_Q_B : int *  

k_max_Q_B : int *  

Q_M : FLT_OR_DBL ***  

l_min_Q_M : int **  

l_max_Q_M : int **  

k_min_Q_M : int *  

k_max_Q_M : int *  

Q_M1 : FLT_OR_DBL ***  

l_min_Q_M1 : int **  

l_max_Q_M1 : int **  

k_min_Q_M1 : int *  

k_max_Q_M1 : int *  

Q_M2 : FLT_OR_DBL ***  

l_min_Q_M2 : int **  

l_max_Q_M2 : int **  

k_min_Q_M2 : int *  

k_max_Q_M2 : int *  

Q_c : FLT_OR_DBL **  

l_min_Q_c : int *  

l_max_Q_c : int *  

k_min_Q_c : int  

k_max_Q_c : int  

Q_cH : FLT_OR_DBL **  

l_min_Q_cH : int *  

l_max_Q_cH : int *  

k_min_Q_cH : int  

k_max_Q_cH : int  

Q_cI : FLT_OR_DBL **  

l_min_Q_cI : int *  

l_max_Q_cI : int *  

k_min_Q_cI : int  

k_max_Q_cI : int  

Q_cM : FLT_OR_DBL **  

l_min_Q_cM : int *  

l_max_Q_cM : int *  

k_min_Q_cM : int  

k_max_Q_cM : int  

Q_rem : FLT_OR_DBL *  

Q_B_rem : FLT_OR_DBL *  

Q_M_rem : FLT_OR_DBL *  

Q_M1_rem : FLT_OR_DBL *  

Q_M2_rem : FLT_OR_DBL *  

Q_c_rem : FLT_OR_DBL  

Q_cH_rem : FLT_OR_DBL  

Q_cI_rem : FLT_OR_DBL  

Q_cM_rem : FLT_OR_DBL  

 : union vrna_mx_pf_s  

C++ includes: ViennaRNA/dp_matrices.h
";

/*
 Common fields for DP matrices 
*/

/*
 Default PF matrices 
*/

/*
Note
----
These data fields are available if  

*/

/*
 Local Folding DP matrices using window approach 
*/

/*
Note
----
These data fields are available if  

*/

/*
 Distance Class DP matrices 
*/

/*
Note
----
These data fields are available if  

*/

// File: structvrna__param__s.xml


%feature("docstring") vrna_param_t "

The datastructure that contains temperature scaled energy parameters.  

Attributes
----------
id : int  

stack : int  

hairpin : int  

bulge : int  

internal_loop : int  

mismatchExt : int  

mismatchI : int  

mismatch1nI : int  

mismatch23I : int  

mismatchH : int  

mismatchM : int  

dangle5 : int  

dangle3 : int  

int11 : int  

int21 : int  

int22 : int  

ninio : int  

lxc : double  

MLbase : int  

MLintern : int  

MLclosing : int  

TerminalAU : int  

DuplexInit : int  

Tetraloop_E : int  

Tetraloops : char  

Triloop_E : int  

Triloops : char  

Hexaloop_E : int  

Hexaloops : char  

TripleC : int  

MultipleCA : int  

MultipleCB : int  

gquad : int  

gquadLayerMismatch : int  

gquadLayerMismatchMax : int  

temperature : double  
    Temperature used for loop contribution scaling.  

model_details : vrna_md_t  
    Model details to be used in the recursions.  

param_file : char  
    The filename the parameters were derived from, or empty string if they represent the default.  

SaltStack : int  

SaltLoop : int  

SaltLoopDbl : double  

SaltMLbase : int  

SaltMLintern : int  

SaltMLclosing : int  

SaltDPXInit : int  

C++ includes: ViennaRNA/params/basic.h
";

// File: structvrna__path__s.xml


%feature("docstring") vrna_path_s "

An element of a refolding path list.  

Usually, one has to deal with an array of RNA.path(), e.g. returned from one of the refolding-path
algorithms.  

Since in most cases the length of the list is not known in advance, such lists have an *end-of-list*
marker, which is either:  

*   a value of *NULL* for RNA.path()::s if RNA.path()::type = RNA.PATH_TYPE_DOT_BRACKET, or  
*   a RNA.path()::move with zero in both fields RNA.move()::pos_5 and RNA.move()::pos_3 if
    RNA.path()::type = RNA.PATH_TYPE_MOVES.  

In the following we show an example for how to cover both cases of iteration:  

See Also
--------
RNA.path_free()  

Attributes
----------
type : unsigned int  
    The type of the path element.  

    A value of RNA.PATH_TYPE_DOT_BRACKET indicates that RNA.path()::s consists of the secondary
    structure in dot-bracket notation, and RNA.path()::en the corresponding free energy.  
     On the other hand, if the value is RNA.PATH_TYPE_MOVES, RNA.path()::s is *NULL* and
    RNA.path()::move is set to the transition move that transforms a previous structure into it's
    neighbor along the path. In this case, the attribute RNA.path()::en states the change in free
    energy with respect to the structure before application of RNA.path()::move.  

en : double  
    Free energy of current structure.  

s : char *  
    Secondary structure in dot-bracket notation.  

move : vrna_move_t  
    Move that transforms the previous structure into it's next neighbor along the path.  

C++ includes: ViennaRNA/landscape/paths.h
";

// File: structvrna__pinfo__s.xml


%feature("docstring") vrna_pinfo_s "

A base pair info structure.  

For each base pair (i,j) with i,j in [0, n-1] the structure lists:  

*   its probability 'p'  
*   an entropy-like measure for its well-definedness 'ent'  
*   the frequency of each type of pair in 'bp[]'
    -   'bp[0]' contains the number of non-compatible sequences  
    -   'bp[1]' the number of CG pairs, etc.  

Attributes
----------
i : unsigned  
    nucleotide position i  

j : unsigned  
    nucleotide position j  

p : float  
    Probability.  

ent : float  
    Pseudo entropy for :math:`p(i,j) = S_{i} + S_{j} - p_{i}j*ln(p_{i}j)`.  

bp : short  
    Frequencies of pair_types.  

comp : char  
    1 iff pair is in mfe structure  

C++ includes: ViennaRNA/utils/alignments.h
";

// File: structvrna__pk__plex__result__s.xml


%feature("docstring") vrna_pk_plex_result_s "

A result of the RNA PKplex interaction prediction.  

See Also
--------
RNA.pk_plex()  

Attributes
----------
structure : char *  
    Secondary Structure in dot-bracket notation.  

energy : double  
    Net free energy in kcal/mol.  

dGpk : double  
    Free energy of PK loop in kcal/mol.  

dGint : double  
    Free energy of PK forming duplex interaction.  

dG1 : double  
    Opening energy for the 5' interaction site used in the heuristic.  

dG2 : double  
    Opening energy for the 3' interaction site used in the heuristic.  

start_5 : unsigned int  
    Start coordinate of the 5' interaction site.  

end_5 : unsigned int  
    End coordinate of the 5' interaction site.  

start_3 : unsigned int  
    Start coordinate of the 3' interaction site.  

end_3 : unsigned int  
    End coordinate of the 3' interaction site.  

C++ includes: ViennaRNA/pk_plex.h
";

// File: structvrna__plot__layout__s.xml


%feature("docstring") vrna_plot_layout_s "

Attributes
----------
length : unsigned int  

x : float *  

y : float *  

arcs : double *  

bbox : int  
";

// File: structvrna__plot__options__puzzler__t.xml


%feature("docstring") vrna_plot_options_puzzler_t "

Options data structure for RNApuzzler algorithm implementation.  

Attributes
----------
drawArcs : short  

paired : double  

unpaired : double  

checkAncestorIntersections : short  

checkSiblingIntersections : short  

checkExteriorIntersections : short  

allowFlipping : short  

optimize : short  

maximumNumberOfConfigChangesAllowed : int  

config : char *  

filename : const char *  

numberOfChangesAppliedToConfig : int  

psNumber : int  

C++ includes: ViennaRNA/plotting/RNApuzzler/RNApuzzler.h
";

// File: structvrna__sc__bp__storage__t.xml


%feature("docstring") vrna_sc_bp_storage_t "

A base pair constraint.  

Attributes
----------
interval_start : unsigned int  

interval_end : unsigned int  

e : int  

C++ includes: ViennaRNA/constraints/soft.h
";

// File: structvrna__sc__mod__param__s.xml


%feature("docstring") vrna_sc_mod_param_s "

Attributes
----------
available : unsigned int  

name : char *  

one_letter_code : char  

unmodified : char  

fallback : char  

pairing_partners : char  

pairing_partners_encoding : unsigned int  

unmodified_encoding : unsigned int  

fallback_encoding : unsigned int  

num_ptypes : size_t  

ptypes : size_t  

stack_dG : int  

stack_dH : int  

dangle5_dG : int  

dangle5_dH : int  

dangle3_dG : int  

dangle3_dH : int  

mismatch_dG : int  

mismatch_dH : int  

terminal_dG : int  

terminal_dH : int  
";

// File: structvrna__sc__motif__s.xml


%feature("docstring") vrna_sc_motif_s "

Attributes
----------
i : int  

j : int  

k : int  

l : int  

number : int  
";

// File: structvrna__sc__s.xml


%feature("docstring") vrna_sc_t "

The soft constraints data structure.  

Attributes
----------
type : const vrna_sc_type_e  

n : unsigned int  

state : unsigned char  

energy_up : int **  
    Energy contribution for stretches of unpaired nucleotides.  

exp_energy_up : FLT_OR_DBL **  
    Boltzmann Factors of the energy contributions for unpaired sequence stretches.  

up_storage : int *  
    Storage container for energy contributions per unpaired nucleotide.  

bp_storage : vrna_sc_bp_storage_t **  
    Storage container for energy contributions per base pair.  

energy_bp : int *  
    Energy contribution for base pairs.  

exp_energy_bp : FLT_OR_DBL *  
    Boltzmann Factors of the energy contribution for base pairs.  

energy_bp_local : int **  
    Energy contribution for base pairs (sliding window approach)  

exp_energy_bp_local : FLT_OR_DBL **  
    Boltzmann Factors of the energy contribution for base pairs (sliding window approach)  

 : union vrna_sc_s  

energy_stack : int *  
    Pseudo Energy contribution per base pair involved in a stack.  

exp_energy_stack : FLT_OR_DBL *  
    Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack.  

f : vrna_sc_f  
    A function pointer used for pseudo energy contribution in MFE calculations.  

    See Also
    --------
    RNA.fold_compound.sc_add()  

bt : vrna_sc_bt_f  
    A function pointer used to obtain backtraced base pairs in loop regions that were altered by
    soft constrained pseudo energy contributions.  

    See Also
    --------
    RNA.fold_compound.sc_add_bt()  

exp_f : vrna_sc_exp_f  
    A function pointer used for pseudo energy contribution boltzmann factors in PF calculations.  

    See Also
    --------
    RNA.fold_compound.sc_add_exp()  

data : void *  
    A pointer to the data object provided for for pseudo energy contribution functions of the
    generic soft constraints feature.  

prepare_data : vrna_auxdata_prepare_f  

free_data : vrna_auxdata_free_f  

C++ includes: ViennaRNA/constraints/soft.h
";

// File: structvrna__sect__s.xml


%feature("docstring") vrna_sect_s "

Stack of partial structures for backtracking.  

Attributes
----------
i : int  

j : int  

ml : int  

C++ includes: ViennaRNA/datastructures/basic.h
";

// File: structvrna__sequence__s.xml


%feature("docstring") vrna_sequence_s "

Data structure representing a nucleotide sequence.  

Attributes
----------
type : vrna_seq_type_e  
    The type of sequence.  

name : char *  

string : char *  
    The string representation of the sequence.  

encoding : short *  
    The integer representation of the sequence.  

encoding5 : short *  

encoding3 : short *  

length : unsigned int  
    The length of the sequence.  

C++ includes: ViennaRNA/sequence.h
";

// File: structvrna__sol__TwoD__pf__t.xml


%feature("docstring") vrna_sol_TwoD_pf_t "

Solution element returned from RNA.pf_TwoD()  

This element contains the partition function for the appropriate kappa (k), lambda (l) neighborhood
The datastructure contains two integer attributes 'k' and 'l' as well as an attribute 'q' of type
FLT_OR_DBL  

A value of INF in k denotes the end of a list  

See Also
--------
RNA.pf_TwoD()  

Attributes
----------
k : int  
    Distance to first reference.  

l : int  
    Distance to second reference.  

q : FLT_OR_DBL  
    partition function  

C++ includes: ViennaRNA/2Dpfold.h
";

// File: structvrna__sol__TwoD__t.xml


%feature("docstring") vrna_sol_TwoD_t "

Solution element returned from RNA.mfe_TwoD()  

This element contains free energy and structure for the appropriate kappa (k), lambda (l)
neighborhood The datastructure contains two integer attributes 'k' and 'l' as well as an attribute
'en' of type float representing the free energy in kcal/mol and an attribute 's' of type char*
containg the secondary structure representative,  

A value of INF in k denotes the end of a list  

See Also
--------
RNA.mfe_TwoD()  

Attributes
----------
k : int  
    Distance to first reference.  

l : int  
    Distance to second reference.  

en : float  
    Free energy in kcal/mol.  

s : char *  
    MFE representative structure in dot-bracket notation.  

C++ includes: ViennaRNA/2Dfold.h
";

// File: structvrna__string__header__s.xml


%feature("docstring") vrna_string_header_s "

The header of an array.  

Attributes
----------
len : size_t  
    The length of the string.  

size : size_t  
    The actual capacity of an array.  

shift_post : size_t  

backup : char  

C++ includes: ViennaRNA/datastructures/string.h
";

// File: structvrna__structured__domains__s.xml


%feature("docstring") vrna_structured_domains_s "

Attributes
----------
__placeholder : char  
";

// File: structvrna__subopt__sol__s.xml


%feature("docstring") vrna_subopt_sol_s "

Solution element from subopt.c.  

Attributes
----------
energy : float  
    Free Energy of structure in kcal/mol.  

structure : char *  
    Structure in dot-bracket notation.  

C++ includes: ViennaRNA/subopt.h
";

// File: structvrna__unstructured__domain__motif__s.xml


%feature("docstring") vrna_unstructured_domain_motif_s "

Attributes
----------
start : int  

number : int  
";

// File: structvrna__unstructured__domain__s.xml


%feature("docstring") vrna_unstructured_domain_s "

Data structure to store all functionality for ligand binding.  

Attributes
----------
uniq_motif_count : int  
    The unique number of motifs of different lengths.  

uniq_motif_size : unsigned int *  
    An array storing a unique list of motif lengths.  

motif_count : int  
    Total number of distinguished motifs.  

motif : char **  
    Motif sequences.  

motif_name : char **  
    Motif identifier/name.  

motif_size : unsigned int *  
    Motif lengths.  

motif_en : double *  
    Ligand binding free energy contribution.  

motif_type : unsigned int *  
    Type of motif, i.e. loop type the ligand binds to.  

prod_cb : vrna_ud_production_f  
    Callback to ligand binding production rule, i.e. create/fill DP free energy matrices.  

    This callback will be executed right before the actual secondary structure decompositions, and,
    therefore, any implementation must not interleave with the regular DP matrices.  

exp_prod_cb : vrna_ud_exp_production_f  
    Callback to ligand binding production rule, i.e. create/fill DP partition function matrices.  

energy_cb : vrna_ud_f  
    Callback to evaluate free energy of ligand binding to a particular unpaired stretch.  

exp_energy_cb : vrna_ud_exp_f  
    Callback to evaluate Boltzmann factor of ligand binding to a particular unpaired stretch.  

data : void *  
    Auxiliary data structure passed to energy evaluation callbacks.  

free_data : vrna_auxdata_free_f  
    Callback to free auxiliary data structure.  

probs_add : vrna_ud_add_probs_f  
    Callback to store/add outside partition function.  

probs_get : vrna_ud_get_probs_f  
    Callback to retrieve outside partition function.  

C++ includes: ViennaRNA/unstructured_domains.h
";

// File: 00-layout_8dox.xml

// File: distance__measures_8dox.xml

// File: main_8dox.xml

// File: doc_2doxygen_2refman_8include_2plotting_8dox.xml

// File: interfaces_2plotting_8dox.xml

// File: aln__utils_8dox.xml

// File: basic__algorithms_8dox.xml

// File: boltzmann__sampling_8dox.xml

// File: combinatorics_8dox.xml

// File: commands_8dox.xml

// File: constraints__hard_8dox.xml

// File: constraints__ligand_8dox.xml

// File: constraints__SHAPE_8dox.xml

// File: constraints__soft_8dox.xml

// File: eval_8dox.xml

// File: file__formats_8dox.xml

// File: fold__compound_8dox.xml

// File: grammar_8dox.xml

// File: heat__capacity_8dox.xml

// File: mfe_8dox.xml

// File: model__details_8dox.xml

// File: neighbor_8dox.xml

// File: params_8dox.xml

// File: part__func_8dox.xml

// File: paths_8dox.xml

// File: sequence_8dox.xml

// File: structure__utils_8dox.xml

// File: subopt_8dox.xml

// File: utils_8dox.xml

// File: walk_8dox.xml

// File: 2Dfold_8h.xml

// File: 2Dpfold_8h.xml

%feature("docstring") TwoDpfold_solution "
";

%feature("docstring") get_TwoDpfold_variables "

Get a datastructure containing all necessary attributes and global folding switches.  

This function prepares all necessary attributes and matrices etc which are needed for a call of
TwoDpfold() . A snapshot of all current global model switches (dangles, temperature and so on) is
done and stored in the returned datastructure. Additionally, all matrices that will hold the
partition function values are prepared.  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.pf_TwoD(), and RNA.fold_compound_free() instead!  

Parameters
----------
seq : const char *
    the RNA sequence in uppercase format with letters from the alphabet {AUCG}  
structure1 : const char *
    the first reference structure in dot-bracket notation  
structure2 : char *
    the second reference structure in dot-bracket notation  
circ : int
    a switch indicating if the sequence is linear (0) or circular (1)  

Returns
-------
TwoDpfold_vars *  
    the datastructure containing all necessary partition function attributes  
";

%feature("docstring") destroy_TwoDpfold_variables "

Free all memory occupied by a TwoDpfold_vars datastructure.  

This function free's all memory occupied by a datastructure obtained from from
get_TwoDpfold_variabless() or get_TwoDpfold_variables_from_MFE()  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.pf_TwoD(), and RNA.fold_compound_free() instead!  

Parameters
----------
vars : TwoDpfold_vars *
    the datastructure to be free'd  

See Also
--------
get_TwoDpfold_variables(), get_TwoDpfold_variables_from_MFE()  
";

%feature("docstring") TwoDpfoldList "

Compute the partition function for all distance classes.  

This function computes the partition functions for all distance classes according the two reference
structures specified in the datastructure 'vars'. Similar to TwoDfold() the arguments maxDistance1
and maxDistance2 specify the maximum distance to both reference structures. A value of '-1' in
either of them makes the appropriate distance restrictionless, i.e. all basepair distancies to the
reference are taken into account during computation. In case there is a restriction, the returned
solution contains an entry where the attribute k=l=-1 contains the partition function for all
structures exceeding the restriction. A values of INF in the attribute 'k' of the returned list
denotes the end of the list  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.pf_TwoD(), and RNA.fold_compound_free() instead!  

Parameters
----------
vars : TwoDpfold_vars *
    the datastructure containing all necessary folding attributes and matrices  
maxDistance1 : int
    the maximum basepair distance to reference1 (may be -1)  
maxDistance2 : int
    the maximum basepair distance to reference2 (may be -1)  

Returns
-------
TwoDpfold_solution *  
    a list of partition funtions for the appropriate distance classes  

See Also
--------
get_TwoDpfold_variables(), destroy_TwoDpfold_variables(), RNA.sol_TwoD_pf()  
";

%feature("docstring") TwoDpfold_pbacktrack "

Sample secondary structure representatives from a set of distance classes according to their
Boltzmann probability.  

If the argument 'd1' is set to '-1', the structure will be backtracked in the distance class where
all structures exceeding the maximum basepair distance to either of the references reside.  

**Precondition**
    The argument 'vars' must contain precalculated partition function matrices, i.e. a call to
    TwoDpfold() preceding this function is mandatory!  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.pf_TwoD(), RNA.pbacktrack_TwoD(), and RNA.fold_compound_free()
    instead!  

Parameters
----------
vars : TwoDpfold_vars *
    the datastructure containing all necessary folding attributes and matrices  
d1 : int
    the distance to reference1 (may be -1)  
d2 : int
    the distance to reference2  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  

See Also
--------
TwoDpfold()  
";

%feature("docstring") TwoDpfold_pbacktrack5 "

Sample secondary structure representatives with a specified length from a set of distance classes
according to their Boltzmann probability.  

This function does essentially the same as TwoDpfold_pbacktrack() with the only difference that
partial structures, i.e. structures beginning from the 5' end with a specified length of the
sequence, are backtracked  

**Precondition**
    The argument 'vars' must contain precalculated partition function matrices, i.e. a call to
    TwoDpfold() preceding this function is mandatory!  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.pf_TwoD(), RNA.pbacktrack5_TwoD(), and
    RNA.fold_compound_free() instead!  

Note
----
This function does not work (since it makes no sense) for circular RNA sequences!  

Parameters
----------
vars : TwoDpfold_vars *
    the datastructure containing all necessary folding attributes and matrices  
d1 : int
    the distance to reference1 (may be -1)  
d2 : int
    the distance to reference2  
length : unsigned int
    the length of the structure beginning from the 5' end  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  

See Also
--------
TwoDpfold_pbacktrack(), TwoDpfold()  
";

%feature("docstring") TwoDpfold "
";

%feature("docstring") TwoDpfold_circ "
";

// File: ali__plex_8h.xml

%feature("docstring") aliLduplexfold "

aliLduplexfold computes the duplexes between two alignments  
";

%feature("docstring") aliLduplexfold_XS "

aliLduplexfold computes the duplexes between two alignments. It also takes the average accessibility
into account  
";

// File: alifold_8h.xml

%feature("docstring") energy_of_alistruct "

Calculate the free energy of a consensus structure given a set of aligned sequences.  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! Use RNA.fold_compound.eval_structure(), and
    RNA.fold_compound.eval_covar_structure() instead!  

Parameters
----------
sequences : const char **
    The NULL terminated array of sequences  
structure : const char *
    The consensus structure  
n_seq : int
    The number of sequences in the alignment  
energy : float *
    A pointer to an array of at least two floats that will hold the free energies (energy[0] will
    contain the free energy, energy[1] will be filled with the covariance energy term)  

Returns
-------
float  
    free energy in kcal/mol  
";

%feature("docstring") energy_of_ali_gquad_structure "
";

%feature("docstring") update_alifold_params "

Update the energy parameters for alifold function.  

Call this to recalculate the pair matrix and energy parameters after a change in folding parameters
like temperature  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! The new API uses RNA.fold_compound() to lump all folding
    related necessities together, including the energy parameters. Use RNA.update_fold_params() to
    update the energy parameters within a RNA.fold_compound().  
";

// File: aln__util_8h.xml

// File: alphabet_8h.xml

%feature("docstring") get_ptypes "
";

// File: boltzmann__sampling_8h.xml

// File: centroid_8h.xml

%feature("docstring") get_centroid_struct_pl "

Get the centroid structure of the ensemble.  

.. deprecated:: 2.6.4
    This function was renamed to RNA.centroid_from_plist()  
";

%feature("docstring") get_centroid_struct_pr "

Get the centroid structure of the ensemble.  

.. deprecated:: 2.6.4
    This function was renamed to RNA.centroid_from_probs()  
";

// File: char__stream_8h.xml

// File: datastructures_2char__stream_8h.xml

// File: cofold_8h.xml

// File: combinatorics_8h.xml

// File: commands_8h.xml

// File: concentrations_8h.xml

%feature("docstring") get_concentrations "

Given two start monomer concentrations a and b, compute the concentrations in thermodynamic
equilibrium of all dimers and the monomers.  

This function takes an array 'startconc' of input concentrations with alternating entries for the
initial concentrations of molecules A and B (terminated by two zeroes), then computes the resulting
equilibrium concentrations from the free energies for the dimers. Dimer free energies should be the
dimer-only free energies, i.e. the FcAB entries from the RNA.dimer_pf() struct.  

.. deprecated:: 2.6.4  

Parameters
----------
FEAB : double
    Free energy of AB dimer (FcAB entry)  
FEAA : double
    Free energy of AA dimer (FcAB entry)  
FEBB : double
    Free energy of BB dimer (FcAB entry)  
FEA : double
    Free energy of monomer A  
FEB : double
    Free energy of monomer B  
startconc : double *
    List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]  

Returns
-------
RNA.dimer_conc() *  
    RNA.dimer_conc() array containing the equilibrium energies and start concentrations  
";

// File: constraints_8h.xml

// File: hard_8h.xml

%feature("docstring") CONSTRAINT_NO_HEADER "

do not print the header information line  

.. deprecated:: 2.6.4
    This mode is not supported anymore!  
";

%feature("docstring") CONSTRAINT_DB_ANG_BRACK "

angle brackets '<', '>' switch for structure constraint (paired downstream/upstream)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_CANONICAL_BP "
";

%feature("docstring") CONSTRAINT_CONTEXT_ENFORCE "

Hard constraint flag to indicate enforcement of constraints.  
";

%feature("docstring") CONSTRAINT_CONTEXT_NO_REMOVE "

Hard constraint flag to indicate not to remove base pairs that conflict with a given constraint.  
";

%feature("docstring") CONSTRAINT_CONTEXT_NONE "

Constraint context flag that forbids a nucleotide or base pair to appear in any loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_CLOSING_LOOPS "

Constraint context flag indicating base pairs that close any loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_ENCLOSED_LOOPS "

Constraint context flag indicating base pairs enclosed by any loop.  
";

%feature("docstring") vrna_hc_init_window "
";

%feature("docstring") vrna_hc_prepare "
";

%feature("docstring") vrna_hc_update "
";

%feature("docstring") vrna_hc_add_up_strand "
";

%feature("docstring") vrna_hc_add_up_strand_batch "
";

%feature("docstring") vrna_hc_add_bp_strand "
";

%feature("docstring") vrna_hc_add_f "

Add a function pointer pointer for the generic hard constraint feature.  
";

%feature("docstring") vrna_hc_add_data "

Add an auxiliary data structure for the generic hard constraints callback function.  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound the generic hard constraint function should be bound to  
data : void *
    A pointer to the data structure that holds required data for function 'f'  
f : RNA.auxdata_free
    A pointer to a function that free's the memory occupied by `data` (Maybe `NULL`)  

See Also
--------
RNA.hc_add()  
";

%feature("docstring") print_tty_constraint "

Print structure constraint characters to stdout. (constraint support is specified by option
parameter)  

.. deprecated:: 2.6.4
    Use RNA.message_constraints() instead!  

Parameters
----------
option : unsigned int
    Option switch that tells which constraint help will be printed  
";

%feature("docstring") print_tty_constraint_full "

Print structure constraint characters to stdout (full constraint support)  

.. deprecated:: 2.6.4
    Use RNA.message_constraint_options_all() instead!  
";

%feature("docstring") constrain_ptypes "

Insert constraining pair types according to constraint structure string.  

.. deprecated:: 2.6.4
    Do not use this function anymore! Structure constraints are now handled through RNA.hc() and
    related functions.  

Parameters
----------
constraint : const char *
    The structure constraint string  
length : unsigned int
    The actual length of the sequence (constraint may be shorter)  
ptype : char *
    A pointer to the basepair type array  
BP : int *
    (not used anymore)  
min_loop_size : int
    The minimal loop size (usually TURN )  
idx_type : unsigned int
    Define the access type for base pair type array (0 = indx, 1 = iindx)  
";

// File: ligand_8h.xml

// File: sc__cb__intern_8h.xml

%feature("docstring") MOD_PARAMS_STACK_dG "
";

%feature("docstring") MOD_PARAMS_STACK_dH "
";

%feature("docstring") MOD_PARAMS_MISMATCH_dG "
";

%feature("docstring") MOD_PARAMS_MISMATCH_dH "
";

%feature("docstring") MOD_PARAMS_TERMINAL_dG "
";

%feature("docstring") MOD_PARAMS_TERMINAL_dH "
";

%feature("docstring") MOD_PARAMS_DANGLES_dG "
";

%feature("docstring") MOD_PARAMS_DANGLES_dH "
";

%feature("docstring") MAX_ALPHABET "
";

%feature("docstring") MAX_PAIRS "
";

// File: SHAPE_8h.xml

%feature("docstring") vrna_sc_SHAPE_parse_method "

Parse a character string and extract the encoded SHAPE reactivity conversion method and possibly the
parameters for conversion into pseudo free energies.  

Parameters
----------
method_string : const char *
    The string that contains the encoded SHAPE reactivity conversion method  
method : char *
    A pointer to the memory location where the method character will be stored  
param_1 : float *
    A pointer to the memory location where the first parameter of the corresponding method will be
    stored  
param_2 : float *
    A pointer to the memory location where the second parameter of the corresponding method will be
    stored  

Returns
-------
int  
    1 on successful extraction of the method, 0 on errors  
";

// File: soft_8h.xml

%feature("docstring") vrna_sc_prepare "
";

%feature("docstring") vrna_sc_update "
";

%feature("docstring") vrna_fold_compound_t::sc_set_stack "
";

%feature("docstring") vrna_sc_set_stack_comparative "
";

%feature("docstring") vrna_fold_compound_t::sc_add_stack "
";

%feature("docstring") vrna_sc_add_stack_comparative "
";

%feature("docstring") vrna_sc_add_auxdata "
";

%feature("docstring") vrna_sc_add_data_comparative "
";

%feature("docstring") vrna_sc_multi_cb_add "
";

%feature("docstring") vrna_sc_add_f_comparative "
";

%feature("docstring") vrna_sc_add_exp_f_comparative "
";

// File: soft__special_8h.xml

// File: constraints__hard_8h.xml

// File: constraints__ligand_8h.xml

// File: constraints__SHAPE_8h.xml

// File: constraints__soft_8h.xml

// File: convert__epars_8h.xml

// File: data__structures_8h.xml

// File: array_8h.xml

// File: hash__tables_8h.xml

/*
 Abstract interface 
*/

/*
 Dot-Bracket / Free Energy entries 
*/

// File: heap_8h.xml

// File: lists_8h.xml

%feature("docstring") LST_USERSPACE "
";

%feature("docstring") LST_HEADER "
";

%feature("docstring") LST_HEAD "
";

%feature("docstring") LST_EMPTY "
";

%feature("docstring") lst_newnode "
";

%feature("docstring") lst_freenode "
";

%feature("docstring") lst_init "
";

%feature("docstring") lst_kill "
";

%feature("docstring") lst_insertafter "
";

%feature("docstring") lst_deletenext "
";

%feature("docstring") lst_first "
";

%feature("docstring") lst_next "
";

%feature("docstring") lst_mergesort "
";

// File: string_8h.xml

// File: dist__vars_8h.xml

// File: dp__matrices_8h.xml

// File: duplex_8h.xml

%feature("docstring") my_duplexfold "
";

%feature("docstring") my_duplex_subopt "
";

%feature("docstring") my_aliduplexfold "
";

%feature("docstring") my_aliduplex_subopt "
";

// File: edit__cost_8h.xml

%feature("docstring") PRIVATE "
";

%feature("docstring") DIST_INF "
";

// File: energy__const_8h.xml

// File: energy__par_8h.xml

// File: equilibrium__probs_8h.xml

/*
 Base pair probabilities and derived computations 
*/

/*
 Multimer probabilities computations 
*/

/*
 Structure probability computations 
*/

// File: eval_8h.xml

/*
 Basic Energy Evaluation Interface with Dot-Bracket Structure String 
*/

/*
 Basic Energy Evaluation Interface with Structure Pair Table 
*/

/*
 Simplified Energy Evaluation with Sequence and Dot-Bracket Strings 
*/

/*
 Simplified Energy Evaluation with Sequence Alignments and Consensus Structure Dot-Bracket String 
*/

/*
 Simplified Energy Evaluation with Sequence String and Structure Pair Table 
*/

/*
 Simplified Energy Evaluation with Sequence Alignment and Consensus Structure Pair Table 
*/

// File: exterior__loops_8h.xml

// File: file__formats_8h.xml

// File: io_2file__formats_8h.xml

// File: file__formats__msa_8h.xml

// File: io_2file__formats__msa_8h.xml

// File: file__utils_8h.xml

// File: findpath_8h.xml

// File: landscape_2findpath_8h.xml

// File: fold_8h.xml

// File: fold__compound_8h.xml

// File: fold__vars_8h.xml

// File: gquad_8h.xml

%feature("docstring") INLINE "
";

// File: grammar_8h.xml

// File: hairpin__loops_8h.xml

// File: heat__capacity_8h.xml

/*
 Basic heat capacity function interface 
*/

/*
 Simplified heat capacity computation 
*/

// File: interior__loops_8h.xml

// File: inverse_8h.xml

// File: move_8h.xml

// File: paths_8h.xml

// File: Lfold_8h.xml

// File: loop__energies_8h.xml

// File: all_8h.xml

// File: external_8h.xml

/*
 Boltzmann weight (partition function) interface 
*/

/*
 Basic free energy interface 
*/

// File: hairpin_8h.xml

/*
 Basic free energy interface 
*/

/*
 Boltzmann weight (partition function) interface 
*/

%feature("docstring") INLINE "
";

// File: internal_8h.xml

/*
 Basic free energy interface 
*/

/*
 Boltzmann weight (partition function) interface 
*/

%feature("docstring") INLINE "
";

// File: multibranch_8h.xml

/*
 Boltzmann weight (partition function) interface 
*/

/*
 Basic free energy interface 
*/

%feature("docstring") INLINE "
";

// File: LPfold_8h.xml

%feature("docstring") putoutpU_prob_par "
";

%feature("docstring") putoutpU_prob_bin_par "
";

%feature("docstring") init_pf_foldLP "

Dunno if this function was ever used by external programs linking to RNAlib, but it was declared
PUBLIC before. Anyway, never use this function as it will be removed soon and does nothing at all  
";

// File: MEA_8h.xml

%feature("docstring") MEA_seq "
";

// File: mfe_8h.xml

/*
 Basic global MFE prediction interface 
*/

/*
 Simplified global MFE prediction using sequence(s) or multiple sequence alignment(s) 
*/

// File: mfe__window_8h.xml

/*
 Basic local (sliding window) MFE prediction interface 
*/

/*
 Simplified local MFE prediction using sequence(s) or multiple sequence alignment(s) 
*/

// File: mm_8h.xml

%feature("docstring") vrna_maximum_matching "

**SWIG Wrapper Notes**
    This function is attached as method `maximum_matching()` to objects of type `fold_compound`. See
    e.g.  :py:meth:`RNA.fold_compound.maximum_matching()` in the :doc:`/api_python`.  
";

%feature("docstring") my_maximum_matching "

**SWIG Wrapper Notes**
    This function is available as global function `maximum_matching()`. See e.g.
    :py:func:`RNA.maximum_matching()` in the :doc:`/api_python`.  
";

%feature("docstring") maximumMatching "
";

%feature("docstring") maximumMatchingConstraint "
";

%feature("docstring") maximumMatching2Constraint "
";

// File: model_8h.xml

// File: move__set_8h.xml

%feature("docstring") print_stren "
";

%feature("docstring") print_str "
";

%feature("docstring") copy_arr "
";

%feature("docstring") allocopy "
";

%feature("docstring") move_gradient "
";

%feature("docstring") move_first "
";

%feature("docstring") move_adaptive "
";

%feature("docstring") move_standard "
";

%feature("docstring") browse_neighs_pt "
";

%feature("docstring") browse_neighs "
";

// File: multibranch__loops_8h.xml

// File: naview_8h.xml

// File: landscape_2neighbor_8h.xml

// File: neighbor_8h.xml

// File: pair__mat_8h.xml

%feature("docstring") NBASES "
";

%feature("docstring") INLINE "
";

%feature("docstring") MAXALPHA "
";

%feature("docstring") ENCODE "
";

%feature("docstring") encode_char "
";

%feature("docstring") make_pair_matrix "
";

%feature("docstring") encode_sequence "
";

// File: params_8h.xml

// File: 1_88_84__epars_8h.xml

%feature("docstring") K0 "
";

%feature("docstring") INF "
";

%feature("docstring") NBPAIRS "
";

%feature("docstring") NST "
";

%feature("docstring") DEF "
";

%feature("docstring") NSM "
";

// File: 1_88_84__intloops_8h.xml

// File: constraints_2basic_8h.xml

%feature("docstring") DECOMP_PAIR_ML_EXT "
";

%feature("docstring") DECOMP_PAIR_ML_OUTSIDE "
";

%feature("docstring") DECOMP_EXT_STEM_EXT1 "
";

%feature("docstring") DECOMP_EXT_L "
";

%feature("docstring") DECOMP_EXT_EXT_L "
";

%feature("docstring") DECOMP_TYPES_MAX "
";

// File: datastructures_2basic_8h.xml

// File: params_2basic_8h.xml

// File: utils_2basic_8h.xml

%feature("docstring") get_indx "
";

%feature("docstring") get_iindx "
";

%feature("docstring") get_line "

Read a line of arbitrary length from a stream.  

Returns a pointer to the resulting string. The necessary memory is allocated and should be released
using *free()* when the string is no longer needed.  

.. deprecated:: 2.6.4
    Use RNA.read_line() as a substitute!  

Parameters
----------
fp : FILE *
    A file pointer to the stream where the function should read from  

Returns
-------
char *  
    A pointer to the resulting string  
";

%feature("docstring") print_tty_input_seq "

Print a line to *stdout* that asks for an input sequence.  

There will also be a ruler (scale line) printed that helps orientation of the sequence positions  

.. deprecated:: 2.6.4
    Use RNA.message_input_seq_simple() instead!  
";

%feature("docstring") print_tty_input_seq_str "

Print a line with a user defined string and a ruler to stdout.  

(usually this is used to ask for user input) There will also be a ruler (scale line) printed that
helps orientation of the sequence positions  

.. deprecated:: 2.6.4
    Use RNA.message_input_seq() instead!  
";

%feature("docstring") warn_user "

Print a warning message.  

Print a warning message to *stderr*  

.. deprecated:: 2.3.0
    Use RNA.message_warning() instead! (since v2.3.0)  
";

%feature("docstring") nrerror "

Die with an error message.  

.. deprecated:: 2.3.0
    Use RNA.message_error() instead! (since v2.3.0)  
";

%feature("docstring") space "

Allocate space safely.  

.. deprecated:: 2.2.0
    Use RNA.alloc() instead! (since v2.2.0)  
";

%feature("docstring") xrealloc "

Reallocate space safely.  

.. deprecated:: 2.2.0
    Use RNA.realloc() instead! (since v2.2.0)  
";

%feature("docstring") init_rand "

Make random number seeds.  

.. deprecated:: 2.6.4
    Use RNA.init_rand() instead!  
";

%feature("docstring") urn "

get a random number from [0..1]  

.. deprecated:: 2.6.4
    Use RNA.urn() instead!  
";

%feature("docstring") int_urn "

Generates a pseudo random integer in a specified range.  

.. deprecated:: 2.6.4
    Use RNA.int_urn() instead!  
";

%feature("docstring") filecopy "

Inefficient `cp`  

.. deprecated:: 2.6.4
    Use RNA.file_copy() instead!  
";

%feature("docstring") time_stamp "

Get a timestamp.  

.. deprecated:: 2.6.4
    Use RNA.time_stamp() instead!  
";

// File: constants_8h.xml

%feature("docstring") GASCONST "

The gas constant  
";

%feature("docstring") K0 "

0 deg Celsius in Kelvin  
";

%feature("docstring") INF "

Infinity as used in minimization routines  
";

%feature("docstring") EMAX "
";

%feature("docstring") FORBIDDEN "

forbidden  
";

%feature("docstring") BONUS "

bonus contribution  
";

%feature("docstring") NBPAIRS "

The number of distinguishable base pairs  
";

%feature("docstring") TURN "

The minimum loop length  
";

%feature("docstring") MAXLOOP "

The maximum loop length  
";

%feature("docstring") UNIT "
";

%feature("docstring") MINPSCORE "
";

// File: convert_8h.xml

// File: default_8h.xml

%feature("docstring") PUBLIC "
";

// File: intl11_8h.xml

// File: intl11dH_8h.xml

// File: intl21_8h.xml

// File: intl21dH_8h.xml

// File: intl22_8h.xml

// File: intl22dH_8h.xml

// File: io_8h.xml

// File: salt_8h.xml

// File: part__func_8h.xml

/*
 Basic global partition function interface 
*/

/*
 Simplified global partition function computation using sequence(s) or multiple sequence alignment(s) 
*/

%feature("docstring") centroid "

.. deprecated:: 2.6.4
    This function is deprecated and should not be used anymore as it is not threadsafe!  

See Also
--------
get_centroid_struct_pl(), get_centroid_struct_pr()  
";

%feature("docstring") get_centroid_struct_gquad_pr "

.. deprecated:: 2.6.4
    This function is deprecated and should not be used anymore as it is not threadsafe!  

See Also
--------
RNA.fold_compound.centroid(), RNA.centroid_from_probs(), RNA.centroid_from_plist()  
";

%feature("docstring") mean_bp_dist "

get the mean pair distance of ensemble  

.. deprecated:: 2.6.4
    This function is not threadsafe and should not be used anymore. Use mean_bp_distance() instead!  
";

%feature("docstring") expLoopEnergy "

.. deprecated:: 2.6.4
    Use exp_E_IntLoop() from loop_energies.h instead  
";

%feature("docstring") expHairpinEnergy "

.. deprecated:: 2.6.4
    Use exp_E_Hairpin() from loop_energies.h instead  
";

%feature("docstring") assign_plist_gquad_from_pr "
";

// File: part__func__co_8h.xml

%feature("docstring") get_plist "

DO NOT USE THIS FUNCTION ANYMORE  

.. deprecated:: 2.6.4
    use assign_plist_from_pr() instead!  
";

// File: part__func__up_8h.xml

%feature("docstring") RNA_UP_MODE_1 "
";

%feature("docstring") RNA_UP_MODE_2 "
";

%feature("docstring") RNA_UP_MODE_3 "
";

// File: part__func__window_8h.xml

/*
 Basic local partition function interface 
*/

/*
 Simplified global partition function computation using sequence(s) or multiple sequence alignment(s) 
*/

// File: perturbation__fold_8h.xml

// File: pf__multifold_8h.xml

%feature("docstring") vrna_pf_multifold_prepare "
";

// File: pk__plex_8h.xml

// File: PKplex_8h.xml

%feature("docstring") PKLduplexfold_XS "
";

// File: plex_8h.xml

%feature("docstring") Lduplexfold "

Lduplexfold Computes duplexes between two single sequences  
";

%feature("docstring") Lduplexfold_XS "

Lduplexfold_XS Computes duplexes between two single sequences with accessibility  
";

%feature("docstring") Lduplexfold_C "

Lduplexfold_C Computes duplexes between two single sequences and takes constraint into account  
";

%feature("docstring") Lduplexfold_CXS "

Lduplexfold_CXS Computes duplexes between two single sequences and takes constraint as well as
accessibility into account  
";

%feature("docstring") arraySize "
";

%feature("docstring") freeDuplexT "
";

// File: plot__aln_8h.xml

// File: plot__layouts_8h.xml

// File: plot__structure_8h.xml

// File: plot__utils_8h.xml

// File: plotting_2alignments_8h.xml

// File: utils_2alignments_8h.xml

// File: layouts_8h.xml

// File: probabilities_8h.xml

// File: RNApuzzler_8h.xml

// File: RNAturtle_8h.xml

// File: plotting_2structures_8h.xml

// File: utils_2structures_8h.xml

// File: ProfileAln_8h.xml

%feature("docstring") profile_aln "
";

%feature("docstring") set_paln_params "
";

// File: profiledist_8h.xml

%feature("docstring") profile_edit_distance "

Align the 2 probability profiles T1, T2  
.  

This is like a Needleman-Wunsch alignment, we should really use affine gap-costs ala Gotoh  
";

%feature("docstring") Make_bp_profile_bppm "

condense pair probability matrix into a vector containing probabilities for unpaired, upstream
paired and downstream paired.  

This resulting probability profile is used as input for profile_edit_distance  

Parameters
----------
bppm : FLT_OR_DBL *
    A pointer to the base pair probability matrix  
length : int
    The length of the sequence  

Returns
-------
float *  
    The bp profile  
";

%feature("docstring") print_bppm "

print string representation of probability profile  
";

%feature("docstring") free_profile "

free space allocated in Make_bp_profile  

Backward compatibility only. You can just use plain free()  
";

%feature("docstring") Make_bp_profile "

.. deprecated:: 2.6.4
    This function is deprecated and will be removed soon! See Make_bp_profile_bppm() for a
    replacement  

See Also
--------
Make_bp_profile_bppm()  

Note
----
This function is NOT threadsafe  
";

// File: PS__dot_8h.xml

// File: read__epars_8h.xml

// File: ribo_8h.xml

// File: RNAstruct_8h.xml

// File: BoyerMoore_8h.xml

// File: sequence_8h.xml

// File: snofold_8h.xml

%feature("docstring") MISMATCH "
";

%feature("docstring") snofold "

snofold is the stem folding array for RNAsnoop  
";

%feature("docstring") snofree_arrays "

Free arrays and structure related to snofold  
";

%feature("docstring") snoinitialize_fold "
";

%feature("docstring") snoupdate_fold_params "
";

%feature("docstring") snoloop_energy "
";

%feature("docstring") snoexport_fold_arrays "
";

%feature("docstring") snobacktrack_fold_from_pair "
";

%feature("docstring") alisnofold "
";

%feature("docstring") alisnofree_arrays "
";

%feature("docstring") alisnobacktrack_fold_from_pair "
";

// File: snoop_8h.xml

%feature("docstring") snoopfold "

computes snoRNA-RNA interactions in RNAduplex manner  
";

%feature("docstring") snoop_subopt "

computes snoRNA-RNA suboptimal interactions in RNAduplex manner  
";

%feature("docstring") Lsnoop_subopt "

computes snoRNA-RNA suboptimal interactions in a RNAplex manner  
";

%feature("docstring") Lsnoop_subopt_list "

computes snoRNA-RNA suboptimal interactions in a RNAplex manner. The stem energy is saved into a
list of struct, leading to a runtime improvement of 20%  
";

%feature("docstring") Lsnoop_subopt_list_XS "

computes snoRNA-RNA suboptimal interactions in a RNAplex manner. The stem energy is saved into a
list of struct, leading to a runtime improvement of 20%. It considers accessibility  
";

%feature("docstring") snoop_subopt_XS "

computes snoRNA-RNA suboptimal interactions in a RNAduplex manner, and considers accessibility  
";

%feature("docstring") alisnoop_subopt "

aliduplex-like alignment version of snoop_subopt  
";

%feature("docstring") aliLsnoop_subopt_list "

RNAplex-like Alignment version of snoop_subopt  
";

%feature("docstring") alisnoopfold "

RNAaliduplex-like version of snoopfold  
";

%feature("docstring") snoopfold_XS "

RNAduplex-like version of snoopfold with accessibility information  
";

// File: special__const_8h.xml

// File: datastructures_2stream__output_8h.xml

// File: stream__output_8h.xml

// File: string__utils_8h.xml

// File: stringdist_8h.xml

%feature("docstring") Make_swString "

Convert a structure into a format suitable for string_edit_distance().  

Parameters
----------
string : char *

Returns
-------
swString *  
";

%feature("docstring") string_edit_distance "

Calculate the string edit distance of T1 and T2.  

Parameters
----------
T1 : swString *
T2 : swString *

Returns
-------
float  
";

// File: structure__utils_8h.xml

// File: structured__domains_8h.xml

// File: subopt_8h.xml

%feature("docstring") MAXDOS "

Maximum density of states discretization for subopt.  
";

%feature("docstring") UNSORTED "
";

%feature("docstring") SORT_BY_ENERGY_LEXICOGRAPHIC_ASC "
";

%feature("docstring") SORT_BY_ENERGY_ASC "
";

// File: subopt__zuker_8h.xml

// File: svm__utils_8h.xml

// File: treedist_8h.xml

%feature("docstring") make_tree "

Constructs a Tree ( essentially the postorder list ) of the structure 'struc', for use in
tree_edit_distance().  

Parameters
----------
struc : char *
    may be any rooted structure representation.  

Returns
-------
Tree *  
";

%feature("docstring") tree_edit_distance "

Calculates the edit distance of the two trees.  

Parameters
----------
T1 : Tree *
T2 : Tree *

Returns
-------
float  
";

%feature("docstring") print_tree "

Print a tree (mainly for debugging)  
";

%feature("docstring") free_tree "

Free the memory allocated for Tree t.  

Parameters
----------
t : Tree *  
";

// File: ugly__bt_8h.xml

%feature("docstring") sanitize_input "
";

// File: unistd__win_8h.xml

%feature("docstring") srandom "
";

%feature("docstring") random "
";

%feature("docstring") R_OK "
";

%feature("docstring") W_OK "
";

%feature("docstring") F_OK "
";

%feature("docstring") access "
";

%feature("docstring") dup2 "
";

%feature("docstring") execve "
";

%feature("docstring") ftruncate "
";

%feature("docstring") unlink "
";

%feature("docstring") fileno "
";

%feature("docstring") getcwd "
";

%feature("docstring") chdir "
";

%feature("docstring") isatty "
";

%feature("docstring") lseek "
";

%feature("docstring") ssize_t "
";

%feature("docstring") STDIN_FILENO "
";

%feature("docstring") STDOUT_FILENO "
";

%feature("docstring") STDERR_FILENO "
";

// File: units_8h.xml

// File: utils_2units_8h.xml

// File: unstructured__domains_8h.xml

%feature("docstring") vrna_ud_get_motif_size_at "

Get a list of unique motif sizes that start at a certain position within the sequence.  
";

%feature("docstring") vrna_ud_get_motifs_at "
";

%feature("docstring") vrna_ud_detect_motifs "
";

%feature("docstring") vrna_fold_compound_t::ud_set_prob_cb "

**SWIG Wrapper Notes**
    This function is attached as method `ud_set_prob_cb()` to objects of type `fold_compound`. See,
    e.g.  :py:meth:`RNA.fold_compound.ud_set_prob_cb()` in the :doc:`/api_python`.  
";

// File: io_2utils_8h.xml

// File: plotting_2utils_8h.xml

// File: utils_8h.xml

// File: cpu_8h.xml

%feature("docstring") CPU_SIMD_NONE "
";

%feature("docstring") CPU_SIMD_SSE2 "
";

%feature("docstring") CPU_SIMD_SSE3 "
";

%feature("docstring") CPU_SIMD_SSE41 "
";

%feature("docstring") CPU_SIMD_SSE42 "
";

%feature("docstring") CPU_SIMD_AVX "
";

%feature("docstring") CPU_SIMD_AVX2 "
";

%feature("docstring") CPU_SIMD_AVX512F "
";

%feature("docstring") vrna_cpu_vendor_string "
";

%feature("docstring") vrna_cpu_simd_capabilities "
";

// File: higher__order__functions_8h.xml

%feature("docstring") vrna_fun_dispatch_disable "
";

%feature("docstring") vrna_fun_dispatch_enable "
";

%feature("docstring") vrna_fun_zip_add_min "
";

// File: strings_8h.xml

%feature("docstring") str_uppercase "

Convert an input sequence to uppercase.  

.. deprecated:: 2.6.4
    Use RNA.seq_toupper() instead!  
";

%feature("docstring") str_DNA2RNA "

Convert a DNA input sequence to RNA alphabet.  

.. deprecated:: 2.6.4
    Use RNA.seq_toRNA() instead!  
";

%feature("docstring") random_string "

Create a random string using characters from a specified symbol set.  

.. deprecated:: 2.6.4
    Use RNA.random_string() instead!  
";

%feature("docstring") hamming "

Calculate hamming distance between two sequences.  

.. deprecated:: 2.6.4
    Use RNA.hamming_distance() instead!  
";

%feature("docstring") hamming_bound "

Calculate hamming distance between two sequences up to a specified length.  

.. deprecated:: 2.6.4
    Use RNA.hamming_distance_bound() instead!  
";

// File: svm_8h.xml

%feature("docstring") get_z "
";

%feature("docstring") avg_regression "
";

%feature("docstring") sd_regression "
";

%feature("docstring") minimal_sd "
";

%feature("docstring") svm_load_model_string "
";

%feature("docstring") get_seq_composition "
";

// File: vrna__config_8h.xml

%feature("docstring") VERSION "
";

%feature("docstring") VERSION_MAJOR "
";

%feature("docstring") VERSION_MINOR "
";

%feature("docstring") VERSION_PATCH "
";

%feature("docstring") WITH_OPENMP "
";

%feature("docstring") WITH_JSON_SUPPORT "
";

%feature("docstring") WITH_SVM "
";

%feature("docstring") WITH_GSL "
";

%feature("docstring") WITH_NAVIEW_LAYOUT "
";

// File: landscape_2walk_8h.xml

// File: walk_8h.xml

// File: wrap__dlib_8h.xml

%feature("docstring") vrna_equilibrium_conc "
";

// File: zscore_8h.xml

%feature("docstring") ZSCORE_OPTIONS_NONE "
";

%feature("docstring") ZSCORE_FILTER_ON "
";

%feature("docstring") ZSCORE_PRE_FILTER "
";

%feature("docstring") ZSCORE_REPORT_SUBSUMED "
";

%feature("docstring") ZSCORE_MODEL_DEFAULT "
";

%feature("docstring") ZSCORE_SETTINGS_DEFAULT "
";

%feature("docstring") vrna_fold_compound_t::zsc_filter_init "
";

%feature("docstring") vrna_fold_compound_t::zsc_filter_update "
";

%feature("docstring") vrna_fold_compound_t::zsc_filter_free "
";

%feature("docstring") vrna_fold_compound_t::zsc_filter_on "
";

%feature("docstring") vrna_fold_compound_t::zsc_filter_threshold "
";

%feature("docstring") vrna_fold_compound_t::zsc_compute "
";

%feature("docstring") vrna_fold_compound_t::zsc_compute_raw "
";

// File: group__eval.xml

%feature("docstring") vrna_fold_compound_t::eval_structure "

Calculate the free energy of an already folded RNA.  

This function allows for energy evaluation of a given pair of structure and sequence (alignment).
Model details, energy parameters, and possibly soft constraints are used as provided via the
parameter 'fc'. The RNA.fold_compound() does not need to contain any DP matrices, but requires all
most basic init values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as method `eval_structure()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.eval_structure()` in the :doc:`/api_python` .  

Parameters
----------
structure : const char *
    Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_pt(), RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt_verbose(),
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.fold_compound.eval_covar_structure()  

Note
----
Accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE and RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_fold_compound_t::eval_covar_structure "

Calculate the pseudo energy derived by the covariance scores of a set of aligned sequences.  

Consensus structure prediction is driven by covariance scores of base pairs in rows of the provided
alignment. This function allows one to retrieve the total amount of this covariance pseudo energy
scores. The RNA.fold_compound() does not need to contain any DP matrices, but requires all most
basic init values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as method `eval_covar_structure()` to objects of type `fold_compound`.
    See, e.g.   :py:meth:`RNA.fold_compound.eval_covar_structure()` in the :doc:`/api_python` .  

Parameters
----------
structure : const char *
    Secondary (consensus) structure in dot-bracket notation  

Returns
-------
float  
    The covariance pseudo energy score of the input structure given the input sequence alignment in
    kcal/mol  

See Also
--------
RNA.fold_compound_comparative(), RNA.fold_compound.eval_structure()  

Note
----
Accepts RNA.fold_compound() of type RNA.FC_TYPE_COMPARATIVE only!  
";

%feature("docstring") vrna_fold_compound_t::eval_structure_verbose "

Calculate the free energy of an already folded RNA and print contributions on a per-loop base.  

This function is a simplyfied version of RNA.eval_structure_v() that uses the *default* verbosity
level.  

**SWIG Wrapper Notes**
    This function is attached as method `eval_structure_verbose()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.eval_structure_verbose()` in the
    :doc:`/api_python` .  

Parameters
----------
structure : const char *
    Secondary structure in dot-bracket notation  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_pt(), RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt_verbose(),  
";

%feature("docstring") vrna_eval_structure_v "

Calculate the free energy of an already folded RNA and print contributions on a per-loop base.  

This function allows for detailed energy evaluation of a given sequence/structure pair. In contrast
to RNA.fold_compound.eval_structure() this function prints detailed energy contributions based on individual
loops to a file handle. If NULL is passed as file handle, this function defaults to print to stdout.
Any positive `verbosity_level` activates potential warning message of the energy evaluting
functions, while values :math:`\\ge 1` allow for detailed control of what data is printed. A
negative parameter `verbosity_level` turns off printing all together.  

Model details, energy parameters, and possibly soft constraints are used as provided via the
parameter 'fc'. The fold_compound does not need to contain any DP matrices, but all the most basic
init values as one would get from a call like this:  

Parameters
----------
fc : RNA.fold_compound() *
    A RNA.fold_compound() containing the energy parameters and model details  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_pt(), RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt_verbose(),  
";

%feature("docstring") vrna_eval_structure_cstr "
";

%feature("docstring") vrna_fold_compound_t::eval_structure_pt "

Calculate the free energy of an already folded RNA.  

This function allows for energy evaluation of a given sequence/structure pair where the structure is
provided in pair_table format as obtained from RNA.ptable(). Model details, energy parameters, and
possibly soft constraints are used as provided via the parameter 'fc'. The fold_compound does not
need to contain any DP matrices, but all the most basic init values as one would get from a call
like this:  

**SWIG Wrapper Notes**
    This function is attached as method `eval_structure_pt()` to objects of type `fold_compound`.
    See, e.g.   :py:meth:`RNA.fold_compound.eval_structure_pt()` in the :doc:`/api_python` .  

Parameters
----------
pt : const short *
    Secondary structure as pair_table  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.ptable(), RNA.fold_compound.eval_structure(), RNA.fold_compound.eval_structure_pt_verbose()  
";

%feature("docstring") vrna_fold_compound_t::eval_structure_pt_verbose "

Calculate the free energy of an already folded RNA.  

This function is a simplyfied version of RNA.eval_structure_simple_v() that uses the *default*
verbosity level.  

**SWIG Wrapper Notes**
    This function is attached as method `eval_structure_pt_verbose()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.eval_structure_pt_verbose()` in the
    :doc:`/api_python` .  

Parameters
----------
pt : const short *
    Secondary structure as pair_table  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.eval_structure_pt_v(), RNA.ptable(), RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_verbose()  
";

%feature("docstring") vrna_eval_structure_pt_v "

Calculate the free energy of an already folded RNA.  

This function allows for energy evaluation of a given sequence/structure pair where the structure is
provided in pair_table format as obtained from RNA.ptable(). Model details, energy parameters, and
possibly soft constraints are used as provided via the parameter 'fc'. The fold_compound does not
need to contain any DP matrices, but all the most basic init values as one would get from a call
like this:  In contrast to RNA.fold_compound.eval_structure_pt() this function prints detailed energy
contributions based on individual loops to a file handle. If NULL is passed as file handle, this
function defaults to print to stdout. Any positive `verbosity_level` activates potential warning
message of the energy evaluting functions, while values :math:`\\ge 1` allow for detailed control of
what data is printed. A negative parameter `verbosity_level` turns off printing all together.  

Parameters
----------
fc : RNA.fold_compound() *
    A RNA.fold_compound() containing the energy parameters and model details  
pt : const short *
    Secondary structure as pair_table  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.ptable(), RNA.fold_compound.eval_structure_pt(), RNA.fold_compound.eval_structure_verbose()  
";

%feature("docstring") vrna_eval_structure_simple "

Calculate the free energy of an already folded RNA.  

This function allows for energy evaluation of a given sequence/structure pair. In contrast to
RNA.fold_compound.eval_structure() this function assumes default model details and default energy parameters in
order to evaluate the free energy of the secondary structure. Therefore, it serves as a simple
interface function for energy evaluation for situations where no changes on the energy model are
required.  

**SWIG Wrapper Notes**
    In the target scripting language, this function serves as a wrapper for
    RNA.eval_structure_simple_v() and, thus, allows for two additional, optional arguments, the
    verbosity level and a file handle which default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively.. See, e.g.   :py:func:`RNA.eval_structure_simple()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure(), RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt_verbose(),  
";

%feature("docstring") vrna_eval_circ_structure "

Evaluate the free energy of a sequence/structure pair where the sequence is circular.  

**SWIG Wrapper Notes**
    In the target scripting language, this function serves as a wrapper for
    RNA.eval_circ_structure_v() and, thus, allows for two additional, optional arguments, the
    verbosity level and a file handle which default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively.. See, e.g.   :py:func:`RNA.eval_circ_structure()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the structure given the circular input sequence in kcal/mol  

See Also
--------
RNA.eval_structure_simple(), RNA.eval_gquad_structure(), RNA.eval_circ_consensus_structure(),
RNA.eval_circ_structure_v(), RNA.fold_compound.eval_structure()  
";

%feature("docstring") vrna_eval_gquad_structure "

Evaluate the free energy of a sequence/structure pair where the structure may contain
G-Quadruplexes.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    In the target scripting language, this function serves as a wrapper for
    RNA.eval_gquad_structure_v() and, thus, allows for two additional, optional arguments, the
    verbosity level and a file handle which default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively.. See, e.g.   :py:func:`RNA.eval_gquad_structure()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the structure including contributions of G-quadruplexes in kcal/mol  

See Also
--------
RNA.eval_structure_simple(), RNA.eval_circ_structure(), RNA.eval_gquad_consensus_structure(),
RNA.eval_gquad_structure_v(), RNA.fold_compound.eval_structure()  
";

%feature("docstring") vrna_eval_circ_gquad_structure "

Evaluate the free energy of a sequence/structure pair where the sequence is circular and the
structure may contain G-Quadruplexes.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    In the target scripting language, this function serves as a wrapper for
    RNA.eval_circ_gquad_structure_v() and, thus, allows for two additional, optional arguments, the
    verbosity level and a file handle which default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively.. See, e.g.   :py:func:`RNA.eval_circ_gquad_structure()` in the :doc:`/api_python`
    .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the structure including contributions of G-quadruplexes in kcal/mol  

See Also
--------
RNA.eval_structure_simple(), RNA.eval_circ_gquad_consensus_structure(),
RNA.eval_circ_gquad_structure_v(), RNA.fold_compound.eval_structure()  
";

%feature("docstring") vrna_eval_structure_simple_verbose "

Calculate the free energy of an already folded RNA and print contributions per loop.  

This function is a simplyfied version of RNA.eval_structure_simple_v() that uses the *default*
verbosity level.  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.eval_structure_simple_v(), RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_verbose(), RNA.fold_compound.eval_structure_pt_verbose()  
";

%feature("docstring") my_eval_structure_simple "

Calculate the free energy of an already folded RNA and print contributions per loop.  

This function allows for detailed energy evaluation of a given sequence/structure pair. In contrast
to RNA.fold_compound.eval_structure() this function prints detailed energy contributions based on individual
loops to a file handle. If NULL is passed as file handle, this function defaults to print to stdout.
Any positive `verbosity_level` activates potential warning message of the energy evaluting
functions, while values :math:`\\ge 1` allow for detailed control of what data is printed. A
negative parameter `verbosity_level` turns off printing all together.  

In contrast to RNA.fold_compound.eval_structure_verbose() this function assumes default model details and default
energy parameters in order to evaluate the free energy of the secondary structure. Threefore, it
serves as a simple interface function for energy evaluation for situations where no changes on the
energy model are required.  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_structure_simple(). The
    last two arguments for this function are optional and default to RNA.VERBOSITY_QUIET and
    `NULL`, respectively. See, e.g.   :py:func:`RNA.eval_structure_simple()` in the
    :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_verbose(), RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_pt_verbose(),  
";

%feature("docstring") my_eval_circ_structure "

Evaluate free energy of a sequence/structure pair, assume sequence to be circular and print
contributions per loop.  

This function is the same as RNA.eval_structure_simple_v() but assumes the input sequence to be
circularized.  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_circ_structure(). The last
    two arguments for this function are optional and default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively. See, e.g.   :py:func:`RNA.eval_circ_structure()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.eval_structure_simple_v(), RNA.eval_circ_structure(), RNA.fold_compound.eval_structure_verbose()  
";

%feature("docstring") my_eval_gquad_structure "

Evaluate free energy of a sequence/structure pair, allow for G-Quadruplexes in the structure and
print contributions per loop.  

This function is the same as RNA.eval_structure_simple_v() but allows for annotated G-Quadruplexes
in the dot-bracket structure input.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_gquad_structure(). The
    last two arguments for this function are optional and default to RNA.VERBOSITY_QUIET and
    `NULL`, respectively. See, e.g.   :py:func:`RNA.eval_gquad_structure()` in the
    :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.eval_structure_simple_v(), RNA.eval_gquad_structure(),
RNA.fold_compound.eval_structure_verbose()  
";

%feature("docstring") my_eval_circ_gquad_structure "

Evaluate free energy of a sequence/structure pair, assume sequence to be circular, allow for
G-Quadruplexes in the structure, and print contributions per loop.  

This function is the same as RNA.eval_structure_simple_v() but assumes the input sequence to be
circular and allows for annotated G-Quadruplexes in the dot-bracket structure input.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_circ_gquad_structure().
    The last two arguments for this function are optional and default to RNA.VERBOSITY_QUIET and
    `NULL`, respectively. See, e.g.   :py:func:`RNA.eval_circ_gquad_structure()` in the
    :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  
";

%feature("docstring") vrna_eval_consensus_structure_simple "

Calculate the free energy of an already folded RNA sequence alignment.  

This function allows for energy evaluation for a given multiple sequence alignment and consensus
structure pair. In contrast to RNA.fold_compound.eval_structure() this function assumes default model details and
default energy parameters in order to evaluate the free energy of the secondary structure.
Therefore, it serves as a simple interface function for energy evaluation for situations where no
changes on the energy model are required.  

**SWIG Wrapper Notes**
    This function is available through an overloadeded version of RNA.eval_structure_simple().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. See, e.g.   :py:func:`RNA.eval_structure_simple()` in the
    :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters and hyphen ('-') to denote gaps  
structure : const char *
    Consensus Secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the consensus structure given the input alignment in kcal/mol  

See Also
--------
RNA.fold_compound.eval_covar_structure(), RNA.fold_compound.eval_structure(),
RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_verbose(), RNA.fold_compound.eval_structure_pt_verbose()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_circ_consensus_structure "

Evaluate the free energy of a multiple sequence alignment/consensus structure pair where the
sequences are circular.  

**SWIG Wrapper Notes**
    This function is available through an overloadeded version of RNA.eval_circ_structure(). Simply
    pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. See, e.g.   :py:func:`RNA.eval_circ_structure()` in the
    :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters  
structure : const char *
    Consensus secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the consensus structure given the circular input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_covar_structure(), RNA.eval_consensus_structure_simple(),
RNA.eval_gquad_consensus_structure(), RNA.eval_circ_structure(),
RNA.eval_circ_consensus_structure_v(), RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_gquad_consensus_structure "

Evaluate the free energy of a multiple sequence alignment/consensus structure pair where the
structure may contain G-Quadruplexes.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloadeded version of RNA.eval_gquad_structure().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. See, e.g.   :py:func:`RNA.eval_gquad_structure()` in the
    :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters  
structure : const char *
    Consensus secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the consensus structure including contributions of G-quadruplexes in kcal/mol  

See Also
--------
RNA.fold_compound.eval_covar_structure(), RNA.eval_consensus_structure_simple(),
RNA.eval_circ_consensus_structure(), RNA.eval_gquad_structure(),
RNA.eval_gquad_consensus_structure_v(), RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_circ_gquad_consensus_structure "

Evaluate the free energy of a multiple sequence alignment/consensus structure pair where the
sequence is circular and the structure may contain G-Quadruplexes.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloadeded version of RNA.eval_circ_gquad_structure().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. See, e.g.   :py:func:`RNA.eval_circ_gquad_structure()` in the
    :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters  
structure : const char *
    Consensus secondary structure in dot-bracket notation  

Returns
-------
float  
    The free energy of the consensus structure including contributions of G-quadruplexes in kcal/mol  

See Also
--------
RNA.fold_compound.eval_covar_structure(), RNA.eval_consensus_structure_simple(),
RNA.eval_circ_consensus_structure(), RNA.eval_gquad_structure(),
RNA.eval_circ_gquad_consensus_structure_v(), RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_consensus_structure_simple_verbose "

Evaluate the free energy of a consensus structure for an RNA sequence alignment and print
contributions per loop.  

This function is a simplyfied version of RNA.eval_consensus_structure_simple_v() that uses the
*default* verbosity level.  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
structure : const char *
    Consensus secondary structure in dot-bracket notation  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the conensus structure given the aligned input sequences in kcal/mol  

See Also
--------
RNA.eval_consensus_structure_simple_v(), RNA.fold_compound.eval_structure_verbose(),
RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure_pt_verbose()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_consensus_structure_simple_v "

Evaluate the free energy of a consensus structure for an RNA sequence alignment and print
contributions per loop.  

This function allows for detailed energy evaluation of a given sequence alignment/consensus
structure pair. In contrast to RNA.eval_consensus_structure_simple() this function prints detailed
energy contributions based on individual loops to a file handle. If NULL is passed as file handle,
this function defaults to print to stdout. Any positive `verbosity_level` activates potential
warning message of the energy evaluting functions, while values :math:`\\ge 1` allow for detailed
control of what data is printed. A negative parameter `verbosity_level` turns off printing all
together.  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_structure_simple(). Simply
    pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. The last two arguments are optional and default to
    RNA.VERBOSITY_QUIET and `NULL`, respectively. See, e.g.
    :py:func:`RNA.eval_structure_simple()` in the :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
structure : const char *
    Consensus secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the consensus structure given the sequence alignment in kcal/mol  

See Also
--------
RNA.eval_consensus_structure(), RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_circ_consensus_structure_v "

Evaluate the free energy of a consensus structure for an alignment of circular RNA sequences and
print contributions per loop.  

This function is identical with RNA.eval_consensus_structure_simple_v() but assumed the aligned
sequences to be circular.  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_circ_structure(). Simply
    pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. The last two arguments are optional and default to
    RNA.VERBOSITY_QUIET and `NULL`, respectively. See, e.g.   :py:func:`RNA.eval_circ_structure()`
    in the :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
structure : const char *
    Consensus secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the consensus structure given the sequence alignment in kcal/mol  

See Also
--------
RNA.eval_consensus_structure_simple_v(), RNA.eval_circ_consensus_structure(),
RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_gquad_consensus_structure_v "

Evaluate the free energy of a consensus structure for an RNA sequence alignment, allow for annotated
G-Quadruplexes in the structure and print contributions per loop.  

This function is identical with RNA.eval_consensus_structure_simple_v() but allows for annotated
G-Quadruplexes in the consensus structure.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_gquad_structure(). Simply
    pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. The last two arguments are optional and default to
    RNA.VERBOSITY_QUIET and `NULL`, respectively. See, e.g.   :py:func:`RNA.eval_gquad_structure()`
    in the :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
structure : const char *
    Consensus secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the consensus structure given the sequence alignment in kcal/mol  

See Also
--------
RNA.eval_consensus_structure_simple_v(), RNA.eval_gquad_consensus_structure(),
RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_circ_gquad_consensus_structure_v "

Evaluate the free energy of a consensus structure for an alignment of circular RNA sequences, allow
for annotated G-Quadruplexes in the structure and print contributions per loop.  

This function is identical with RNA.eval_consensus_structure_simple_v() but assumes the sequences
in the alignment to be circular and allows for annotated G-Quadruplexes in the consensus structure.  

G-Quadruplexes are annotated as plus signs ('+') for each G involved in the motif. Linker sequences
must be denoted by dots ('.') as they are considered unpaired. Below is an example of a 2-layer
G-quadruplex:  

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_circ_gquad_structure().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. The last two arguments are optional and default to
    RNA.VERBOSITY_QUIET and `NULL`, respectively. See, e.g.
    :py:func:`RNA.eval_circ_gquad_structure()` in the :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
structure : const char *
    Consensus secondary structure in dot-bracket notation  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
float  
    The free energy of the consensus structure given the sequence alignment in kcal/mol  

See Also
--------
RNA.eval_consensus_structure_simple_v(), RNA.eval_circ_gquad_consensus_structure(),
RNA.fold_compound.eval_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_structure_pt_simple "

Calculate the free energy of an already folded RNA.  

In contrast to RNA.fold_compound.eval_structure_pt() this function assumes default model details and default
energy parameters in order to evaluate the free energy of the secondary structure. Threefore, it
serves as a simple interface function for energy evaluation for situations where no changes on the
energy model are required.  

**SWIG Wrapper Notes**
    In the target scripting language, this function serves as a wrapper for
    RNA.eval_structure_pt_v() and, thus, allows for two additional, optional arguments, the
    verbosity level and a file handle which default to RNA.VERBOSITY_QUIET and `NULL`,
    respectively. See, e.g.   :py:func:`RNA.eval_structure_pt_simple()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
pt : const short *
    Secondary structure as pair_table  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.ptable(), RNA.eval_structure_simple(), RNA.fold_compound.eval_structure_pt()  
";

%feature("docstring") vrna_eval_structure_pt_simple_verbose "

Calculate the free energy of an already folded RNA.  

This function is a simplyfied version of RNA.eval_structure_pt_simple_v() that uses the *default*
verbosity level.  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
pt : const short *
    Secondary structure as pair_table  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.eval_structure_pt_simple_v(), RNA.ptable(), RNA.fold_compound.eval_structure_pt_verbose(),
RNA.eval_structure_simple()  
";

%feature("docstring") my_eval_structure_pt_simple "

Calculate the free energy of an already folded RNA.  

This function allows for energy evaluation of a given sequence/structure pair where the structure is
provided in pair_table format as obtained from RNA.ptable(). Model details, energy parameters, and
possibly soft constraints are used as provided via the parameter 'fc'. The fold_compound does not
need to contain any DP matrices, but all the most basic init values as one would get from a call
like this:  In contrast to RNA.fold_compound.eval_structure_pt_verbose() this function assumes default model
details and default energy parameters in order to evaluate the free energy of the secondary
structure. Threefore, it serves as a simple interface function for energy evaluation for situations
where no changes on the energy model are required.  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
pt : const short *
    Secondary structure as pair_table  
verbosity_level : int
    The level of verbosity of this function  
file : FILE *
    A file handle where this function should print to (may be NULL).  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10cal/mol  

See Also
--------
RNA.ptable(), RNA.eval_structure_pt_v(), RNA.eval_structure_simple()  
";

%feature("docstring") vrna_eval_consensus_structure_pt_simple "

Evaluate the Free Energy of a Consensus Secondary Structure given a Sequence Alignment.  

**SWIG Wrapper Notes**
    This function is available through an overloadeded version of RNA.eval_structure_pt_simple().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. See, e.g.   :py:func:`RNA.eval_structure_pt_simple()` in the
    :doc:`/api_python` .  

Parameters
----------
alignment : const char **
    RNA sequence alignment in uppercase letters. Gaps are denoted by hyphens ('-')  
pt : const short *
    Secondary structure in pair table format  

Returns
-------
int  
    Free energy of the consensus structure in 10cal/mol  

See Also
--------
RNA.eval_consensus_structure_simple(), RNA.fold_compound.eval_structure_pt(),
RNA.fold_compound.eval_structure(),
RNA.fold_compound.eval_covar_structure()  

Note
----
The free energy returned from this function already includes the covariation pseudo energies that is
used fir comparative structure prediction within this library.  
";

%feature("docstring") vrna_eval_consensus_structure_pt_simple_verbose "
";

%feature("docstring") vrna_eval_consensus_structure_pt_simple_v "

**SWIG Wrapper Notes**
    This function is available through an overloaded version of RNA.eval_structure_pt_simple().
    Simply pass a sequence alignment as list of strings (including gaps) as first, and the consensus
    structure as second argument. The last two arguments are optional and default to
    RNA.VERBOSITY_QUIET and `NULL`, respectively. See, e.g.
    :py:func:`RNA.eval_structure_pt_simple()` in the :doc:`/api_python` .  
";

%feature("docstring") VERBOSITY_QUIET "

Quiet level verbosity setting.  
";

%feature("docstring") VERBOSITY_DEFAULT "

Default level verbosity setting.  
";

// File: group__eval__loops.xml

%feature("docstring") vrna_fold_compound_t::eval_loop_pt "

Calculate energy of a loop.  

**SWIG Wrapper Notes**
    This function is attached as method `eval_loop_pt()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.eval_loop_pt()` in the :doc:`/api_python` .  

Parameters
----------
i : int
    position of covering base pair  
pt : const short *
    the pair table of the secondary structure  

Returns
-------
int  
    free energy of the loop in 10cal/mol  
";

%feature("docstring") vrna_eval_loop_pt_v "

Calculate energy of a loop.  

Parameters
----------
fc : RNA.fold_compound() *
    A RNA.fold_compound() containing the energy parameters and model details  
i : int
    position of covering base pair  
pt : const short *
    the pair table of the secondary structure  
verbosity_level : int
    The level of verbosity of this function  

Returns
-------
int  
    free energy of the loop in 10cal/mol  
";

// File: group__eval__move.xml

%feature("docstring") vrna_fold_compound_t::eval_move "

Calculate energy of a move (closing or opening of a base pair)  

If the parameters m1 and m2 are negative, it is deletion (opening) of a base pair, otherwise it is
insertion (opening).  

**SWIG Wrapper Notes**
    This function is attached as method `eval_move()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.eval_move()` in the :doc:`/api_python` .  

Parameters
----------
structure : const char *
    secondary structure in dot-bracket notation  
m1 : int
    first coordinate of base pair  
m2 : int
    second coordinate of base pair  

Returns
-------
float  
    energy change of the move in kcal/mol (INF / 100. upon any error)  

See Also
--------
RNA.fold_compound.eval_move_pt()  
";

%feature("docstring") vrna_fold_compound_t::eval_move_pt "

Calculate energy of a move (closing or opening of a base pair)  

If the parameters m1 and m2 are negative, it is deletion (opening) of a base pair, otherwise it is
insertion (opening).  

**SWIG Wrapper Notes**
    This function is attached as method `eval_move_pt()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.eval_move_pt()` in the :doc:`/api_python` .  

Parameters
----------
pt : short *
    the pair table of the secondary structure  
m1 : int
    first coordinate of base pair  
m2 : int
    second coordinate of base pair  

Returns
-------
int  
    energy change of the move in 10cal/mol  

See Also
--------
RNA.fold_compound.eval_move()  
";

%feature("docstring") vrna_eval_move_pt_simple "
";

%feature("docstring") vrna_eval_move_shift_pt "
";

// File: group__eval__deprecated.xml

%feature("docstring") energy_of_structure "

Calculate the free energy of an already folded RNA using global model detail settings.  

If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure() or RNA.fold_compound.eval_structure_verbose() instead!  

Note
----
OpenMP: This function relies on several global model settings variables and thus is not to be
considered threadsafe. See energy_of_struct_par() for a completely threadsafe implementation.  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    secondary structure in dot-bracket notation  
verbosity_level : int
    a flag to turn verbose output on/off  

Returns
-------
float  
    the free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure()  
";

%feature("docstring") energy_of_struct_par "

Calculate the free energy of an already folded RNA.  

If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure() or RNA.fold_compound.eval_structure_verbose() instead!  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
structure : const char *
    Secondary structure in dot-bracket notation  
parameters : RNA.param() *
    A data structure containing the prescaled energy contributions and the model details.  
verbosity_level : int
    A flag to turn verbose output on/off  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure()  
";

%feature("docstring") energy_of_circ_structure "

Calculate the free energy of an already folded circular RNA.  


If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

Note
----
OpenMP: This function relies on several global model settings variables and thus is not to be
considered threadsafe. See energy_of_circ_struct_par() for a completely threadsafe implementation.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure() or RNA.fold_compound.eval_structure_verbose() instead!  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    Secondary structure in dot-bracket notation  
verbosity_level : int
    A flag to turn verbose output on/off  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure()  
";

%feature("docstring") energy_of_circ_struct_par "

Calculate the free energy of an already folded circular RNA.  

If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure() or RNA.fold_compound.eval_structure_verbose() instead!  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    Secondary structure in dot-bracket notation  
parameters : RNA.param() *
    A data structure containing the prescaled energy contributions and the model details.  
verbosity_level : int
    A flag to turn verbose output on/off  

Returns
-------
float  
    The free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure()  
";

%feature("docstring") energy_of_gquad_structure "
";

%feature("docstring") energy_of_gquad_struct_par "
";

%feature("docstring") energy_of_structure_pt "

Calculate the free energy of an already folded RNA.  

If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure_pt() or RNA.fold_compound.eval_structure_pt_verbose()
instead!  

Note
----
OpenMP: This function relies on several global model settings variables and thus is not to be
considered threadsafe. See energy_of_struct_pt_par() for a completely threadsafe implementation.  

Parameters
----------
string : const char *
    RNA sequence  
ptable : short *
    the pair table of the secondary structure  
s : short *
    encoded RNA sequence  
s1 : short *
    encoded RNA sequence  
verbosity_level : int
    a flag to turn verbose output on/off  

Returns
-------
int  
    the free energy of the input structure given the input sequence in 10kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_pt()  
";

%feature("docstring") energy_of_struct_pt_par "

Calculate the free energy of an already folded RNA.  

If verbosity level is set to a value >0, energies of structure elements are printed to stdout  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_structure_pt() or RNA.fold_compound.eval_structure_pt_verbose()
instead!  

Parameters
----------
string : const char *
    RNA sequence in uppercase letters  
ptable : short *
    The pair table of the secondary structure  
s : short *
    Encoded RNA sequence  
s1 : short *
    Encoded RNA sequence  
parameters : RNA.param() *
    A data structure containing the prescaled energy contributions and the model details.  
verbosity_level : int
    A flag to turn verbose output on/off  

Returns
-------
int  
    The free energy of the input structure given the input sequence in 10kcal/mol  

See Also
--------
RNA.fold_compound.eval_structure_pt()  
";

%feature("docstring") energy_of_move "

Calculate energy of a move (closing or opening of a base pair)  

If the parameters m1 and m2 are negative, it is deletion (opening) of a base pair, otherwise it is
insertion (opening).  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_move() instead!  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    secondary structure in dot-bracket notation  
m1 : int
    first coordinate of base pair  
m2 : int
    second coordinate of base pair  

Returns
-------
float  
    energy change of the move in kcal/mol  

See Also
--------
RNA.fold_compound.eval_move()  
";

%feature("docstring") energy_of_move_pt "

Calculate energy of a move (closing or opening of a base pair)  

If the parameters m1 and m2 are negative, it is deletion (opening) of a base pair, otherwise it is
insertion (opening).  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_move_pt() instead!  

Parameters
----------
pt : short *
    the pair table of the secondary structure  
s : short *
    encoded RNA sequence  
s1 : short *
    encoded RNA sequence  
m1 : int
    first coordinate of base pair  
m2 : int
    second coordinate of base pair  

Returns
-------
int  
    energy change of the move in 10cal/mol  

See Also
--------
RNA.fold_compound.eval_move_pt()  
";

%feature("docstring") loop_energy "

Calculate energy of a loop.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.eval_loop_pt() instead!  

Parameters
----------
ptable : short *
    the pair table of the secondary structure  
s : short *
    encoded RNA sequence  
s1 : short *
    encoded RNA sequence  
i : int
    position of covering base pair  

Returns
-------
int  
    free energy of the loop in 10cal/mol  

See Also
--------
RNA.fold_compound.eval_loop_pt()  
";

%feature("docstring") energy_of_struct "

Calculate the free energy of an already folded RNA  

.. deprecated:: 2.6.4
    This function is deprecated and should not be used in future programs! Use energy_of_structure()
    instead!  

Note
----
This function is not entirely threadsafe! Depending on the state of the global variable eos_debug it
prints energy information to stdout or not...  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    secondary structure in dot-bracket notation  

Returns
-------
float  
    the free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
energy_of_structure, energy_of_circ_struct(), energy_of_struct_pt()  
";

%feature("docstring") energy_of_struct_pt "

Calculate the free energy of an already folded RNA  

.. deprecated:: 2.6.4
    This function is deprecated and should not be used in future programs! Use
    energy_of_structure_pt() instead!  

Note
----
This function is not entirely threadsafe! Depending on the state of the global variable eos_debug it
prints energy information to stdout or not...  

Parameters
----------
string : const char *
    RNA sequence  
ptable : short *
    the pair table of the secondary structure  
s : short *
    encoded RNA sequence  
s1 : short *
    encoded RNA sequence  

Returns
-------
int  
    the free energy of the input structure given the input sequence in 10kcal/mol  

See Also
--------
make_pair_table(), energy_of_structure()  
";

%feature("docstring") energy_of_circ_struct "

Calculate the free energy of an already folded circular RNA  

.. deprecated:: 2.6.4
    This function is deprecated and should not be used in future programs Use
    energy_of_circ_structure() instead!  

Note
----
This function is not entirely threadsafe! Depending on the state of the global variable eos_debug it
prints energy information to stdout or not...  

Parameters
----------
string : const char *
    RNA sequence  
structure : const char *
    secondary structure in dot-bracket notation  

Returns
-------
float  
    the free energy of the input structure given the input sequence in kcal/mol  

See Also
--------
energy_of_circ_structure(), energy_of_struct(), energy_of_struct_pt()  
";

%feature("docstring") E_Stem "

Compute the energy contribution of a stem branching off a loop-region.  

This function computes the energy contribution of a stem that branches off a loop region. This can
be the case in multiloops, when a stem branching off increases the degree of the loop but also
*immediately interior base pairs* of an exterior loop contribute free energy. To switch the behavior
of the function according to the evaluation of a multiloop- or exterior-loop-stem, you pass the flag
'extLoop'. The returned energy contribution consists of a TerminalAU penalty if the pair type is
greater than 2, dangling end contributions of mismatching nucleotides adjacent to the stem if only
one of the si1, sj1 parameters is greater than 0 and mismatch energies if both mismatching
nucleotides are positive values. Thus, to avoid incorporating dangling end or mismatch energies just
pass a negative number, e.g. -1 to the mismatch argument.  

This is an illustration of how the energy contribution is assembled:
      3'  5'
      |   |
      X - Y
5'-si1     sj1-3'  

Here, (X,Y) is the base pair that closes the stem that branches off a loop region. The nucleotides
si1 and sj1 are the 5'- and 3'- mismatches, respectively. If the base pair type of (X,Y) is greater
than 2 (i.e. an A-U or G-U pair, the TerminalAU penalty will be included in the energy contribution
returned. If si1 and sj1 are both nonnegative numbers, mismatch energies will also be included. If
one of si1 or sj1 is a negative value, only 5' or 3' dangling end contributions are taken into
account. To prohibit any of these mismatch contributions to be incorporated, just pass a negative
number to both, si1 and sj1. In case the argument extLoop is 0, the returned energy contribution
also includes the *internal-loop-penalty* of a multiloop stem with closing pair type.  

.. deprecated:: 2.6.4
    Please use one of the functions RNA.E_ext_stem() and E_MLstem() instead! Use the former for
    cases where `extLoop` != 0 and the latter otherwise.  

See Also
--------
E_MLstem(), _ExtLoop()  

Note
----
This function is threadsafe  

Parameters
----------
type : int
    The pair type of the first base pair un the stem  
si1 : int
    The 5'-mismatching nucleotide  
sj1 : int
    The 3'-mismatching nucleotide  
extLoop : int
    A flag that indicates whether the contribution reflects the one of an exterior loop or not  
P : RNA.param() *
    The data structure containing scaled energy parameters  

Returns
-------
int  
    The Free energy of the branch off the loop in dcal/mol  
";

%feature("docstring") E_ExtLoop "
";

%feature("docstring") exp_E_ExtLoop "

This is the partition function variant of E_ExtLoop()  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.exp_E_ext_stem() instead!  

Returns
-------
FLT_OR_DBL  
    The Boltzmann weighted energy contribution of the introduced exterior-loop stem  

See Also
--------
E_ExtLoop()  
";

%feature("docstring") exp_E_Stem "

Compute the Boltzmann weighted energy contribution of a stem branching off a loop-region
----------------------------------------------------------------------------------------  
This is the partition function variant of E_Stem()  

Returns
-------
FLT_OR_DBL  
    The Boltzmann weighted energy contribution of the branch off the loop  

See Also
--------
E_Stem()  

Note
----
This function is threadsafe  
";

%feature("docstring") E_IntLoop "

Compute the Energy of an interior-loop
--------------------------------------  
This function computes the free energy :math:`\\Delta G` of an interior-loop with the following
structure:  
      3'  5'
      |   |
      U - V
  a_n       b_1
   .        .
   .        .
   .        .
  a_1       b_m
      X - Y
      |   |
      5'  3'
 This general structure depicts an interior-loop that is closed by the base pair (X,Y). The enclosed
base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n that constitute the loop. In
this example, the length of the interior-loop is :math:`(n+m)` where n or m may be 0 resulting in a
bulge-loop or base pair stack. The mismatching nucleotides for the closing pair (X,Y) are:  
 5'-mismatch: a_1  
 3'-mismatch: b_m  
 and for the enclosed base pair (V,U):  
 5'-mismatch: b_1  
 3'-mismatch: a_n  

Parameters
----------
n1 : int
    The size of the 'left'-loop (number of unpaired nucleotides)  
n2 : int
    The size of the 'right'-loop (number of unpaired nucleotides)  
type : int
    The pair type of the base pair closing the interior loop  
type_2 : int
    The pair type of the enclosed base pair  
si1 : int
    The 5'-mismatching nucleotide of the closing pair  
sj1 : int
    The 3'-mismatching nucleotide of the closing pair  
sp1 : int
    The 3'-mismatching nucleotide of the enclosed pair  
sq1 : int
    The 5'-mismatching nucleotide of the enclosed pair  
P : RNA.param() *
    The datastructure containing scaled energy parameters  

Returns
-------
int  
    The Free energy of the Interior-loop in dcal/mol  

See Also
--------
scale_parameters(), RNA.param()  

Note
----
Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair must be 'turned
arround' when evaluating the free energy of the interior-loop  
 This function is threadsafe  
";

%feature("docstring") exp_E_IntLoop "

Compute Boltzmann weight :math:`e^{-\\Delta G/kT}` of interior loop
-------------------------------------------------------------------  
multiply by scale[u1+u2+2] for scaling  

Parameters
----------
u1 : int
    The size of the 'left'-loop (number of unpaired nucleotides)  
u2 : int
    The size of the 'right'-loop (number of unpaired nucleotides)  
type : int
    The pair type of the base pair closing the interior loop  
type2 : int
    The pair type of the enclosed base pair  
si1 : short
    The 5'-mismatching nucleotide of the closing pair  
sj1 : short
    The 3'-mismatching nucleotide of the closing pair  
sp1 : short
    The 3'-mismatching nucleotide of the enclosed pair  
sq1 : short
    The 5'-mismatching nucleotide of the enclosed pair  
P : RNA.exp_param() *
    The datastructure containing scaled Boltzmann weights of the energy parameters  

Returns
-------
FLT_OR_DBL  
    The Boltzmann weight of the Interior-loop  

See Also
--------
get_scaled_pf_parameters(), RNA.exp_param(), E_IntLoop()  

Note
----
This function is threadsafe  
";

%feature("docstring") E_IntLoop_Co "
";

%feature("docstring") ubf_eval_int_loop "
";

%feature("docstring") ubf_eval_int_loop2 "
";

%feature("docstring") ubf_eval_ext_int_loop "
";

%feature("docstring") E_MLstem "
";

%feature("docstring") exp_E_MLstem "
";

%feature("docstring") ON_SAME_STRAND "
";

// File: group__grammar.xml

%feature("docstring") vrna_gr_set_aux_f "
";

%feature("docstring") vrna_gr_set_aux_exp_f "
";

%feature("docstring") vrna_gr_set_aux_c "
";

%feature("docstring") vrna_gr_set_aux_exp_c "
";

%feature("docstring") vrna_gr_set_aux_m "
";

%feature("docstring") vrna_gr_set_aux_exp_m "
";

%feature("docstring") vrna_gr_set_aux_m1 "
";

%feature("docstring") vrna_gr_set_aux_exp_m1 "
";

%feature("docstring") vrna_gr_set_aux "
";

%feature("docstring") vrna_gr_set_aux_exp "
";

%feature("docstring") vrna_gr_set_data "
";

%feature("docstring") vrna_gr_set_cond "
";

%feature("docstring") vrna_gr_reset "
";

// File: group__model__details.xml

%feature("docstring") vrna_md_t::reset "

Apply default model details to a provided RNA.md() data structure.  

Use this function to initialize a RNA.md() data structure with its default values  

Parameters
----------  
";

%feature("docstring") vrna_md_update "

Update the model details data structure.  

This function should be called after changing the RNA.md().energy_set attribute since it re-
initializes base pairing related arrays within the RNA.md() data structure. In particular,
RNA.md().pair, RNA.md().alias, and RNA.md().rtype are set to the values that correspond to the
specified RNA.md().energy_set option  

See Also
--------
RNA.md(), RNA.md().energy_set, RNA.md().pair, RNA.md().rtype, RNA.md().alias,
RNA.md.reset()  
";

%feature("docstring") vrna_md_copy "

Copy/Clone a RNA.md() model.  

Use this function to clone a given model either inplace (target container `md_to` given) or create a
copy by cloning the source model and returning it (`md_to` == NULL).  

Parameters
----------
md_to : RNA.md() *
    The model to be overwritten (if non-NULL and `md_to` != `md_from`)  
md_from : const RNA.md() *
    The model to copy (if non-NULL)  

Returns
-------
RNA.md() *  
    A pointer to the copy model (or NULL if `md_from` == NULL)  
";

%feature("docstring") vrna_md_t::option_string "

Get a corresponding commandline parameter string of the options in a RNA.md().  

Note
----
This function is not threadsafe!  
";

%feature("docstring") vrna_md_set_nonstandards "
";

%feature("docstring") vrna_md_defaults_reset "

Reset the global default model details to a specific set of parameters, or their initial values.  

This function resets the global default model details to their initial values, i.e. as specified by
the ViennaRNA Package release, upon passing NULL as argument. Alternatively it resets them according
to a set of provided parameters.  

Parameters
----------
md_p : RNA.md() *
    A set of model details to use as global default (if NULL is passed, factory defaults are
    restored)  

Warnings
--------
This function first resets the global default settings to factory defaults, and only then applies
user provided settings (if any). User settings that do not meet specifications are skipped.  

See Also
--------
RNA.md.reset(), RNA.md()  

Note
----
The global default parameters affect all function calls of RNAlib where model details are not
explicitly provided. Hence, any change of them is not considered threadsafe  
";

%feature("docstring") vrna_md_defaults_temperature "

Set default temperature for energy evaluation of loops.  

Parameters
----------
T : double
    Temperature in centigrade  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_TEMPERATURE  
";

%feature("docstring") vrna_md_defaults_temperature_get "

Get default temperature for energy evaluation of loops.  

Returns
-------
double  
    The global default settings for temperature in centigrade  

See Also
--------
RNA.md_defaults_temperature(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_TEMPERATURE  
";

%feature("docstring") vrna_md_defaults_betaScale "

Set default scaling factor of thermodynamic temperature in Boltzmann factors.  

Bolzmann factors are then computed as :math:`exp(-E / (b \\cdot kT))`.  

Parameters
----------
b : double
    The scaling factor, default is 1.0  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_BETA_SCALE  
";

%feature("docstring") vrna_md_defaults_betaScale_get "

Get default scaling factor of thermodynamic temperature in Boltzmann factors.  

Returns
-------
double  
    The global default thermodynamic temperature scaling factor  

See Also
--------
RNA.md_defaults_betaScale(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_BETA_SCALE  
";

%feature("docstring") vrna_md_defaults_pf_smooth "
";

%feature("docstring") vrna_md_defaults_pf_smooth_get "
";

%feature("docstring") vrna_md_defaults_dangles "

Set default dangle model for structure prediction.  

Parameters
----------
d : int
    The dangle model  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_DANGLES  
";

%feature("docstring") vrna_md_defaults_dangles_get "

Get default dangle model for structure prediction.  

Returns
-------
int  
    The global default settings for the dangle model  

See Also
--------
RNA.md_defaults_dangles(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_DANGLES  
";

%feature("docstring") vrna_md_defaults_special_hp "

Set default behavior for lookup of tabulated free energies for special hairpin loops, such as Tri-,
Tetra-, or Hexa-loops.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_SPECIAL_HP  
";

%feature("docstring") vrna_md_defaults_special_hp_get "

Get default behavior for lookup of tabulated free energies for special hairpin loops, such as Tri-,
Tetra-, or Hexa-loops.  

Returns
-------
int  
    The global default settings for the treatment of special hairpin loops  

See Also
--------
RNA.md_defaults_special_hp(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_SPECIAL_HP  
";

%feature("docstring") vrna_md_defaults_noLP "

Set default behavior for prediction of canonical secondary structures.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_NO_LP  
";

%feature("docstring") vrna_md_defaults_noLP_get "

Get default behavior for prediction of canonical secondary structures.  

Returns
-------
int  
    The global default settings for predicting canonical secondary structures  

See Also
--------
RNA.md_defaults_noLP(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_NO_LP  
";

%feature("docstring") vrna_md_defaults_noGU "

Set default behavior for treatment of G-U wobble pairs.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_NO_GU  
";

%feature("docstring") vrna_md_defaults_noGU_get "

Get default behavior for treatment of G-U wobble pairs.  

Returns
-------
int  
    The global default settings for treatment of G-U wobble pairs  

See Also
--------
RNA.md_defaults_noGU(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_NO_GU  
";

%feature("docstring") vrna_md_defaults_noGUclosure "

Set default behavior for G-U pairs as closing pair for loops.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_NO_GU_CLOSURE  
";

%feature("docstring") vrna_md_defaults_noGUclosure_get "

Get default behavior for G-U pairs as closing pair for loops.  

Returns
-------
int  
    The global default settings for treatment of G-U pairs closing a loop  

See Also
--------
RNA.md_defaults_noGUclosure(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_NO_GU_CLOSURE  
";

%feature("docstring") vrna_md_defaults_logML "

Set default behavior recomputing free energies of multi-branch loops using a logarithmic model.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_LOG_ML  
";

%feature("docstring") vrna_md_defaults_logML_get "

Get default behavior recomputing free energies of multi-branch loops using a logarithmic model.  

Returns
-------
int  
    The global default settings for logarithmic model in multi-branch loop free energy evaluation  

See Also
--------
RNA.md_defaults_logML(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_LOG_ML  
";

%feature("docstring") vrna_md_defaults_circ "

Set default behavior whether input sequences are circularized.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_CIRC  
";

%feature("docstring") vrna_md_defaults_circ_get "

Get default behavior whether input sequences are circularized.  

Returns
-------
int  
    The global default settings for treating input sequences as circular  

See Also
--------
RNA.md_defaults_circ(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_CIRC  
";

%feature("docstring") vrna_md_defaults_gquad "

Set default behavior for treatment of G-Quadruplexes.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_GQUAD  
";

%feature("docstring") vrna_md_defaults_gquad_get "

Get default behavior for treatment of G-Quadruplexes.  

Returns
-------
int  
    The global default settings for treatment of G-Quadruplexes  

See Also
--------
RNA.md_defaults_gquad(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_GQUAD  
";

%feature("docstring") vrna_md_defaults_uniq_ML "

Set default behavior for creating additional matrix for unique multi-branch loop prediction.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_UNIQ_ML  

Note
----
Activating this option usually results in higher memory consumption!  
";

%feature("docstring") vrna_md_defaults_uniq_ML_get "

Get default behavior for creating additional matrix for unique multi-branch loop prediction.  

Returns
-------
int  
    The global default settings for creating additional matrices for unique multi-branch loop
    prediction  

See Also
--------
RNA.md_defaults_uniq_ML(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_UNIQ_ML  
";

%feature("docstring") vrna_md_defaults_energy_set "

Set default energy set.  

Parameters
----------
e : int
    Energy set (0, 1, 2, 3)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_ENERGY_SET  
";

%feature("docstring") vrna_md_defaults_energy_set_get "

Get default energy set.  

Returns
-------
int  
    The global default settings for the energy set  

See Also
--------
RNA.md_defaults_energy_set(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_ENERGY_SET  
";

%feature("docstring") vrna_md_defaults_backtrack "

Set default behavior for whether to backtrack secondary structures.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_BACKTRACK  
";

%feature("docstring") vrna_md_defaults_backtrack_get "

Get default behavior for whether to backtrack secondary structures.  

Returns
-------
int  
    The global default settings for backtracking structures  

See Also
--------
RNA.md_defaults_backtrack(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_BACKTRACK  
";

%feature("docstring") vrna_md_defaults_backtrack_type "

Set default backtrack type, i.e. which DP matrix is used.  

Parameters
----------
t : char
    The type ('F', 'C', or 'M')  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_BACKTRACK_TYPE  
";

%feature("docstring") vrna_md_defaults_backtrack_type_get "

Get default backtrack type, i.e. which DP matrix is used.  

Returns
-------
char  
    The global default settings that specify which DP matrix is used for backtracking  

See Also
--------
RNA.md_defaults_backtrack_type(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_BACKTRACK_TYPE  
";

%feature("docstring") vrna_md_defaults_compute_bpp "

Set the default behavior for whether to compute base pair probabilities after partition function
computation.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_COMPUTE_BPP  
";

%feature("docstring") vrna_md_defaults_compute_bpp_get "

Get the default behavior for whether to compute base pair probabilities after partition function
computation.  

Returns
-------
int  
    The global default settings that specify whether base pair probabilities are computed together
    with partition function  

See Also
--------
RNA.md_defaults_compute_bpp(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_COMPUTE_BPP  
";

%feature("docstring") vrna_md_defaults_max_bp_span "

Set default maximal base pair span.  

Parameters
----------
span : int
    Maximal base pair span  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_MAX_BP_SPAN  
";

%feature("docstring") vrna_md_defaults_max_bp_span_get "

Get default maximal base pair span.  

Returns
-------
int  
    The global default settings for maximum base pair span  

See Also
--------
RNA.md_defaults_max_bp_span(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_MAX_BP_SPAN  
";

%feature("docstring") vrna_md_defaults_min_loop_size "

Set default minimal loop size.  

Parameters
----------
size : int
    Minimal size, i.e. number of unpaired nucleotides for a hairpin loop  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), TURN  
";

%feature("docstring") vrna_md_defaults_min_loop_size_get "

Get default minimal loop size.  

Returns
-------
int  
    The global default settings for minimal size of hairpin loops  

See Also
--------
RNA.md_defaults_min_loop_size(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), TURN  
";

%feature("docstring") vrna_md_defaults_window_size "

Set default window size for sliding window structure prediction approaches.  

Parameters
----------
size : int
    The size of the sliding window  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_WINDOW_SIZE  
";

%feature("docstring") vrna_md_defaults_window_size_get "

Get default window size for sliding window structure prediction approaches.  

Returns
-------
int  
    The global default settings for the size of the sliding window  

See Also
--------
RNA.md_defaults_window_size(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_WINDOW_SIZE  
";

%feature("docstring") vrna_md_defaults_oldAliEn "

Set default behavior for whether to use old energy model for comparative structure prediction.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_ALI_OLD_EN  

Note
----
This option is outdated. Activating the old energy model usually results in worse consensus
structure predictions.  
";

%feature("docstring") vrna_md_defaults_oldAliEn_get "

Get default behavior for whether to use old energy model for comparative structure prediction.  

Returns
-------
int  
    The global default settings for using old energy model for comparative structure prediction  

See Also
--------
RNA.md_defaults_oldAliEn(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_ALI_OLD_EN  
";

%feature("docstring") vrna_md_defaults_ribo "

Set default behavior for whether to use Ribosum Scoring in comparative structure prediction.  

Parameters
----------
flag : int
    On/Off switch (0 = OFF, else = ON)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_ALI_RIBO  
";

%feature("docstring") vrna_md_defaults_ribo_get "

Get default behavior for whether to use Ribosum Scoring in comparative structure prediction.  

Returns
-------
int  
    The global default settings for using Ribosum scoring in comparative structure prediction  

See Also
--------
RNA.md_defaults_ribo(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_ALI_RIBO  
";

%feature("docstring") vrna_md_defaults_cv_fact "

Set the default co-variance scaling factor used in comparative structure prediction.  

Parameters
----------
factor : double
    The co-variance factor  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_ALI_CV_FACT  
";

%feature("docstring") vrna_md_defaults_cv_fact_get "

Get the default co-variance scaling factor used in comparative structure prediction.  

Returns
-------
double  
    The global default settings for the co-variance factor  

See Also
--------
RNA.md_defaults_cv_fact(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_ALI_CV_FACT  
";

%feature("docstring") vrna_md_defaults_nc_fact "

Parameters
----------
factor : double

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(), RNA.MODEL_DEFAULT_ALI_NC_FACT  
";

%feature("docstring") vrna_md_defaults_nc_fact_get "

Returns
-------
double  

See Also
--------
RNA.md_defaults_nc_fact(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md(),
RNA.MODEL_DEFAULT_ALI_NC_FACT  
";

%feature("docstring") vrna_md_defaults_sfact "

Set the default scaling factor used to avoid under-/overflows in partition function computation.  

Parameters
----------
factor : double
    The scaling factor (default: 1.07)  

See Also
--------
RNA.md_defaults_reset(), RNA.md.reset(), RNA.md()  
";

%feature("docstring") vrna_md_defaults_sfact_get "

Get the default scaling factor used to avoid under-/overflows in partition function computation.  

Returns
-------
double  
    The global default settings of the scaling factor  

See Also
--------
RNA.md_defaults_sfact(), RNA.md_defaults_reset(), RNA.md.reset(), RNA.md()  
";

%feature("docstring") vrna_md_defaults_salt "

Set the default salt concentration.  

Parameters
----------
salt : double
    The sodium concentration in M (default: 1.021)  
";

%feature("docstring") vrna_md_defaults_salt_get "

Get the default salt concentration.  

Returns
-------
double  
    The default salt concentration  
";

%feature("docstring") vrna_md_defaults_saltMLLower "

Set the default multiloop size lower bound for loop salt correciton linear fitting.  

Parameters
----------
lower : int
    Size lower bound (number of backbone in loop)  
";

%feature("docstring") vrna_md_defaults_saltMLLower_get "

Get the default multiloop size lower bound for loop salt correciton linear fitting.  

Returns
-------
int  
    The default lower bound  
";

%feature("docstring") vrna_md_defaults_saltMLUpper "

Set the default multiloop size upper bound for loop salt correciton linear fitting.  

Parameters
----------
upper : int
    Size Upper bound (number of backbone in loop)  
";

%feature("docstring") vrna_md_defaults_saltMLUpper_get "

Get the default multiloop size upper bound for loop salt correciton linear fitting.  

Returns
-------
int  
    The default upper bound  
";

%feature("docstring") vrna_md_defaults_saltDPXInit "

Set user-provided salt correciton for duplex initialization If value is 99999 the default value from
fitting is used.  

Parameters
----------
value : int
    The value of salt correction for duplex initialization (in dcal/mol)  
";

%feature("docstring") vrna_md_defaults_saltDPXInit_get "

Get user-provided salt correciton for duplex initialization If value is 99999 the default value from
fitting is used.  

Returns
-------
int  
    The user-provided salt correction for duplex initialization  
";

%feature("docstring") vrna_md_defaults_saltDPXInitFact "
";

%feature("docstring") vrna_md_defaults_saltDPXInitFact_get "
";

%feature("docstring") vrna_md_defaults_helical_rise "
";

%feature("docstring") vrna_md_defaults_helical_rise_get "
";

%feature("docstring") vrna_md_defaults_backbone_length "
";

%feature("docstring") vrna_md_defaults_backbone_length_get "
";

%feature("docstring") set_model_details "

Set default model details.  

Use this function if you wish to initialize a RNA.md() data structure with its default values, i.e.
the global model settings as provided by the deprecated global variables.  

.. deprecated:: 2.6.4
    This function will vanish as soon as backward compatibility of RNAlib is dropped (expected in
    version 3). Use RNA.md.reset() instead!  

Parameters
----------
md : RNA.md() *
    A pointer to the data structure that is about to be initialized  
";

%feature("docstring") option_string "
";

%feature("docstring") NBASES "
";

%feature("docstring") MODEL_DEFAULT_TEMPERATURE "

 Default temperature for structure prediction and free energy evaluation in &#176C  Default
temperature for structure prediction and free energy evaluation in $^\\circ C$  

See Also
--------
RNA.md().temperature, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_PF_SCALE "

Default scaling factor for partition function computations.  

See Also
--------
RNA.exp_param().pf_scale, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_BETA_SCALE "

Default scaling factor for absolute thermodynamic temperature in Boltzmann factors.  

See Also
--------
RNA.exp_param().alpha, RNA.md().betaScale, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_DANGLES "

Default dangling end model.  

See Also
--------
RNA.md().dangles, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_SPECIAL_HP "

Default model behavior for lookup of special tri-, tetra-, and hexa-loops.  

See Also
--------
RNA.md().special_hp, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_NO_LP "

Default model behavior for so-called 'lonely pairs'.  

See Also
--------
RNA.md().noLP, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_NO_GU "

Default model behavior for G-U base pairs.  

See Also
--------
RNA.md().noGU, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_NO_GU_CLOSURE "

Default model behavior for G-U base pairs closing a loop.  

See Also
--------
RNA.md().noGUclosure, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_CIRC "

Default model behavior to treat a molecule as a circular RNA (DNA)  

See Also
--------
RNA.md().circ, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_GQUAD "

Default model behavior regarding the treatment of G-Quadruplexes.  

See Also
--------
RNA.md().gquad, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_UNIQ_ML "

Default behavior of the model regarding unique multi-branch loop decomposition.  

See Also
--------
RNA.md().uniq_ML, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_ENERGY_SET "

Default model behavior on which energy set to use.  

See Also
--------
RNA.md().energy_set, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_BACKTRACK "

Default model behavior with regards to backtracking of structures.  

See Also
--------
RNA.md().backtrack, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_BACKTRACK_TYPE "

Default model behavior on what type of backtracking to perform.  

See Also
--------
RNA.md().backtrack_type, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_COMPUTE_BPP "

Default model behavior with regards to computing base pair probabilities.  

See Also
--------
RNA.md().compute_bpp, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_MAX_BP_SPAN "

Default model behavior for the allowed maximum base pair span.  

See Also
--------
RNA.md().max_bp_span, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_WINDOW_SIZE "

Default model behavior for the sliding window approach.  

See Also
--------
RNA.md().window_size, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_LOG_ML "

Default model behavior on how to evaluate the energy contribution of multi-branch loops.  

See Also
--------
RNA.md().logML, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_ALI_OLD_EN "

Default model behavior for consensus structure energy evaluation.  

See Also
--------
RNA.md().oldAliEn, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_ALI_RIBO "

Default model behavior for consensus structure co-variance contribution assessment.  

See Also
--------
RNA.md().ribo, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_ALI_CV_FACT "

Default model behavior for weighting the co-variance score in consensus structure prediction.  

See Also
--------
RNA.md().cv_fact, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_ALI_NC_FACT "

Default model behavior for weighting the nucleotide conservation? in consensus structure prediction.  

See Also
--------
RNA.md().nc_fact, RNA.md_defaults_reset(), RNA.md.reset()  
";

%feature("docstring") MODEL_DEFAULT_PF_SMOOTH "
";

%feature("docstring") MODEL_DEFAULT_SALT "

Default model salt concentration (M)  
";

%feature("docstring") MODEL_DEFAULT_SALT_MLLOWER "

Default model lower bound of multiloop size for salt correction fiting.  
";

%feature("docstring") MODEL_DEFAULT_SALT_MLUPPER "

Default model upper bound of multiloop size for salt correction fiting.  
";

%feature("docstring") MODEL_DEFAULT_SALT_DPXINIT "

Default model value to turn off user-provided salt correction for duplex initializtion.  
";

%feature("docstring") MODEL_SALT_DPXINIT_FACT_RNA "
";

%feature("docstring") MODEL_SALT_DPXINIT_FACT_DNA "
";

%feature("docstring") MODEL_DEFAULT_SALT_DPXINIT_FACT "
";

%feature("docstring") MODEL_HELICAL_RISE_RNA "
";

%feature("docstring") MODEL_HELICAL_RISE_DNA "
";

%feature("docstring") MODEL_DEFAULT_HELICAL_RISE "

Default helical rise.  
";

%feature("docstring") MODEL_BACKBONE_LENGTH_RNA "
";

%feature("docstring") MODEL_BACKBONE_LENGTH_DNA "
";

%feature("docstring") MODEL_DEFAULT_BACKBONE_LENGTH "

Default backbone length.  
";

%feature("docstring") MAXALPHA "

Maximal length of alphabet.  
";

%feature("docstring") model_detailsT "
";

// File: group__energy__parameters.xml

%feature("docstring") vrna_params "

Get a data structure containing prescaled free energy parameters.  

If a NULL pointer is passed for the model details parameter, the default model parameters are stored
within the requested RNA.param() structure.  

Parameters
----------
md : RNA.md() *
    A pointer to the model details to store inside the structure (Maybe NULL)  

Returns
-------
RNA.param() *  
    A pointer to the memory location where the requested parameters are stored  

See Also
--------
RNA.md(), RNA.md.reset(), RNA.exp_params()  
";

%feature("docstring") vrna_params_copy "

Get a copy of the provided free energy parameters.  

If NULL is passed as parameter, a default set of energy parameters is created and returned.  

Parameters
----------
par : RNA.param() *
    The free energy parameters that are to be copied (Maybe NULL)  

Returns
-------
RNA.param() *  
    A copy or a default set of the (provided) parameters  

See Also
--------
RNA.params(), RNA.param()  
";

%feature("docstring") vrna_exp_params "

Get a data structure containing prescaled free energy parameters already transformed to Boltzmann
factors.  

This function returns a data structure that contains all necessary precomputed energy contributions
for each type of loop.  

In contrast to RNA.params(), the free energies within this data structure are stored as their
Boltzmann factors, i.e.  

:math:`exp(-E / kT)`  

where :math:`E` is the free energy.  

If a NULL pointer is passed for the model details parameter, the default model parameters are stored
within the requested RNA.exp_param() structure.  

Parameters
----------
md : RNA.md() *
    A pointer to the model details to store inside the structure (Maybe NULL)  

Returns
-------
RNA.exp_param() *  
    A pointer to the memory location where the requested parameters are stored  

See Also
--------
RNA.md(), RNA.md.reset(), RNA.params(), RNA.rescale_pf_params()  
";

%feature("docstring") vrna_exp_params_comparative "

Get a data structure containing prescaled free energy parameters already transformed to Boltzmann
factors (alifold version)  

If a NULL pointer is passed for the model details parameter, the default model parameters are stored
within the requested RNA.exp_param() structure.  

Parameters
----------
n_seq : unsigned int
    The number of sequences in the alignment  
md : RNA.md() *
    A pointer to the model details to store inside the structure (Maybe NULL)  

Returns
-------
RNA.exp_param() *  
    A pointer to the memory location where the requested parameters are stored  

See Also
--------
RNA.md(), RNA.md.reset(), RNA.exp_params(), RNA.params()  
";

%feature("docstring") vrna_exp_params_copy "

Get a copy of the provided free energy parameters (provided as Boltzmann factors)  

If NULL is passed as parameter, a default set of energy parameters is created and returned.  

Parameters
----------
par : RNA.exp_param() *
    The free energy parameters that are to be copied (Maybe NULL)  

Returns
-------
RNA.exp_param() *  
    A copy or a default set of the (provided) parameters  

See Also
--------
RNA.exp_params(), RNA.exp_param()  
";

%feature("docstring") vrna_fold_compound_t::params_subst "

Update/Reset energy parameters data structure within a RNA.fold_compound().  

Passing NULL as second argument leads to a reset of the energy parameters within fc to their default
values. Otherwise, the energy parameters provided will be copied over into fc.  

**SWIG Wrapper Notes**
    This function is attached to RNA.fc() objects as overloaded `params_subst()` method.  

    When no parameter is passed, the resulting action is the same as passing `NULL` as second
    parameter to RNA.fold_compound.params_subst(), i.e. resetting the parameters to the global
    defaults. See,
    e.g.  :py:meth:`RNA.fold_compound.params_subst()` in the :doc:`/api_python`.  

Parameters
----------
par : RNA.param() *
    The energy parameters used to substitute those within fc (Maybe NULL)  

See Also
--------
RNA.fold_compound.params_reset(), RNA.param(), RNA.md(), RNA.params()  
";

%feature("docstring") vrna_fold_compound_t::exp_params_subst "

Update the energy parameters for subsequent partition function computations.  

This function can be used to properly assign new energy parameters for partition function
computations to a RNA.fold_compound(). For this purpose, the data of the provided pointer `params`
will be copied into `fc` and a recomputation of the partition function scaling factor is issued, if
the `pf_scale` attribute of `params` is less than `1.0`.  

Passing NULL as second argument leads to a reset of the energy parameters within fc to their default
values  

**SWIG Wrapper Notes**
    This function is attached to RNA.fc() objects as overloaded `exp_params_subst()` method.  

    When no parameter is passed, the resulting action is the same as passing `NULL` as second
    parameter to RNA.fold_compound.exp_params_subst(), i.e. resetting the parameters to the global
    defaults. See,
    e.g.  :py:meth:`RNA.fold_compound.exp_params_subst()` in the :doc:`/api_python`.  

Parameters
----------
params : RNA.exp_param() *
    A pointer to the new energy parameters  

See Also
--------
RNA.fold_compound.exp_params_reset(), RNA.fold_compound.exp_params_rescale(), RNA.exp_param(),
RNA.md(), RNA.exp_params()  
";

%feature("docstring") vrna_fold_compound_t::exp_params_rescale "

Rescale Boltzmann factors for partition function computations.  

This function may be used to (automatically) rescale the Boltzmann factors used in partition
function computations. Since partition functions over subsequences can easily become extremely
large, the RNAlib internally rescales them to avoid numerical over- and/or underflow. Therefore, a
proper scaling factor :math:`s` needs to be chosen that in turn is then used to normalize the
corresponding partition functions :math:`\\hat{q}[i,j] = q[i,j] / s^{(j-i+1)}`.  

This function provides two ways to automatically adjust the scaling factor.  

1.  Automatic guess  
2.  Automatic adjustment according to MFE  

Passing `NULL` as second parameter activates the *automatic guess mode*. Here, the scaling factor is
recomputed according to a mean free energy of `184.3*length` cal for random sequences.
On the other hand, if the MFE for a sequence is known, it can be used to recompute a more robust
scaling factor, since it represents the lowest free energy of the entire ensemble of structures,
i.e. the highest Boltzmann factor. To activate this second mode of *automatic adjustment according
to MFE*, a pointer to the MFE value needs to be passed as second argument. This value is then taken
to compute the scaling factor as :math:`s = exp((sfact * MFE) / kT / length )`, where sfact is an
additional scaling weight located in the RNA.md() data structure of `exp_params` in `fc`.  

Note
----
This recomputation only takes place if the `pf_scale` attribute of the `exp_params` data structure
contained in `fc` has a value below `1.0`.  

The computed scaling factor :math:`s` will be stored as `pf_scale` attribute of the `exp_params`
data structure in `fc`.  

**SWIG Wrapper Notes**
    This function is attached to RNA.fc() objects as overloaded `exp_params_rescale()` method.  

    When no parameter is passed to this method, the resulting action is the same as passing `NULL`
    as second parameter to RNA.fold_compound.exp_params_rescale(), i.e. default scaling of the
    partition
    function. Passing an energy in kcal/mol, e.g. as retrieved by a previous call to the `mfe()`
    method, instructs all subsequent calls to scale the partition function accordingly. See, e.g.
    :py:meth:`RNA.fold_compound.exp_params_rescale()` in the :doc:`/api_python`.  

Parameters
----------
mfe : double *
    A pointer to the MFE (in kcal/mol) or NULL  

See Also
--------
RNA.fold_compound.exp_params_subst(), RNA.md(), RNA.exp_param(), RNA.fold_compound()  
";

%feature("docstring") vrna_fold_compound_t::params_reset "

Reset free energy parameters within a RNA.fold_compound() according to provided, or default model
details.  

This function allows one to rescale free energy parameters for subsequent structure prediction or
evaluation according to a set of model details, e.g. temperature values. To do so, the caller
provides either a pointer to a set of model details to be used for rescaling, or NULL if global
default setting should be used.  

**SWIG Wrapper Notes**
    This function is attached to RNA.fc() objects as overloaded `params_reset()` method.  

    When no parameter is passed to this method, the resulting action is the same as passing `NULL`
    as second parameter to RNA.fold_compound.params_reset(), i.e. global default model settings are
    used. Passing
    an object of type RNA.md() resets the fold compound according to the specifications stored
    within the RNA.md() object. See, e.g.  :py:meth:`RNA.fold_compound.params_reset()` in the
    :doc:`/api_python`.  

Parameters
----------
md : RNA.md() *
    A pointer to the new model details (or NULL for reset to defaults)  

See Also
--------
RNA.fold_compound.exp_params_reset(), RNA.params_subs()  
";

%feature("docstring") vrna_fold_compound_t::exp_params_reset "

Reset Boltzmann factors for partition function computations within a RNA.fold_compound() according
to provided, or default model details.  

This function allows one to rescale Boltzmann factors for subsequent partition function computations
according to a set of model details, e.g. temperature values. To do so, the caller provides either a
pointer to a set of model details to be used for rescaling, or NULL if global default setting should
be used.  

**SWIG Wrapper Notes**
    This function is attached to RNA.fc() objects as overloaded `exp_params_reset()` method.  

    When no parameter is passed to this method, the resulting action is the same as passing `NULL`
    as second parameter to RNA.fold_compound.exp_params_reset(), i.e. global default model settings
    are used.
    Passing an object of type RNA.md() resets the fold compound according to the specifications
    stored within the RNA.md() object. See, e.g.  :py:meth:`RNA.fold_compound.exp_params_reset()`
    in the :doc:`/api_python`.  

Parameters
----------
md : RNA.md() *
    A pointer to the new model details (or NULL for reset to defaults)  

See Also
--------
RNA.fold_compound.params_reset(), RNA.fold_compound.exp_params_subst(),
RNA.fold_compound.exp_params_rescale()  
";

%feature("docstring") vrna_params_prepare "
";

%feature("docstring") get_parameter_copy "
";

%feature("docstring") get_scaled_pf_parameters "

get a data structure of type RNA.exp_param() which contains the Boltzmann weights of several energy
parameters scaled according to the current temperature  

.. deprecated:: 2.6.4
    Use RNA.exp_params() instead!  

Returns
-------
RNA.exp_param() *  
    The data structure containing Boltzmann weights for use in partition function calculations  
";

%feature("docstring") get_boltzmann_factors "

Get precomputed Boltzmann factors of the loop type dependent energy contributions with independent
thermodynamic temperature.  

This function returns a data structure that contains all necessary precalculated Boltzmann factors
for each loop type contribution.  
 In contrast to get_scaled_pf_parameters(), this function enables setting of independent
temperatures for both, the individual energy contributions as well as the thermodynamic temperature
used in :math:`exp(-\\Delta G / kT)`  

.. deprecated:: 2.6.4
    Use RNA.exp_params() instead!  

Parameters
----------
temperature : double
    The temperature in degrees Celcius used for (re-)scaling the energy contributions  
betaScale : double
    A scaling value that is used as a multiplication factor for the absolute temperature of the
    system  
md : RNA.md()
    The model details to be used  
pf_scale : double
    The scaling factor for the Boltzmann factors  

Returns
-------
RNA.exp_param() *  
    A set of precomputed Boltzmann factors  

See Also
--------
get_scaled_pf_parameters(), get_boltzmann_factor_copy()  
";

%feature("docstring") get_boltzmann_factor_copy "

Get a copy of already precomputed Boltzmann factors.  

.. deprecated:: 2.6.4
    Use RNA.exp_params_copy() instead!  

Parameters
----------
parameters : RNA.exp_param() *
    The input data structure that shall be copied  

Returns
-------
RNA.exp_param() *  
    A copy of the provided Boltzmann factor data set  

See Also
--------
get_boltzmann_factors(), get_scaled_pf_parameters()  
";

%feature("docstring") get_scaled_alipf_parameters "

Get precomputed Boltzmann factors of the loop type dependent energy contributions (alifold variant)  

.. deprecated:: 2.6.4
    Use RNA.exp_params_comparative() instead!  
";

%feature("docstring") get_boltzmann_factors_ali "

Get precomputed Boltzmann factors of the loop type dependent energy contributions (alifold variant)
with independent thermodynamic temperature.  

.. deprecated:: 2.6.4
    Use RNA.exp_params_comparative() instead!  
";

%feature("docstring") scale_parameters "

Get precomputed energy contributions for all the known loop types.  

.. deprecated:: 2.6.4
    Use RNA.params() instead!  

Note
----
OpenMP: This function relies on several global model settings variables and thus is not to be
considered threadsafe. See get_scaled_parameters() for a completely threadsafe implementation.  

Returns
-------
RNA.param() *  
    A set of precomputed energy contributions  
";

%feature("docstring") get_scaled_parameters "

Get precomputed energy contributions for all the known loop types.  

Call this function to retrieve precomputed energy contributions, i.e. scaled according to the
temperature passed. Furthermore, this function assumes a data structure that contains the model
details as well, such that subsequent folding recursions are able to retrieve the correct model
settings  

.. deprecated:: 2.6.4
    Use RNA.params() instead!  

Parameters
----------
temperature : double
    The temperature in degrees Celcius  
md : RNA.md()
    The model details  

Returns
-------
RNA.param() *  
    precomputed energy contributions and model settings  

See Also
--------
RNA.md(), set_model_details()  
";

%feature("docstring") copy_parameters "
";

%feature("docstring") set_parameters "
";

%feature("docstring") scale_pf_parameters "
";

%feature("docstring") copy_pf_param "
";

%feature("docstring") set_pf_param "
";

%feature("docstring") GQUAD_MAX_STACK_SIZE "
";

%feature("docstring") GQUAD_MIN_STACK_SIZE "
";

%feature("docstring") GQUAD_MAX_LINKER_LENGTH "
";

%feature("docstring") GQUAD_MIN_LINKER_LENGTH "
";

%feature("docstring") GQUAD_MIN_BOX_SIZE "
";

%feature("docstring") GQUAD_MAX_BOX_SIZE "
";

// File: group__domains.xml

// File: group__domains__up.xml

%feature("docstring") vrna_ud_motifs_centroid "

Detect unstructured domains in centroid structure.  

Given a centroid structure and a set of unstructured domains compute the list of unstructured domain
motifs present in the centroid. Since we do not explicitly annotate unstructured domain motifs in
dot-bracket strings, this function can be used to check for the presence and location of
unstructured domain motifs under the assumption that the dot-bracket string is the centroid
structure of the equiibrium ensemble.  

Parameters
----------
fc : RNA.fold_compound() *
    The fold_compound data structure with pre-computed equilibrium probabilities and model settings  
structure : const char *
    The centroid structure in dot-bracket notation  

Returns
-------
RNA.ud_motif() *  
    A list of unstructured domain motifs (possibly NULL). The last element terminates the list with
    `start=0`, `number=-1`  

See Also
--------
RNA.fold_compound.centroid()  
";

%feature("docstring") vrna_ud_motifs_MEA "

Detect unstructured domains in MEA structure.  

Given an MEA structure and a set of unstructured domains compute the list of unstructured domain
motifs present in the MEA structure. Since we do not explicitly annotate unstructured domain motifs
in dot-bracket strings, this function can be used to check for the presence and location of
unstructured domain motifs under the assumption that the dot-bracket string is the MEA structure of
the equiibrium ensemble.  

Parameters
----------
fc : RNA.fold_compound() *
    The fold_compound data structure with pre-computed equilibrium probabilities and model settings  
structure : const char *
    The MEA structure in dot-bracket notation  
probability_list : RNA.ep() *
    The list of probabilities to extract the MEA structure from  

Returns
-------
RNA.ud_motif() *  
    A list of unstructured domain motifs (possibly NULL). The last element terminates the list with
    `start=0`, `number=-1`  

See Also
--------
MEA()  
";

%feature("docstring") vrna_ud_motifs_MFE "

Detect unstructured domains in MFE structure.  

Given an MFE structure and a set of unstructured domains compute the list of unstructured domain
motifs present in the MFE structure. Since we do not explicitly annotate unstructured domain motifs
in dot-bracket strings, this function can be used to check for the presence and location of
unstructured domain motifs under the assumption that the dot-bracket string is the MFE structure of
the equiibrium ensemble.  

Parameters
----------
fc : RNA.fold_compound() *
    The fold_compound data structure with model settings  
structure : const char *
    The MFE structure in dot-bracket notation  

Returns
-------
RNA.ud_motif() *  
    A list of unstructured domain motifs (possibly NULL). The last element terminates the list with
    `start=0`, `number=-1`  

See Also
--------
RNA.fold_compound.mfe()  
";

%feature("docstring") vrna_fold_compound_t::ud_add_motif "

Add an unstructured domain motif, e.g. for ligand binding.  

This function adds a ligand binding motif and the associated binding free energy to the RNA.ud()
attribute of a RNA.fold_compound(). The motif data will then be used in subsequent secondary
structure predictions. Multiple calls to this function with different motifs append all additional
data to a list of ligands, which all will be evaluated. Ligand motif data can be removed from the
RNA.fold_compound() again using the RNA.fold_compound.ud_remove() function. The loop type parameter allows one
to limit the ligand binding to particular loop type, such as the exterior loop, hairpin loops,
interior loops, or multibranch loops.  

Parameters
----------
motif : const char *
    The sequence motif the ligand binds to  
motif_en : double
    The binding free energy of the ligand in kcal/mol  
motif_name : const char *
    The name/id of the motif (may be `NULL`)  
loop_type : unsigned int
    The loop type the ligand binds to  

See Also
--------
RNA.UNSTRUCTURED_DOMAIN_EXT_LOOP, RNA.UNSTRUCTURED_DOMAIN_HP_LOOP,
RNA.UNSTRUCTURED_DOMAIN_INT_LOOP, RNA.UNSTRUCTURED_DOMAIN_MB_LOOP,
RNA.UNSTRUCTURED_DOMAIN_ALL_LOOPS, RNA.fold_compound.ud_remove()  
";

%feature("docstring") vrna_fold_compound_t::ud_remove "

Remove ligand binding to unpaired stretches.  

This function removes all ligand motifs that were bound to a RNA.fold_compound() using the
RNA.fold_compound.ud_add_motif() function.  

**SWIG Wrapper Notes**
    This function is attached as method `ud_remove()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.ud_remove()` in the :doc:`/api_python`.  

Parameters
----------  
";

%feature("docstring") vrna_fold_compound_t::ud_set_data "

Attach an auxiliary data structure.  

This function binds an arbitrary, auxiliary data structure for user-implemented ligand binding. The
optional callback `free_cb` will be passed the bound data structure whenever the
RNA.fold_compound() is removed from memory to avoid memory leaks.  

**SWIG Wrapper Notes**
    This function is attached as method `ud_set_data()` to objects of type `fold_compound`. See,
    e.g.  :py:meth:`RNA.fold_compound.ud_set_data()` in the :doc:`/api_python`.  

Parameters
----------
data : void *
    A pointer to the auxiliary data structure  
free_cb : RNA.auxdata_free
    A pointer to a callback function that free's memory occupied by `data`  

See Also
--------
RNA.fold_compound.ud_set_prod_rule_cb(), RNA.fold_compound.ud_set_exp_prod_rule_cb(),
RNA.fold_compound.ud_remove()  
";

%feature("docstring") vrna_fold_compound_t::ud_set_prod_rule_cb "

Attach production rule callbacks for free energies computations.  

Use this function to bind a user-implemented grammar extension for unstructured domains.  

The callback `e_cb` needs to evaluate the free energy contribution :math:`f(i,j)` of the unpaired
segment :math:`[i,j]`. It will be executed in each of the regular secondary structure production
rules. Whenever the callback is passed the RNA.UNSTRUCTURED_DOMAIN_MOTIF flag via its `loop_type`
parameter the contribution of any ligand that consecutively binds from position :math:`i` to
:math:`j` (the white box) is requested. Otherwise, the callback usually performs a lookup in the
precomputed `B` matrices. Which `B` matrix is addressed will be indicated by the flags
RNA.UNSTRUCTURED_DOMAIN_EXT_LOOP,
RNA.UNSTRUCTURED_DOMAIN_HP_LOOPRNA.UNSTRUCTURED_DOMAIN_INT_LOOP, and
RNA.UNSTRUCTURED_DOMAIN_MB_LOOP. As their names already imply, they specify exterior loops (`F`
production rule), hairpin loops and interior loops (`C` production rule), and multibranch loops (`M`
and `M1` production rule).  


The `pre_cb` callback will be executed as a pre-processing step right before the regular secondary
structure rules. Usually one would use this callback to fill the dynamic programming matrices `U`
and preparations of the auxiliary data structure RNA.unstructured_domain().data  


**SWIG Wrapper Notes**
    This function is attached as method `ud_set_prod_rule_cb()` to objects of type `fold_compound`.
    See, e.g.  :py:meth:`RNA.fold_compound.ud_set_prod_rule_cb()` in the :doc:`/api_python`.  

Parameters
----------
pre_cb : RNA.ud_production
    A pointer to a callback function for the `B` production rule  
e_cb : RNA.ud
    A pointer to a callback function for free energy evaluation  
";

%feature("docstring") vrna_fold_compound_t::ud_set_exp_prod_rule_cb "

Attach production rule for partition function.  

This function is the partition function companion of RNA.fold_compound.ud_set_prod_rule_cb().  

Use it to bind callbacks to (i) fill the `U` production rule dynamic programming matrices and/or
prepare the RNA.unstructured_domain().data, and (ii) provide a callback to retrieve partition
functions for subsegments :math:`[i,j]`.  



**SWIG Wrapper Notes**
    This function is attached as method `ud_set_exp_prod_rule_cb()` to objects of type
    `fold_compound`. See, e.g.  :py:meth:`RNA.fold_compound.ud_set_exp_prod_rule_cb()` in the
    :doc:`/api_python`.  

Parameters
----------
pre_cb : RNA.ud_exp_production
    A pointer to a callback function for the `B` production rule  
exp_e_cb : RNA.ud_exp
    A pointer to a callback function that retrieves the partition function for a segment
    :math:`[i,j]` that may be bound by one or more ligands.  

See Also
--------
RNA.fold_compound.ud_set_prod_rule_cb()  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_EXT_LOOP "

Flag to indicate ligand bound to unpiared stretch in the exterior loop.  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_HP_LOOP "

Flag to indicate ligand bound to unpaired stretch in a hairpin loop.  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_INT_LOOP "

Flag to indicate ligand bound to unpiared stretch in an interior loop.  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_MB_LOOP "

Flag to indicate ligand bound to unpiared stretch in a multibranch loop.  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_MOTIF "

Flag to indicate ligand binding without additional unbound nucleotides (motif-only)  
";

%feature("docstring") UNSTRUCTURED_DOMAIN_ALL_LOOPS "

Flag to indicate ligand bound to unpiared stretch in any loop (convenience macro)  
";

// File: group__domains__struc.xml

// File: group__constraints.xml

%feature("docstring") vrna_message_constraint_options "

Print a help message for pseudo dot-bracket structure constraint characters to stdout. (constraint
support is specified by option parameter)  

Currently available options are:  RNA.CONSTRAINT_DB_PIPE (paired with another base)
RNA.CONSTRAINT_DB_DOT (no constraint at all)  RNA.CONSTRAINT_DB_X (base must not pair)
RNA.CONSTRAINT_DB_ANG_BRACK (paired downstream/upstream)  RNA.CONSTRAINT_DB_RND_BRACK (base i
pairs base j)  
 pass a collection of options as one value like this:  

    RNA.message_constraints(option_1 | option_2 | option_n)  

Parameters
----------
option : unsigned int
    Option switch that tells which constraint help will be printed  

See Also
--------
RNA.message_constraint_options_all(), RNA.fold_compound.constraints_add(), RNA.CONSTRAINT_DB,
RNA.CONSTRAINT_DB_PIPE, RNA.CONSTRAINT_DB_DOT, RNA.CONSTRAINT_DB_X, RNA.CONSTRAINT_DB_ANG_BRACK,
RNA.CONSTRAINT_DB_RND_BRACK, RNA.CONSTRAINT_DB_INTERMOL, RNA.CONSTRAINT_DB_INTRAMOL  
";

%feature("docstring") vrna_message_constraint_options_all "

Print structure constraint characters to stdout (full constraint support)  

See Also
--------
RNA.message_constraint_options(), RNA.fold_compound.constraints_add(), RNA.CONSTRAINT_DB,
RNA.CONSTRAINT_DB_PIPE, RNA.CONSTRAINT_DB_DOT, RNA.CONSTRAINT_DB_X, RNA.CONSTRAINT_DB_ANG_BRACK,
RNA.CONSTRAINT_DB_RND_BRACK, RNA.CONSTRAINT_DB_INTERMOL, RNA.CONSTRAINT_DB_INTRAMOL  
";

%feature("docstring") CONSTRAINT_FILE "

Flag for RNA.fold_compound.constraints_add() to indicate that constraints are present in a text file.  

See Also
--------
RNA.fold_compound.constraints_add()  

.. deprecated:: 2.6.4
    Use 0 instead!  
";

%feature("docstring") CONSTRAINT_SOFT_MFE "

Indicate generation of constraints for MFE folding.  

.. deprecated:: 2.2.6
    This flag has no meaning anymore, since constraints are now always stored! (since v2.2.6)  
";

%feature("docstring") CONSTRAINT_SOFT_PF "

Indicate generation of constraints for partition function computation.  

.. deprecated:: 2.6.4
    Use RNA.OPTION_PF instead!  
";

%feature("docstring") DECOMP_PAIR_HP "

Flag passed to generic softt constraints callback to indicate hairpin loop decomposition step.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a hairpin loop enclosed by the base pair :math:`(i,j)`.  

";

%feature("docstring") DECOMP_PAIR_IL "

Indicator for interior loop decomposition step.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an interior loop enclosed by the base pair :math:`(i,j)`, and enclosing the base pair
:math:`(k,l)`.  

";

%feature("docstring") DECOMP_PAIR_ML "

Indicator for multibranch loop decomposition step.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop enclosed by the base pair :math:`(i,j)`, and consisting of some
enclosed multi loop content from k to l.  

";

%feature("docstring") DECOMP_ML_ML_ML "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, which will be decomposed into two
multibranch loop parts :math:`[i:k]`, and :math:`[l:j]`.  

";

%feature("docstring") DECOMP_ML_STEM "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, which will be considered a single
stem branching off with base pair :math:`(k,l)`.  

";

%feature("docstring") DECOMP_ML_ML "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, which will be decomposed into a
(usually) smaller multibranch loop part :math:`[k:l]`.  

";

%feature("docstring") DECOMP_ML_UP "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, which will be considered a
multibranch loop part that only consists of unpaired nucleotides.  

";

%feature("docstring") DECOMP_ML_ML_STEM "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, which will decomposed into a
multibranch loop part :math:`[i:k]`, and a stem with enclosing base pair :math:`(l,j)`.  

";

%feature("docstring") DECOMP_ML_COAXIAL "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, where two stems with enclosing
pairs :math:`(i,k)` and :math:`(l,j)` are coaxially stacking onto each other.  

";

%feature("docstring") DECOMP_ML_COAXIAL_ENC "

Indicator for decomposition of multibranch loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates a multibranch loop part in the interval :math:`[i:j]`, where two stems with enclosing
pairs :math:`(i,k)` and :math:`(l,j)` are coaxially stacking onto each other.  

";

%feature("docstring") DECOMP_EXT_EXT "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be decomposed into a
(usually) smaller exterior loop part :math:`[k:l]`.  

";

%feature("docstring") DECOMP_EXT_UP "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be considered as an
exterior loop component consisting of only unpaired nucleotides.  

";

%feature("docstring") DECOMP_EXT_STEM "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be considered a stem with
enclosing pair :math:`(k,l)`.  

";

%feature("docstring") DECOMP_EXT_EXT_EXT "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be decomposed into two
exterior loop parts :math:`[i:k]` and :math:`[l:j]`.  

";

%feature("docstring") DECOMP_EXT_STEM_EXT "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be decomposed into a stem
branching off with base pair :math:`(i,k)`, and an exterior loop part :math:`[l:j]`.  

";

%feature("docstring") DECOMP_EXT_STEM_OUTSIDE "

Indicator for decomposition of exterior loop part.  
";

%feature("docstring") DECOMP_EXT_EXT_STEM "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be decomposed into an
exterior loop part :math:`[i:k]`, and a stem branching off with base pair :math:`(l,j)`.  

";

%feature("docstring") DECOMP_EXT_EXT_STEM1 "

Indicator for decomposition of exterior loop part.  

This flag notifies the soft or hard constraint callback function that the current decomposition step
evaluates an exterior loop part in the interval :math:`[i:j]`, which will be decomposed into an
exterior loop part :math:`[i:k]`, and a stem branching off with base pair :math:`(l,j-1)`.  

";

// File: group__hard__constraints.xml

%feature("docstring") vrna_fold_compound_t::constraints_add "

Add constraints to a RNA.fold_compound() data structure.  

Use this function to add/update the hard/soft constraints The function allows for passing a string
'constraint' that can either be a filename that points to a constraints definition file or it may be
a pseudo dot-bracket notation indicating hard constraints. For the latter, the user has to pass the
RNA.CONSTRAINT_DB option. Also, the user has to specify, which characters are allowed to be
interpreted as constraints by passing the corresponding options via the third parameter.  

The following is an example for adding hard constraints given in pseudo dot-bracket notation. Here,
`fc` is the RNA.fold_compound() object, `structure` is a char array with the hard constraint in
dot-bracket notation, and `enforceConstraints` is a flag indicating whether or not constraints for
base pairs should be enforced instead of just doing a removal of base pair that conflict with the
constraint.  


In constrat to the above, constraints may also be read from file:  


Parameters
----------
constraint : const char *
    A string with either the filename of the constraint definitions or a pseudo dot-bracket notation
    of the hard constraint. May be NULL.  
options : unsigned int
    The option flags  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.hc_add_up(), RNA.hc_add_up_batch()
RNA.hc_add_bp_unspecific(),
RNA.fold_compound.hc_add_bp(), RNA.fold_compound.hc_init(), RNA.fold_compound.sc_set_up(),
RNA.fold_compound.sc_set_bp(), RNA.fold_compound.sc_add_SHAPE_deigan(),
RNA.fold_compound.sc_add_SHAPE_zarringhalam(), RNA.hc_free(), RNA.sc_free(), RNA.CONSTRAINT_DB,
RNA.CONSTRAINT_DB_DEFAULT, RNA.CONSTRAINT_DB_PIPE, RNA.CONSTRAINT_DB_DOT, RNA.CONSTRAINT_DB_X,
RNA.CONSTRAINT_DB_ANG_BRACK, RNA.CONSTRAINT_DB_RND_BRACK, RNA.CONSTRAINT_DB_INTRAMOL,
RNA.CONSTRAINT_DB_INTERMOL, RNA.CONSTRAINT_DB_GQUAD  
";

%feature("docstring") vrna_fold_compound_t::hc_init "

Initialize/Reset hard constraints to default values.  

This function resets the hard constraints to their default values, i.e. all positions may be
unpaired in all contexts, and base pairs are allowed in all contexts, if they resemble canonical
pairs. Previously set hard constraints will be removed before initialization.  

**SWIG Wrapper Notes**
    This function is attached as method `hc_init()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.hc_init()` in the :doc:`/api_python` .  

Parameters
----------

See Also
--------
RNA.fold_compound.hc_add_bp(), RNA.fold_compound.hc_add_bp_nonspecific(),
RNA.fold_compound.hc_add_up()  
";

%feature("docstring") vrna_fold_compound_t::hc_add_up "

Make a certain nucleotide unpaired.  

Parameters
----------
i : int
    The position that needs to stay unpaired (1-based)  
option : unsigned char
    The options flag indicating how/where to store the hard constraints  

See Also
--------
RNA.fold_compound.hc_add_bp(), RNA.fold_compound.hc_add_bp_nonspecific(),
RNA.fold_compound.hc_init(), RNA.CONSTRAINT_CONTEXT_EXT_LOOP,
RNA.CONSTRAINT_CONTEXT_HP_LOOP, RNA.CONSTRAINT_CONTEXT_INT_LOOP, RNA.CONSTRAINT_CONTEXT_MB_LOOP,
RNA.CONSTRAINT_CONTEXT_ALL_LOOPS  
";

%feature("docstring") vrna_hc_add_up_batch "

Apply a list of hard constraints for single nucleotides.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() the hard constraints are associated with  
constraints : RNA.hc_up() *
    The list off constraints to apply, last entry must have position attribute set to 0  
";

%feature("docstring") vrna_fold_compound_t::hc_add_bp "

Favorize/Enforce a certain base pair (i,j)  

Parameters
----------
i : int
    The 5' located nucleotide position of the base pair (1-based)  
j : int
    The 3' located nucleotide position of the base pair (1-based)  
option : unsigned char
    The options flag indicating how/where to store the hard constraints  

See Also
--------
RNA.fold_compound.hc_add_bp_nonspecific(), RNA.fold_compound.hc_add_up(),
RNA.fold_compound.hc_init(), RNA.CONSTRAINT_CONTEXT_EXT_LOOP,
RNA.CONSTRAINT_CONTEXT_HP_LOOP, RNA.CONSTRAINT_CONTEXT_INT_LOOP,
RNA.CONSTRAINT_CONTEXT_INT_LOOP_ENC, RNA.CONSTRAINT_CONTEXT_MB_LOOP,
RNA.CONSTRAINT_CONTEXT_MB_LOOP_ENC, RNA.CONSTRAINT_CONTEXT_ENFORCE,
RNA.CONSTRAINT_CONTEXT_ALL_LOOPS  
";

%feature("docstring") vrna_fold_compound_t::hc_add_bp_nonspecific "

Enforce a nucleotide to be paired (upstream/downstream)  

Parameters
----------
i : int
    The position that needs to stay unpaired (1-based)  
d : int
    The direction of base pairing ( :math:`d < 0`: pairs upstream, :math:`d > 0`: pairs downstream,
    :math:`d == 0`: no direction)  
option : unsigned char
    The options flag indicating in which loop type context the pairs may appear  

See Also
--------
RNA.fold_compound.hc_add_bp(), RNA.fold_compound.hc_add_up(), RNA.fold_compound.hc_init(),
RNA.CONSTRAINT_CONTEXT_EXT_LOOP,
RNA.CONSTRAINT_CONTEXT_HP_LOOP, RNA.CONSTRAINT_CONTEXT_INT_LOOP,
RNA.CONSTRAINT_CONTEXT_INT_LOOP_ENC, RNA.CONSTRAINT_CONTEXT_MB_LOOP,
RNA.CONSTRAINT_CONTEXT_MB_LOOP_ENC, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS  
";

%feature("docstring") vrna_hc_free "

Free the memory allocated by a RNA.hc() data structure.  

Use this function to free all memory that was allocated for a data structure of type RNA.hc() .  

See Also
--------
get_hard_constraints(), RNA.hc()  
";

%feature("docstring") vrna_fold_compound_t::hc_add_from_db "

Add hard constraints from pseudo dot-bracket notation.  

This function allows one to apply hard constraints from a pseudo dot-bracket notation. The `options`
parameter controls, which characters are recognized by the parser. Use the
RNA.CONSTRAINT_DB_DEFAULT convenience macro, if you want to allow all known characters  

**SWIG Wrapper Notes**
    This function is attached as method `hc_add_from_db()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.hc_add_from_db()` in the :doc:`/api_python` .  

Parameters
----------
constraint : const char *
    A pseudo dot-bracket notation of the hard constraint.  
options : unsigned int
    The option flags  

See Also
--------
RNA.CONSTRAINT_DB_PIPE, RNA.CONSTRAINT_DB_DOT, RNA.CONSTRAINT_DB_X, RNA.CONSTRAINT_DB_ANG_BRACK,
RNA.CONSTRAINT_DB_RND_BRACK, RNA.CONSTRAINT_DB_INTRAMOL, RNA.CONSTRAINT_DB_INTERMOL,
RNA.CONSTRAINT_DB_GQUAD  
";

%feature("docstring") CONSTRAINT_DB "

Flag for RNA.fold_compound.constraints_add() to indicate that constraint is passed in pseudo dot-bracket
notation.  

See Also
--------
RNA.fold_compound.constraints_add(), RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_ENFORCE_BP "

Switch for dot-bracket structure constraint to enforce base pairs.  

This flag should be used to really enforce base pairs given in dot-bracket constraint rather than
just weakly-enforcing them.  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_PIPE "

Flag that is used to indicate the pipe '|' sign in pseudo dot-bracket notation of hard constraints.  

Use this definition to indicate the pipe sign '|' (paired with another base)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_DOT "

dot '.' switch for structure constraints (no constraint at all)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_X "

'x' switch for structure constraint (base must not pair)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_RND_BRACK "

round brackets '(',')' switch for structure constraint (base i pairs base j)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_INTRAMOL "

Flag that is used to indicate the character 'l' in pseudo dot-bracket notation of hard constraints.  

Use this definition to indicate the usage of 'l' character (intramolecular pairs only)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_INTERMOL "

Flag that is used to indicate the character 'e' in pseudo dot-bracket notation of hard constraints.  

Use this definition to indicate the usage of 'e' character (intermolecular pairs only)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_DB_GQUAD "

'+' switch for structure constraint (base is involved in a gquad)  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  

Warnings
--------
This flag is for future purposes only! No implementation recognizes it yet.  
";

%feature("docstring") CONSTRAINT_DB_WUSS "

Flag to indicate Washington University Secondary Structure (WUSS) notation of the hard constraint
string.  

This secondary structure notation for RNAs is usually used as consensus secondary structure
(SS_cons) entry in Stockholm formatted files  
";

%feature("docstring") CONSTRAINT_DB_DEFAULT "

Switch for dot-bracket structure constraint with default symbols.  

This flag conveniently combines all possible symbols in dot-bracket notation for hard constraints
and RNA.CONSTRAINT_DB  

See Also
--------
RNA.fold_compound.hc_add_from_db(), RNA.fold_compound.constraints_add(),
RNA.message_constraint_options(),
RNA.message_constraint_options_all()  
";

%feature("docstring") CONSTRAINT_CONTEXT_EXT_LOOP "

Hard constraints flag, base pair in the exterior loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_HP_LOOP "

Hard constraints flag, base pair encloses hairpin loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_INT_LOOP "

Hard constraints flag, base pair encloses an interior loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_INT_LOOP_ENC "

Hard constraints flag, base pair encloses a multi branch loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_MB_LOOP "

Hard constraints flag, base pair is enclosed in an interior loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_MB_LOOP_ENC "

Hard constraints flag, base pair is enclosed in a multi branch loop.  
";

%feature("docstring") CONSTRAINT_CONTEXT_ALL_LOOPS "

Constraint context flag indicating any loop context.  
";

// File: group__soft__constraints.xml

%feature("docstring") vrna_fold_compound_t::sc_init "

Initialize an empty soft constraints data structure within a RNA.fold_compound().  

This function adds a proper soft constraints data structure to the RNA.fold_compound() data
structure. If soft constraints already exist within the fold compound, they are removed.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_init()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_init()` in the :doc:`/api_python` .  

Parameters
----------

See Also
--------
RNA.fold_compound.sc_set_bp(), RNA.fold_compound.sc_set_up(),
RNA.fold_compound.sc_add_SHAPE_deigan(), RNA.fold_compound.sc_add_SHAPE_zarringhalam(),
RNA.fold_compound.sc_remove(), RNA.fold_compound.sc_add(), RNA.fold_compound.sc_add_exp(),
RNA.sc_add_pre(), RNA.sc_add_post()  

Note
----
Accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE and RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_fold_compound_t::sc_set_bp "

Set soft constraints for paired nucleotides.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_set_bp()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_set_bp()` in the :doc:`/api_python` .  

Parameters
----------
constraints : const FLT_OR_DBL **
    A two-dimensional array of pseudo free energies in :math:`kcal / mol`  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    Non-zero on successful application of the constraint, 0 otherwise.  

See Also
--------
RNA.fold_compound.sc_add_bp(), RNA.fold_compound.sc_set_up(), RNA.fold_compound.sc_add_up()  

Note
----
This function replaces any pre-exisitng soft constraints with the ones supplied in `constraints`.  
";

%feature("docstring") vrna_fold_compound_t::sc_add_bp "

Add soft constraints for paired nucleotides.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `sc_add_bp()` to objects of type
    `fold_compound`. The method either takes arguments for a single base pair (i,j) with the
    corresponding energy value:  
 or an entire 2-dimensional matrix with dimensions n x n that stores free energy contributions for
any base pair (i,j) with :math:`1 \\leq i < j \\leq n`:  In both variants, the optional argument
`options` defaults to RNA.OPTION_DEFAULT. See, e.g.   :py:meth:`RNA.fold_compound.sc_add_bp()` in
the :doc:`/api_python` .  

Parameters
----------
i : int
    The 5' position of the base pair the soft constraint is added for  
j : int
    The 3' position of the base pair the soft constraint is added for  
energy : FLT_OR_DBL
    The free energy (soft-constraint) in :math:`kcal / mol`  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    Non-zero on successful application of the constraint, 0 otherwise.  

See Also
--------
RNA.fold_compound.sc_set_bp(), RNA.fold_compound.sc_set_up(), RNA.fold_compound.sc_add_up()  
";

%feature("docstring") vrna_fold_compound_t::sc_set_up "

Set soft constraints for unpaired nucleotides.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_set_up()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_set_up()` in the :doc:`/api_python` .  

Parameters
----------
constraints : const FLT_OR_DBL *
    A vector of pseudo free energies in :math:`kcal / mol`  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    Non-zero on successful application of the constraint, 0 otherwise.  

See Also
--------
RNA.fold_compound.sc_add_up(), RNA.fold_compound.sc_set_bp(), RNA.fold_compound.sc_add_bp()  

Note
----
This function replaces any pre-exisitng soft constraints with the ones supplied in `constraints`.  
";

%feature("docstring") vrna_fold_compound_t::sc_add_up "

Add soft constraints for unpaired nucleotides.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `sc_add_up()` to objects of type
    `fold_compound`. The method either takes arguments for a single nucleotide :math:`i` with the
    corresponding energy value:  
 or an entire vector that stores free energy contributions for each nucleotide :math:`i` with
:math:`1 \\leq i \\leq n`:  In both variants, the optional argument `options` defaults to
RNA.OPTION_DEFAULT. See, e.g.   :py:meth:`RNA.fold_compound.sc_add_up()` in the :doc:`/api_python`
.  

Parameters
----------
i : int
    The nucleotide position the soft constraint is added for  
energy : FLT_OR_DBL
    The free energy (soft-constraint) in :math:`kcal / mol`  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    Non-zero on successful application of the constraint, 0 otherwise.  

See Also
--------
RNA.fold_compound.sc_set_up(), RNA.fold_compound.sc_add_bp(), RNA.fold_compound.sc_set_bp()  
";

%feature("docstring") vrna_fold_compound_t::sc_remove "

Remove soft constraints from RNA.fold_compound().  

**SWIG Wrapper Notes**
    This function is attached as method `sc_remove()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_remove()` in the :doc:`/api_python` .  

Parameters
----------

Note
----
Accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE and RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_sc_free "

Free memory occupied by a RNA.sc() data structure.  

Parameters
----------
sc : RNA.sc() *
    The data structure to free from memory  
";

%feature("docstring") vrna_fold_compound_t::sc_add_data "

Add an auxiliary data structure for the generic soft constraints callback function.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_data()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.sc_add_data()` in the :doc:`/api_python` .  

Parameters
----------
data : void *
    A pointer to the data structure that holds required data for function 'f'  
free_data : RNA.auxdata_free
    A pointer to a function that free's the memory occupied by `data` (Maybe NULL)  

Returns
-------
int  
    Non-zero on successful binding the data (and free-function), 0 otherwise  

See Also
--------
RNA.fold_compound.sc_add(), RNA.fold_compound.sc_add_exp(), RNA.fold_compound.sc_add_bt()  
";

%feature("docstring") vrna_fold_compound_t::sc_add_f "

Bind a function pointer for generic soft constraint feature (MFE version)  

This function allows one to easily bind a function pointer and corresponding data structure to the
soft constraint part RNA.sc() of the RNA.fold_compound(). The function for evaluating the generic
soft constraint feature has to return a pseudo free energy :math:`\\hat{E}` in :math:`dacal/mol`,
where :math:`1 dacal/mol = 10 cal/mol`.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_add()` in the :doc:`/api_python` .  

Parameters
----------
f : RNA.sc
    A pointer to the function that evaluates the generic soft constraint feature  

Returns
-------
int  
    Non-zero on successful binding the callback function, 0 otherwise  

See Also
--------
RNA.fold_compound.sc_add_data(), RNA.fold_compound.sc_add_bt(), RNA.fold_compound.sc_add_exp()  
";

%feature("docstring") vrna_fold_compound_t::sc_add_bt "

Bind a backtracking function pointer for generic soft constraint feature.  

This function allows one to easily bind a function pointer to the soft constraint part RNA.sc() of
the RNA.fold_compound(). The provided function should be used for backtracking purposes in loop
regions that were altered via the generic soft constraint feature. It has to return an array of
RNA.basepair() data structures, were the last element in the list is indicated by a value of -1 in
it's i position.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_bt()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.sc_add_bt()` in the :doc:`/api_python` .  

Parameters
----------
f : RNA.sc_bt
    A pointer to the function that returns additional base pairs  

Returns
-------
int  
    Non-zero on successful binding the callback function, 0 otherwise  

See Also
--------
RNA.fold_compound.sc_add_data(), RNA.fold_compound.sc_add(), RNA.fold_compound.sc_add_exp()  
";

%feature("docstring") vrna_fold_compound_t::sc_add_exp_f "

Bind a function pointer for generic soft constraint feature (PF version)  

This function allows one to easily bind a function pointer and corresponding data structure to the
soft constraint part RNA.sc() of the RNA.fold_compound(). The function for evaluating the generic
soft constraint feature has to return a pseudo free energy :math:`\\hat{E}` as Boltzmann factor,
i.e. :math:`exp(- \\hat{E} / kT)`. The required unit for :math:`E` is :math:`cal/mol`.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_exp()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.sc_add_exp()` in the :doc:`/api_python` .  

Parameters
----------
exp : RNA.sc_exp
    A pointer to the function that evaluates the generic soft constraint feature  

Returns
-------
int  
    Non-zero on successful binding the callback function, 0 otherwise  

See Also
--------
RNA.fold_compound.sc_add_bt(), RNA.fold_compound.sc_add(), RNA.fold_compound.sc_add_data()  
";

// File: group__landscape.xml

// File: group__mfe.xml

// File: group__pf__fold.xml

%feature("docstring") vrna_pf_float_precision "

Find out whether partition function computations are using single precision floating points.  

Returns
-------
int  
    1 if single precision is used, 0 otherwise  

See Also
--------
FLT_OR_DBL  
";

// File: group__mfe__global.xml

%feature("docstring") vrna_fold_compound_t::mfe "

Compute minimum free energy and an appropriate secondary structure of an RNA sequence, or RNA
sequence alignment.  

Depending on the type of the provided RNA.fold_compound(), this function predicts the MFE for a
single sequence (or connected component of multiple sequences), or an averaged MFE for a sequence
alignment. If backtracking is activated, it also constructs the corresponding secondary structure,
or consensus structure. Therefore, the second parameter, *structure*, has to point to an allocated
block of memory with a size of at least :math:`\\mathrm{strlen}(\\mathrm{sequence})+1` to store the
backtracked MFE structure. (For consensus structures, this is the length of the alignment + 1. If
`NULL` is passed, no backtracking will be performed.  

**SWIG Wrapper Notes**
    This function is attached as method `mfe()` to objects of type `fold_compound`. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:meth:`RNA.fold_compound.mfe()` in the :doc:`/api_python`.  

Parameters
----------
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to (Maybe NULL)  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.fold_compound(), RNA.fold_compound(), RNA.fold(), RNA.circfold(),
RNA.fold_compound_comparative(), RNA.alifold(), RNA.circalifold()  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_fold_compound_t::mfe_dimer "

Compute the minimum free energy of two interacting RNA molecules.  

The code is analog to the RNA.fold_compound.mfe() function.  

.. deprecated:: 2.6.4
    This function is obsolete since RNA.fold_compound.mfe() can handle complexes multiple sequences
since v2.5.0.
    Use RNA.fold_compound.mfe() for connected component MFE instead and compute MFEs of unconnected
states
    separately.  

**SWIG Wrapper Notes**
    This function is attached as method `mfe_dimer()` to objects of type `fold_compound`. The
    parameter `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:meth:`RNA.fold_compound.mfe_dimer()` in the :doc:`/api_python`.  

Parameters
----------
structure : char *
    Will hold the barcket dot structure of the dimer molecule  

Returns
-------
float  
    minimum free energy of the structure  

See Also
--------
RNA.fold_compound.mfe()  
";

%feature("docstring") my_fold "

Compute Minimum Free Energy (MFE), and a corresponding secondary structure for an RNA sequence.  

This simplified interface to RNA.fold_compound.mfe() computes the MFE and, if required, a secondary structure for
an RNA sequence using default options. Memory required for dynamic programming (DP) matrices will be
allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
matrices are not available any more for any post-processing, e.g. suboptimal backtracking, etc.  

**SWIG Wrapper Notes**
    This function is available as function `fold()` in the global namespace. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:func:`RNA.fold()` in the :doc:`/api_python`.  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.circfold(), RNA.fold_compound.mfe()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.mfe(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") my_circfold "

Compute Minimum Free Energy (MFE), and a corresponding secondary structure for a circular RNA
sequence.  

This simplified interface to RNA.fold_compound.mfe() computes the MFE and, if required, a secondary structure for
a circular RNA sequence using default options. Memory required for dynamic programming (DP) matrices
will be allocated and free'd on-the-fly. Hence, after return of this function, the recursively
filled matrices are not available any more for any post-processing, e.g. suboptimal backtracking,
etc.  

Folding of circular RNA sequences is handled as a post-processing step of the forward recursions.
See   :cite:t:`hofacker:2006`  for further details.  

**SWIG Wrapper Notes**
    This function is available as function `circfold()` in the global namespace. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:func:`RNA.circfold()` in the :doc:`/api_python`.  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.fold(), RNA.fold_compound.mfe()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.mfe(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") my_alifold "

Compute Minimum Free Energy (MFE), and a corresponding consensus secondary structure for an RNA
sequence alignment using a comparative method.  

This simplified interface to RNA.fold_compound.mfe() computes the MFE and, if required, a consensus secondary
structure for an RNA sequence alignment using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing, e.g.
suboptimal backtracking, etc.  

**SWIG Wrapper Notes**
    This function is available as function `alifold()` in the global namespace. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:func:`RNA.alifold()` in the :doc:`/api_python`.  

Parameters
----------
sequences : const char **
    RNA sequence alignment  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.circalifold(), RNA.fold_compound.mfe()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.mfe(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") my_circalifold "

Compute Minimum Free Energy (MFE), and a corresponding consensus secondary structure for a sequence
alignment of circular RNAs using a comparative method.  

This simplified interface to RNA.fold_compound.mfe() computes the MFE and, if required, a consensus secondary
structure for an RNA sequence alignment using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing, e.g.
suboptimal backtracking, etc.  

Folding of circular RNA sequences is handled as a post-processing step of the forward recursions.
See   :cite:t:`hofacker:2006`  for further details.  

**SWIG Wrapper Notes**
    This function is available as function `circalifold()` in the global namespace. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:func:`RNA.circalifold()` in the :doc:`/api_python`.  

Parameters
----------
sequences : const char **
    Sequence alignment of circular RNAs  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.alifold(), RNA.fold_compound.mfe()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.mfe(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") my_cofold "

Compute Minimum Free Energy (MFE), and a corresponding secondary structure for two dimerized RNA
sequences.  

This simplified interface to RNA.fold_compound.mfe() computes the MFE and, if required, a secondary structure for
two RNA sequences upon dimerization using default options. Memory required for dynamic programming
(DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this function, the
recursively filled matrices are not available any more for any post-processing, e.g. suboptimal
backtracking, etc.  

.. deprecated:: 2.6.4
    This function is obsolete since RNA.mfe()/RNA.fold() can handle complexes multiple sequences
    since v2.5.0. Use RNA.mfe()/RNA.fold() for connected component MFE instead and compute MFEs of
    unconnected states separately.  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.mfe(),
and the data
structure RNA.fold_compound() instead.  

**SWIG Wrapper Notes**
    This function is available as function `cofold()` in the global namespace. The parameter
    `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:func:`RNA.cofold()` in the :doc:`/api_python`.  

Parameters
----------
sequence : const char *
    two RNA sequences separated by the '&' character  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.fold(), RNA.fold_compound.mfe(), RNA.fold_compound(), RNA.fold_compound(),
RNA.cut_point_insert()  
";

// File: group__mfe__window.xml

%feature("docstring") vrna_fold_compound_t::mfe_window "

Local MFE prediction using a sliding window approach.  

Computes minimum free energy structures using a sliding window approach, where base pairs may not
span outside the window. In contrast to RNA.fold_compound.mfe(), where a maximum base pair span may be set using
the RNA.md().max_bp_span attribute and one globally optimal structure is predicted, this function
uses a sliding window to retrieve all locally optimal structures within each window. The size of the
sliding window is set in the RNA.md().window_size attribute, prior to the retrieval of the
RNA.fold_compound() using RNA.fold_compound() with option RNA.OPTION_WINDOW  

The predicted structures are written on-the-fly, either to stdout, if a NULL pointer is passed as
file parameter, or to the corresponding filehandle.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `mfe_window()` to objects of type
    `fold_compound`. The parameter `FILE` has default value of `NULL` and can be omitted. See e.g.
    :py:meth:`RNA.fold_compound.mfe_window()` in the :doc:`/api_python`.  

Parameters
----------
file : FILE *
    The output file handle where predictions are written to (maybe NULL)  

See Also
--------
RNA.fold_compound(), RNA.fold_compound.mfe_window_zscore(), RNA.fold_compound.mfe(), RNA.Lfold(),
RNA.Lfoldz(),
RNA.OPTION_WINDOW, RNA.md().max_bp_span, RNA.md().window_size  
";

%feature("docstring") vrna_fold_compound_t::mfe_window_cb "

**SWIG Wrapper Notes**
    This function is attached as overloaded method `mfe_window_cb()` to objects of type
    `fold_compound`. The parameter `data` has default value of `NULL` and can be omitted. See e.g.
    :py:meth:`RNA.fold_compound.mfe_window_cb()` in the :doc:`/api_python`.  
";

%feature("docstring") vrna_fold_compound_t::mfe_window_zscore "

Local MFE prediction using a sliding window approach (with z-score cut-off)  

Computes minimum free energy structures using a sliding window approach, where base pairs may not
span outside the window. This function is the z-score version of RNA.fold_compound.mfe_window(), i.e. only
predictions above a certain z-score cut-off value are printed. As for RNA.fold_compound.mfe_window(), the size of
the sliding window is set in the RNA.md().window_size attribute, prior to the retrieval of the
RNA.fold_compound() using RNA.fold_compound() with option RNA.OPTION_WINDOW.  

The predicted structures are written on-the-fly, either to stdout, if a NULL pointer is passed as
file parameter, or to the corresponding filehandle.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `mfe_window_zscore()` to objects of type
    `fold_compound`. The parameter `FILE` has default value of `NULL` and can be omitted. See e.g.
    :py:meth:`RNA.fold_compound.mfe_window_zscore()` in the :doc:`/api_python`.  

Parameters
----------
min_z : double
    The minimal z-score for a predicted structure to appear in the output  
file : FILE *
    The output file handle where predictions are written to (maybe NULL)  

See Also
--------
RNA.fold_compound(), RNA.fold_compound.mfe_window_zscore(), RNA.fold_compound.mfe(), RNA.Lfold(),
RNA.Lfoldz(),
RNA.OPTION_WINDOW, RNA.md().max_bp_span, RNA.md().window_size  
";

%feature("docstring") vrna_mfe_window_zscore_cb "

**SWIG Wrapper Notes**
    This function is attached as overloaded method `mfe_window_zscore_cb()` to objects of type
    `fold_compound`. The parameter `data` has default value of `NULL` and can be omitted. See e.g.
    :py:meth:`RNA.fold_compound.mfe_window_zscore()` in the :doc:`/api_python`.  
";

%feature("docstring") my_Lfold "

Local MFE prediction using a sliding window approach (simplified interface)  

This simplified interface to RNA.fold_compound.mfe_window() computes the MFE and locally optimal secondary
structure using default options. Structures are predicted using a sliding window approach, where
base pairs may not span outside the window. Memory required for dynamic programming (DP) matrices
will be allocated and free'd on-the-fly. Hence, after return of this function, the recursively
filled matrices are not available any more for any post-processing.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `Lfold()` in the global namespace. The
    parameter `file` defaults to `NULL` and may be omitted. See e.g.  :py:func:`RNA.Lfold()` in the
    :doc:`/api_python`.  

Parameters
----------
string : const char *
    The nucleic acid sequence  
window_size : int
    The window size for locally optimal structures  
file : FILE *
    The output file handle where predictions are written to (if NULL, output is written to stdout)  

See Also
--------
RNA.fold_compound.mfe_window(), RNA.Lfoldz(), RNA.fold_compound.mfe_window_zscore()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use
RNA.fold_compound.mfe_window(), and the
data structure RNA.fold_compound() instead.  
";

%feature("docstring") vrna_Lfold_cb "

**SWIG Wrapper Notes**
    This function is available as overloaded function `Lfold_cb()` in the global namespace. The
    parameter `data` defaults to `NULL` and may be omitted. See e.g.  :py:func:`RNA.Lfold_cb()` in
    the :doc:`/api_python`.  
";

%feature("docstring") my_Lfoldz "

Local MFE prediction using a sliding window approach with z-score cut-off (simplified interface)  

This simplified interface to RNA.fold_compound.mfe_window_zscore() computes the MFE and locally optimal secondary
structure using default options. Structures are predicted using a sliding window approach, where
base pairs may not span outside the window. Memory required for dynamic programming (DP) matrices
will be allocated and free'd on-the-fly. Hence, after return of this function, the recursively
filled matrices are not available any more for any post-processing. This function is the z-score
version of RNA.Lfold(), i.e. only predictions above a certain z-score cut-off value are printed.  

Parameters
----------
string : const char *
    The nucleic acid sequence  
window_size : int
    The window size for locally optimal structures  
min_z : double
    The minimal z-score for a predicted structure to appear in the output  
file : FILE *
    The output file handle where predictions are written to (if NULL, output is written to stdout)  

See Also
--------
RNA.fold_compound.mfe_window_zscore(), RNA.Lfold(), RNA.fold_compound.mfe_window()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use
RNA.fold_compound.mfe_window(), and the
data structure RNA.fold_compound() instead.  
";

%feature("docstring") vrna_Lfoldz_cb "

**SWIG Wrapper Notes**
    This function is available as overloaded function `Lfoldz_cb()` in the global namespace. The
    parameter `data` defaults to `NULL` and may be omitted. See e.g.  :py:func:`RNA.Lfoldz_cb()` in
    the :doc:`/api_python`.  
";

%feature("docstring") my_aliLfold "

**SWIG Wrapper Notes**
    This function is available as overloaded function `aliLfold()` in the global namespace. The
    parameter `fp` defaults to `NULL` and may be omitted. See e.g.  :py:func:`RNA.aliLfold()` in the
    :doc:`/api_python`.  
";

%feature("docstring") vrna_aliLfold_cb "

**SWIG Wrapper Notes**
    This function is available as overloaded function `aliLfold_cb()` in the global namespace. The
    parameter `data` defaults to `NULL` and may be omitted. See e.g.  :py:func:`RNA.aliLfold_cb()`
    in the :doc:`/api_python`.  
";

// File: group__mfe__backtracking.xml

%feature("docstring") vrna_backtrack_from_intervals "
";

%feature("docstring") vrna_fold_compound_t::backtrack "

Backtrack an MFE (sub)structure.  

This function allows one to backtrack the MFE structure for a (sub)sequence  

**Precondition**
    Requires pre-filled MFE dynamic programming matrices, i.e. one has to call
RNA.fold_compound.mfe() prior to
    calling this function  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `backtrack()` to objects of type `fold_compound`.
    The parameter `length` defaults to the total length of the RNA sequence and may be omitted. The
    parameter `structure` is returned along with the MFE und must not be provided. See e.g.
    :py:meth:`RNA.fold_compound.backtrack()` in the :doc:`/api_python`.  

Parameters
----------
length : unsigned int
    The length of the subsequence, starting from the 5' end  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to. (Must have size of at least $p length + 1)  

Returns
-------
float  
    The minimum free energy (MFE) for the specified `length` in kcal/mol and a corresponding
    secondary structure in dot-bracket notation (stored in `structure`)  

See Also
--------
RNA.fold_compound.mfe(), RNA.fold_compound.pbacktrack5()  

Note
----
On error, the function returns INF / 100. and stores the empty string in `structure`.  
";

%feature("docstring") vrna_backtrack_window "
";

%feature("docstring") vrna_BT_ext_loop_f5 "
";

%feature("docstring") vrna_BT_ext_loop_f3 "
";

%feature("docstring") vrna_BT_ext_loop_f3_pp "
";

%feature("docstring") vrna_BT_hp_loop "

Backtrack a hairpin loop closed by :math:`(i,j)`.  

Note
----
This function is polymorphic! The provided RNA.fold_compound() may be of type RNA.FC_TYPE_SINGLE
or RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_BT_stack "

Backtrack a stacked pair closed by :math:`(i,j)`.  
";

%feature("docstring") vrna_BT_int_loop "

Backtrack an interior loop closed by :math:`(i,j)`.  
";

%feature("docstring") vrna_BT_mb_loop "

Backtrack the decomposition of a multi branch loop closed by :math:`(i,j)`.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() filled with all relevant data for backtracking  
i : int *
    5' position of base pair closing the loop (will be set to 5' position of leftmost decomposed
    block upon successful backtracking)  
j : int *
    3' position of base pair closing the loop (will be set to 3' position of rightmost decomposed
    block upon successful backtracking)  
k : int *
    Split position that delimits leftmost from rightmost block, [i,k] and [k+1, j], respectively.
    (Will be set upon successful backtracking)  
en : int
    The energy contribution of the substructure enclosed by :math:`(i,j)`  
component1 : int *
    Type of leftmost block (1 = ML, 2 = C)  
component2 : int *
    Type of rightmost block (1 = ML, 2 = C)  

Returns
-------
int  
    1, if backtracking succeeded, 0 otherwise.  
";

%feature("docstring") vrna_BT_mb_loop_split "
";

// File: group__part__func__global.xml

%feature("docstring") vrna_fold_compound_t::pf "

Compute the partition function :math:`Q` for a given RNA sequence, or sequence alignment.  

If *structure* is not a NULL pointer on input, it contains on return a string consisting of the
letters \" . , | { } ( ) \" denoting bases that are essentially unpaired, weakly paired, strongly
paired without preference, weakly upstream (downstream) paired, or strongly up- (down-)stream paired
bases, respectively. If the model's compute_bpp is set to 0 base pairing probabilities will not be
computed (saving CPU time), otherwise after calculations took place pr will contain the probability
that bases *i* and *j* pair.  

**SWIG Wrapper Notes**
    This function is attached as method `pf()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.pf()` in the :doc:`/api_python`.  

Parameters
----------
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  

Returns
-------
FLT_OR_DBL  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.fold_compound(), RNA.fold_compound(), RNA.pf_fold(), RNA.pf_circfold(),
RNA.fold_compound_comparative(), RNA.pf_alifold(), RNA.pf_circalifold(), RNA.db_from_probs(),
RNA.exp_params(), RNA.aln_pinfo()  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE. Also, this function may return INF / 100. in case of contradicting
constraints or numerical over-/underflow. In the latter case, a corresponding warning will be issued
to `stdout`.  
";

%feature("docstring") vrna_fold_compound_t::pf_dimer "

Calculate partition function and base pair probabilities of nucleic acid/nucleic acid dimers.  

This is the cofold partition function folding.  

**SWIG Wrapper Notes**
    This function is attached as method `pf_dimer()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.pf_dimer()` in the :doc:`/api_python`.  

Parameters
----------
structure : char *
    Will hold the structure or constraints  

Returns
-------
RNA.dimer_pf()  
    RNA.dimer_pf() structure containing a set of energies needed for concentration computations.  

See Also
--------
RNA.fold_compound() for how to retrieve the necessary data structure  

Note
----
This function may return INF / 100. for the `FA`, `FB`, `FAB`, `F0AB` members of the output data
structure in case of contradicting constraints or numerical over-/underflow. In the latter case, a
corresponding warning will be issued to `stdout`.  
";

%feature("docstring") vrna_pf_substrands "
";

%feature("docstring") my_pf_add "
";

%feature("docstring") vrna_pf_fold "

Compute Partition function :math:`Q` (and base pair probabilities) for an RNA sequence using a
comparative method.  

This simplified interface to RNA.fold_compound.pf() computes the partition function and, if required, base pair
probabilities for an RNA sequence using default options. Memory required for dynamic programming
(DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this function, the
recursively filled matrices are not available any more for any post-processing.  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.pf_circfold(), RNA.fold_compound.pf(), RNA.fold_compound(), RNA.fold_compound()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.pf(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") vrna_pf_circfold "

Compute Partition function :math:`Q` (and base pair probabilities) for a circular RNA sequences
using a comparative method.  

This simplified interface to RNA.fold_compound.pf() computes the partition function and, if required, base pair
probabilities for a circular RNA sequence using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing.  


Folding of circular RNA sequences is handled as a post-processing step of the forward recursions.
See   :cite:t:`hofacker:2006`  for further details.  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.pf(),
and the data
structure RNA.fold_compound() instead.  

Parameters
----------
sequence : const char *
    A circular RNA sequence  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.pf_fold(), RNA.fold_compound.pf(), RNA.fold_compound(), RNA.fold_compound()  
";

%feature("docstring") vrna_pf_alifold "

Compute Partition function :math:`Q` (and base pair probabilities) for an RNA sequence alignment
using a comparative method.  

This simplified interface to RNA.fold_compound.pf() computes the partition function and, if required, base pair
probabilities for an RNA sequence alignment using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing.  

Parameters
----------
sequences : const char **
    RNA sequence alignment  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.pf_circalifold(), RNA.fold_compound.pf(), RNA.fold_compound_comparative(), RNA.fold_compound()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.pf(),
and the data
structure RNA.fold_compound() instead.  
";

%feature("docstring") vrna_pf_circalifold "

Compute Partition function :math:`Q` (and base pair probabilities) for an alignment of circular RNA
sequences using a comparative method.  

This simplified interface to RNA.fold_compound.pf() computes the partition function and, if required, base pair
probabilities for an RNA sequence alignment using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing.  


Folding of circular RNA sequences is handled as a post-processing step of the forward recursions.
See   :cite:t:`hofacker:2006`  for further details.  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use RNA.fold_compound.pf(),
and the data
structure RNA.fold_compound() instead.  

Parameters
----------
sequences : const char **
    Sequence alignment of circular RNAs  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.pf_alifold(), RNA.fold_compound.pf(), RNA.fold_compound_comparative(), RNA.fold_compound()  
";

%feature("docstring") vrna_pf_co_fold "

Calculate partition function and base pair probabilities of nucleic acid/nucleic acid dimers.  

This simplified interface to RNA.fold_compound.pf_dimer() computes the partition function and, if required, base
pair probabilities for an RNA-RNA interaction using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing.  

Parameters
----------
seq : const char *
    Two concatenated RNA sequences with a delimiting '&' in between  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
RNA.dimer_pf()  
    RNA.dimer_pf() structure containing a set of energies needed for concentration computations.  

See Also
--------
RNA.fold_compound.pf_dimer()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use
RNA.fold_compound.pf_dimer(), and the
data structure RNA.fold_compound() instead.  
";

%feature("docstring") vrna_fold_compound_t::plist_from_probs "

Create a RNA.ep() from base pair probability matrix.  

The probability matrix provided via the RNA.fold_compound() is parsed and all pair probabilities
above the given threshold are used to create an entry in the plist  

The end of the plist is marked by sequence positions i as well as j equal to 0. This condition
should be used to stop looping over its entries  

Parameters
----------
cut_off : double
    The cutoff value  

Returns
-------
RNA.ep() *  
    A pointer to the plist that is to be created  
";

// File: group__part__func__window.xml

%feature("docstring") vrna_fold_compound_t::probs_window "

Compute various equilibrium probabilities under a sliding window approach.  

This function applies a sliding window scan for the sequence provided with the argument `fc` and
reports back equilibrium probabilities through the callback function `cb`. The data reported to the
callback depends on the `options` flag.  

#### Options:  
Note
----
The parameter `ulength` only affects computation and resulting data if unpaired probability
computations are requested through the `options` flag.  

*   RNA.PROBS_WINDOW_BPP -  Trigger base pairing probabilities.  
*   RNA.PROBS_WINDOW_UP -  Trigger unpaired probabilities.  
*   RNA.PROBS_WINDOW_UP_SPLIT -  Trigger detailed unpaired probabilities split up into different
    loop type contexts.  

Options may be OR-ed together  

Parameters
----------
ulength : int
    The maximal length of an unpaired segment (only for unpaired probability computations)  
cb : RNA.probs_window
    The callback function which collects the pair probability data for further processing  
data : void *
    Some arbitrary data structure that is passed to the callback `cb`  
options : unsigned int
    Option flags to control the behavior of this function  

Returns
-------
int  
    0 on failure, non-zero on success  

See Also
--------
RNA.pfl_fold_cb(), RNA.pfl_fold_up_cb()  
";

%feature("docstring") my_pfl_fold "

Compute base pair probabilities using a sliding-window approach.  

This is a simplified wrapper to RNA.fold_compound.probs_window() that given a nucleid acid sequence, a window
size, a maximum base pair span, and a cutoff value computes the pair probabilities for any base pair
in any window. The pair probabilities are returned as a list and the user has to take care to free()
the memory occupied by the list.  

Parameters
----------
sequence : const char *
    The nucleic acid input sequence  
window_size : int
    The size of the sliding window  
max_bp_span : int
    The maximum distance along the backbone between two nucleotides that form a base pairs  
cutoff : float
    A cutoff value that omits all pairs with lower probability  

Returns
-------
RNA.ep() *  
    A list of base pair probabilities, terminated by an entry with RNA.ep().i and RNA.ep().j set
    to 0  

See Also
--------
RNA.fold_compound.probs_window(), RNA.pfl_fold_cb(), RNA.pfl_fold_up()  

Note
----
This function uses default model settings! For custom model settings, we refer to the function
RNA.fold_compound.probs_window().  
 In case of any computation errors, this function returns `NULL`  
";

%feature("docstring") vrna_pfl_fold_cb "

Compute base pair probabilities using a sliding-window approach (callback version)  

This is a simplified wrapper to RNA.fold_compound.probs_window() that given a nucleid acid sequence, a window
size, a maximum base pair span, and a cutoff value computes the pair probabilities for any base pair
in any window. It is similar to RNA.pfl_fold() but uses a callback mechanism to return the pair
probabilities.  

Read the details for RNA.fold_compound.probs_window() for details on the callback implementation!  

Parameters
----------
sequence : const char *
    The nucleic acid input sequence  
window_size : int
    The size of the sliding window  
max_bp_span : int
    The maximum distance along the backbone between two nucleotides that form a base pairs  
cb : RNA.probs_window
    The callback function which collects the pair probability data for further processing  
data : void *
    Some arbitrary data structure that is passed to the callback `cb`  

Returns
-------
int  
    0 on failure, non-zero on success  

See Also
--------
RNA.fold_compound.probs_window(), RNA.pfl_fold(), RNA.pfl_fold_up_cb()  

Note
----
This function uses default model settings! For custom model settings, we refer to the function
RNA.fold_compound.probs_window().  
";

%feature("docstring") pfl_fold_up "

Compute probability of contiguous unpaired segments.  

This is a simplified wrapper to RNA.fold_compound.probs_window() that given a nucleic acid sequence, a maximum
length of unpaired segments (`ulength`), a window size, and a maximum base pair span computes the
equilibrium probability of any segment not exceeding `ulength`. The probabilities to be unpaired are
returned as a 1-based, 2-dimensional matrix with dimensions :math:`N \\times M`, where :math:`N` is
the length of the sequence and :math:`M` is the maximum segment length. As an example, the
probability of a segment of size 5 starting at position 100 is stored in the matrix entry
:math:`X[100][5]`.  

It is the users responsibility to free the memory occupied by this matrix.  

Parameters
----------
sequence : const char *
    The nucleic acid input sequence  
ulength : int
    The maximal length of an unpaired segment  
window_size : int
    The size of the sliding window  
max_bp_span : int
    The maximum distance along the backbone between two nucleotides that form a base pairs  

Returns
-------
double **  
    The probabilities to be unpaired for any segment not exceeding `ulength`  

Note
----
This function uses default model settings! For custom model settings, we refer to the function
RNA.fold_compound.probs_window().  
";

%feature("docstring") vrna_pfl_fold_up_cb "

Compute probability of contiguous unpaired segments.  

This is a simplified wrapper to RNA.fold_compound.probs_window() that given a nucleic acid sequence, a maximum
length of unpaired segments (`ulength`), a window size, and a maximum base pair span computes the
equilibrium probability of any segment not exceeding `ulength`. It is similar to RNA.pfl_fold_up()
but uses a callback mechanism to return the unpaired probabilities.  

Read the details for RNA.fold_compound.probs_window() for details on the callback implementation!  

Parameters
----------
sequence : const char *
    The nucleic acid input sequence  
ulength : int
    The maximal length of an unpaired segment  
window_size : int
    The size of the sliding window  
max_bp_span : int
    The maximum distance along the backbone between two nucleotides that form a base pairs  
cb : RNA.probs_window
    The callback function which collects the pair probability data for further processing  
data : void *
    Some arbitrary data structure that is passed to the callback `cb`  

Returns
-------
int  
    0 on failure, non-zero on success  

Note
----
This function uses default model settings! For custom model settings, we refer to the function
RNA.fold_compound.probs_window().  
";

%feature("docstring") EXT_LOOP "

Exterior loop.  
";

%feature("docstring") HP_LOOP "

Hairpin loop.  
";

%feature("docstring") INT_LOOP "

Internal loop.  
";

%feature("docstring") MB_LOOP "

Multibranch loop.  
";

%feature("docstring") ANY_LOOP "

Any loop.  
";

%feature("docstring") PROBS_WINDOW_BPP "

Trigger base pairing probabilities.  

Passing this flag to RNA.fold_compound.probs_window() activates callback execution for base pairing
probabilities. In turn, the corresponding callback receives this flag through the `type` argument
whenever base pairing probabilities are provided.  

Detailed information for the algorithm to compute unpaired probabilities can be taken from
:cite:t:`bernhart:2005` .  

See Also
--------
RNA.fold_compound.probs_window()  
";

%feature("docstring") PROBS_WINDOW_UP "

Trigger unpaired probabilities.  

Passing this flag to RNA.fold_compound.probs_window() activates callback execution for unpaired probabilities. In
turn, the corresponding callback receives this flag through the `type` argument whenever unpaired
probabilities are provided.  

Detailed information for the algorithm to compute unpaired probabilities can be taken from
:cite:t:`bernhart:2011` .  

See Also
--------
RNA.fold_compound.probs_window()  
";

%feature("docstring") PROBS_WINDOW_STACKP "

Trigger base pair stack probabilities.  

Passing this flag to RNA.fold_compound.probs_window() activates callback execution for stacking probabilities. In
turn, the corresponding callback receives this flag through the `type` argument whenever stack
probabilities are provided.  

**Bug**
    Currently, this flag is a placeholder doing nothing as the corresponding implementation for
    stack probability computation is missing.  

See Also
--------
RNA.fold_compound.probs_window()  
";

%feature("docstring") PROBS_WINDOW_UP_SPLIT "

Trigger detailed unpaired probabilities split up into different loop type contexts.  

Passing this flag to RNA.fold_compound.probs_window() activates callback execution for unpaired probabilities. In
contrast to RNA.PROBS_WINDOW_UP this flag requests unpaired probabilities to be split up into
different loop type contexts. In turn, the corresponding callback receives the RNA.PROBS_WINDOW_UP
flag OR-ed together with the corresponding loop type, i.e.:  

*   RNA.EXT_LOOP -  Exterior loop.  
*   RNA.HP_LOOP -  Hairpin loop.  
*   RNA.INT_LOOP -  Internal loop.  
*   RNA.MB_LOOP -  Multibranch loop.  
*   RNA.ANY_LOOP -  Any loop.  

See Also
--------
RNA.fold_compound.probs_window(), RNA.PROBS_WINDOW_UP  
";

%feature("docstring") PROBS_WINDOW_PF "

Trigger partition function.  

Passing this flag to RNA.fold_compound.probs_window() activates callback execution for partition function. In
turn, the corresponding callback receives this flag through it's `type` argument whenever partition
function data is provided.  

Note
----
Instead of actually providing the partition function :math:`Z`, the callback is always provided with
the corresponding enemble free energy :math:`\\Delta G = - RT \\ln Z`.  

See Also
--------
RNA.fold_compound.probs_window()  
";

// File: group__subopt__and__representatives.xml

// File: group__subopt__zuker.xml

%feature("docstring") zukersubopt "

Compute Zuker type suboptimal structures.  

Compute Suboptimal structures according to M. Zuker, i.e. for every possible base pair the minimum
energy structure containing the resp. base pair. Returns a list of these structures and their
energies.  

.. deprecated:: 2.6.4
    use RNA.zukersubopt() instead  

Parameters
----------
string : const char *
    RNA sequence  

Returns
-------
SOLUTION *  
    List of zuker suboptimal structures  
";

%feature("docstring") zukersubopt_par "

Compute Zuker type suboptimal structures.  

.. deprecated:: 2.6.4
    use RNA.zukersubopt() instead  
";

%feature("docstring") vrna_fold_compound_t::subopt_zuker "

Compute Zuker type suboptimal structures.  

Compute Suboptimal structures according to   :cite:t:`zuker:1989`  , i.e. for every possible base
pair the minimum energy structure containing the resp. base pair. Returns a list of these structures
and their energies.  

**SWIG Wrapper Notes**
    This function is attached as method **subopt_zuker()** to objects of type `fold_compound`. See,
    e.g.  :py:meth:`RNA.fold_compound.subopt_zuker()` in the :doc:`/api_python`.  

Parameters
----------

Returns
-------
RNA.subopt_solution() *  
    List of zuker suboptimal structures  

See Also
--------
RNA.fold_compound.subopt(), zukersubopt(), zukersubopt_par()  
";

// File: group__subopt__wuchty.xml

%feature("docstring") vrna_fold_compound_t::subopt "

Returns list of subopt structures or writes to fp.  

This function produces **all** suboptimal secondary structures within 'delta' * 0.01 kcal/mol of the
optimum, see   :cite:t:`wuchty:1999` . The results are either directly written to a 'fp' (if 'fp' is
not NULL), or (fp==NULL) returned in a RNA.subopt_solution() * list terminated by an entry were the
'structure' member is NULL.  

**SWIG Wrapper Notes**
    This function is attached as method **subopt()** to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.subopt()` in the :doc:`/api_python`.  

Parameters
----------
delta : int
sorted : int
    Sort results by energy in ascending order  
fp : FILE *

Returns
-------
RNA.subopt_solution() *  

See Also
--------
RNA.fold_compound.subopt_cb(), RNA.fold_compound.subopt_zuker()  

Note
----
This function requires all multibranch loop DP matrices for unique multibranch loop backtracing.
Therefore, the supplied RNA.fold_compound()`fc` (argument 1) must be initialized with
RNA.md().uniq_ML = 1, for instance like this:  
";

%feature("docstring") vrna_fold_compound_t::subopt_cb "

Generate suboptimal structures within an energy band arround the MFE.  

This is the most generic implementation of the suboptimal structure generator according to
:cite:t:`wuchty:1999` . Identical to RNA.fold_compound.subopt(), it computes all secondary structures within an
energy band `delta` arround the MFE. However, this function does not print the resulting structures
and their corresponding free energies to a file pointer, or returns them as a list. Instead, it
calls a user-provided callback function which it passes the structure in dot-bracket format, the
corresponding free energy in kcal/mol, and a user-provided data structure each time a structure was
backtracked successfully. This function indicates the final output, i.e. the end of the backtracking
procedure by passing NULL instead of an actual dot-bracket string to the callback.  

**SWIG Wrapper Notes**
    This function is attached as method **subopt_cb()** to objects of type `fold_compound`. See,
    e.g.  :py:meth:`RNA.fold_compound.subopt_cb()` in the :doc:`/api_python`.  

Parameters
----------
delta : int
    Energy band arround the MFE in 10cal/mol, i.e. deka-calories  
cb : RNA.subopt_result
    Pointer to a callback function that handles the backtracked structure and its free energy in
    kcal/mol  
data : void *
    Pointer to some data structure that is passed along to the callback  

See Also
--------
RNA.subopt_result, RNA.fold_compound.subopt(), RNA.fold_compound.subopt_zuker()  

Note
----
This function requires all multibranch loop DP matrices for unique multibranch loop backtracing.
Therefore, the supplied RNA.fold_compound()`fc` (argument 1) must be initialized with
RNA.md().uniq_ML = 1, for instance like this:  
";

%feature("docstring") subopt "

Returns list of subopt structures or writes to fp.  

This function produces **all** suboptimal secondary structures within 'delta' * 0.01 kcal/mol of the
optimum. The results are either directly written to a 'fp' (if 'fp' is not NULL), or (fp==NULL)
returned in a SOLUTION * list terminated by an entry were the 'structure' pointer is NULL.  

Parameters
----------
seq : char *
structure : char *
delta : int
fp : FILE *

Returns
-------
SOLUTION *  
";

%feature("docstring") subopt_par "

Returns list of subopt structures or writes to fp.  
";

%feature("docstring") subopt_circ "

Returns list of circular subopt structures or writes to fp.  

This function is similar to subopt() but calculates secondary structures assuming the RNA sequence
to be circular instead of linear  

Parameters
----------
seq : char *
sequence : char *
delta : int
fp : FILE *

Returns
-------
SOLUTION *  
";

// File: group__subopt__stochbt.xml

%feature("docstring") vrna_fold_compound_t::pbacktrack5 "

Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a
secondary structure. The parameter `length` specifies the length of the substructure starting from
the 5' end.  

The structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack5()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack5()` in the
    :doc:`/api_python`  and the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
length : unsigned int
    The length of the subsequence to consider (starting with 5' end)  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation (or NULL on error)  

See Also
--------
RNA.pbacktrack5_num(), RNA.pbacktrack5_cb(), RNA.fold_compound.pbacktrack()  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack5_num "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack5()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack5()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
length : unsigned int
    The length of the subsequence to consider (starting with 5' end)  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack5(), RNA.pbacktrack5_cb(), RNA.pbacktrack_num(), RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack5_cb "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.fold_compound.pbacktrack5() and RNA.pbacktrack5_num() this function yields the structure
samples through a callback mechanism.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack5()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack5()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
length : unsigned int
    The length of the subsequence to consider (starting with 5' end)  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack5(), RNA.pbacktrack5_num(), RNA.pbacktrack_cb(), RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack5_resume "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack5_cb() this function allows for resuming a previous sampling round in
specialized Boltzmann sampling, such as non-redundant backtracking. For that purpose, the user
passes the address of a Boltzmann sampling data structure (RNA.pbacktrack_mem()) which will be re-
used in each round of sampling, i.e. each successive call to RNA.pbacktrack5_resume_cb() or
RNA.pbacktrack5_resume().  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack5()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to
    the list of structures, this function also returns the `nr_mem` data structure as first return
    value. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack5()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
length : unsigned int
    The length of the subsequence to consider (starting with 5' end)  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack5_resume_cb(), RNA.pbacktrack5_cb(), RNA.pbacktrack_resume(),
RNA.pbacktrack_mem(), RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT,
RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack5_resume_cb "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack5_resume() this function yields the structure samples through a
callback mechanism.  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack5()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to
    the number of structures backtraced, this function also returns the `nr_mem` data structure as
    first return value. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack5()` in the
    :doc:`/api_python`  and the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
length : unsigned int
    The length of the subsequence to consider (starting with 5' end)  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack5_resume(), RNA.pbacktrack5_cb(), RNA.pbacktrack_resume_cb(),
RNA.pbacktrack_mem(), RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT,
RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_fold_compound_t::pbacktrack "

Sample a secondary structure from the Boltzmann ensemble according its probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a
secondary structure.  

The structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`
    and the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation (or NULL on error)  

See Also
--------
RNA.fold_compound.pbacktrack5(), RNA.pbacktrack_num, RNA.pbacktrack_cb()  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_num "

Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type `fold_compound`
    with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack(), RNA.pbacktrack_cb(), RNA.pbacktrack5_num(), RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_cb "

Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.fold_compound.pbacktrack() and RNA.pbacktrack_num() this function yields the structure
samples through a callback mechanism.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type `fold_compound`
    with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack(), RNA.pbacktrack_num(), RNA.pbacktrack5_cb(), RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_resume "

Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack_cb() this function allows for resuming a previous sampling round in
specialized Boltzmann sampling, such as non-redundant backtracking. For that purpose, the user
passes the address of a Boltzmann sampling data structure (RNA.pbacktrack_mem()) which will be re-
used in each round of sampling, i.e. each successive call to RNA.pbacktrack_resume_cb() or
RNA.pbacktrack_resume().  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type `fold_compound`
    with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to the list of
    structures, this function also returns the `nr_mem` data structure as first return value. See,
    e.g.   :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack_resume_cb(), RNA.pbacktrack_cb(), RNA.pbacktrack5_resume(), RNA.pbacktrack_mem(),
RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT, RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_resume_cb "

Obtain a set of secondary structure samples from the Boltzmann ensemble according their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack5_resume() this function yields the structure samples through a
callback mechanism.  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type `fold_compound`
    with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to the number of
    structures backtraced, this function also returns the `nr_mem` data structure as first return
    value. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack_resume(), RNA.pbacktrack_cb(), RNA.pbacktrack5_resume_cb(), RNA.pbacktrack_mem(),
RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT, RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_fold_compound_t::pbacktrack_sub "

Sample a secondary structure of a subsequence from the Boltzmann ensemble according its probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a
secondary structure. The parameters `start` and `end` specify the interval :math:`[start:end]` of
the subsequence with :math:`1 \\leq start < end \\leq n` for sequence length :math:`n`, the
structure :math:`s_{start,end}` should be drawn from.  

The resulting substructure :math:`s_{start,end}` with free energy :math:`E(s_{start, end})` is
picked from the Boltzmann distributed sub ensemble of all structures within the interval
:math:`[start:end]` according to its probability  

.. math::

  p(s_{start,end}) = \\frac{exp(-E(s_{start,end}) / kT)}{Z_{start,end}}  

with partition function :math:`Z_{start,end} = \\sum_{s_{start,end}} exp(-E(s_{start,end}) / kT)`,
Boltzmann constant :math:`k` and thermodynamic temperature :math:`T`.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack_sub()` to objects of type
    *fold_compound*. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack_sub()` in the
    :doc:`/api_python`  and the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
start : unsigned int
    The start of the subsequence to consider, i.e. 5'-end position(1-based)  
end : unsigned int
    The end of the subsequence to consider, i.e. 3'-end position (1-based)  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation (or NULL on error)  

See Also
--------
RNA.pbacktrack_sub_num(), RNA.pbacktrack_sub_cb(), RNA.fold_compound.pbacktrack()  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_sub_num "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack_sub()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack_sub()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
start : unsigned int
    The start of the subsequence to consider, i.e. 5'-end position(1-based)  
end : unsigned int
    The end of the subsequence to consider, i.e. 3'-end position (1-based)  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack_sub(), RNA.pbacktrack_sub_cb(), RNA.pbacktrack_num(),
RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_sub_cb "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.fold_compound.pbacktrack5() and RNA.pbacktrack5_num() this function yields the structure
samples through a callback mechanism.  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack()` to objects of type `fold_compound`
    with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.pbacktrack()` in the :doc:`/api_python`  and the
    :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
start : unsigned int
    The start of the subsequence to consider, i.e. 5'-end position(1-based)  
end : unsigned int
    The end of the subsequence to consider, i.e. 3'-end position (1-based)  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.fold_compound.pbacktrack5(), RNA.pbacktrack5_num(), RNA.pbacktrack_cb(), RNA.PBACKTRACK_DEFAULT,
RNA.PBACKTRACK_NON_REDUNDANT  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_sub_resume "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack5_cb() this function allows for resuming a previous sampling round in
specialized Boltzmann sampling, such as non-redundant backtracking. For that purpose, the user
passes the address of a Boltzmann sampling data structure (RNA.pbacktrack_mem()) which will be re-
used in each round of sampling, i.e. each successive call to RNA.pbacktrack5_resume_cb() or
RNA.pbacktrack5_resume().  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack_sub()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to
    the list of structures, this function also returns the `nr_mem` data structure as first return
    value. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack_sub()` in the :doc:`/api_python`  and
    the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
start : unsigned int
    The start of the subsequence to consider, i.e. 5'-end position(1-based)  
end : unsigned int
    The end of the subsequence to consider, i.e. 3'-end position (1-based)  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
char **  
    A set of secondary structure samples in dot-bracket notation terminated by NULL (or NULL on
    error)  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack5_resume_cb(), RNA.pbacktrack5_cb(), RNA.pbacktrack_resume(),
RNA.pbacktrack_mem(), RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT,
RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_sub_resume_cb "

Obtain a set of secondary structure samples for a subsequence from the Boltzmann ensemble according
their probability.  

Perform a probabilistic (stochastic) backtracing in the partition function DP arrays to obtain a set
of `num_samples` secondary structures. The parameter `length` specifies the length of the
substructure starting from the 5' end.  

Any structure :math:`s` with free energy :math:`E(s)` is picked from the Boltzmann distributed
ensemble according to its probability  

.. math::

  p(s) = \\frac{exp(-E(s) / kT)}{Z}  

with partition function :math:`Z = \\sum_{s} exp(-E(s) / kT)`, Boltzmann constant :math:`k` and
thermodynamic temperature :math:`T`.  

Using the `options` flag one can switch between regular (RNA.PBACKTRACK_DEFAULT) backtracing mode,
and non-redundant sampling (RNA.PBACKTRACK_NON_REDUNDANT) along the lines of
:cite:t:`michalik:2017` .  

In contrast to RNA.pbacktrack5_resume() this function yields the structure samples through a
callback mechanism.  

A successive sample call to this function may look like:  

**Precondition**
    Unique multiloop decomposition has to be active upon creation of `fc` with RNA.fold_compound()
    or similar. This can be done easily by passing RNA.fold_compound() a model details parameter
    with RNA.md().uniq_ML = 1.  RNA.fold_compound.pf() has to be called first to fill the partition
function
    matrices  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `pbacktrack_sub()` to objects of type
    `fold_compound` with optional last argument `options` = RNA.PBACKTRACK_DEFAULT. In addition to
    the number of structures backtraced, this function also returns the `nr_mem` data structure as
    first return value. See, e.g.   :py:meth:`RNA.fold_compound.pbacktrack_sub()` in the
    :doc:`/api_python`  and the   :ref:`examples/python:boltzmann sampling` Python examples .  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure  
num_samples : unsigned int
    The size of the sample set, i.e. number of structures  
start : unsigned int
    The start of the subsequence to consider, i.e. 5'-end position(1-based)  
end : unsigned int
    The end of the subsequence to consider, i.e. 3'-end position (1-based)  
cb : RNA.bs_result
    The callback that receives the sampled structure  
data : void *
    A data structure passed through to the callback `cb`  
nr_mem : RNA.pbacktrack_mem() *
    The address of the Boltzmann sampling memory data structure  
options : unsigned int
    A bitwise OR-flag indicating the backtracing mode.  

Returns
-------
unsigned int  
    The number of structures actually backtraced  

Warnings
--------
In non-redundant sampling mode (RNA.PBACKTRACK_NON_REDUNDANT), this function may not yield the full
number of requested samples. This may happen if a) the number of requested structures is larger than
the total number of structuresin the ensemble, b) numeric instabilities prevent the backtracking
function to enumerate structures with high free energies, or c) any other error occurs.  

See Also
--------
RNA.pbacktrack5_resume(), RNA.pbacktrack5_cb(), RNA.pbacktrack_resume_cb(),
RNA.pbacktrack_mem(), RNA.PBACKTRACK_DEFAULT, RNA.PBACKTRACK_NON_REDUNDANT,
RNA.pbacktrack_mem_free  

Note
----
This function is polymorphic. It accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE, and
RNA.FC_TYPE_COMPARATIVE.  
";

%feature("docstring") vrna_pbacktrack_mem_free "

Release memory occupied by a Boltzmann sampling memory data structure.  

Parameters
----------
s : RNA.pbacktrack_mem()
    The non-redundancy memory data structure  

See Also
--------
RNA.pbacktrack_mem(), RNA.pbacktrack5_resume(), RNA.pbacktrack5_resume_cb(),
RNA.pbacktrack_resume(), RNA.pbacktrack_resume_cb()  
";

%feature("docstring") PBACKTRACK_DEFAULT "

Boltzmann sampling flag indicating default backtracing mode.  

See Also
--------
RNA.pbacktrack5_num(), RNA.pbacktrack5_cb(), RNA.pbacktrack5_resume(),
RNA.pbacktrack5_resume_cb(), RNA.pbacktrack_num(), RNA.pbacktrack_cb(), RNA.pbacktrack_resume(),
RNA.pbacktrack_resume_cb()  
";

%feature("docstring") PBACKTRACK_NON_REDUNDANT "

Boltzmann sampling flag indicating non-redundant backtracing mode.  

This flag will turn the Boltzmann sampling into non-redundant backtracing mode along the lines of
:cite:t:`michalik:2017`  

See Also
--------
RNA.pbacktrack5_num(), RNA.pbacktrack5_cb(), RNA.pbacktrack5_resume(),
RNA.pbacktrack5_resume_cb(), RNA.pbacktrack_num(), RNA.pbacktrack_cb(), RNA.pbacktrack_resume(),
RNA.pbacktrack_resume_cb()  
";

// File: group__mea__fold.xml

%feature("docstring") vrna_fold_compound_t::MEA "

Compute a MEA (maximum expected accuracy) structure.  

The algorithm maximizes the expected accuracy  

.. math::

  A(S) = \\sum_{(i,j) \\in S} 2 \\gamma p_{ij} + \\sum_{i \\notin S} p^u_{i}  

Higher values of :math:`\\gamma` result in more base pairs of lower probability and thus higher
sensitivity. Low values of :math:`\\gamma` result in structures containing only highly likely pairs
(high specificity). The code of the MEA function also demonstrates the use of sparse dynamic
programming scheme to reduce the time and memory complexity of folding.  

**Precondition**
    RNA.fold_compound.pf() must be executed on input parameter `fc`  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `MEA`(gamma = 1.) to objects of type
    `fold_compound`. Note, that it returns the MEA structure and MEA value as a tuple
    (MEA_structure, MEA). See, e.g.  :py:meth:`RNA.fold_compound.MEA()` in the :doc:`/api_python`.  

Parameters
----------
gamma : double
    The weighting factor for base pairs vs. unpaired nucleotides  
mea : float *
    A pointer to a variable where the MEA value will be written to  

Returns
-------
char *  
    An MEA structure (or NULL on any error)  
";

%feature("docstring") my_MEA_from_plist "

Compute a MEA (maximum expected accuracy) structure from a list of probabilities.  

The algorithm maximizes the expected accuracy  

.. math::

  A(S) = \\sum_{(i,j) \\in S} 2 \\gamma p_{ij} + \\sum_{i \\notin S} p^u_{i}  

Higher values of :math:`\\gamma` result in more base pairs of lower probability and thus higher
sensitivity. Low values of :math:`\\gamma` result in structures containing only highly likely pairs
(high specificity). The code of the MEA function also demonstrates the use of sparse dynamic
programming scheme to reduce the time and memory complexity of folding.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `MEA_from_plist`(gamma = 1., md = NULL). Note,
    that it returns the MEA structure and MEA value as a tuple (MEA_structure, MEA). See, e.g.
    :py:func:`RNA.MEA_from_plist()` in the :doc:`/api_python`.  

Parameters
----------
plist : RNA.ep() *
    A list of base pair probabilities the MEA structure is computed from  
sequence : const char *
    The RNA sequence that corresponds to the list of probability values  
gamma : double
    The weighting factor for base pairs vs. unpaired nucleotides  
md : RNA.md() *
    A model details data structure (maybe NULL)  
mea : float *
    A pointer to a variable where the MEA value will be written to  

Returns
-------
char *  
    An MEA structure (or NULL on any error)  

Note
----
The unpaired probabilities :math:`p^u_{i} = 1 - \\sum_{j \\neq i} p_{ij}` are usually computed from
the supplied pairing probabilities :math:`p_{ij}` as stored in `plist` entries of type
RNA.PLIST_TYPE_BASEPAIR. To overwrite individual :math:`p^u_{o}` values simply add entries with
type RNA.PLIST_TYPE_UNPAIRED  
 To include G-Quadruplex support, the corresponding field in `md` must be set.  
";

%feature("docstring") MEA "

Computes a MEA (maximum expected accuracy) structure.  

The algorithm maximizes the expected accuracy  

.. math::

  A(S) = \\sum_{(i,j) \\in S} 2 \\gamma p_{ij} + \\sum_{i \\notin S} p^u_{i}  

Higher values of :math:`\\gamma` result in more base pairs of lower probability and thus higher
sensitivity. Low values of :math:`\\gamma` result in structures containing only highly likely pairs
(high specificity). The code of the MEA function also demonstrates the use of sparse dynamic
programming scheme to reduce the time and memory complexity of folding.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.MEA() or RNA.MEA_from_plist() instead!  
";

// File: group__centroid__fold.xml

%feature("docstring") vrna_fold_compound_t::centroid "

Get the centroid structure of the ensemble.  

The centroid is the structure with the minimal average distance to all other structures
:math:`<d(S)> = \\sum_{(i,j) \\in S} (1-p_{ij}) + \\sum_{(i,j) \\notin S} p_{ij}`  
Thus, the centroid is simply the structure containing all pairs with :math:`p_{i}j>0.5` The distance
of the centroid to the ensemble is written to the memory adressed by *dist*.  

Parameters
----------
dist : double *
    A pointer to the distance variable where the centroid distance will be written to  

Returns
-------
char *  
    The centroid structure of the ensemble in dot-bracket notation (`NULL` on error)  
";

%feature("docstring") vrna_centroid_from_plist "

Get the centroid structure of the ensemble.  

This function is a threadsafe replacement for centroid() with a RNA.ep() input  

The centroid is the structure with the minimal average distance to all other structures
:math:`<d(S)> = \\sum_{(i,j) \\in S} (1-p_{ij}) + \\sum_{(i,j) \\notin S} p_{ij}`  
Thus, the centroid is simply the structure containing all pairs with :math:`p_{i}j>0.5` The distance
of the centroid to the ensemble is written to the memory adressed by *dist*.  

Parameters
----------
length : int
    The length of the sequence  
dist : double *
    A pointer to the distance variable where the centroid distance will be written to  
pl : RNA.ep() *
    A pair list containing base pair probability information about the ensemble  

Returns
-------
char *  
    The centroid structure of the ensemble in dot-bracket notation (`NULL` on error)  
";

%feature("docstring") vrna_centroid_from_probs "

Get the centroid structure of the ensemble.  

This function is a threadsafe replacement for centroid() with a probability array input  

The centroid is the structure with the minimal average distance to all other structures
:math:`<d(S)> = \\sum_{(i,j) \\in S} (1-p_{ij}) + \\sum_{(i,j) \\notin S} p_{ij}`  
Thus, the centroid is simply the structure containing all pairs with :math:`p_{i}j>0.5` The distance
of the centroid to the ensemble is written to the memory adressed by *dist*.  

Parameters
----------
length : int
    The length of the sequence  
dist : double *
    A pointer to the distance variable where the centroid distance will be written to  
probs : FLT_OR_DBL *
    An upper triangular matrix containing base pair probabilities (access via iindx
    RNA.idx_row_wise() )  

Returns
-------
char *  
    The centroid structure of the ensemble in dot-bracket notation (`NULL` on error)  
";

// File: group__cofold.xml

// File: group__class__fold.xml

// File: group__kl__neighborhood.xml

// File: group__kl__neighborhood__mfe.xml

%feature("docstring") vrna_mfe_TwoD "

Compute MFE's and representative for distance partitioning.  

This function computes the minimum free energies and a representative secondary structure for each
distance class according to the two references specified in the datastructure 'vars'. The maximum
basepair distance to each of both references may be set by the arguments 'distance1' and
'distance2', respectively. If both distance arguments are set to '-1', no restriction is assumed and
the calculation is performed for each distance class possible.  

The returned list contains an entry for each distance class. If a maximum basepair distance to
either of the references was passed, an entry with k=l=-1 will be appended in the list, denoting the
class where all structures exceeding the maximum will be thrown into. The end of the list is denoted
by an attribute value of INF in the k-attribute of the list entry.  

Parameters
----------
fc : RNA.fold_compound() *
    The datastructure containing all precomputed folding attributes  
distance1 : int
    maximum distance to reference1 (-1 means no restriction)  
distance2 : int
    maximum distance to reference2 (-1 means no restriction)  

Returns
-------
RNA.sol_TwoD() *  
    A list of minimum free energies (and corresponding structures) for each distance class  

See Also
--------
RNA.fold_compound_TwoD(), RNA.fold_compound_free(), RNA.pf_TwoD()RNA.backtrack5_TwoD(),
RNA.sol_TwoD(), RNA.fold_compound()  
";

%feature("docstring") vrna_backtrack5_TwoD "

Backtrack a minimum free energy structure from a 5' section of specified length.  

This function allows one to backtrack a secondary structure beginning at the 5' end, a specified
length and residing in a specific distance class. If the argument 'k' gets a value of -1, the
structure that is backtracked is assumed to reside in the distance class where all structures
exceeding the maximum basepair distance specified in RNA.mfe_TwoD() belong to.  

Parameters
----------
fc : RNA.fold_compound() *
    The datastructure containing all precomputed folding attributes  
j : unsigned int
    The length in nucleotides beginning from the 5' end  
k : int
    distance to reference1 (may be -1)  
l : int
    distance to reference2  

See Also
--------
RNA.mfe_TwoD()  

Note
----
The argument 'vars' must contain precalculated energy values in the energy matrices, i.e. a call to
RNA.mfe_TwoD() preceding this function is mandatory!  
";

%feature("docstring") get_TwoDfold_variables "

Get a structure of type TwoDfold_vars prefilled with current global settings.  

This function returns a datastructure of type TwoDfold_vars. The data fields inside the
TwoDfold_vars are prefilled by global settings and all memory allocations necessary to start a
computation are already done for the convenience of the user  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.mfe_TwoD(), and RNA.fold_compound_free() instead!  

Note
----
Make sure that the reference structures are compatible with the sequence according to Watson-Crick-
and Wobble-base pairing  

Parameters
----------
seq : const char *
    The RNA sequence  
structure1 : const char *
    The first reference structure in dot-bracket notation  
structure2 : const char *
    The second reference structure in dot-bracket notation  
circ : int
    A switch to indicate the assumption to fold a circular instead of linear RNA (0=OFF, 1=ON)  

Returns
-------
TwoDfold_vars *  
    A datastructure prefilled with folding options and allocated memory  
";

%feature("docstring") destroy_TwoDfold_variables "

Destroy a TwoDfold_vars datastructure without memory loss.  

This function free's all allocated memory that depends on the datastructure given.  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.mfe_TwoD(), and RNA.fold_compound_free() instead!  

Parameters
----------
our_variables : TwoDfold_vars *
    A pointer to the datastructure to be destroyed  
";

%feature("docstring") TwoDfoldList "

Compute MFE's and representative for distance partitioning.  

This function computes the minimum free energies and a representative secondary structure for each
distance class according to the two references specified in the datastructure 'vars'. The maximum
basepair distance to each of both references may be set by the arguments 'distance1' and
'distance2', respectively. If both distance arguments are set to '-1', no restriction is assumed and
the calculation is performed for each distance class possible.  

The returned list contains an entry for each distance class. If a maximum basepair distance to
either of the references was passed, an entry with k=l=-1 will be appended in the list, denoting the
class where all structures exceeding the maximum will be thrown into. The end of the list is denoted
by an attribute value of INF in the k-attribute of the list entry.  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.mfe_TwoD(), and RNA.fold_compound_free() instead!  

Parameters
----------
vars : TwoDfold_vars *
    the datastructure containing all predefined folding attributes  
distance1 : int
    maximum distance to reference1 (-1 means no restriction)  
distance2 : int
    maximum distance to reference2 (-1 means no restriction)  
";

%feature("docstring") TwoDfold_backtrack_f5 "

Backtrack a minimum free energy structure from a 5' section of specified length.  

This function allows one to backtrack a secondary structure beginning at the 5' end, a specified
length and residing in a specific distance class. If the argument 'k' gets a value of -1, the
structure that is backtracked is assumed to reside in the distance class where all structures
exceeding the maximum basepair distance specified in TwoDfold() belong to.  

.. deprecated:: 2.6.4
    Use the new API that relies on RNA.fold_compound() and the corresponding functions
    RNA.fold_compound_TwoD(), RNA.mfe_TwoD(), RNA.backtrack5_TwoD(), and
    RNA.fold_compound_free() instead!  

Note
----
The argument 'vars' must contain precalculated energy values in the energy matrices, i.e. a call to
TwoDfold() preceding this function is mandatory!  

Parameters
----------
j : unsigned int
    The length in nucleotides beginning from the 5' end  
k : int
    distance to reference1 (may be -1)  
l : int
    distance to reference2  
vars : TwoDfold_vars *
    the datastructure containing all predefined folding attributes  
";

%feature("docstring") TwoDfold "
";

%feature("docstring") TwoDfold_solution "
";

// File: group__kl__neighborhood__pf.xml

%feature("docstring") vrna_pf_TwoD "

Compute the partition function for all distance classes.  

This function computes the partition functions for all distance classes according the two reference
structures specified in the datastructure 'vars'. Similar to RNA.mfe_TwoD() the arguments
maxDistance1 and maxDistance2 specify the maximum distance to both reference structures. A value of
'-1' in either of them makes the appropriate distance restrictionless, i.e. all basepair distancies
to the reference are taken into account during computation. In case there is a restriction, the
returned solution contains an entry where the attribute k=l=-1 contains the partition function for
all structures exceeding the restriction. A value of INF in the attribute 'k' of the returned list
denotes the end of the list  

Parameters
----------
fc : RNA.fold_compound() *
    The datastructure containing all necessary folding attributes and matrices  
maxDistance1 : int
    The maximum basepair distance to reference1 (may be -1)  
maxDistance2 : int
    The maximum basepair distance to reference2 (may be -1)  

Returns
-------
RNA.sol_TwoD_pf() *  
    A list of partition funtions for the corresponding distance classes  

See Also
--------
RNA.fold_compound_TwoD(), RNA.fold_compound_free(), RNA.fold_compoundRNA.sol_TwoD_pf()  
";

// File: group__kl__neighborhood__stochbt.xml

%feature("docstring") vrna_pbacktrack_TwoD "

Sample secondary structure representatives from a set of distance classes according to their
Boltzmann probability.  

If the argument 'd1' is set to '-1', the structure will be backtracked in the distance class where
all structures exceeding the maximum basepair distance to either of the references reside.  

**Precondition**
    The argument 'vars' must contain precalculated partition function matrices, i.e. a call to
    RNA.pf_TwoD() preceding this function is mandatory!  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() datastructure containing all necessary folding attributes and matrices  
d1 : int
    The distance to reference1 (may be -1)  
d2 : int
    The distance to reference2  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  

See Also
--------
RNA.pf_TwoD()  
";

%feature("docstring") vrna_pbacktrack5_TwoD "

Sample secondary structure representatives with a specified length from a set of distance classes
according to their Boltzmann probability.  

This function does essentially the same as RNA.pbacktrack_TwoD() with the only difference that
partial structures, i.e. structures beginning from the 5' end with a specified length of the
sequence, are backtracked  

**Precondition**
    The argument 'vars' must contain precalculated partition function matrices, i.e. a call to
    RNA.pf_TwoD() preceding this function is mandatory!  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() datastructure containing all necessary folding attributes and matrices  
d1 : int
    The distance to reference1 (may be -1)  
d2 : int
    The distance to reference2  
length : unsigned int
    The length of the structure beginning from the 5' end  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  

See Also
--------
RNA.pbacktrack_TwoD(), RNA.pf_TwoD()  

Note
----
This function does not work (since it makes no sense) for circular RNA sequences!  
";

// File: group__thermodynamics.xml

%feature("docstring") vrna_pairing_probs "
";

%feature("docstring") vrna_mean_bp_distance_pr "

Get the mean base pair distance in the thermodynamic ensemble from a probability matrix.  

.. math::

  <d> = \\sum_{a,b} p_{a} p_{b} d(S_{a},S_{b})  

this can be computed from the pair probs :math:`p_{ij}` as  

.. math::

  <d> = \\sum_{ij} p_{ij}(1-p_{ij})  

Parameters
----------
length : int
    The length of the sequence  
pr : FLT_OR_DBL *
    The matrix containing the base pair probabilities  

Returns
-------
double  
    The mean pair distance of the structure ensemble  
";

%feature("docstring") vrna_fold_compound_t::mean_bp_distance "

Get the mean base pair distance in the thermodynamic ensemble.  

.. math::

  <d> = \\sum_{a,b} p_{a} p_{b} d(S_{a},S_{b})  

this can be computed from the pair probs :math:`p_{ij}` as  

.. math::

  <d> = \\sum_{ij} p_{ij}(1-p_{ij})  

**SWIG Wrapper Notes**
    This function is attached as method `mean_bp_distance()` to objects of type `fold_compound`.
    See, e.g.  :py:meth:`RNA.fold_compound.mean_bp_distance()` in the :doc:`/api_python`.  

Parameters
----------

Returns
-------
double  
    The mean pair distance of the structure ensemble  
";

%feature("docstring") vrna_ensemble_defect_pt "

Compute the Ensemble Defect for a given target structure provided as a **RNA.ptable**.  

Given a target structure :math:`s`, compute the average dissimilarity of a randomly drawn structure
from the ensemble, i.e.:  

.. math::

  ED(s) = 1 - \\frac{1}{n} \\sum_{ij, (i,j) \\in s} p_{ij} - \\frac{1}{n} \\sum_{i}(1 - s_{i})q_{i}  

with sequence length :math:`n`, the probability :math:`p_{ij}` of a base pair :math:`(i,j)`, the
probability :math:`q_{i} = 1 - \\sum_{j} p_{ij}` of nucleotide :math:`i` being unpaired, and the
indicator variable :math:`s_{i} = 1` if :math:`\\exists (i,j) \\in s`, and :math:`s_{i} = 0`
otherwise.  

**Precondition**
    The RNA.fold_compound() input parameter `fc` must contain a valid base pair probability matrix.
    This means that partition function and base pair probabilities must have been computed using
    `fc` before execution of this function!  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `ensemble_defect()` to objects of type
    `fold_compound`. See, e.g.  :py:meth:`RNA.fold_compound.ensemble_defect()` in the
    :doc:`/api_python`.  

Parameters
----------
fc : RNA.fold_compound() *
    A fold_compound with pre-computed base pair probabilities  
pt : const short *
    A pair table representing a target structure  

Returns
-------
double  
    The ensemble defect with respect to the target structure, or -1. upon failure, e.g. pre-
    conditions are not met  

See Also
--------
RNA.fold_compound.pf(), RNA.pairing_probs(), RNA.fold_compound.ensemble_defect()  
";

%feature("docstring") vrna_fold_compound_t::ensemble_defect "

Compute the Ensemble Defect for a given target structure.  

This is a wrapper around **RNA.ensemble_defect_pt()**. Given a target structure :math:`s`, compute
the average dissimilarity of a randomly drawn structure from the ensemble, i.e.:  

.. math::

  ED(s) = 1 - \\frac{1}{n} \\sum_{ij, (i,j) \\in s} p_{ij} - \\frac{1}{n} \\sum_{i}(1 - s_{i})q_{i}  

with sequence length :math:`n`, the probability :math:`p_{ij}` of a base pair :math:`(i,j)`, the
probability :math:`q_{i} = 1 - \\sum_{j} p_{ij}` of nucleotide :math:`i` being unpaired, and the
indicator variable :math:`s_{i} = 1` if :math:`\\exists (i,j) \\in s`, and :math:`s_{i} = 0`
otherwise.  

**Precondition**
    The RNA.fold_compound() input parameter `fc` must contain a valid base pair probability matrix.
    This means that partition function and base pair probabilities must have been computed using
    `fc` before execution of this function!  

**SWIG Wrapper Notes**
    This function is attached as method `ensemble_defect()` to objects of type `fold_compound`. Note
    that the SWIG wrapper takes a structure in dot-bracket notation and converts it into a pair
    table using RNA.ptable_from_string(). The resulting pair table is then internally passed to
    RNA.ensemble_defect_pt(). To control which kind of matching brackets will be used during
    conversion, the optional argument `options` can be used. See also the description of
    RNA.ptable_from_string() for available options. (default: `RNA.BRACKETS_RND`). See, e.g.
    :py:meth:`RNA.fold_compound.ensemble_defect()` in the :doc:`/api_python`.  

Parameters
----------
structure : const char *
    A target structure in dot-bracket notation  

Returns
-------
double  
    The ensemble defect with respect to the target structure, or -1. upon failure, e.g. pre-
    conditions are not met  

See Also
--------
RNA.fold_compound.pf(), RNA.pairing_probs(), RNA.ensemble_defect_pt()  
";

%feature("docstring") vrna_fold_compound_t::positional_entropy "

Compute a vector of positional entropies.  

This function computes the positional entropies from base pair probabilities as  

.. math::

  S(i) = - \\sum_{j} p_{ij} \\log(p_{ij}) - q_{i} \\log(q_{i})  

with unpaired probabilities :math:`q_{i} = 1 - \\sum_{j} p_{ij}`.  

Low entropy regions have little structural flexibility and the reliability of the predicted
structure is high. High entropy implies many structural alternatives. While these alternatives may
be functionally important, they make structure prediction more difficult and thus less reliable.  

**Precondition**
    This function requires pre-computed base pair probabilities! Thus, RNA.fold_compound.pf() must
be called
    beforehand.  

**SWIG Wrapper Notes**
    This function is attached as method `positional_entropy()` to objects of type `fold_compound`.
    See, e.g.  :py:meth:`RNA.fold_compound.positional_entropy()` in the :doc:`/api_python`.  

Parameters
----------

Returns
-------
double *  
    A 1-based vector of positional entropies :math:`S(i)`. (position 0 contains the sequence length)  
";

%feature("docstring") vrna_stack_prob "

Compute stacking probabilities.  

For each possible base pair :math:`(i,j)`, compute the probability of a stack :math:`(i,j)`,
:math:`(i+1, j-1)`.  

Parameters
----------
fc : RNA.fold_compound() *
    The fold compound data structure with precomputed base pair probabilities  
cutoff : double
    A cutoff value that limits the output to stacks with :math:`p > \\textrm{cutoff}`.  

Returns
-------
RNA.ep() *  
    A list of stacks with enclosing base pair :math:`(i,j)` and probabiltiy :math:`p`  
";

%feature("docstring") vrna_pf_dimer_probs "

Compute Boltzmann probabilities of dimerization without homodimers.  

Given the pair probabilities and free energies (in the null model) for a dimer AB and the two
constituent monomers A and B, compute the conditional pair probabilities given that a dimer AB
actually forms. Null model pair probabilities are given as a list as produced by
RNA.fold_compound.plist_from_probs(), the dimer probabilities 'prAB' are modified in place.  

Parameters
----------
FAB : double
    free energy of dimer AB  
FA : double
    free energy of monomer A  
FB : double
    free energy of monomer B  
prAB : RNA.ep() *
    pair probabilities for dimer  
prA : const RNA.ep() *
    pair probabilities monomer  
prB : const RNA.ep() *
    pair probabilities monomer  
Alength : int
    Length of molecule A  
exp_params : const RNA.exp_param() *
    The precomputed Boltzmann factors  
";

%feature("docstring") vrna_fold_compound_t::pr_structure "

Compute the equilibrium probability of a particular secondary structure.  

The probability :math:`p(s)` of a particular secondary structure :math:`s` can be computed as  

.. math::

  p(s) = \\frac{exp(-\\beta E(s)}{Z}  

from the structures free energy :math:`E(s)` and the partition function  

.. math::

  Z = \\sum_{s} exp(-\\beta E(s)),\\quad\\mathrm{with}\\quad\\beta = \\frac{1}{RT}  

where :math:`R` is the gas constant and :math:`T` the thermodynamic temperature.  

**Precondition**
    The fold compound `fc` must have went through a call to RNA.fold_compound.pf() to fill the
dynamic
    programming matrices with the corresponding partition function.  

**SWIG Wrapper Notes**
    This function is attached as method `pr_structure()` to objects of type `fold_compound`. See,
    e.g.  :py:meth:`RNA.fold_compound.pr_structure()` in the :doc:`/api_python`.  

Parameters
----------
structure : const char *
    The secondary structure to compute the probability for in dot-bracket notation  

Returns
-------
double  
    The probability of the input structure (range :math:`[0:1]`)  
";

%feature("docstring") vrna_fold_compound_t::pr_energy "

**SWIG Wrapper Notes**
    This function is attached as method `pr_energy()` to objects of type `fold_compound`. See, e.g.
    :py:meth:`RNA.fold_compound.pr_energy()` in the :doc:`/api_python`.  
";

%feature("docstring") vrna_fold_compound_t::heat_capacity "

Compute the specific heat for an RNA.  

This function computes an RNAs specific heat in a given temperature range from the partition
function by numeric differentiation. The result is returned as a list of pairs of temperature in C
and specific heat in Kcal/(Mol*K).  

Users can specify the temperature range for the computation from `T_min` to `T_max`, as well as the
increment step size `T_increment`. The latter also determines how many times the partition function
is computed. Finally, the parameter `mpoints` determines how smooth the curve should be. The
algorithm itself fits a parabola to :math:`2 \\cdot mpoints + 1` data points to calculate 2nd
derivatives. Increasing this parameter produces a smoother curve.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `heat_capacity()` to objects of type
    `fold_compound`. If the optional function arguments `T_min`, `T_max`, `T_increment`, and
    `mpoints` are omitted, they default to 0.0, 100.0, 1.0 and 2, respectively. See, e.g.
    :py:meth:`RNA.fold_compound.heat_capacity()` in the :doc:`/api_python`.  

Parameters
----------
T_min : float
    Lowest temperature in C  
T_max : float
    Highest temperature in C  
T_increment : float
    Stepsize for temperature incrementation in C (a reasonable choice might be 1C)  
mpoints : unsigned int
    The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2,
    min: 1, max: 100)  

Returns
-------
RNA.heat_capacity() *  
    A list of pairs of temperatures and corresponding heat capacity or *NULL* upon any failure. The
    last entry of the list is indicated by a **temperature** field set to a value smaller than
    `T_min`  

See Also
--------
RNA.fold_compound.heat_capacity_cb(), RNA.heat_capacity(), RNA.heat_capacity()  
";

%feature("docstring") vrna_fold_compound_t::heat_capacity_cb "

Compute the specific heat for an RNA (callback variant)  

Similar to RNA.fold_compound.heat_capacity(), this function computes an RNAs specific heat in a given temperature
range from the partition function by numeric differentiation. Instead of returning a list of
temperature/specific heat pairs, however, this function returns the individual results through a
callback mechanism. The provided function will be called for each result and passed the
corresponding temperature and specific heat values along with the arbitrary data as provided through
the `data` pointer argument.  

Users can specify the temperature range for the computation from `T_min` to `T_max`, as well as the
increment step size `T_increment`. The latter also determines how many times the partition function
is computed. Finally, the parameter `mpoints` determines how smooth the curve should be. The
algorithm itself fits a parabola to :math:`2 \\cdot mpoints + 1` data points to calculate 2nd
derivatives. Increasing this parameter produces a smoother curve.  

**SWIG Wrapper Notes**
    This function is attached as method `heat_capacity_cb()` to objects of type `fold_compound`.
    See, e.g.  :py:meth:`RNA.fold_compound.heat_capacity_cb()` in the :doc:`/api_python`.  

Parameters
----------
T_min : float
    Lowest temperature in C  
T_max : float
    Highest temperature in C  
T_increment : float
    Stepsize for temperature incrementation in C (a reasonable choice might be 1C)  
mpoints : unsigned int
    The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2,
    min: 1, max: 100)  
cb : RNA.heat_capacity
    The user-defined callback function that receives the individual results  
data : void *
    An arbitrary data structure that will be passed to the callback in conjunction with the results  

Returns
-------
int  
    Returns 0 upon failure, and non-zero otherwise  

See Also
--------
RNA.fold_compound.heat_capacity(), RNA.heat_capacity  
";

%feature("docstring") my_heat_capacity "

Compute the specific heat for an RNA (simplified variant)  

Similar to RNA.fold_compound.heat_capacity(), this function computes an RNAs specific heat in a given temperature
range from the partition function by numeric differentiation. This simplified version, however, only
requires the RNA sequence as input instead of a RNA.fold_compound() data structure. The result is
returned as a list of pairs of temperature in C and specific heat in Kcal/(Mol*K).  

Users can specify the temperature range for the computation from `T_min` to `T_max`, as well as the
increment step size `T_increment`. The latter also determines how many times the partition function
is computed. Finally, the parameter `mpoints` determines how smooth the curve should be. The
algorithm itself fits a parabola to :math:`2 \\cdot mpoints + 1` data points to calculate 2nd
derivatives. Increasing this parameter produces a smoother curve.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `heat_capacity()`. If the optional function
    arguments `T_min`, `T_max`, `T_increment`, and `mpoints` are omitted, they default to 0.0,
    100.0, 1.0 and 2, respectively. See, e.g.  :py:func:`RNA.head_capacity()` in the
    :doc:`/api_python`.  

Parameters
----------
sequence : const char *
    The RNA sequence input (must be uppercase)  
T_min : float
    Lowest temperature in C  
T_max : float
    Highest temperature in C  
T_increment : float
    Stepsize for temperature incrementation in C (a reasonable choice might be 1C)  
mpoints : unsigned int
    The number of interpolation points to calculate 2nd derivative (a reasonable choice might be 2,
    min: 1, max: 100)  

Returns
-------
RNA.heat_capacity() *  
    A list of pairs of temperatures and corresponding heat capacity or *NULL* upon any failure. The
    last entry of the list is indicated by a **temperature** field set to a value smaller than
    `T_min`  

See Also
--------
RNA.fold_compound.heat_capacity(), RNA.fold_compound.heat_capacity_cb(), RNA.heat_capacity(),
RNA.heat_capacity()  
";

// File: group__dos.xml

// File: group__inverse__fold.xml

%feature("docstring") my_inverse_fold "

Find sequences with predefined structure.  

This function searches for a sequence with minimum free energy structure provided in the parameter
'target', starting with sequence 'start'. It returns 0 if the search was successful, otherwise a
structure distance in terms of the energy difference between the search result and the actual target
'target' is returned. The found sequence is returned in 'start'. If give_up is set to 1, the
function will return as soon as it is clear that the search will be unsuccessful, this speeds up the
algorithm if you are only interested in exact solutions.  

Parameters
----------
start : char *
    The start sequence  
target : const char *
    The target secondary structure in dot-bracket notation  

Returns
-------
float  
    The distance to the target in case a search was unsuccessful, 0 otherwise  
";

%feature("docstring") my_inverse_pf_fold "

Find sequence that maximizes probability of a predefined structure.  

This function searches for a sequence with maximum probability to fold into the provided structure
'target' using the partition function algorithm. It returns :math:`-kT \\cdot \\log(p)` where
:math:`p` is the frequency of 'target' in the ensemble of possible structures. This is usually much
slower than inverse_fold().  

Parameters
----------
start : char *
    The start sequence  
target : const char *
    The target secondary structure in dot-bracket notation  

Returns
-------
float  
    The distance to the target in case a search was unsuccessful, 0 otherwise  
";

// File: group__neighbors.xml

%feature("docstring") vrna_move_init "

Create an atomic move.  

Parameters
----------
pos_5 : int
    The 5' position of the move (positive for insertions, negative for removal, any value for shift
    moves)  
pos_3 : int
    The 3' position of the move (positive for insertions, negative for removal, any value for shift
    moves)  

Returns
-------
RNA.move()  
    An atomic move as specified by `pos_5` and `pos_3`  

See Also
--------
RNA.move()  
";

%feature("docstring") vrna_move_list_free "

delete all moves in a zero terminated list.  
";

%feature("docstring") vrna_move_apply "

Apply a particular move / transition to a secondary structure, i.e. transform a structure.  

Parameters
----------
pt : short *
    The pair table representation of the secondary structure  
m : const RNA.move() *
    The move to apply  
";

%feature("docstring") vrna_move_apply_db "
";

%feature("docstring") vrna_move_t::is_removal "

Test whether a move is a base pair removal.  

Parameters
----------

Returns
-------
int  
    Non-zero if the move is a base pair removal, 0 otherwise  
";

%feature("docstring") vrna_move_t::is_insertion "

Test whether a move is a base pair insertion.  

Parameters
----------

Returns
-------
int  
    Non-zero if the move is a base pair insertion, 0 otherwise  
";

%feature("docstring") vrna_move_t::is_shift "

Test whether a move is a base pair shift.  

Parameters
----------

Returns
-------
int  
    Non-zero if the move is a base pair shift, 0 otherwise  
";

%feature("docstring") vrna_move_t::compare "

Compare two moves.  

The function compares two moves `m` and `b` and returns whether move `m` is lexicographically
smaller (-1), larger (1) or equal to move `b`.  

If any of the moves `m` or `b` is a shift move, this comparison only makes sense in a structure
context. Thus, the third argument with the current structure must be provided.  

Parameters
----------
b : const RNA.move() *
    The second move of the comparison  
pt : const short *
    The pair table of the current structure that is compatible with both moves (maybe NULL if moves
    are guaranteed to be no shifts)  

Returns
-------
int  
    -1 if `m` < `b`, 1 if `m` > `b`, 0 otherwise  

Warnings
--------
Currently, shift moves are not supported!  

Note
----
This function returns 0 (equality) upon any error, e.g. missing input  
";

%feature("docstring") vrna_loopidx_update "

Alters the loopIndices array that was constructed with RNA.loopidx_from_ptable().  

The loopIndex of the current move will be inserted. The correctness of the input will not be checked
because the speed should be optimized.  

Parameters
----------
loopidx : int *
    The loop index data structure that needs an update  
pt : const short *
    A pair table on which the move will be executed  
length : int
    The length of the structure  
m : const RNA.move() *
    The move that is applied to the current structure  
";

%feature("docstring") vrna_fold_compound_t::neighbors "

Generate neighbors of a secondary structure.  

This function allows one to generate all structural neighbors (according to a particular move set)
of an RNA secondary structure. The neighborhood is then returned as a list of transitions / moves
required to transform the current structure into the actual neighbor.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `neighbors()` to objects of type
    `fold_compound`. The optional parameter `options` defaults to RNA.MOVESET_DEFAULT if it is
    omitted. See, e.g.  :py:meth:`RNA.fold_compound.neighbors()` in the :doc:`/api_python`.  

Parameters
----------
pt : const short *
    The pair table representation of the structure  
options : unsigned int
    Options to modify the behavior of this function, e.g. available move set  

Returns
-------
RNA.move() *  
    Neighbors as a list of moves / transitions (the last element in the list has both of its fields
    set to 0)  

See Also
--------
RNA.neighbors_successive(), RNA.move_apply(), RNA.MOVESET_INSERTION, RNA.MOVESET_DELETION,
RNA.MOVESET_SHIFT, RNA.MOVESET_DEFAULT  
";

%feature("docstring") vrna_neighbors_successive "

Generate neighbors of a secondary structure (the fast way)  

This function implements a fast way to generate all neighbors of a secondary structure that results
from successive applications of individual moves. The speed-up results from updating an already
known list of valid neighbors before the individual move towards the current structure took place.
In essence, this function removes neighbors that are not accessible anymore and inserts neighbors
emerging after a move took place.  

Parameters
----------
fc : const RNA.fold_compound() *
    A RNA.fold_compound() containing the energy parameters and model details  
curr_move : const RNA.move() *
    The move that was/will be applied to `prev_pt`  
prev_pt : const short *
    A pair table representation of the structure before `curr_move` is/was applied  
prev_neighbors : const RNA.move() *
    The list of neighbors of `prev_pt`  
size_prev_neighbors : int
    The size of `prev_neighbors`, i.e. the lists length  
size_neighbors : int *
    A pointer to store the size / length of the new neighbor list  
options : unsigned int
    Options to modify the behavior of this function, e.g. available move set  

Returns
-------
RNA.move() *  
    Neighbors as a list of moves / transitions (the last element in the list has both of its fields
    set to 0)  

See Also
--------
RNA.fold_compound.neighbors(), RNA.move_apply(), RNA.MOVESET_INSERTION, RNA.MOVESET_DELETION,
RNA.MOVESET_SHIFT, RNA.MOVESET_DEFAULT  
";

%feature("docstring") vrna_move_neighbor_diff_cb "

Apply a move to a secondary structure and indicate which neighbors have changed consequentially.  

This function applies a move to a secondary structure and explores the local neighborhood of the
affected loop. Any changes to previously compatible neighbors that have been affected by this loop
will be reported through a callback function. In particular, any of the three cases might appear:  

*   A previously available neighbor move has changed, usually the free energy change of the move
    (RNA.NEIGHBOR_CHANGE)  
*   A previously available neighbor move became invalid (RNA.NEIGHBOR_INVALID)  
*   A new neighbor move becomes available (RNA.NEIGHBOR_NEW)  

Parameters
----------
fc : RNA.fold_compound() *
    A fold compound for the RNA sequence(s) that this function operates on  
ptable : short *
    The current structure as pair table  
move : RNA.move()
    The move to apply  
cb : RNA.move_update
    The address of the callback function that is passed the neighborhood changes  
data : void *
    An arbitrary data pointer that will be passed through to the callback function `cb`  
options : unsigned int
    Options to modify the behavior of this function, .e.g available move set  

Returns
-------
int  
    Non-zero on success, 0 otherwise  

See Also
--------
RNA.fold_compound.move_neighbor_diff(), RNA.NEIGHBOR_CHANGE, RNA.NEIGHBOR_INVALID, RNA.NEIGHBOR_NEW,
RNA.move_update, RNA.MOVE_NO_APPLY  
";

%feature("docstring") vrna_fold_compound_t::move_neighbor_diff "

Apply a move to a secondary structure and indicate which neighbors have changed consequentially.  

Similar to RNA.move_neighbor_diff_cb(), this function applies a move to a secondary structure and
reports back the neighbors of the current structure become affected by this move. Instead of
executing a callback for each of the affected neighbors, this function compiles two lists of
neighbor moves, one that is returned and consists of all moves that are novel or may have changed in
energy, and a second, `invalid_moves`, that consists of all the neighbor moves that become invalid,
respectively.  

Parameters
----------
ptable : short *
    The current structure as pair table  
move : RNA.move()
    The move to apply  
invalid_moves : RNA.move() **
    The address of a move list where the function stores those moves that become invalid  
options : unsigned int
    Options to modify the behavior of this function, .e.g available move set  

Returns
-------
RNA.move() *  
    A list of moves that might have changed in energy or are novel compared to the structure before
    application of the move  
";

%feature("docstring") MOVESET_INSERTION "

Option flag indicating insertion move.  

See Also
--------
RNA.fold_compound.neighbors(), RNA.neighbors_successive, RNA.fold_compound.path()  
";

%feature("docstring") MOVESET_DELETION "

Option flag indicating deletion move.  

See Also
--------
RNA.fold_compound.neighbors(), RNA.neighbors_successive, RNA.fold_compound.path()  
";

%feature("docstring") MOVESET_SHIFT "

Option flag indicating shift move.  

See Also
--------
RNA.fold_compound.neighbors(), RNA.neighbors_successive, RNA.fold_compound.path()  
";

%feature("docstring") MOVESET_NO_LP "

Option flag indicating moves without lonely base pairs.  

See Also
--------
RNA.fold_compound.neighbors(), RNA.neighbors_successive, RNA.fold_compound.path()  
";

%feature("docstring") MOVESET_DEFAULT "

Option flag indicating default move set, i.e. insertions/deletion of a base pair.  

See Also
--------
RNA.fold_compound.neighbors(), RNA.neighbors_successive, RNA.fold_compound.path()  
";

%feature("docstring") MOVE_NO_APPLY "
";

%feature("docstring") NEIGHBOR_CHANGE "

State indicator for a neighbor that has been changed.  

See Also
--------
RNA.move_neighbor_diff_cb()  
";

%feature("docstring") NEIGHBOR_INVALID "

State indicator for a neighbor that has been invalidated.  

See Also
--------
RNA.move_neighbor_diff_cb()  
";

%feature("docstring") NEIGHBOR_NEW "

State indicator for a neighbor that has become newly available.  

See Also
--------
RNA.move_neighbor_diff_cb()  
";

// File: group__paths.xml

%feature("docstring") vrna_path_free "

Release (free) memory occupied by a (re-)folding path.  

Parameters
----------
path : RNA.path() *
    The refolding path to be free'd  

See Also
--------
RNA.path_direct(), RNA.fold_compound.path_direct(), RNA.path_findpath(), RNA.path_findpath_ub()  
";

%feature("docstring") vrna_path_options_free "

Release (free) memory occupied by an options data structure for (re-)folding path implementations.  

Parameters
----------
options : RNA.path_options()
    The options data structure to be free'd  

See Also
--------
RNA.path_options_findpath(), RNA.path_direct(), RNA.fold_compound.path_direct()  
";

%feature("docstring") PATH_TYPE_DOT_BRACKET "

Flag to indicate producing a (re-)folding path as list of dot-bracket structures.  

See Also
--------
RNA.path(), RNA.path_options_findpath(), RNA.path_direct(), RNA.fold_compound.path_direct()  
";

%feature("docstring") PATH_TYPE_MOVES "

Flag to indicate producing a (re-)folding path as list of transition moves.  

See Also
--------
RNA.path(), RNA.path_options_findpath(), RNA.path_direct(), RNA.fold_compound.path_direct()  
";

// File: group__paths__direct.xml

%feature("docstring") vrna_path_findpath_saddle "

Find energy of a saddle point between 2 structures (search only direct path)  

This function uses an inplementation of the *findpath* algorithm   :cite:p:`flamm:2001`  for near-
optimal direct refolding path prediction.  

Model details, and energy parameters are used as provided via the parameter 'fc'. The
RNA.fold_compound() does not require memory for any DP matrices, but requires all most basic init
values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_findpath_saddle()` to objects of type
    `fold_compound`. The optional parameter `width` defaults to 1 if it is omitted. See, e.g.
    :py:meth:`RNA.fold_compound.path_findpath_saddle()` in the :doc:`/api_python`.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() with precomputed sequence encoding and model details  
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
width : int
    A number specifying how many strutures are being kept at each step during the search  

Returns
-------
int  
    The saddle energy in 10cal/mol  

See Also
--------
RNA.fold_compound.path_findpath_saddle(), RNA.fold_compound(), RNA.fold_compound(),
RNA.path_findpath()  
";

%feature("docstring") vrna_fold_compound_t::path_findpath_saddle "

Find energy of a saddle point between 2 structures (search only direct path)  

This function uses an inplementation of the *findpath* algorithm   :cite:p:`flamm:2001`  for near-
optimal direct refolding path prediction.  

Model details, and energy parameters are used as provided via the parameter 'fc'. The
RNA.fold_compound() does not require memory for any DP matrices, but requires all most basic init
values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_findpath_saddle()` to objects of type
    `fold_compound`. The optional parameter `width` defaults to 1 if it is omitted, while the
    optional parameter `maxE` defaults to INF. In case the function did not find a path with
    :math:`E_{saddle} < E_{max}` the function returns a `NULL` object, i.e. `undef` for Perl and
    `None` for Python. See, e.g.  :py:meth:`RNA.fold_compound.path_findpath_saddle()` in the
    :doc:`/api_python`.  

Parameters
----------
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
width : int
    A number specifying how many strutures are being kept at each step during the search  
maxE : int
    An upper bound for the saddle point energy in 10cal/mol  

Returns
-------
int  
    The saddle energy in 10cal/mol  

Warnings
--------
The argument `maxE` ( :math:`E_{max}`) enables one to specify an upper bound, or maximum free energy
for the saddle point between the two input structures. If no path with :math:`E_{saddle} < E_{max}`
is found, the function simply returns `maxE`  

See Also
--------
RNA.path_findpath_saddle(), RNA.fold_compound(), RNA.fold_compound(), RNA.path_findpath()  
";

%feature("docstring") vrna_path_findpath "

Find refolding path between 2 structures (search only direct path)  

This function uses an inplementation of the *findpath* algorithm   :cite:p:`flamm:2001`  for near-
optimal direct refolding path prediction.  

Model details, and energy parameters are used as provided via the parameter 'fc'. The
RNA.fold_compound() does not require memory for any DP matrices, but requires all most basic init
values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_findpath()` to objects of type
    `fold_compound`. The optional parameter `width` defaults to 1 if it is omitted. See, e.g.
    :py:meth:`RNA.fold_compound.path_findpath()` in the :doc:`/api_python`.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() with precomputed sequence encoding and model details  
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
width : int
    A number specifying how many strutures are being kept at each step during the search  

Returns
-------
RNA.path() *  
    The saddle energy in 10cal/mol  

See Also
--------
RNA.path_findpath_ub(), RNA.fold_compound(), RNA.fold_compound(), RNA.path_findpath_saddle()  
";

%feature("docstring") vrna_path_findpath_ub "

Find refolding path between 2 structures (search only direct path)  

This function uses an inplementation of the *findpath* algorithm   :cite:p:`flamm:2001`  for near-
optimal direct refolding path prediction.  

Model details, and energy parameters are used as provided via the parameter 'fc'. The
RNA.fold_compound() does not require memory for any DP matrices, but requires all most basic init
values as one would get from a call like this:  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_findpath()` to objects of type
    `fold_compound`. The optional parameter `width` defaults to 1 if it is omitted, while the
    optional parameter `maxE` defaults to INF. In case the function did not find a path with
    :math:`E_{saddle} < E_{max}` the function returns an empty list. See, e.g.
    :py:meth:`RNA.fold_compound.path_findpath()` in the :doc:`/api_python`.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() with precomputed sequence encoding and model details  
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
width : int
    A number specifying how many strutures are being kept at each step during the search  
maxE : int
    An upper bound for the saddle point energy in 10cal/mol  

Returns
-------
RNA.path() *  
    The saddle energy in 10cal/mol  

Warnings
--------
The argument `maxE` enables one to specify an upper bound, or maximum free energy for the saddle
point between the two input structures. If no path with :math:`E_{saddle} < E_{max}` is found, the
function simply returns *NULL*  

See Also
--------
RNA.path_findpath(), RNA.fold_compound(), RNA.fold_compound(), RNA.path_findpath_saddle()  
";

%feature("docstring") my_path_options_findpath "

Create options data structure for findpath direct (re-)folding path heuristic.  

This function returns an options data structure that switches the RNA.path_direct() and
RNA.fold_compound.path_direct() API functions to use the *findpath* :cite:p:`flamm:2001`  heuristic. The
parameter `width` specifies the width of the breadth-first search while the second parameter `type`
allows one to set the type of the returned (re-)folding path.  

Currently, the following return types are available:  

*   A list of dot-bracket structures and corresponding free energy (flag:
    RNA.PATH_TYPE_DOT_BRACKET)  
*   A list of transition moves and corresponding free energy changes (flag: RNA.PATH_TYPE_MOVES)  

**SWIG Wrapper Notes**
    This function is available as overloaded function `path_options_findpath()`. The optional
    parameter `width` defaults to 10 if omitted, while the optional parameter `type` defaults to
    RNA.PATH_TYPE_DOT_BRACKET. See, e.g.  :py:func:`RNA.path_options_findpath()` in the
    :doc:`/api_python`.  

Parameters
----------
width : int
    Width of the breath-first search strategy  
type : unsigned int
    Setting that specifies how the return (re-)folding path should be encoded  

Returns
-------
RNA.path_options()  
    An options data structure with settings for the findpath direct path heuristic  

See Also
--------
RNA.PATH_TYPE_DOT_BRACKET, RNA.PATH_TYPE_MOVES, RNA.path_options_free(), RNA.path_direct(),
RNA.fold_compound.path_direct()  
";

%feature("docstring") vrna_path_direct "

Determine an optimal direct (re-)folding path between two secondary structures.  

This is the generic wrapper function to retrieve (an optimal) (re-)folding path between two
secondary structures `s1` and `s2`. The actual algorithm that is used to generate the (re-)folding
path is determined by the settings specified in the `options` data structure. This data structure
also determines the return type, which might be either:  

*   a list of dot-bracket structures with corresponding free energy, or  
*   a list of transition moves with corresponding free energy change  

If the `options` parameter is passed a *NULL* pointer, this function defaults to the *findpath
heuristic* :cite:p:`flamm:2001`  with a breadth-first search width of :math:`10`, and the returned
path consists of dot-bracket structures with corresponding free energies.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_direct()` to objects of type
    `fold_compound`. The optional parameter `options` defaults to `NULL` if it is omitted. See, e.g.
    :py:meth:`RNA.fold_compound.path_direct()` in the :doc:`/api_python`.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() with precomputed sequence encoding and model details  
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
options : RNA.path_options()
    An options data structure that specifies the path heuristic and corresponding settings (maybe
    *NULL*)  

Returns
-------
RNA.path() *  
    An optimal (re-)folding path between the two input structures  

See Also
--------
RNA.fold_compound.path_direct(), RNA.path_options_findpath(), RNA.path_options_free(),
RNA.path_free()  
";

%feature("docstring") vrna_fold_compound_t::path_direct "

Determine an optimal direct (re-)folding path between two secondary structures.  

This function is similar to RNA.path_direct(), but allows to specify an *upper-bound* for the
saddle point energy. The underlying algorithms will stop determining an (optimal) (re-)folding path,
if none can be found that has a saddle point below the specified upper-bound threshold `maxE`.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_direct()` to objects of type
    `fold_compound`. The optional parameter `maxE` defaults to #INT_MAX - 1 if it is omitted, while
    the optional parameter `options` defaults to `NULL`. In case the function did not find a path
    with :math:`E_{saddle} < E_{max}` it returns an empty list. See, e.g.
    :py:meth:`RNA.fold_compound.path_direct()` in the :doc:`/api_python`.  

Parameters
----------
s1 : const char *
    The start structure in dot-bracket notation  
s2 : const char *
    The target structure in dot-bracket notation  
maxE : int
    Upper bound for the saddle point along the (re-)folding path  
options : RNA.path_options()
    An options data structure that specifies the path heuristic and corresponding settings (maybe
    *NULL*)  

Returns
-------
RNA.path() *  
    An optimal (re-)folding path between the two input structures  

Warnings
--------
The argument `maxE` enables one to specify an upper bound, or maximum free energy for the saddle
point between the two input structures. If no path with :math:`E_{saddle} < E_{max}` is found, the
function simply returns *NULL*  

See Also
--------
RNA.fold_compound.path_direct(), RNA.path_options_findpath(), RNA.path_options_free(),
RNA.path_free()  
";

// File: group__paths__walk.xml

%feature("docstring") vrna_fold_compound_t::path "

Compute a path, store the final structure, and return a list of transition moves from the start to
the final structure.  

This function computes, given a start structure in pair table format, a transition path, updates the
pair table to the final structure of the path. Finally, if not requested otherwise by using the
RNA.PATH_NO_TRANSITION_OUTPUT flag in the `options` field, this function returns a list of
individual transitions that lead from the start to the final structure if requested.  

The currently available transition paths are  

*   Steepest Descent / Gradient walk (flag: RNA.PATH_STEEPEST_DESCENT)  
*   Random walk (flag: RNA.PATH_RANDOM)  

The type of transitions must be set through the `options` parameter  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path()` to objects of type `fold_compound`.
    The optional parameter `options` defaults to RNA.PATH_DEFAULT if it is omitted. See, e.g.
    :py:meth:`RNA.fold_compound.path()` in the :doc:`/api_python`.  

Parameters
----------
pt : short *
    The pair table containing the start structure. Used to update to the final structure after
    execution of this function  
options : unsigned int
    Options to modify the behavior of this function  

Returns
-------
RNA.move() *  
    A list of transition moves (default), or NULL (if options & RNA.PATH_NO_TRANSITION_OUTPUT)  

See Also
--------
RNA.fold_compound.path_gradient(), RNA.fold_compound.path_random(), RNA.ptable(), RNA.ptable_copy(),
RNA.fold_compound()RNA.PATH_STEEPEST_DESCENT, RNA.PATH_RANDOM, RNA.MOVESET_DEFAULT,
RNA.MOVESET_SHIFT, RNA.PATH_NO_TRANSITION_OUTPUT  

Note
----
Since the result is written to the input structure you may want to use RNA.ptable_copy() before
calling this function to keep the initial structure  
";

%feature("docstring") vrna_fold_compound_t::path_gradient "

Compute a steepest descent / gradient path, store the final structure, and return a list of
transition moves from the start to the final structure.  

This function computes, given a start structure in pair table format, a steepest descent path,
updates the pair table to the final structure of the path. Finally, if not requested otherwise by
using the RNA.PATH_NO_TRANSITION_OUTPUT flag in the `options` field, this function returns a list
of individual transitions that lead from the start to the final structure if requested.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_gradient()` to objects of type
    `fold_compound`. The optional parameter `options` defaults to RNA.PATH_DEFAULT if it is
    omitted. See, e.g.  :py:meth:`RNA.fold_compound.path_gradient()` in the :doc:`/api_python`.  

Parameters
----------
pt : short *
    The pair table containing the start structure. Used to update to the final structure after
    execution of this function  
options : unsigned int
    Options to modify the behavior of this function  

Returns
-------
RNA.move() *  
    A list of transition moves (default), or NULL (if options & RNA.PATH_NO_TRANSITION_OUTPUT)  

See Also
--------
RNA.fold_compound.path_random(), RNA.fold_compound.path(), RNA.ptable(), RNA.ptable_copy(),
RNA.fold_compound()RNA.MOVESET_DEFAULT, RNA.MOVESET_SHIFT, RNA.PATH_NO_TRANSITION_OUTPUT  

Note
----
Since the result is written to the input structure you may want to use RNA.ptable_copy() before
calling this function to keep the initial structure  
";

%feature("docstring") vrna_fold_compound_t::path_random "

Generate a random walk / path of a given length, store the final structure, and return a list of
transition moves from the start to the final structure.  

This function generates, given a start structure in pair table format, a random walk / path, updates
the pair table to the final structure of the path. Finally, if not requested otherwise by using the
RNA.PATH_NO_TRANSITION_OUTPUT flag in the `options` field, this function returns a list of
individual transitions that lead from the start to the final structure if requested.  

**SWIG Wrapper Notes**
    This function is attached as an overloaded method `path_gradient()` to objects of type
    `fold_compound`. The optional parameter `options` defaults to RNA.PATH_DEFAULT if it is
    omitted. See, e.g.  :py:meth:`RNA.fold_compound.path_random()` in the :doc:`/api_python`.  

Parameters
----------
pt : short *
    The pair table containing the start structure. Used to update to the final structure after
    execution of this function  
steps : unsigned int
    The length of the path, i.e. the total number of transitions / moves  
options : unsigned int
    Options to modify the behavior of this function  

Returns
-------
RNA.move() *  
    A list of transition moves (default), or NULL (if options & RNA.PATH_NO_TRANSITION_OUTPUT)  

See Also
--------
RNA.fold_compound.path_gradient(), RNA.fold_compound.path(), RNA.ptable(), RNA.ptable_copy(),
RNA.fold_compound()RNA.MOVESET_DEFAULT, RNA.MOVESET_SHIFT, RNA.PATH_NO_TRANSITION_OUTPUT  

Note
----
Since the result is written to the input structure you may want to use RNA.ptable_copy() before
calling this function to keep the initial structure  
";

%feature("docstring") PATH_STEEPEST_DESCENT "

Option flag to request a steepest descent / gradient path.  

See Also
--------
RNA.fold_compound.path()  
";

%feature("docstring") PATH_RANDOM "

Option flag to request a random walk path.  

See Also
--------
RNA.fold_compound.path()  
";

%feature("docstring") PATH_NO_TRANSITION_OUTPUT "

Option flag to omit returning the transition path.  

See Also
--------
RNA.fold_compound.path(), RNA.fold_compound.path_gradient(), RNA.fold_compound.path_random()  
";

%feature("docstring") PATH_DEFAULT "

Option flag to request defaults (steepest descent / default move set)  

See Also
--------
RNA.fold_compound.path(), RNA.PATH_STEEPEST_DESCENT, RNA.MOVESET_DEFAULT  
";

// File: group__probing__data.xml

// File: group__SHAPE__reactivities.xml

%feature("docstring") vrna_constraints_add_SHAPE "
";

%feature("docstring") vrna_constraints_add_SHAPE_ali "
";

%feature("docstring") vrna_fold_compound_t::sc_add_SHAPE_deigan "

Add SHAPE reactivity data as soft constraints (Deigan et al. method)  

This approach of SHAPE directed RNA folding uses the simple linear ansatz  

.. math::

  \\Delta G_{\\text{SHAPE}}(i) = m \\ln(\\text{SHAPE reactivity}(i)+1)+ b  

to convert SHAPE reactivity values to pseudo energies whenever a nucleotide :math:`i` contributes to
a stacked pair. A positive slope :math:`m` penalizes high reactivities in paired regions, while a
negative intercept :math:`b` results in a confirmatory `bonus' free energy for correctly predicted
base pairs. Since the energy evaluation of a base pair stack involves two pairs, the pseudo energies
are added for all four contributing nucleotides. Consequently, the energy term is applied twice for
pairs inside a helix and only once for pairs adjacent to other structures. For all other loop types
the energy model remains unchanged even when the experimental data highly disagrees with a certain
motif.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_SHAPE_deigan()` to objects of type `fold_compound`.
    See, e.g.   :py:meth:`RNA.fold_compound.sc_add_SHAPE_deigan()` in the :doc:`/api_python` .  

Parameters
----------
reactivities : const double *
    A vector of normalized SHAPE reactivities  
m : double
    The slope of the conversion function  
b : double
    The intercept of the conversion function  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    1 on successful extraction of the method, 0 on errors  

See Also
--------
RNA.fold_compound.sc_remove(), RNA.fold_compound.sc_add_SHAPE_zarringhalam(),
RNA.sc_minimize_pertubation()  

Note
----
For further details, we refer to   :cite:t:`deigan:2009` .  
";

%feature("docstring") vrna_fold_compound_t::sc_add_SHAPE_deigan_ali "

Add SHAPE reactivity data from files as soft constraints for consensus structure prediction (Deigan
et al. method)  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_SHAPE_deigan_ali()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.sc_add_SHAPE_deigan_ali()` in the
    :doc:`/api_python` .  

Parameters
----------
shape_files : const char **
    A set of filenames that contain normalized SHAPE reactivity data  
shape_file_association : const int *
    An array of integers that associate the files with sequences in the alignment  
m : double
    The slope of the conversion function  
b : double
    The intercept of the conversion function  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    1 on successful extraction of the method, 0 on errors  
";

%feature("docstring") vrna_fold_compound_t::sc_add_SHAPE_zarringhalam "

Add SHAPE reactivity data as soft constraints (Zarringhalam et al. method)  

This method first converts the observed SHAPE reactivity of nucleotide :math:`i` into a probability
:math:`q_{i}` that position :math:`i` is unpaired by means of a non-linear map. Then pseudo-energies
of the form  

.. math::

  \\Delta G_{\\text{SHAPE}}(x,i) = \\beta\\ |x_{i} - q_{i}|  

are computed, where :math:`x_{i}=0` if position :math:`i` is unpaired and :math:`x_{i}=1` if
:math:`i` is paired in a given secondary structure. The parameter :math:`\\beta` serves as scaling
factor. The magnitude of discrepancy between prediction and experimental observation is represented
by :math:`|x_{i} - q_{i}|`.  

**SWIG Wrapper Notes**
    This function is attached as method `sc_add_SHAPE_zarringhalam()` to objects of type
    `fold_compound`. See, e.g.   :py:meth:`RNA.fold_compound.sc_add_SHAPE_zarringhalam()` in the
    :doc:`/api_python` .  

Parameters
----------
reactivities : const double *
    A vector of normalized SHAPE reactivities  
b : double
    The scaling factor :math:`\\beta` of the conversion function  
default_value : double
    The default value for a nucleotide where reactivity data is missing for  
shape_conversion : const char *
    A flag that specifies how to convert reactivities to probabilities  
options : unsigned int
    The options flag indicating how/where to store the soft constraints  

Returns
-------
int  
    1 on successful extraction of the method, 0 on errors  

See Also
--------
RNA.fold_compound.sc_remove(), RNA.fold_compound.sc_add_SHAPE_deigan(),
RNA.sc_minimize_pertubation()  

Note
----
For further details, we refer to   :cite:t:`zarringhalam:2012`  
";

%feature("docstring") vrna_sc_SHAPE_to_pr "

Convert SHAPE reactivity values to probabilities for being unpaired.  

This function parses the informations from a given file and stores the result in the preallocated
string sequence and the FLT_OR_DBL array values.  

Parameters
----------
shape_conversion : const char *
    String definining the method used for the conversion process  
values : double *
    Pointer to an array of SHAPE reactivities  
length : int
    Length of the array of SHAPE reactivities  
default_value : double
    Result used for position with invalid/missing reactivity values  

See Also
--------
RNA.file_SHAPE_read()  
";

// File: group__perturbation.xml

%feature("docstring") vrna_sc_minimize_pertubation "

Find a vector of perturbation energies that minimizes the discripancies between predicted and
observed pairing probabilities and the amount of neccessary adjustments.  

Use an iterative minimization algorithm to find a vector of perturbation energies whose
incorporation as soft constraints shifts the predicted pairing probabilities closer to the
experimentally observed probabilities. The algorithm aims to minimize an objective function that
penalizes discripancies between predicted and observed pairing probabilities and energy model
adjustments, i.e. an appropriate vector of perturbation energies satisfies  

.. math::

  F(\\vec\\epsilon) = \\sum_{\\mu}{ \\frac{\\epsilon_{\\mu}^2}{\\tau^2} } + \\sum_{i =
  1}^n{ \\frac{(p_{i}(\\vec\\epsilon) - q_{i})^2}{\\sigma^2} } \\to \\min.  

An initialized fold compound and an array containing the observed probability for each nucleotide to
be unbound are required as input data. The parameters objective_function, sigma_squared and
tau_squared are responsible for adjusting the aim of the objective function. Dependend on which type
of objective function is selected, either squared or absolute aberrations are contributing to the
objective function. The ratio of the parameters sigma_squared and tau_squared can be used to adjust
the algorithm to find a solution either close to the thermodynamic prediction (sigma_squared >>
tau_squared) or close to the experimental data (tau_squared >> sigma_squared). The minimization can
be performed by makeing use of a custom gradient descent implementation or using one of the
minimizing algorithms provided by the GNU Scientific Library. All algorithms require the evaluation
of the gradient of the objective function, which includes the evaluation of conditional pairing
probabilites. Since an exact evaluation is expensive, the probabilities can also be estimated from
sampling by setting an appropriate sample size. The found vector of perturbation energies will be
stored in the array epsilon. The progress of the minimization process can be tracked by implementing
and passing a callback function.  

Parameters
----------
fc : RNA.fold_compound() *
    Pointer to a fold compound  
q_prob_unpaired : const double *
    Pointer to an array containing the probability to be unpaired for each nucleotide  
objective_function : int
    The type of objective function to be used (RNA.OBJECTIVE_FUNCTION_QUADRATIC /
    RNA.OBJECTIVE_FUNCTION_LINEAR)  
sigma_squared : double
    A factor used for weighting the objective function. More weight on this factor will lead to a
    solution close to the null vector.  
tau_squared : double
    A factor used for weighting the objective function. More weight on this factor will lead to a
    solution close to the data provided in q_prob_unpaired.  
algorithm : int
    The minimization algorithm (RNA.MINIMIZER_*)  
sample_size : int
    The number of sampled sequences used for estimating the pairing probabilities. A value <= 0 will
    lead to an exact evaluation.  
epsilon : double *
    A pointer to an array used for storing the calculated vector of perturbation energies  
callback : progress_callback
    A pointer to a callback function used for reporting the current minimization progress  

See Also
--------
For further details we refer to   :cite:t:`washietl:2012` .  
";

%feature("docstring") OBJECTIVE_FUNCTION_QUADRATIC "

Use the sum of squared aberrations as objective function.  

:math:`F(\\vec\\epsilon) = \\sum_{i = 1}^n{ \\frac{\\epsilon_{i}^2}{\\tau^2} } + \\sum_{i = 1}^n{
\\frac{(p_{i}(\\vec\\epsilon) - q_{i})^2}{\\sigma^2} } \\to min`  
";

%feature("docstring") OBJECTIVE_FUNCTION_ABSOLUTE "

Use the sum of absolute aberrations as objective function.  

:math:`F(\\vec\\epsilon) = \\sum_{i = 1}^n{ \\frac{|\\epsilon_{i}|}{\\tau^2} } + \\sum_{i = 1}^n{
\\frac{|p_{i}(\\vec\\epsilon) - q_{i}|}{\\sigma^2} } \\to min`  
";

%feature("docstring") MINIMIZER_DEFAULT "

Use a custom implementation of the gradient descent algorithm to minimize the objective function.  
";

%feature("docstring") MINIMIZER_CONJUGATE_FR "

Use the GNU Scientific Library implementation of the Fletcher-Reeves conjugate gradient algorithm to
minimize the objective function.  

Please note that this algorithm can only be used when the GNU Scientific Library is available on
your system  
";

%feature("docstring") MINIMIZER_CONJUGATE_PR "

Use the GNU Scientific Library implementation of the Polak-Ribiere conjugate gradient algorithm to
minimize the objective function.  

Please note that this algorithm can only be used when the GNU Scientific Library is available on
your system  
";

%feature("docstring") MINIMIZER_VECTOR_BFGS "

Use the GNU Scientific Library implementation of the vector Broyden-Fletcher-Goldfarb-Shanno
algorithm to minimize the objective function.  

Please note that this algorithm can only be used when the GNU Scientific Library is available on
your system  
";

%feature("docstring") MINIMIZER_VECTOR_BFGS2 "

Use the GNU Scientific Library implementation of the vector Broyden-Fletcher-Goldfarb-Shanno
algorithm to minimize the objective function.  

Please note that this algorithm can only be used when the GNU Scientific Library is available on
your system  
";

%feature("docstring") MINIMIZER_STEEPEST_DESCENT "

Use the GNU Scientific Library implementation of the steepest descent algorithm to minimize the
objective function.  

Please note that this algorithm can only be used when the GNU Scientific Library is available on
your system  
";

// File: group__ligand__binding.xml

// File: group__ligands__up.xml

// File: group__constraints__ligand.xml

%feature("docstring") vrna_fold_compound_t::sc_add_hi_motif "

Add soft constraints for hairpin or interior loop binding motif.  

Here is an example that adds a theophylline binding motif. Free energy contribution is derived from
:math:`k_{d} = 0.1 \\mu M`, taken from Jenison et al. 1994. At :math:`1M` concentration the
corresponding binding free energy amounts to :math:`-9.93~kcal/mol`.  



**SWIG Wrapper Notes**
    This function is attached as method `sc_add_hi_motif()` to objects of type `fold_compound`. The
    last parameter is optional an defaults to `options` = RNA.OPTION_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_add_hi_motif()` in the :doc:`/api_python` .  

Parameters
----------
seq : const char *
    The sequence motif (may be interspaced by '&' character  
structure : const char *
    The structure motif (may be interspaced by '&' character  
energy : FLT_OR_DBL
    The free energy of the motif (e.g. binding free energy)  
options : unsigned int
    Options  

Returns
-------
int  
    non-zero value if application of the motif using soft constraints was successful  
";

%feature("docstring") vrna_sc_ligand_detect_motifs "
";

%feature("docstring") vrna_sc_ligand_get_all_motifs "
";

// File: group__paired__modules.xml

// File: group__pseudoknots.xml

%feature("docstring") vrna_pk_plex "

Predict Pseudoknot interactions in terms of a two-step folding process.  

Computes simple pseudoknot interactions according to the PKplex algorithm. This simple heuristic
first compiles a list of potential interaction sites that may form a pseudoknot. The resulting
candidate interactions are then fixed and an PK-free MFE structure for the remainder of the sequence
is computed.  

The `accessibility` argument is a list of opening energies for potential interaction sites. It is
used in the first step of the algorithm to identify potential interactions. Upon passing *NULL*, the
opening energies are determined automatically based on the current model settings.  

Depending on the `options`, the function can return the MFE (incl. PK loops) or suboptimal
structures within an energy band around the MFE. The PK loop is internally scored by a scoring
function that in the simplest cases assigns a constant value for each PK loop. More complicated
scoring functions can be passed as well, see RNA.pk_plex_score and RNA.pk_plex_opt_fun().  

The function returns *NULL* on any error. Otherwise, a list of structures and interaction
coordinates with corresponding energy contributions is returned. If no PK-interaction that satisfies
the options is found, the list only consists of the PK-free MFE structure.  

Parameters
----------
fc : RNA.fold_compound() *
    fold compound with the input sequence and model settings  
accessibility : const int **
    An array of opening energies for the implemented heuristic (maybe *NULL*)  
options : RNA.pk_plex_opt()
    An RNA.pk_plex_opt() options data structure that determines the algorithm parameters  

Returns
-------
RNA.pk_plex() *  
    A list of potentially pseudoknotted structures (Last element in the list indicated by *NULL*
    value in RNA.pk_plex_result().structure)  
";

%feature("docstring") vrna_pk_plex_accessibility "

Obtain a list of opening energies suitable for PKplex computations.  

Parameters
----------
sequence : const char *
    The RNA sequence  
unpaired : unsigned int
    The maximum number of unpaired nucleotides, i.e. length of interaction  
cutoff : double
    A cutoff value for unpaired probabilities  

Returns
-------
int **  
    Opening energies as required for RNA.pk_plex()  

See Also
--------
RNA.pk_plex()  
";

%feature("docstring") vrna_pk_plex_opt_defaults "

Default options for PKplex algorithm.  

Returns
-------
RNA.pk_plex_opt()  
    An options data structure suitabe for PKplex computations  

See Also
--------
RNA.pk_plex(), RNA.pk_plex_opt(), RNA.pk_plex_opt_fun()  
";

%feature("docstring") vrna_pk_plex_opt "

Simple options for PKplex algorithm.  

Parameters
----------
delta : unsigned int
    Size of energy band around MFE for suboptimal results in dekacal/mol  
max_interaction_length : unsigned int
    Maximum length of interaction  
pk_penalty : int
    Energy constant to score the PK forming loop  

Returns
-------
RNA.pk_plex_opt()  
    An options data structure suitabe for PKplex computations  

See Also
--------
RNA.pk_plex(), RNA.pk_plex_opt_defaults(), RNA.pk_plex_opt_fun()  
";

%feature("docstring") vrna_pk_plex_opt_fun "

Simple options for PKplex algorithm.  

Parameters
----------
delta : unsigned int
    Size of energy band around MFE for suboptimal results in dekacal/mol  
max_interaction_length : unsigned int
    Maximum length of interaction  
scoring_function : RNA.pk_plex_score
    Energy evaluating function to score the PK forming loop  
scoring_data : void *
    An arbitrary data structure passed to the scoring function (maybe *NUL*)  

Returns
-------
RNA.pk_plex_opt()  
    An options data structure suitabe for PKplex computations  

See Also
--------
RNA.pk_plex(), RNA.pk_plex_opt_defaults(), RNA.pk_plex_opt(), RNA.pk_plex_score  
";

// File: group__gquads.xml

%feature("docstring") E_gquad "
";

%feature("docstring") exp_E_gquad "
";

%feature("docstring") E_gquad_ali_en "
";

%feature("docstring") get_gquad_matrix "

Get a triangular matrix prefilled with minimum free energy contributions of G-quadruplexes.  

At each position ij in the matrix, the minimum free energy of any G-quadruplex delimited by i and j
is stored. If no G-quadruplex formation is possible, the matrix element is set to INF. Access the
elements in the matrix via matrix[indx[j]+i]. To get the integer array indx see get_jindx().  

Parameters
----------
S : short *
    The encoded sequence  
P : RNA.param() *
    A pointer to the data structure containing the precomputed energy contributions  

Returns
-------
int *  
    A pointer to the G-quadruplex contribution matrix  

See Also
--------
get_jindx(), encode_sequence()  
";

%feature("docstring") get_gquad_ali_matrix "
";

%feature("docstring") get_gquad_pf_matrix "
";

%feature("docstring") get_gquad_pf_matrix_comparative "
";

%feature("docstring") get_gquad_L_matrix "
";

%feature("docstring") vrna_gquad_mx_local_update "
";

%feature("docstring") get_gquad_pattern_mfe "
";

%feature("docstring") get_gquad_pattern_exhaustive "
";

%feature("docstring") get_gquad_pattern_pf "
";

%feature("docstring") get_plist_gquad_from_pr "
";

%feature("docstring") get_plist_gquad_from_pr_max "
";

%feature("docstring") get_plist_gquad_from_db "
";

%feature("docstring") vrna_get_plist_gquad_from_pr "
";

%feature("docstring") vrna_get_plist_gquad_from_pr_max "
";

%feature("docstring") get_gquad_count "
";

%feature("docstring") get_gquad_layer_count "
";

%feature("docstring") get_gquad_pattern_mfe_ali "
";

%feature("docstring") parse_gquad "

given a dot-bracket structure (possibly) containing gquads encoded by '+' signs, find first gquad,
return end position or 0 if none found Upon return L and l[] contain the number of stacked layers,
as well as the lengths of the linker regions. To parse a string with many gquads, call parse_gquad
repeatedly e.g. end1 = parse_gquad(struc, &L, l); ... ; end2 = parse_gquad(struc+end1, &L, l);
end2+=end1; ... ; end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ;  
";

%feature("docstring") backtrack_GQuad_IntLoop "

backtrack an interior loop like enclosed g-quadruplex with closing pair (i,j)  

Parameters
----------
c : int
    The total contribution the loop should resemble  
i : int
    position i of enclosing pair  
j : int
    position j of enclosing pair  
type : int
    base pair type of enclosing pair (must be reverse type)  
S : short *
    integer encoded sequence  
ggg : int *
    triangular matrix containing g-quadruplex contributions  
index : int *
    the index for accessing the triangular matrix  
p : int *
    here the 5' position of the gquad is stored  
q : int *
    here the 3' position of the gquad is stored  
P : RNA.param() *
    the datastructure containing the precalculated contibutions  

Returns
-------
int  
    1 on success, 0 if no gquad found  
";

%feature("docstring") backtrack_GQuad_IntLoop_comparative "
";

%feature("docstring") backtrack_GQuad_IntLoop_L "

backtrack an interior loop like enclosed g-quadruplex with closing pair (i,j) with underlying Lfold
matrix  

Parameters
----------
c : int
    The total contribution the loop should resemble  
i : int
    position i of enclosing pair  
j : int
    position j of enclosing pair  
type : int
    base pair type of enclosing pair (must be reverse type)  
S : short *
    integer encoded sequence  
ggg : int **
    triangular matrix containing g-quadruplex contributions  
p : int *
    here the 5' position of the gquad is stored  
q : int *
    here the 3' position of the gquad is stored  
P : RNA.param() *
    the datastructure containing the precalculated contibutions  

Returns
-------
int  
    1 on success, 0 if no gquad found  
";

%feature("docstring") vrna_BT_gquad_int "
";

%feature("docstring") vrna_BT_gquad_mfe "
";

%feature("docstring") backtrack_GQuad_IntLoop_L_comparative "
";

%feature("docstring") E_GQuad_IntLoop "
";

%feature("docstring") E_GQuad_IntLoop_comparative "
";

%feature("docstring") E_GQuad_IntLoop_L_comparative "
";

%feature("docstring") E_GQuad_IntLoop_exhaustive "
";

%feature("docstring") E_GQuad_IntLoop_L "
";

%feature("docstring") exp_E_GQuad_IntLoop "
";

%feature("docstring") exp_E_GQuad_IntLoop_comparative "
";

// File: group__modified__bases.xml

%feature("docstring") my_sc_mod_read_from_jsonfile "

Parse and extract energy parameters for a modified base from a JSON file.  

**SWIG Wrapper Notes**
    This function is available as an overloaded function `sc_mod_read_from_jsonfile()` where the
    `md` parameter may be omitted and defaults to `NULL`. See, e.g.
    :py:func:`RNA.sc_mod_read_from_jsonfile()` in the :doc:`/api_python` .  

Parameters
----------
filename : const char *
    The JSON file containing the specifications of the modified base  
md : RNA.md() *
    A model-details data structure (for look-up of canonical base pairs)  

Returns
-------
RNA.sc_mod_param()  
    Parameters of the modified base  

See Also
--------
RNA.sc_mod_read_from_json(), RNA.sc_mod_parameters_free(), RNA.fold_compound.sc_mod(), modified-
bases-params  
";

%feature("docstring") my_sc_mod_read_from_json "

Parse and extract energy parameters for a modified base from a JSON string.  

**SWIG Wrapper Notes**
    This function is available as an overloaded function `sc_mod_read_from_json()` where the `md`
    parameter may be omitted and defaults to `NULL`. See, e.g.
    :py:func:`RNA.sc_mod_read_from_json()` in the :doc:`/api_python` .  

Parameters
----------
filename :
    The JSON file containing the specifications of the modified base  
md : RNA.md() *
    A model-details data structure (for look-up of canonical base pairs)  

Returns
-------
RNA.sc_mod_param()  
    Parameters of the modified base  

See Also
--------
RNA.sc_mod_read_from_jsonfile(), RNA.sc_mod_parameters_free(), RNA.fold_compound.sc_mod(), modified-
bases-
params  
";

%feature("docstring") vrna_sc_mod_parameters_free "

Release memory occupied by a modified base parameter data structure.  

Properly free a RNA.sc_mod_param() data structure  

Parameters
----------
params : RNA.sc_mod_param()
    The data structure to free  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_json "

Prepare soft constraint callbacks for modified base as specified in JSON string.  

This function prepares all requirements to acknowledge modified bases as specified in the provided
`json` string. All subsequent predictions will treat each modification site special and adjust
energy contributions if necessary.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_json()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_json()` in the :doc:`/api_python` .  

Parameters
----------
json : const char *
    The JSON formatted string with the modified base parameters  
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod(), RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(),
RNA.fold_compound.sc_mod_inosine(), RNA.fold_compound.sc_mod_7DA(),
RNA.fold_compound.sc_mod_purine(), RNA.fold_compound.sc_mod_dihydrouridine(),
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT,
modified-bases-params  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_jsonfile "

Prepare soft constraint callbacks for modified base as specified in JSON string.  

Similar to RNA.fold_compound.sc_mod_json(), this function prepares all requirements to acknowledge modified bases
as specified in the provided `json` file. All subsequent predictions will treat each modification
site special and adjust energy contributions if necessary.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_jsonfile()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_jsonfile()` in the :doc:`/api_python` .  

Parameters
----------
json :
    The JSON formatted string with the modified base parameters  
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.fold_compound.sc_mod_json(), RNA.fold_compound.sc_mod(), RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(),
RNA.fold_compound.sc_mod_inosine(), RNA.fold_compound.sc_mod_7DA(),
RNA.fold_compound.sc_mod_purine(), RNA.fold_compound.sc_mod_dihydrouridine(),
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT,
modified-bases-params  
";

%feature("docstring") vrna_fold_compound_t::sc_mod "

Prepare soft constraint callbacks for modified base as specified in JSON string.  

This function takes a RNA.sc_mod_param() data structure as obtained from
RNA.sc_mod_read_from_json() or RNA.sc_mod_read_from_jsonfile() and prepares all requirements to
acknowledge modified bases as specified in the provided `params` data structure. All subsequent
predictions will treat each modification site special and adjust energy contributions if necessary.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod()` to objects of type `fold_compound`
    with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.   :py:meth:`RNA.fold_compound.sc_mod()`
    in the :doc:`/api_python` .  

Parameters
----------
json :
    The JSON formatted string with the modified base parameters  
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.sc_mod_read_from_json(), RNA.sc_mod_read_from_jsonfile(), RNA.fold_compound.sc_mod_json(),
RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(), RNA.fold_compound.sc_mod_inosine(),
RNA.fold_compound.sc_mod_7DA(), RNA.fold_compound.sc_mod_purine(),
RNA.sc_mod_dihydrouridine()RNA.SC_MOD_CHECK_FALLBACK,
RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_m6A "

Add soft constraint callbacks for N6-methyl-adenosine (m6A)  

This is a convenience wrapper to add support for m6A using the soft constraint callback mechanism.
Modification sites are provided as a list of sequence positions (1-based). Energy parameter
corrections are derived from   :cite:t:`kierzek:2022` .  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_m6A()` to objects of type `fold_compound`
    with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_m6A()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_pseudouridine "

Add soft constraint callbacks for Pseudouridine.  

This is a convenience wrapper to add support for pseudouridine using the soft constraint callback
mechanism. Modification sites are provided as a list of sequence positions (1-based). Energy
parameter corrections are derived from   :cite:t:`hudson:2013` .  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_pseudouridine()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_pseudouridine()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_inosine "

Add soft constraint callbacks for Inosine.  

This is a convenience wrapper to add support for inosine using the soft constraint callback
mechanism. Modification sites are provided as a list of sequence positions (1-based). Energy
parameter corrections are derived from   :cite:t:`wright:2007`  and   :cite:t:`wright:2018` .  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_inosine()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_inosine()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_7DA "

Add soft constraint callbacks for 7-deaza-adenosine (7DA)  

This is a convenience wrapper to add support for 7-deaza-adenosine using the soft constraint
callback mechanism. Modification sites are provided as a list of sequence positions (1-based).
Energy parameter corrections are derived from   :cite:t:`richardson:2016` .  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_7DA()` to objects of type `fold_compound`
    with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_7DA()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_purine "

Add soft constraint callbacks for Purine (a.k.a. nebularine)  

This is a convenience wrapper to add support for Purine using the soft constraint callback
mechanism. Modification sites are provided as a list of sequence positions (1-based). Energy
parameter corrections are derived from   :cite:t:`jolley:2017` .  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_purine()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_purine()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") vrna_fold_compound_t::sc_mod_dihydrouridine "

Add soft constraint callbacks for dihydrouridine.  

This is a convenience wrapper to add support for dihydrouridine using the soft constraint callback
mechanism. Modification sites are provided as a list of sequence positions (1-based). Energy
parameter corrections are derived from Rosetta/RECESS predictions.  

**SWIG Wrapper Notes**
    This function is attached as overloaded method `sc_mod_dihydrouridine()` to objects of type
    `fold_compound` with default `options` = RNA.SC_MOD_DEFAULT. See, e.g.
    :py:meth:`RNA.fold_compound.sc_mod_dihydrouridine()` in the :doc:`/api_python` .  

Parameters
----------
modification_sites : const unsigned int *
    A list of modification site, i.e. positions that contain the modified base (1-based, last
    element in the list indicated by 0)  
options : unsigned int
    A bitvector of options how to handle the input, e.g. RNA.SC_MOD_DEFAULT  

Returns
-------
int  
    Number of sequence positions modified base parameters will be used for  

See Also
--------
RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_SILENT, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") SC_MOD_CHECK_FALLBACK "

Check for sequence positions whether they resemble the fallback base.  

This flag can be used to enable a sanity check within the RNA.sc_mod*() functions to see whether a
supposedly modified position actually resembles the fallback base as specified in the modification
parameters  

See Also
--------
RNA.fold_compound.sc_mod_json(), RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod(),
RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(), RNA.fold_compound.sc_mod_inosine(),
RNA.fold_compound.sc_mod_7DA(), RNA.fold_compound.sc_mod_purine(),
RNA.fold_compound.sc_mod_dihydrouridine(), RNA.SC_MOD_CHECK_UNMOD, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") SC_MOD_CHECK_UNMOD "

Check for sequence positions whether they resemble the unmodified base.  

This flag can be used to enable a sanity check within the RNA.sc_mod*() functions to see whether a
supposedly modified position actually resembles the unmodified base as specified in the modification
parameters  

See Also
--------
RNA.fold_compound.sc_mod_json(), RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod(),
RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(), RNA.fold_compound.sc_mod_inosine(),
RNA.fold_compound.sc_mod_7DA(), RNA.fold_compound.sc_mod_purine(),
RNA.fold_compound.sc_mod_dihydrouridine(), RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_DEFAULT  
";

%feature("docstring") SC_MOD_SILENT "

Do not produce any warnings within the RNA.sc_mod*() functions.  

See Also
--------
RNA.fold_compound.sc_mod_json(), RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod(),
RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(), RNA.fold_compound.sc_mod_inosine(),
RNA.fold_compound.sc_mod_7DA(), RNA.fold_compound.sc_mod_purine(),
RNA.fold_compound.sc_mod_dihydrouridine()  
";

%feature("docstring") SC_MOD_DEFAULT "

Default settings for the RNA.sc_mod*() functions.  

See Also
--------
RNA.fold_compound.sc_mod_json(), RNA.fold_compound.sc_mod_jsonfile(), RNA.fold_compound.sc_mod(),
RNA.fold_compound.sc_mod_m6A(),
RNA.fold_compound.sc_mod_pseudouridine(), RNA.fold_compound.sc_mod_inosine(),
RNA.fold_compound.sc_mod_7DA(), RNA.fold_compound.sc_mod_purine(),
RNA.fold_compound.sc_mod_dihydrouridine(), RNA.SC_MOD_CHECK_FALLBACK, RNA.SC_MOD_CHECK_UNMOD,
RNA.SC_MOD_SILENT  
";

// File: group__utils.xml

%feature("docstring") vrna_alloc "

Allocate space safely.  

Parameters
----------
size : unsigned
    The size of the memory to be allocated in bytes  

Returns
-------
void *  
    A pointer to the allocated memory  
";

%feature("docstring") vrna_realloc "

Reallocate space safely.  

Parameters
----------
p : void *
    A pointer to the memory region to be reallocated  
size : unsigned
    The size of the memory to be allocated in bytes  

Returns
-------
void *  
    A pointer to the newly allocated memory  
";

%feature("docstring") vrna_init_rand "

Initialize seed for random number generator.  

See Also
--------
RNA.init_rand_seed(), RNA.urn()  
";

%feature("docstring") vrna_init_rand_seed "

Initialize the random number generator with a pre-defined seed.  

**SWIG Wrapper Notes**
    This function is available as an overloaded function **init_rand()** where the argument `seed`
    is optional. See, e.g.  :py:func:`RNA.init_rand()` in the :doc:`/api_python`.  

Parameters
----------
seed : unsigned int
    The seed for the random number generator  

See Also
--------
RNA.init_rand(), RNA.urn()  
";

%feature("docstring") vrna_urn "

get a random number from [0..1]  

Returns
-------
double  
    A random number in range [0..1]  

See Also
--------
RNA.int_urn(), RNA.init_rand(), RNA.init_rand_seed()  

Note
----
Usually implemented by calling *erand48()*.  
";

%feature("docstring") vrna_int_urn "

Generates a pseudo random integer in a specified range.  

Parameters
----------
from : int
    The first number in range  
to : int
    The last number in range  

Returns
-------
int  
    A pseudo random number in range [from, to]  

See Also
--------
RNA.urn(), RNA.init_rand()  
";

%feature("docstring") vrna_time_stamp "

Get a timestamp.  

Returns a string containing the current date in the format  

    Fri Mar 19 21:10:57 1993  

Returns
-------
char *  
    A string containing the timestamp  
";

%feature("docstring") get_input_line "

Retrieve a line from 'stdin' savely while skipping comment characters and other features This
function returns the type of input it has read if recognized. An option argument allows one to
switch between different reading modes.  
Currently available options are:  RNA.INPUT_COMMENT, RNA.INPUT_NOSKIP_COMMENTS,
RNA.INPUT_NO_TRUNCATION  

pass a collection of options as one value like this:  

    get_input_line(string, option_1 | option_2 | option_n)  

If the function recognizes the type of input, it will report it in the return value. It also reports
if a user defined 'quit' command (-sign on 'stdin') was given. Possible return values are:
RNA.INPUT_FASTA_HEADER, RNA.INPUT_ERROR, RNA.INPUT_MISC, RNA.INPUT_QUIT  

Parameters
----------
string : char **
    A pointer to the character array that contains the line read  
options : unsigned int
    A collection of options for switching the functions behavior  

Returns
-------
unsigned int  
    A flag with information about what has been read  
";

%feature("docstring") vrna_idx_row_wise "

Get an index mapper array (iindx) for accessing the energy matrices, e.g. in partition function
related functions.  

Access of a position \"(i,j)\" is then accomplished by using  

    (i,j) ~ iindx[i]-j  This function is necessary as most of the two-dimensional energy matrices
are actually one-dimensional arrays throughout the ViennaRNA Package  

Consult the implemented code to find out about the mapping formula ;)  

Parameters
----------
length : unsigned int
    The length of the RNA sequence  

Returns
-------
int *  
    The mapper array  

See Also
--------
RNA.idx_col_wise()  
";

%feature("docstring") vrna_idx_col_wise "

Get an index mapper array (indx) for accessing the energy matrices, e.g. in MFE related functions.  

Access of a position \"(i,j)\" is then accomplished by using  

    (i,j) ~ indx[j]+i  This function is necessary as most of the two-dimensional energy matrices are
actually one-dimensional arrays throughout the ViennaRNAPackage  

Consult the implemented code to find out about the mapping formula ;)  

Parameters
----------
length : unsigned int
    The length of the RNA sequence  

Returns
-------
int *  
    The mapper array  

See Also
--------
RNA.idx_row_wise()  
";

%feature("docstring") PUBLIC "
";

%feature("docstring") PRIVATE "
";

%feature("docstring") INPUT_ERROR "

Output flag of get_input_line(): *\"An ERROR has occured, maybe EOF\"*.  
";

%feature("docstring") INPUT_QUIT "

Output flag of get_input_line(): *\"the user requested quitting the program\"*.  
";

%feature("docstring") INPUT_MISC "

Output flag of get_input_line(): *\"something was read\"*.  
";

%feature("docstring") INPUT_FASTA_HEADER "

Input/Output flag of get_input_line():  
if used as input option this tells get_input_line() that the data to be read should comply with the
FASTA format.  

the function will return this flag if a fasta header was read  
";

%feature("docstring") INPUT_SEQUENCE "
";

%feature("docstring") INPUT_CONSTRAINT "

Input flag for get_input_line():  
Tell get_input_line() that we assume to read a structure constraint.  
";

%feature("docstring") INPUT_NO_TRUNCATION "

Input switch for get_input_line(): *\"do not trunkate the line by eliminating white spaces at end of
line\"*.  
";

%feature("docstring") INPUT_NO_REST "

Input switch for RNA.file_fasta_read_record(): *\"do fill rest array\"*.  
";

%feature("docstring") INPUT_NO_SPAN "

Input switch for RNA.file_fasta_read_record(): *\"never allow data to span more than one line\"*.  
";

%feature("docstring") INPUT_NOSKIP_BLANK_LINES "

Input switch for RNA.file_fasta_read_record(): *\"do not skip empty lines\"*.  
";

%feature("docstring") INPUT_BLANK_LINE "

Output flag for RNA.file_fasta_read_record(): *\"read an empty line\"*.  
";

%feature("docstring") INPUT_NOSKIP_COMMENTS "

Input switch for get_input_line(): *\"do not skip comment lines\"*.  
";

%feature("docstring") INPUT_COMMENT "

Output flag for RNA.file_fasta_read_record(): *\"read a comment\"*.  
";

%feature("docstring") MIN2 "

Get the minimum of two comparable values.  
";

%feature("docstring") MAX2 "

Get the maximum of two comparable values.  
";

%feature("docstring") MIN3 "

Get the minimum of three comparable values.  
";

%feature("docstring") MAX3 "

Get the maximum of three comparable values.  
";

// File: group__eval__loops__ext.xml

%feature("docstring") vrna_fold_compound_t::exp_E_ext_stem "

Evaluate a stem branching off the exterior loop (Boltzmann factor version)  

Given a base pair :math:`(i,j)` encoded by *type*, compute the energy contribution including
dangling-end/terminal-mismatch contributions. Instead of returning the energy contribution per-se,
this function returns the corresponding Boltzmann factor. If either of the adjacent nucleotides
:math:`(i - 1)` and :math:`(j+1)` must not contribute stacking energy, the corresponding encoding
must be :math:`-1`.  

Parameters
----------
type : unsigned int
    The base pair encoding  
n5d : int
    The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)  
n3d : int
    The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)  
p : RNA.exp_param() *
    The pre-computed energy parameters (Boltzmann factor version)  

Returns
-------
FLT_OR_DBL  
    The Boltzmann weighted energy contribution of the introduced exterior-loop stem  

See Also
--------
RNA.E_ext_stem()  
";

%feature("docstring") vrna_exp_E_ext_fast_init "
";

%feature("docstring") vrna_exp_E_ext_fast_rotate "
";

%feature("docstring") vrna_exp_E_ext_fast_free "
";

%feature("docstring") vrna_exp_E_ext_fast "
";

%feature("docstring") vrna_exp_E_ext_fast_update "
";

%feature("docstring") vrna_E_ext_stem "

Evaluate a stem branching off the exterior loop.  

Given a base pair :math:`(i,j)` encoded by *type*, compute the energy contribution including
dangling-end/terminal-mismatch contributions. Instead of returning the energy contribution per-se,
this function returns the corresponding Boltzmann factor. If either of the adjacent nucleotides
:math:`(i - 1)` and :math:`(j+1)` must not contribute stacking energy, the corresponding encoding
must be :math:`-1`.  

Parameters
----------
type : unsigned int
    The base pair encoding  
n5d : int
    The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)  
n3d : int
    The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)  
p : RNA.param() *
    The pre-computed energy parameters  

Returns
-------
int  
    The energy contribution of the introduced exterior-loop stem  

See Also
--------
RNA.E_exp_stem()  
";

%feature("docstring") vrna_fold_compound_t::eval_ext_stem "

Evaluate the free energy of a base pair in the exterior loop.  

Evalue the free energy of a base pair connecting two nucleotides in the exterior loop and take hard
constraints into account.  

Typically, this is simply dangling end contributions of the adjacent nucleotides, potentially a
terminal A-U mismatch penalty, and maybe some generic soft constraint contribution for that
decomposition.  

Parameters
----------
i : int
    5' position of the base pair  
j : int
    3' position of the base pair  

Returns
-------
int  
    Free energy contribution that arises when this pair is formed in the exterior loop  

Note
----
For dangles == 1 || 3 this function also evaluates the three additional pairs (i + 1, j), (i, j -
1), and (i + 1, j - 1) and returns the minimum for all four possibilities in total.  
";

%feature("docstring") vrna_E_ext_loop_5 "
";

%feature("docstring") vrna_E_ext_loop_3 "
";

// File: group__eval__loops__hp.xml

%feature("docstring") vrna_fold_compound_t::E_hp_loop "

Evaluate the free energy of a hairpin loop and consider hard constraints if they apply.  

This function evaluates the free energy of a hairpin loop  

In case the base pair is not allowed due to a constraint conflict, this function returns INF.  

Parameters
----------
i : int
    The 5' nucleotide of the base pair (3' to evaluate the pair as exterior hairpin loop)  
j : int
    The 3' nucleotide of the base pair (5' to evaluate the pair as exterior hairpin loop)  

Returns
-------
int  
    The free energy of the hairpin loop in 10cal/mol  

Note
----
This function is polymorphic! The provided RNA.fold_compound() may be of type RNA.FC_TYPE_SINGLE
or RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_fold_compound_t::E_ext_hp_loop "

Evaluate the free energy of an exterior hairpin loop and consider possible hard constraints.  

Note
----
This function is polymorphic! The provided RNA.fold_compound() may be of type RNA.FC_TYPE_SINGLE
or RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_fold_compound_t::eval_ext_hp_loop "

Evaluate free energy of an exterior hairpin loop.  
";

%feature("docstring") vrna_fold_compound_t::eval_hp_loop "

Evaluate free energy of a hairpin loop.  

**SWIG Wrapper Notes**
    This function is attached as method `eval_hp_loop()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.eval_hp_loop()` in the :doc:`/api_python` .  

Parameters
----------
i : int
    5'-position of the base pair  
j : int
    3'-position of the base pair  

Returns
-------
int  
    Free energy of the hairpin loop closed by :math:`(i,j)` in deka-kal/mol  

Note
----
This function is polymorphic! The provided RNA.fold_compound() may be of type RNA.FC_TYPE_SINGLE
or RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") E_Hairpin "

Compute the Energy of a hairpin-loop.  

To evaluate the free energy of a hairpin-loop, several parameters have to be known. A general
hairpin-loop has this structure:  
      a3 a4
    a2     a5
    a1     a6
      X - Y
      |   |
      5'  3'
 where X-Y marks the closing pair [e.g. a **(G,C)** pair]. The length of this loop is 6 as there are
six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is a1 while the 3'
mismatch is a6. The nucleotide sequence of this loop is \"a1.a2.a3.a4.a5.a6\"  

Parameters
----------
size : int
    The size of the loop (number of unpaired nucleotides)  
type : int
    The pair type of the base pair closing the hairpin  
si1 : int
    The 5'-mismatching nucleotide  
sj1 : int
    The 3'-mismatching nucleotide  
string : const char *
    The sequence of the loop (May be `NULL`, otherwise mst be at least :math:`size + 2` long)  
P : RNA.param() *
    The datastructure containing scaled energy parameters  

Returns
-------
int  
    The Free energy of the Hairpin-loop in dcal/mol  

Warnings
--------
Not (really) thread safe! A threadsafe implementation will replace this function in a future
release!  
 Energy evaluation may change due to updates in global variable \"tetra_loop\"  

See Also
--------
scale_parameters(), RNA.param()  

Note
----
The parameter sequence should contain the sequence of the loop in capital letters of the nucleic
acid alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and
hexa-loops which are treated differently (based on experimental data) if they are tabulated.  
";

%feature("docstring") exp_E_Hairpin "

Compute Boltzmann weight :math:`e^{-\\Delta G/kT}` of a hairpin loop.  

Parameters
----------
u : int
    The size of the loop (number of unpaired nucleotides)  
type : int
    The pair type of the base pair closing the hairpin  
si1 : short
    The 5'-mismatching nucleotide  
sj1 : short
    The 3'-mismatching nucleotide  
string : const char *
    The sequence of the loop (May be `NULL`, otherwise mst be at least :math:`size + 2` long)  
P : RNA.exp_param() *
    The datastructure containing scaled Boltzmann weights of the energy parameters  

Returns
-------
FLT_OR_DBL  
    The Boltzmann weight of the Hairpin-loop  

Warnings
--------
Not (really) thread safe! A threadsafe implementation will replace this function in a future
release!  
 Energy evaluation may change due to updates in global variable \"tetra_loop\"  

See Also
--------
get_scaled_pf_parameters(), RNA.exp_param(), E_Hairpin()  

Note
----
multiply by scale[u+2]  
";

%feature("docstring") vrna_fold_compound_t::exp_E_hp_loop "

High-Level function for hairpin loop energy evaluation (partition function variant)  

See Also
--------
RNA.fold_compound.E_hp_loop() for it's free energy counterpart  

Note
----
This function is polymorphic! The provided RNA.fold_compound() may be of type RNA.FC_TYPE_SINGLE
or RNA.FC_TYPE_COMPARATIVE  
";

// File: group__eval__loops__int.xml

%feature("docstring") vrna_fold_compound_t::E_int_loop "
";

%feature("docstring") vrna_fold_compound_t::eval_int_loop "

Evaluate the free energy contribution of an interior loop with delimiting base pairs :math:`(i,j)`
and :math:`(k,l)`.  

**SWIG Wrapper Notes**
    This function is attached as method `eval_int_loop()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.eval_int_loop()` in the :doc:`/api_python` .  

Note
----
This function is polymorphic, i.e. it accepts RNA.fold_compound() of type RNA.FC_TYPE_SINGLE as
well as RNA.FC_TYPE_COMPARATIVE  
";

%feature("docstring") vrna_fold_compound_t::E_ext_int_loop "
";

%feature("docstring") vrna_fold_compound_t::E_stack "
";

%feature("docstring") vrna_fold_compound_t::exp_E_int_loop "
";

%feature("docstring") vrna_fold_compound_t::exp_E_interior_loop "
";

// File: group__eval__loops__mb.xml

%feature("docstring") vrna_exp_E_mb_loop_fast "
";

%feature("docstring") vrna_exp_E_ml_fast_init "
";

%feature("docstring") vrna_exp_E_ml_fast_rotate "
";

%feature("docstring") vrna_exp_E_ml_fast_free "
";

%feature("docstring") vrna_exp_E_ml_fast_qqm "
";

%feature("docstring") vrna_exp_E_ml_fast_qqm1 "
";

%feature("docstring") vrna_exp_E_ml_fast "
";

%feature("docstring") vrna_E_mb_loop_stack "

Evaluate energy of a multi branch helices stacking onto closing pair (i,j)  

Computes total free energy for coaxial stacking of (i.j) with (i+1.k) or (k+1.j-1)  
";

%feature("docstring") vrna_E_mb_loop_fast "
";

%feature("docstring") E_ml_rightmost_stem "
";

%feature("docstring") vrna_E_ml_stems_fast "
";

// File: group__pf__cofold.xml

%feature("docstring") vrna_pf_dimer_concentrations "

Given two start monomer concentrations a and b, compute the concentrations in thermodynamic
equilibrium of all dimers and the monomers.  

This function takes an array 'startconc' of input concentrations with alternating entries for the
initial concentrations of molecules A and B (terminated by two zeroes), then computes the resulting
equilibrium concentrations from the free energies for the dimers. Dimer free energies should be the
dimer-only free energies, i.e. the FcAB entries from the RNA.dimer_pf() struct.  

Parameters
----------
FcAB : double
    Free energy of AB dimer (FcAB entry)  
FcAA : double
    Free energy of AA dimer (FcAB entry)  
FcBB : double
    Free energy of BB dimer (FcAB entry)  
FEA : double
    Free energy of monomer A  
FEB : double
    Free energy of monomer B  
startconc : const double *
    List of start concentrations [a0],[b0],[a1],[b1],...,[an][bn],[0],[0]  
exp_params : const RNA.exp_param() *
    The precomputed Boltzmann factors  

Returns
-------
RNA.dimer_conc() *  
    RNA.dimer_conc() array containing the equilibrium energies and start concentrations  
";

%feature("docstring") vrna_equilibrium_constants "
";

%feature("docstring") vrna_pf_co_fold "

Calculate partition function and base pair probabilities of nucleic acid/nucleic acid dimers.  

This simplified interface to RNA.fold_compound.pf_dimer() computes the partition function and, if required, base
pair probabilities for an RNA-RNA interaction using default options. Memory required for dynamic
programming (DP) matrices will be allocated and free'd on-the-fly. Hence, after return of this
function, the recursively filled matrices are not available any more for any post-processing.  

Parameters
----------
seq : const char *
    Two concatenated RNA sequences with a delimiting '&' in between  
structure : char *
    A pointer to the character array where position-wise pairing propensity will be stored. (Maybe
    NULL)  
pl : RNA.ep() **
    A pointer to a list of RNA.ep() to store pairing probabilities (Maybe NULL)  

Returns
-------
RNA.dimer_pf()  
    RNA.dimer_pf() structure containing a set of energies needed for concentration computations.  

See Also
--------
RNA.fold_compound.pf_dimer()  

Note
----
In case you want to use the filled DP matrices for any subsequent post-processing step, or you
require other conditions than specified by the default model details, use
RNA.fold_compound.pf_dimer(), and the
data structure RNA.fold_compound() instead.  
";

// File: group__up__cofold.xml

%feature("docstring") pf_unstru "

Calculate the partition function over all unpaired regions of a maximal length.  

You have to call function pf_fold() providing the same sequence before calling pf_unstru(). If you
want to calculate unpaired regions for a constrained structure, set variable 'structure' in function
'pf_fold()' to the constrain string. It returns a pu_contrib struct containing four arrays of
dimension [i = 1 to length(sequence)][j = 0 to u-1] containing all possible contributions to the
probabilities of unpaired regions of maximum length u. Each array in pu_contrib contains one of the
contributions to the total probability of being unpaired: The probability of being unpaired within
an exterior loop is in array pu_contrib->E, the probability of being unpaired within a hairpin loop
is in array pu_contrib->H, the probability of being unpaired within an interior loop is in array
pu_contrib->I and probability of being unpaired within a multi-loop is in array pu_contrib->M. The
total probability of being unpaired is the sum of the four arrays of pu_contrib.  

This function frees everything allocated automatically. To free the output structure call
free_pu_contrib().  

Parameters
----------
sequence : char *
max_w : int

Returns
-------
pu_contrib *  
";

%feature("docstring") pf_interact "

Calculates the probability of a local interaction between two sequences.  

The function considers the probability that the region of interaction is unpaired within 's1' and
's2'. The longer sequence has to be given as 's1'. The shorter sequence has to be given as 's2'.
Function pf_unstru() has to be called for 's1' and 's2', where the probabilities of being unpaired
have to be given in 'p_c' and 'p_c2', respectively. If you do not want to include the probabilities
of being unpaired for 's2' set 'p_c2' to NULL. If variable 'cstruc' is not NULL, constrained folding
is done: The available constrains for intermolecular interaction are: '.' (no constrain), 'x' (the
base has no intermolecular interaction) and '|' (the corresponding base has to be paired
intermolecularily).  
The parameter 'w' determines the maximal length of the interaction. The parameters 'incr5' and
'incr3' allows inclusion of unpaired residues left ('incr5') and right ('incr3') of the region of
interaction in 's1'. If the 'incr' options are used, function pf_unstru() has to be called with
w=w+incr5+incr3 for the longer sequence 's1'.  

It returns a structure of type interact which contains the probability of the best local interaction
including residue i in Pi and the minimum free energy in Gi, where i is the position in sequence
's1'. The member Gikjl of structure interact is the best interaction between region [k,i] k<i in
longer sequence 's1' and region [j,l] j<l in 's2'. Gikjl_wo is Gikjl without the probability of
beeing unpaired.  
Use free_interact() to free the returned structure, all other stuff is freed inside pf_interact().  

Parameters
----------
s1 : const char *
s2 : const char *
p_c : pu_contrib *
p_c2 : pu_contrib *
max_w : int
cstruc : char *
incr3 : int
incr5 : int

Returns
-------
interact *  
";

%feature("docstring") free_interact "

Frees the output of function pf_interact().  
";

%feature("docstring") Up_plot "
";

%feature("docstring") get_pu_contrib_struct "
";

%feature("docstring") free_pu_contrib_struct "

Frees the output of function pf_unstru().  
";

%feature("docstring") free_pu_contrib "
";

// File: group__energy__parameters__rw.xml

%feature("docstring") my_params_load "

Load energy parameters from a file.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `params_load`(fname=\"\",
    options=RNA.PARAMETER_FORMAT_DEFAULT). Here, the empty filename string indicates to load
    default RNA parameters, i.e. this is equivalent to calling RNA.params_load_defaults(). See,
    e.g.  :py:func:`RNA.fold_compound.params_load()` in the :doc:`/api_python`.  

Parameters
----------
fname : const char
    The path to the file containing the energy parameters  
options : unsigned int
    File format bit-mask (usually RNA.PARAMETER_FORMAT_DEFAULT)  

Returns
-------
int  
    Non-zero on success, 0 on failure  

See Also
--------
RNA.params_load_from_string(), RNA.params_save(), RNA.params_load_defaults(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") my_params_save "

Save energy parameters to a file.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `params_save`(fname,
    options=RNA.PARAMETER_FORMAT_DEFAULT). See, e.g.  :py:func:`RNA.params_save()` in the
    :doc:`/api_python`.  

Parameters
----------
fname : const char
    A filename (path) for the file where the current energy parameters will be written to  
options : unsigned int
    File format bit-mask (usually RNA.PARAMETER_FORMAT_DEFAULT)  

Returns
-------
int  
    Non-zero on success, 0 on failure  

See Also
--------
RNA.params_load()  
";

%feature("docstring") my_params_load_from_string "

Load energy paramters from string.  

The string must follow the default energy parameter file convention! The optional `name` argument
allows one to specify a name for the parameter set which is stored internally.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `params_load_from_string`(string, name=\"\",
    options=RNA.PARAMETER_FORMAT_DEFAULT). See, e.g.  :py:func:`RNA.params_load_from_string()` in
    the :doc:`/api_python`.  

Parameters
----------
string : const char *
    A 0-terminated string containing energy parameters  
name : const char *
    A name for the parameter set in `string` (Maybe `NULL`)  
options : unsigned int
    File format bit-mask (usually RNA.PARAMETER_FORMAT_DEFAULT)  

Returns
-------
int  
    Non-zero on success, 0 on failure  

See Also
--------
RNA.params_load(), RNA.params_save(), RNA.params_load_defaults(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_defaults "

Load default RNA energy parameter set.  

This is a convenience function to load the Turner 2004 RNA free energy parameters. It's the same as
calling RNA.params_load_RNA_Turner2004()  

**SWIG Wrapper Notes**
    This function is available as overloaded function `params_load()`. See, e.g.
    :py:func:`RNA.params_load()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_RNA_Turner2004 "

Load Turner 2004 RNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_RNA_Turner2004()`. See, e.g.
    :py:func:`RNA.params_load_RNA_Turner2004()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of RNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(), RNA.params_load_defaults(),
RNA.params_load_RNA_Turner1999(), RNA.params_load_RNA_Andronescu2007(),
RNA.params_load_RNA_Langdon2018(), RNA.params_load_RNA_misc_special_hairpins(),
RNA.params_load_DNA_Mathews2004(), RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_RNA_Turner1999 "

Load Turner 1999 RNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_RNA_Turner1999()`. See, e.g.
    :py:func:`RNA.params_load_RNA_Turner1999()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of RNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_defaults(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_RNA_Andronescu2007 "

Load Andronsecu 2007 RNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_RNA_Andronescu2007()`. See, e.g.
    :py:func:`RNA.params_load_RNA_Andronescu2007()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of RNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(), RNA.params_load_defaults(),
RNA.params_load_RNA_Langdon2018(), RNA.params_load_RNA_misc_special_hairpins(),
RNA.params_load_DNA_Mathews2004(), RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_RNA_Langdon2018 "

Load Langdon 2018 RNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_RNA_Langdon2018()`. See, e.g.
    :py:func:`RNA.params_load_RNA_Langdon2018()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of RNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_defaults(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_RNA_misc_special_hairpins "

Load Misc Special Hairpin RNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_RNA_misc_special_hairpins()`. See, e.g.
    :py:func:`RNA.params_load_RNA_misc_special_hairpins()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of RNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_defaults(), RNA.params_load_DNA_Mathews2004(), RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_DNA_Mathews2004 "

Load Mathews 2004 DNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_DNA_Mathews2004()`. See, e.g.
    :py:func:`RNA.params_load_DNA_Mathews2004()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of DNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_defaults(),
RNA.params_load_DNA_Mathews1999()  
";

%feature("docstring") vrna_params_load_DNA_Mathews1999 "

Load Mathews 1999 DNA energy parameter set.  

**SWIG Wrapper Notes**
    This function is available as function `params_load_DNA_Mathews1999()`. See, e.g.
    :py:func:`RNA.params_load_DNA_Mathews1999()` in the :doc:`/api_python`.  

Returns
-------
int  
    Non-zero on success, 0 on failure  

Warnings
--------
This function also resets the default geometric parameters as stored in RNA.md() to those of DNA.
Only subsequently initialized RNA.md() structures will be affected by this change.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save(),
RNA.params_load_RNA_Turner2004(), RNA.params_load_RNA_Turner1999(),
RNA.params_load_RNA_Andronescu2007(), RNA.params_load_RNA_Langdon2018(),
RNA.params_load_RNA_misc_special_hairpins(), RNA.params_load_DNA_Mathews2004(),
RNA.params_load_defaults()  
";

%feature("docstring") last_parameter_file "

Get the file name of the parameter file that was most recently loaded.  

Returns
-------
const char *  
    The file name of the last parameter file, or NULL if parameters are still at defaults  
";

%feature("docstring") read_parameter_file "

Read energy parameters from a file.  

.. deprecated:: 2.6.4
    Use RNA.params_load() instead!  

Parameters
----------
fname : const char
    The path to the file containing the energy parameters  
";

%feature("docstring") write_parameter_file "

Write energy parameters to a file.  

.. deprecated:: 2.6.4
    Use RNA.params_save() instead!  

Parameters
----------
fname : const char
    A filename (path) for the file where the current energy parameters will be written to  
";

%feature("docstring") gettype "
";

%feature("docstring") settype "
";

%feature("docstring") PARAMETER_FORMAT_DEFAULT "

Default Energy Parameter File format.  

See Also
--------
RNA.params_load(), RNA.params_load_from_string(), RNA.params_save()  
";

// File: group__energy__parameters__salt.xml

%feature("docstring") vrna_salt_loop "

Get salt correction for a loop at a given salt concentration and temperature.  

Parameters
----------
L : int
    backbone number in loop  
salt : double
    salt concentration (M)  
T : double
    absolute temperature (K)  
backbonelen : double
    Backbone Length, phosphate-to-phosphate distance (typically 6 for RNA, 6.76 for DNA)  

Returns
-------
double  
    Salt correction for loop in dcal/mol  
";

%feature("docstring") vrna_salt_loop_int "

Get salt correction for a loop at a given salt concentration and temperature.  

This functions is same as RNA.salt_loop but returns rounded salt correction in integer  

Parameters
----------
L : int
    backbone number in loop  
salt : double
    salt concentration (M)  
T : double
    absolute temperature (K)  
backbonelen : double
    Backbone Length, phosphate-to-phosphate distance (typically 6 for RNA, 6.76 for DNA)  

Returns
-------
int  
    Rounded salt correction for loop in dcal/mol  

See Also
--------
RNA.salt_loop  
";

%feature("docstring") vrna_salt_stack "

Get salt correction for a stack at a given salt concentration and temperature.  

Parameters
----------
salt : double
    salt concentration (M)  
T : double
    absolute temperature (K)  
hrise : double
    Helical Rise (typically 2.8 for RNA, 3.4 for DNA)  

Returns
-------
int  
    Rounded salt correction for stack in dcal/mol  
";

%feature("docstring") vrna_salt_ml "

Fit linear function to loop salt correction.  

For a given range of loop size (backbone number), we perform a linear fitting on loop salt
correction  

.. math::

  \\text{Loop correction} \\approx m \\cdot L + b.  

Parameters
----------
saltLoop : double
    List of loop salt correction of size from 1  
lower : int
    Define the size lower bound for fitting  
upper : int
    Define the size upper bound for fitting  
m : int *
    pointer to store the parameter m in fitting result  
b : int *
    pointer to store the parameter b in fitting result  

See Also
--------
RNA.salt_loop()  
";

%feature("docstring") vrna_salt_duplex_init "

Get salt correction for duplex initialization at a given salt concentration.  

Parameters
----------
md : RNA.md() *
    Model details data structure that specfifies salt concentration in buffer (M)  

Returns
-------
int  
    Rounded correction for duplex initialization in dcal/mol  
";

// File: group__energy__parameters__convert.xml

%feature("docstring") convert_parameter_file "

Convert/dump a Vienna 1.8.4 formatted energy parameter file  

The options argument allows one to control the different output modes.  
Currently available options are:  RNA.CONVERT_OUTPUT_ALL, RNA.CONVERT_OUTPUT_HP,
RNA.CONVERT_OUTPUT_STACK  RNA.CONVERT_OUTPUT_MM_HP, RNA.CONVERT_OUTPUT_MM_INT,
RNA.CONVERT_OUTPUT_MM_INT_1N  RNA.CONVERT_OUTPUT_MM_INT_23, RNA.CONVERT_OUTPUT_MM_MULTI,
RNA.CONVERT_OUTPUT_MM_EXT  RNA.CONVERT_OUTPUT_DANGLE5, RNA.CONVERT_OUTPUT_DANGLE3,
RNA.CONVERT_OUTPUT_INT_11  RNA.CONVERT_OUTPUT_INT_21, RNA.CONVERT_OUTPUT_INT_22,
RNA.CONVERT_OUTPUT_BULGE  RNA.CONVERT_OUTPUT_INT, RNA.CONVERT_OUTPUT_ML, RNA.CONVERT_OUTPUT_MISC
RNA.CONVERT_OUTPUT_SPECIAL_HP, RNA.CONVERT_OUTPUT_VANILLA, RNA.CONVERT_OUTPUT_NINIO
RNA.CONVERT_OUTPUT_DUMP  

The defined options are fine for bitwise compare- and assignment-operations, e. g.: pass a
collection of options as a single value like this:  

    convert_parameter_file(ifile, ofile, option_1 | option_2 | option_n)  

Parameters
----------
iname : const char *
    The input file name (If NULL input is read from stdin)  
oname : const char *
    The output file name (If NULL output is written to stdout)  
options : unsigned int
    The options (as described above)  
";

%feature("docstring") CONVERT_OUTPUT_ALL "

Flag to indicate printing of a complete parameter set  
";

%feature("docstring") CONVERT_OUTPUT_HP "

Flag to indicate printing of hairpin contributions  
";

%feature("docstring") CONVERT_OUTPUT_STACK "

Flag to indicate printing of base pair stack contributions  
";

%feature("docstring") CONVERT_OUTPUT_MM_HP "

Flag to indicate printing of hairpin mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_MM_INT "

Flag to indicate printing of interior loop mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_MM_INT_1N "

Flag to indicate printing of 1:n interior loop mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_MM_INT_23 "

Flag to indicate printing of 2:3 interior loop mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_MM_MULTI "

Flag to indicate printing of multi loop mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_MM_EXT "

Flag to indicate printing of exterior loop mismatch contribution  
";

%feature("docstring") CONVERT_OUTPUT_DANGLE5 "

Flag to indicate printing of 5' dangle conctribution  
";

%feature("docstring") CONVERT_OUTPUT_DANGLE3 "

Flag to indicate printing of 3' dangle contribution  
";

%feature("docstring") CONVERT_OUTPUT_INT_11 "

Flag to indicate printing of 1:1 interior loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_INT_21 "

Flag to indicate printing of 2:1 interior loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_INT_22 "

Flag to indicate printing of 2:2 interior loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_BULGE "

Flag to indicate printing of bulge loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_INT "

Flag to indicate printing of interior loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_ML "

Flag to indicate printing of multi loop contribution  
";

%feature("docstring") CONVERT_OUTPUT_MISC "

Flag to indicate printing of misc contributions (such as terminalAU)  
";

%feature("docstring") CONVERT_OUTPUT_SPECIAL_HP "

Flag to indicate printing of special hairpin contributions (tri-, tetra-, hexa-loops)  
";

%feature("docstring") CONVERT_OUTPUT_VANILLA "

Flag to indicate printing of given parameters only  

Note
----
This option overrides all other output options, except RNA.CONVERT_OUTPUT_DUMP !  
";

%feature("docstring") CONVERT_OUTPUT_NINIO "

Flag to indicate printing of interior loop asymmetry contribution  
";

%feature("docstring") CONVERT_OUTPUT_DUMP "

Flag to indicate dumping the energy contributions from the library instead of an input file  
";

// File: group__alphabet__utils.xml

%feature("docstring") vrna_sequence_length_max "
";

%feature("docstring") vrna_nucleotide_IUPAC_identity "
";

%feature("docstring") vrna_ptypes_prepare "
";

%feature("docstring") vrna_ptypes "

Get an array of the numerical encoding for each possible base pair (i,j)  

See Also
--------
RNA.idx_col_wise(), RNA.fold_compound()  

Note
----
This array is always indexed in column-wise order, in contrast to previously different indexing
between mfe and pf variants!  
";

%feature("docstring") my_seq_encode "

Get a numerical representation of the nucleotide sequence.  

**SWIG Wrapper Notes**
    In the target scripting language, this function is wrapped as overloaded function `seq_encode()`
    where the last parameter, the `model_details` data structure, is optional. If it is omitted,
    default model settings are applied, i.e. default nucleotide letter conversion. The wrapped
    function returns a list/tuple of integer representations of the input sequence. See, e.g.
    :py:func:`RNA.seq_encode()` in the :doc:`/api_python`.  

Parameters
----------
sequence : const char *
    The input sequence in upper-case letters  
md : RNA.md() *
    A pointer to a RNA.md() data structure that specifies the conversion type  

Returns
-------
short *  
    A list of integer encodings for each sequence letter (1-based). Position 0 denotes the length of
    the list  
";

%feature("docstring") vrna_seq_encode_simple "

Get a numerical representation of the nucleotide sequence (simple version)  
";

%feature("docstring") vrna_nucleotide_encode "

Encode a nucleotide character to numerical value.  

This function encodes a nucleotide character to its numerical representation as required by many
functions in RNAlib.  

Parameters
----------
c : char
    The nucleotide character to encode  
md : RNA.md() *
    The model details that determine the kind of encoding  

Returns
-------
int  
    The encoded nucleotide  

See Also
--------
RNA.nucleotide_decode(), RNA.seq_encode()  
";

%feature("docstring") vrna_nucleotide_decode "

Decode a numerical representation of a nucleotide back into nucleotide alphabet.  

This function decodes a numerical representation of a nucleotide character back into nucleotide
alphabet  

Parameters
----------
enc : int
    The encoded nucleotide  
md : RNA.md() *
    The model details that determine the kind of decoding  

Returns
-------
char  
    The decoded nucleotide character  

See Also
--------
RNA.nucleotide_encode(), RNA.seq_encode()  
";

%feature("docstring") vrna_aln_encode "
";

%feature("docstring") vrna_get_ptype_md "
";

%feature("docstring") vrna_get_ptype "
";

%feature("docstring") vrna_get_ptype_window "
";

%feature("docstring") vrna_sequence "
";

%feature("docstring") vrna_fold_compound_t::sequence_add "
";

%feature("docstring") vrna_fold_compound_t::sequence_remove "
";

%feature("docstring") vrna_fold_compound_t::sequence_remove_all "
";

%feature("docstring") vrna_fold_compound_t::sequence_prepare "
";

%feature("docstring") vrna_sequence_order_update "
";

%feature("docstring") vrna_msa_add "
";

%feature("docstring") SEQUENCE_RNA "
";

%feature("docstring") SEQUENCE_DNA "
";

// File: group__string__utils.xml

%feature("docstring") vrna_strdup_printf "

Safely create a formatted string.  

This function is a safe implementation for creating a formatted character array, similar to
*sprintf*. Internally, it uses the *asprintf* function if available to dynamically allocate a large
enough character array to store the supplied content. If *asprintf* is not available, mimic it's
behavior using *vsnprintf*.  

Parameters
----------
format : const char *
    The format string (See also asprintf)  
... :
    The list of variables used to fill the format string  

Returns
-------
char *  
    The formatted, null-terminated string, or NULL if something has gone wrong  

See Also
--------
RNA.strdup_vprintf(), RNA.strcat_printf()  

Note
----
The returned pointer of this function should always be passed to *free()* to release the allocated
memory  
";

%feature("docstring") vrna_strdup_vprintf "

Safely create a formatted string.  

This function is the *va_list* version of RNA.strdup_printf()  

Parameters
----------
format : const char *
    The format string (See also asprintf)  
argp : va_list
    The list of arguments to fill the format string  

Returns
-------
char *  
    The formatted, null-terminated string, or NULL if something has gone wrong  

See Also
--------
RNA.strdup_printf(), RNA.strcat_printf(), RNA.strcat_vprintf()  

Note
----
The returned pointer of this function should always be passed to *free()* to release the allocated
memory  
";

%feature("docstring") vrna_strcat_printf "

Safely append a formatted string to another string.  

This function is a safe implementation for appending a formatted character array, similar to a
cobination of *strcat* and *sprintf*. The function automatically allocates enough memory to store
both, the previous content stored at `dest` and the appended format string. If the `dest` pointer is
NULL, the function allocate memory only for the format string. The function returns the number of
characters in the resulting string or -1 in case of an error.  

Parameters
----------
dest : char **
    The address of a char *pointer where the formatted string is to be appended  
format : const char *
    The format string (See also sprintf)  
... :
    The list of variables used to fill the format string  

Returns
-------
int  
    The number of characters in the final string, or -1 on error  

See Also
--------
RNA.strcat_vprintf(), RNA.strdup_printf(), RNA.strdup_vprintf()  
";

%feature("docstring") vrna_strcat_vprintf "

Safely append a formatted string to another string.  

This function is the *va_list* version of RNA.strcat_printf()  

Parameters
----------
dest : char **
    The address of a char *pointer where the formatted string is to be appended  
format : const char *
    The format string (See also sprintf)  
args : va_list
    The list of argument to fill the format string  

Returns
-------
int  
    The number of characters in the final string, or -1 on error  

See Also
--------
RNA.strcat_printf(), RNA.strdup_printf(), RNA.strdup_vprintf()  
";

%feature("docstring") my_strtrim "

Trim a string by removing (multiple) occurences of a particular character.  

This function removes (multiple) consecutive occurences of a set of characters (`delimiters`) within
an input string. It may be used to remove leading and/or trailing whitespaces or to restrict the
maximum number of consecutive occurences of the delimiting characters `delimiters`. Setting `keep=0`
removes all occurences, while other values reduce multiple consecutive occurences to at most `keep`
delimiters. This might be useful if one would like to reduce multiple whitespaces to a single one,
or to remove empty fields within a comma-separated value string.  

The parameter `delimiters` may be a pointer to a 0-terminated char string containing a set of any
ASCII character. If *NULL* is passed as delimiter set or an empty char string, all whitespace
characters are trimmed. The `options` parameter is a bit vector that specifies which part of the
string should undergo trimming. The implementation distinguishes the leading (RNA.TRIM_LEADING),
trailing (RNA.TRIM_TRAILING), and in-between (RNA.TRIM_IN_BETWEEN) part with respect to the
delimiter set. Combinations of these parts can be specified by using logical-or operator.  

The following example code removes all leading and trailing whitespace characters from the input
string:  

**SWIG Wrapper Notes**
    Since many scripting languages treat strings as immutable objects, this function does not modify
    the input string directly. Instead, it returns the modified string as second return value,
    together with the number of removed delimiters.  

    The scripting language interface provides an overloaded version of this function, with default
    parameters `delimiters=NULL`, `keep=0`, and `options=RNA.TRIM_DEFAULT`. See, e.g.
    :py:func:`RNA.strtrim()` in the :doc:`/api_python`.  

Parameters
----------
string : char *
    The '\\0'-terminated input string to trim  
delimiters : const char *
    The delimiter characters as 0-terminated char array (or *NULL*)  
keep : unsigned int
    The maximum number of consecutive occurences of the delimiter in the output string  
options : unsigned int
    The option bit vector specifying the mode of operation  

Returns
-------
unsigned int  
    The number of delimiters removed from the string  

See Also
--------
RNA.TRIM_LEADING, RNA.TRIM_TRAILING, RNA.TRIM_IN_BETWEEN, RNA.TRIM_SUBST_BY_FIRST,
RNA.TRIM_DEFAULT, RNA.TRIM_ALL  

Note
----
The delimiter always consists of a single character from the set of characters provided. In case of
alternative delimiters and non-null `keep` parameter, the first `keep` delimiters are preserved
within the string. Use RNA.TRIM_SUBST_BY_FIRST to substitute all remaining delimiting characters
with the first from the `delimiters` list.  
";

%feature("docstring") vrna_strsplit "

Split a string into tokens using a delimiting character.  

This function splits a string into an array of strings using a single character that delimits the
elements within the string. The default delimiter is the ampersand `'&'` and will be used when
`NULL` is passed as a second argument. The returned list is NULL terminated, i.e. the last element
is `NULL`. If the delimiter is not found, the returned list contains exactly one element: the input
string.  

For instance, the following code:  

 produces this output:  

    * GGGG
    * CCCC
    * AAAAA
    *  and properly free's the memory occupied by the returned element array.  

Parameters
----------
string : const char *
    The input string that should be split into elements  
delimiter : const char *
    The delimiting character. If `NULL`, the delimiter is `\"&\"`  

Returns
-------
char **  
    A `NULL` terminated list of the elements in the string  

See Also
--------
RNA.strtrim()  

Note
----
This function internally uses *strtok_r()* and is therefore considered to be thread-safe. Also note,
that it is the users responsibility to free the memory of the array and that of the individual
element strings!  
 In case the input string consists of consecutive delimiters, starts or ends with one or multiple
delimiters, empty strings are produced in the output list, indicating the empty fields of data
resulting from the split. Use RNA.strtrim() prior to a call to this function to remove any leading,
trailing, or in-between empty fields.  
";

%feature("docstring") vrna_strjoin "
";

%feature("docstring") vrna_random_string "

Create a random string using characters from a specified symbol set.  

Parameters
----------
l : int
    The length of the sequence  
symbols : const char
    The symbol set  

Returns
-------
char *  
    A random string of length 'l' containing characters from the symbolset  
";

%feature("docstring") my_hamming "

Calculate hamming distance between two sequences.  

Parameters
----------
s1 : const char *
    The first sequence  
s2 : const char *
    The second sequence  

Returns
-------
int  
    The hamming distance between s1 and s2  
";

%feature("docstring") my_hamming_bound "

Calculate hamming distance between two sequences up to a specified length.  

This function is similar to RNA.hamming_distance() but instead of comparing both sequences up to
their actual length only the first 'n' characters are taken into account  

Parameters
----------
s1 : const char *
    The first sequence  
s2 : const char *
    The second sequence  
n : int
    The length of the subsequences to consider (starting from the 5' end)  

Returns
-------
int  
    The hamming distance between s1 and s2  
";

%feature("docstring") vrna_seq_toRNA "

Convert an input sequence (possibly containing DNA alphabet characters) to RNA alphabet.  

This function substitudes *T* and *t* with *U* and *u*, respectively  

Parameters
----------
sequence : char *
    The sequence to be converted  
";

%feature("docstring") vrna_seq_toupper "

Convert an input sequence to uppercase.  

Parameters
----------
sequence : char *
    The sequence to be converted  
";

%feature("docstring") vrna_seq_reverse "

Reverse a string in-place.  

This function reverses a character string in the form of an array of characters in-place, i.e. it
changes the input parameter.  

**Postcondition**
    After execution, the input `sequence` consists of the reverse string prior to the execution.  

Parameters
----------
sequence : char *
    The string to reverse  

See Also
--------
RNA.DNA_complement()  
";

%feature("docstring") vrna_DNA_complement "

Retrieve a DNA sequence which resembles the complement of the input sequence.  

This function returns a mew DNA string which is the complement of the input, i.e. the nucleotide
letters `A`,`C`,`G`, and `T` are substituted by their complements `T`,`G`,`C`, and `A`,
respectively.  

Any characters not belonging to the alphabet of the 4 canonical bases of DNA are not altered.  

Parameters
----------
sequence : const char *
    the input DNA sequence  

Returns
-------
char *  
    The complement of the input DNA sequence  

See Also
--------
RNA.seq_reverse()  

Note
----
This function also handles lower-case input sequences and treats `U` of the RNA alphabet equally to
`T`  
";

%feature("docstring") vrna_seq_ungapped "

Remove gap characters from a nucleotide sequence.  

Parameters
----------
sequence : const char *
    The original, null-terminated nucleotide sequence  

Returns
-------
char *  
    A copy of the input sequence with all gap characters removed  
";

%feature("docstring") vrna_cut_point_insert "

Add a separating '&' character into a string according to cut-point position.  

If the cut-point position is less or equal to zero, this function just returns a copy of the
provided string. Otherwise, the cut-point character is set at the corresponding position  

Parameters
----------
string : const char *
    The original string  
cp : int
    The cut-point position  

Returns
-------
char *  
    A copy of the provided string including the cut-point character  
";

%feature("docstring") vrna_cut_point_remove "

Remove a separating '&' character from a string.  

This function removes the cut-point indicating '&' character from a string and memorizes its
position in a provided integer variable. If not '&' is found in the input, the integer variable is
set to -1. The function returns a copy of the input string with the '&' being sliced out.  

Parameters
----------
string : const char *
    The original string  
cp : int *
    The cut-point position  

Returns
-------
char *  
    A copy of the input string with the '&' being sliced out  
";

%feature("docstring") vrna_strchr "

Find (all) occurrences of a character within a string.  

string The C string to be scanned  

c The character to be searched for  

n The maximum number of occurences to search for (or 0 for all occurrences)  

Returns
-------
size() *  
    An 1-based array of positions(0-based) or **NULL** on error. Position 0 specifies the number of
    occurrences found.  
";

%feature("docstring") XSTR "

Stringify a macro after expansion.  
";

%feature("docstring") STR "

Stringify a macro argument.  
";

%feature("docstring") FILENAME_MAX_LENGTH "

Maximum length of filenames that are generated by our programs.  

This definition should be used throughout the complete ViennaRNA package wherever a static array
holding filenames of output files is declared.  
";

%feature("docstring") FILENAME_ID_LENGTH "

Maximum length of id taken from fasta header for filename generation.  

this has to be smaller than FILENAME_MAX_LENGTH since in most cases, some suffix will be appended to
the ID  
";

%feature("docstring") TRIM_LEADING "

Trim only characters leading the string.  

See Also
--------
RNA.strtrim()  
";

%feature("docstring") TRIM_TRAILING "

Trim only characters trailing the string.  

See Also
--------
RNA.strtrim()  
";

%feature("docstring") TRIM_IN_BETWEEN "

Trim only characters within the string.  

See Also
--------
RNA.strtrim()  
";

%feature("docstring") TRIM_SUBST_BY_FIRST "

Replace remaining characters after trimming with the first delimiter in list.  

See Also
--------
RNA.strtrim()  
";

%feature("docstring") TRIM_DEFAULT "

Default settings for trimming, i.e. trim leading and trailing.  

See Also
--------
RNA.strtrim()  
";

%feature("docstring") TRIM_ALL "

Trim characters anywhere in the string.  

See Also
--------
RNA.strtrim()  
";

// File: group__struct__utils.xml

%feature("docstring") my_loopidx_from_ptable "

Get a loop index representation of a structure.  
";

%feature("docstring") vrna_refBPcnt_matrix "

Make a reference base pair count matrix.  

Get an upper triangular matrix containing the number of basepairs of a reference structure for each
interval [i,j] with i<j. Access it via iindx!!!  
";

%feature("docstring") vrna_refBPdist_matrix "

Make a reference base pair distance matrix.  

Get an upper triangular matrix containing the base pair distance of two reference structures for
each interval [i,j] with i<j. Access it via iindx!!!  
";

%feature("docstring") vrna_db_from_probs "

Create a dot-bracket like structure string from base pair probability matrix.  

**SWIG Wrapper Notes**
    This function is available as parameter-less method **db_from_probs()** bound to objects of type
    *fold_compound*. Parameters `pr` and `length` are implicitely taken from the *fold_compound*
    object the method is bound to. Upon missing base pair probabilities, this method returns an
    empty string. See, e.g.  :py:func:`RNA.db_from_probs()` in the :doc:`/api_python`.  
";

%feature("docstring") vrna_bpp_symbol "

Get a pseudo dot bracket notation for a given probability information.  
";

%feature("docstring") vrna_db_from_bp_stack "

Create a dot-backet/parenthesis structure from backtracking stack.  

This function is capable to create dot-bracket structures from suboptimal structure prediction sensu
M. Zuker  

Parameters
----------
bp : RNA.bp_stack() *
    Base pair stack containing the traced base pairs  
length : unsigned int
    The length of the structure  

Returns
-------
char *  
    The secondary structure in dot-bracket notation as provided in the input  
";

%feature("docstring") vrna_letter_structure "
";

%feature("docstring") make_pair_table_pk "
";

%feature("docstring") make_loop_index_pt "
";

%feature("docstring") letter_structure "
";

// File: group__struct__utils__dot__bracket.xml

%feature("docstring") vrna_db_pack "

Pack secondary secondary structure, 5:1 compression using base 3 encoding.  

Returns a binary string encoding of the secondary structure using a 5:1 compression scheme. The
string is NULL terminated and can therefore be used with standard string functions such as strcmp().
Useful for programs that need to keep many structures in memory.  

Parameters
----------
struc : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
char *  
    The binary encoded structure  

See Also
--------
RNA.db_unpack()  
";

%feature("docstring") vrna_db_unpack "

Unpack secondary structure previously packed with RNA.db_pack()  

Translate a compressed binary string produced by RNA.db_pack() back into the familiar dot-bracket
notation.  

Parameters
----------
packed : const char *
    The binary encoded packed secondary structure  

Returns
-------
char *  
    The unpacked secondary structure in dot-bracket notation  

See Also
--------
RNA.db_pack()  
";

%feature("docstring") db_flatten "

Substitute pairs of brackets in a string with parenthesis.  

This function can be used to replace brackets of unusual types, such as angular brackets `<>` , to
dot-bracket format. The `options` parameter is used tpo specify which types of brackets will be
replaced by round parenthesis ``() .  

**SWIG Wrapper Notes**
    This function flattens an input structure string in-place! The second parameter is optional and
    defaults to RNA.BRACKETS_DEFAULT.  

    An overloaded version of this function exists, where an additional second parameter can be
    passed to specify the target brackets, i.e. the type of matching pair characters all brackets
    will be flattened to. Therefore, in the scripting language interface this function is a
    replacement for RNA.db_flatten_to(). See, e.g.  :py:func:`RNA.db_flatten()` in the
    :doc:`/api_python`.  

Parameters
----------
structure : char *
    The structure string where brackets are flattened in-place  
options : unsigned int
    A bitmask to specify which types of brackets should be flattened out  

See Also
--------
RNA.db_flatten_to(), RNA.BRACKETS_RND, RNA.BRACKETS_ANG, RNA.BRACKETS_CLY, RNA.BRACKETS_SQR,
RNA.BRACKETS_DEFAULT  
";

%feature("docstring") vrna_db_flatten_to "

Substitute pairs of brackets in a string with another type of pair characters.  

This function can be used to replace brackets in a structure annotation string, such as square
brackets ``[] , to another type of pair characters, e.g. angular brackets `<>` .  

The `target` array must contain a character for the 'pair open' annotation at position 0, and one
for 'pair close' at position 1. T`options` parameter is used to specify which types of brackets will
be replaced by the new pairs.  

**SWIG Wrapper Notes**
    This function is available as an overloaded version of RNA.db_flatten(). See, e.g.
    :py:func:`RNA.db_flatten()` in the :doc:`/api_python`.  

Parameters
----------
string : char *
    The structure string where brackets are flattened in-place  
target : const char
    The new pair characters the string will be flattened to  
options : unsigned int
    A bitmask to specify which types of brackets should be flattened out  

See Also
--------
RNA.db_flatten(), RNA.BRACKETS_RND, RNA.BRACKETS_ANG, RNA.BRACKETS_CLY, RNA.BRACKETS_SQR,
RNA.BRACKETS_DEFAULT  
";

%feature("docstring") my_db_from_ptable "

Convert a pair table into dot-parenthesis notation.  

This function also converts pair table formatted structures that contain pseudoknots. Non-nested
base pairs result in additional pairs of parenthesis and brackets within the resulting dot-
parenthesis string. The following pairs are awailable: (), []. {}. <>, as well as pairs of matching
upper-/lower-case characters from the alphabet A-Z.  

Parameters
----------
pt : const short *
    The pair table to be copied  

Returns
-------
char *  
    A char pointer to the dot-bracket string  

Note
----
In cases where the level of non-nested base pairs exceeds the maximum number of 30 different base
pair indicators (4 parenthesis/brackets, 26 matching characters), a warning is printed and the
remaining base pairs are left out from the conversion.  
";

%feature("docstring") db_from_plist "

Convert a list of base pairs into dot-bracket notation.  

Parameters
----------
pairs : RNA.ep() *
    A RNA.ep() containing the pairs to be included in the dot-bracket string  
n : unsigned int
    The length of the structure (number of nucleotides)  

Returns
-------
char *  
    The dot-bracket string containing the provided base pairs  

See Also
--------
RNA.plist()  
";

%feature("docstring") vrna_db_to_element_string "

Convert a secondary structure in dot-bracket notation to a nucleotide annotation of loop contexts.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
char *  
    A string annotating each nucleotide according to it's structural context  
";

%feature("docstring") db_pk_remove "

Remove pseudo-knots from an input structure.  

This function removes pseudo-knots from an input structure by determining the minimum number of base
pairs that need to be removed to make the structure pseudo-knot free.  

To accomplish that, we use a dynamic programming algorithm similar to the Nussinov maxmimum matching
approach.  

The input structure must be in a dot-bracket string like form where crossing base pairs are denoted
by the use of additional types of matching brackets, e.g. `<>`, `{}`, ``[], `{}`. Furthermore,
crossing pairs may be annotated by matching uppercase/lowercase letters from the alphabet `A-Z`. For
the latter, the uppercase letter must be the 5' and the lowercase letter the 3' nucleotide of the
base pair. The actual type of brackets to be recognized by this function must be specifed through
the `options` parameter.  

**SWIG Wrapper Notes**
    This function is available as an overloaded function `db_pk_remove()` where the optional second
    parameter `options` defaults to RNA.BRACKETS_ANY. See, e.g.  :py:func:`RNA.db_pk_remove()` in
    the :doc:`/api_python`.  

Parameters
----------
structure : const char *
    Input structure in dot-bracket format that may include pseudo-knots  
options : unsigned int
    A bitmask to specify which types of brackets should be processed  

Returns
-------
char *  
    The input structure devoid of pseudo-knots in dot-bracket notation  

See Also
--------
RNA.pt_pk_remove(), RNA.db_flatten(), RNA.BRACKETS_RND, RNA.BRACKETS_ANG, RNA.BRACKETS_CLY,
RNA.BRACKETS_SQR, RNA.BRACKETS_ALPHA, RNA.BRACKETS_DEFAULT, RNA.BRACKETS_ANY  

Note
----
Brackets in the input structure string that are not covered by the `options` bitmask will be
silently ignored!  
";

%feature("docstring") BRACKETS_ALPHA "

Bitflag to indicate secondary structure notations using uppercase/lowercase letters from the latin
alphabet.  

See Also
--------
RNA.ptable_from_string()  
";

%feature("docstring") BRACKETS_RND "

Bitflag to indicate secondary structure notations using round brackets (parenthesis), `()`  

See Also
--------
RNA.ptable_from_string(), RNA.db_flatten(), RNA.db_flatten_to()  
";

%feature("docstring") BRACKETS_CLY "

Bitflag to indicate secondary structure notations using curly brackets, `{}`  

See Also
--------
RNA.ptable_from_string(), RNA.db_flatten(), RNA.db_flatten_to()  
";

%feature("docstring") BRACKETS_ANG "

Bitflag to indicate secondary structure notations using angular brackets, `<>`  

See Also
--------
RNA.ptable_from_string(), RNA.db_flatten(), RNA.db_flatten_to()  
";

%feature("docstring") BRACKETS_SQR "

Bitflag to indicate secondary structure notations using square brackets, `[]`  

See Also
--------
RNA.ptable_from_string(), RNA.db_flatten(), RNA.db_flatten_to()  
";

%feature("docstring") BRACKETS_DEFAULT "

Default bitmask to indicate secondary structure notation using any pair of brackets.  

This set of matching brackets/parenthesis is always nested, i.e. pseudo-knot free, in WUSS format.
However, in general different kinds of brackets are mostly used for annotating pseudo-knots. Thus
special care has to be taken to remove pseudo-knots if this bitmask is used in functions that return
secondary structures without pseudo-knots!  

See Also
--------
RNA.ptable_from_string(), RNA.db_flatten(), RNA.db_flatten_to(),
RNA.db_pk_remove()RNA.pt_pk_remove()  
";

%feature("docstring") BRACKETS_ANY "

Bitmask to indicate secondary structure notation using any pair of brackets or uppercase/lowercase
alphabet letters.  

See Also
--------
RNA.ptable_from_string(), RNA.db_pk_remove(), RNA.db_flatten(), RNA.db_flatten_to()  
";

// File: group__struct__utils__wuss.xml

%feature("docstring") db_from_WUSS "

Convert a WUSS annotation string to dot-bracket format.  

Parameters
----------
wuss : const char *
    The input string in WUSS notation  

Returns
-------
char *  
    A dot-bracket notation of the input secondary structure  

Note
----
This function flattens all brackets, and treats pseudo-knots annotated by matching pairs of
upper/lowercase letters as unpaired nucleotides  
";

// File: group__struct__utils__pair__table.xml

%feature("docstring") vrna_ptable "

Create a pair table from a dot-bracket notation of a secondary structure.  

Returns a newly allocated table, such that table[i]=j if (i.j) pair or 0 if i is unpaired, table[0]
contains the length of the structure.  

**SWIG Wrapper Notes**
    This functions is wrapped as overloaded function `ptable()` that takes an optional argument
    `options` to specify which type of matching brackets should be considered during conversion. The
    default set is round brackets, i.e. RNA.BRACKETS_RND. See, e.g.  :py:func:`RNA.ptable()` in the
    :doc:`/api_python`.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
short *  
    A pointer to the created pair_table  

See Also
--------
RNA.ptable_from_string(), RNA.db_from_ptable()  
";

%feature("docstring") my_ptable "

Create a pair table for a secondary structure string.  

This function takes an input string of a secondary structure annotation in dot-bracket-notation or
dot-bracket-ext-notation, and converts it into a pair table representation.  

**SWIG Wrapper Notes**
    This functions is wrapped as overloaded function `ptable()` that takes an optional argument
    `options` to specify which type of matching brackets should be considered during conversion. The
    default set is round brackets, i.e. RNA.BRACKETS_RND. See, e.g.  :py:func:`RNA.ptable()` in the
    :doc:`/api_python`.  

Parameters
----------
structure : const char *
    Secondary structure in dot-bracket-ext-notation  
options : unsigned int
    A bitmask to specify which brackets are recognized during conversion to pair table  

Returns
-------
short *  
    A pointer to a new pair table of the provided secondary structure  

See Also
--------
RNA.ptable(), RNA.db_from_ptable(), RNA.db_flatten_to(), RNA.pt_pk_remove()RNA.BRACKETS_RND,
RNA.BRACKETS_ANG, RNA.BRACKETS_CLY, RNA.BRACKETS_SQR, RNA.BRACKETS_ALPHA, RNA.BRACKETS_DEFAULT,
RNA.BRACKETS_ANY  

Note
----
This function also extracts crossing base pairs, i.e. pseudo-knots if more than a single matching
bracket type is allowed through the bitmask `options`.  
";

%feature("docstring") my_ptable_pk "

Create a pair table of a secondary structure (pseudo-knot version)  

Returns a newly allocated table, such that table[i]=j if (i.j) pair or 0 if i is unpaired, table[0]
contains the length of the structure.  

In contrast to RNA.ptable() this function also recognizes the base pairs denoted by '[' and ']'
brackets. Thus, this function behaves like  

Parameters
----------
structure : const char *
    The secondary structure in (extended) dot-bracket notation  

Returns
-------
short *  
    A pointer to the created pair_table  

See Also
--------
RNA.ptable_from_string()  
";

%feature("docstring") vrna_ptable_copy "

Get an exact copy of a pair table.  

Parameters
----------
pt : const short *
    The pair table to be copied  

Returns
-------
short *  
    A pointer to the copy of 'pt'  
";

%feature("docstring") vrna_pt_ali_get "

Create a pair table of a secondary structure (snoop align version)  
";

%feature("docstring") vrna_pt_snoop_get "

Create a pair table of a secondary structure (snoop version)  

returns a newly allocated table, such that: table[i]=j if (i.j) pair or 0 if i is unpaired, table[0]
contains the length of the structure. The special pseudoknotted H/ACA-mRNA structure is taken into
account.  
";

%feature("docstring") my_pt_pk_remove "

Remove pseudo-knots from a pair table.  

This function removes pseudo-knots from an input structure by determining the minimum number of base
pairs that need to be removed to make the structure pseudo-knot free.  

To accomplish that, we use a dynamic programming algorithm similar to the Nussinov maxmimum matching
approach.  

Parameters
----------
ptable : const short *
    Input structure that may include pseudo-knots  
options : unsigned int

Returns
-------
short *  
    The input structure devoid of pseudo-knots  

See Also
--------
RNA.db_pk_remove()  
";

// File: group__struct__utils__plist.xml

%feature("docstring") my_plist "

Create a RNA.ep() from a dot-bracket string.  

The dot-bracket string is parsed and for each base pair an entry in the plist is created. The
probability of each pair in the list is set by a function parameter.  

The end of the plist is marked by sequence positions i as well as j equal to 0. This condition
should be used to stop looping over its entries  

Parameters
----------
struc : const char *
    The secondary structure in dot-bracket notation  
pr : float
    The probability for each base pair used in the plist  

Returns
-------
RNA.ep() *  
    The plist array  
";

%feature("docstring") PLIST_TYPE_BASEPAIR "

A Base Pair element.  
";

%feature("docstring") PLIST_TYPE_GQUAD "

A G-Quadruplex element.  
";

%feature("docstring") PLIST_TYPE_H_MOTIF "

A Hairpin loop motif element.  
";

%feature("docstring") PLIST_TYPE_I_MOTIF "

An Internal loop motif element.  
";

%feature("docstring") PLIST_TYPE_UD_MOTIF "

An Unstructured Domain motif element.  
";

%feature("docstring") PLIST_TYPE_STACK "

A Base Pair stack element.  
";

%feature("docstring") PLIST_TYPE_UNPAIRED "

An unpaired base.  
";

%feature("docstring") PLIST_TYPE_TRIPLE "

One pair of a base triplet.  
";

// File: group__struct__utils__abstract__shapes.xml

%feature("docstring") abstract_shapes "

Convert a secondary structure in dot-bracket notation to its abstract shapes representation.  

This function converts a secondary structure into its abstract shapes representation as presented by
:cite:t:`giegerich:2004` .  

**SWIG Wrapper Notes**
    This function is available as an overloaded function `abstract_shapes()` where the optional
    second parameter `level` defaults to 5. See, e.g.  :py:func:`RNA.abstract_shapes()` in the
    :doc:`/api_python`.  

Parameters
----------
structure : const char *
    A secondary structure in dot-bracket notation  
level : unsigned int
    The abstraction level (integer in the range of 0 to 5)  

Returns
-------
char *  
    The secondary structure in abstract shapes notation  

See Also
--------
RNA.abstract_shapes_pt()  
";

%feature("docstring") vrna_abstract_shapes_pt "

Convert a secondary structure to its abstract shapes representation.  

This function converts a secondary structure into its abstract shapes representation as presented by
:cite:t:`giegerich:2004` . This function is equivalent to RNA.db_to_shapes(), but requires a pair
table input instead of a dot-bracket structure.  

**SWIG Wrapper Notes**
    This function is available as an overloaded function `abstract_shapes()` where the optional
    second parameter `level` defaults to 5. See, e.g.  :py:func:`RNA.abstract_shapes()` in the
    :doc:`/api_python`.  

Parameters
----------
pt : const short *
    A secondary structure in pair table format  
level : unsigned int
    The abstraction level (integer in the range of 0 to 5)  

Returns
-------
char *  
    The secondary structure in abstract shapes notation  

See Also
--------
RNA.abstract_shapes()  

Note
----
The length of the structure must be present at `pt`[0]!  
";

// File: group__struct__utils__helix__list.xml

%feature("docstring") my_hx_from_ptable "

Convert a pair table representation of a secondary structure into a helix list.  

Parameters
----------
pt : short *
    The secondary structure in pair table representation  

Returns
-------
RNA.hx() *  
    The secondary structure represented as a helix list  
";

%feature("docstring") vrna_hx_merge "

Create a merged helix list from another helix list.  
";

// File: group__struct__utils__tree.xml

%feature("docstring") db_to_tree_string "

Convert a Dot-Bracket structure string into tree string representation.  

This function allows one to convert a secondary structure in dot-bracket notation into one of the
various tree representations for secondary structures. The resulting tree is then represented as a
string of parenthesis and node symbols, similar to to the Newick format.  

Currently we support conversion into the following formats, denoted by the value of parameter
`type:`  

*   RNA.STRUCTURE_TREE_HIT -  Homeomorphically Irreducible Tree (HIT) representation of a secondary
    structure.   (See also   :cite:t:`fontana:1993b` )  
*   RNA.STRUCTURE_TREE_SHAPIRO_SHORT -  (short) Coarse Grained representation of a secondary
    structure   (same as   :cite:t:`shapiro:1988` , but with root node `R` and without `S` nodes for
    the stems)  
*   RNA.STRUCTURE_TREE_SHAPIRO -  (full) Coarse Grained representation of a secondary structure
    (See also   :cite:t:`shapiro:1988` )  
*   RNA.STRUCTURE_TREE_SHAPIRO_EXT -  (extended) Coarse Grained representation of a secondary
    structure   (same as   :cite:t:`shapiro:1988` , but external nodes denoted as `E` )  
*   RNA.STRUCTURE_TREE_SHAPIRO_WEIGHT -  (weighted) Coarse Grained representation of a secondary
    structure   (same as RNA.STRUCTURE_TREE_SHAPIRO_EXT but with additional weights for number of
    unpaired nucleotides in loop, and number of pairs in stems)  
*   RNA.STRUCTURE_TREE_EXPANDED -  Expanded Tree representation of a secondary structure.  

Parameters
----------
structure : const char *
    The null-terminated dot-bracket structure string  
type : unsigned int
    A switch to determine the type of tree string representation  

Returns
-------
char *  
    A tree representation of the input `structure`  

See Also
--------
sec_structure_representations_tree  
";

%feature("docstring") tree_string_unweight "

Remove weights from a linear string tree representation of a secondary structure.  

This function strips the weights of a linear string tree representation such as `HIT`, or Coarse
Grained Tree sensu   :cite:t:`shapiro:1988`  

Parameters
----------
structure : const char *
    A linear string tree representation of a secondary structure with weights  

Returns
-------
char *  
    A linear string tree representation of a secondary structure without weights  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") tree_string_to_db "

Convert a linear tree string representation of a secondary structure back to Dot-Bracket notation.  

Parameters
----------
tree : const char *
    A linear tree string representation of a secondary structure  

Returns
-------
char *  
    A dot-bracket notation of the secondary structure provided in `tree`  

Warnings
--------
This function only accepts *Expanded* and *HIT* tree representations!  

See Also
--------
RNA.db_to_tree_string(), RNA.STRUCTURE_TREE_EXPANDED, RNA.STRUCTURE_TREE_HIT,
sec_structure_representations_tree  
";

%feature("docstring") STRUCTURE_TREE_HIT "

Homeomorphically Irreducible Tree (HIT) representation of a secondary structure.  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") STRUCTURE_TREE_SHAPIRO_SHORT "

(short) Coarse Grained representation of a secondary structure  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") STRUCTURE_TREE_SHAPIRO "

(full) Coarse Grained representation of a secondary structure  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") STRUCTURE_TREE_SHAPIRO_EXT "

(extended) Coarse Grained representation of a secondary structure  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") STRUCTURE_TREE_SHAPIRO_WEIGHT "

(weighted) Coarse Grained representation of a secondary structure  

See Also
--------
RNA.db_to_tree_string()  
";

%feature("docstring") STRUCTURE_TREE_EXPANDED "

Expanded Tree representation of a secondary structure.  

See Also
--------
RNA.db_to_tree_string()  
";

// File: group__struct__utils__metrics.xml

%feature("docstring") vrna_bp_distance_pt "

Compute the \"base pair\" distance between two pair tables pt1 and pt2 of secondary structures.  

The pair tables should have the same length. dist = number of base pairs in one structure but not in
the other same as edit distance with open-pair close-pair as move-set  

**SWIG Wrapper Notes**
    This function is available as an overloaded method **bp_distance()**. See, e.g.
    :py:func:`RNA.bp_distance()` in the :doc:`/api_python`.  

Parameters
----------
pt1 : const short *
    First structure in dot-bracket notation  
pt2 : const short *
    Second structure in dot-bracket notation  

Returns
-------
int  
    The base pair distance between pt1 and pt2  

See Also
--------
RNA.bp_distance()  
";

%feature("docstring") my_bp_distance "

Compute the \"base pair\" distance between two secondary structures s1 and s2.  

This is a wrapper around **RNA.bp_distance_pt()**. The sequences should have the same length. dist
= number of base pairs in one structure but not in the other same as edit distance with open-pair
close-pair as move-set  

**SWIG Wrapper Notes**
    This function is available as an overloaded method **bp_distance()**. Note that the SWIG wrapper
    takes two structure in dot-bracket notation and converts them into pair tables using
    RNA.ptable_from_string(). The resulting pair tables are then internally passed to
    RNA.bp_distance_pt(). To control which kind of matching brackets will be used during
    conversion, the optional argument `options` can be used. See also the description of
    RNA.ptable_from_string() for available options. (default: **RNA.BRACKETS_RND**). See, e.g.
    :py:func:`RNA.bp_distance()` in the :doc:`/api_python`.  

Parameters
----------
str1 : const char *
    First structure in dot-bracket notation  
str2 : const char *
    Second structure in dot-bracket notation  

Returns
-------
int  
    The base pair distance between str1 and str2  

See Also
--------
RNA.bp_distance_pt()  
";

%feature("docstring") my_dist_mountain "
";

// File: group__aln__utils.xml

%feature("docstring") my_aln_mpi "

Get the mean pairwise identity in steps from ?to?(ident)  

**SWIG Wrapper Notes**
    This function is available as function `aln_mpi()`. See e.g.  :py:func:`RNA.aln_mpi()` in the
    :doc:`/api_python`.  

Parameters
----------
alignment : const char **
    Aligned sequences  

Returns
-------
int  
    The mean pairwise identity  
";

%feature("docstring") vrna_aln_pinfo "

Retrieve an array of RNA.pinfo() structures from precomputed pair probabilities.  

This array of structures contains information about positionwise pair probabilies, base pair entropy
and more  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() of type RNA.FC_TYPE_COMPARATIVE with precomputed partition function
    matrices  
structure : const char *
    An optional structure in dot-bracket notation (Maybe NULL)  
threshold : double
    Do not include results with pair probabilities below threshold  

Returns
-------
RNA.pinfo() *  
    The RNA.pinfo() array  

See Also
--------
RNA.pinfo(), and RNA.fold_compound.pf()  
";

%feature("docstring") my_aln_pscore "

**SWIG Wrapper Notes**
    This function is available as overloaded function `aln_pscore()` where the last parameter may be
    omitted, indicating `md` = `NULL`. See e.g.  :py:func:`RNA.aln_pscore()` in the
    :doc:`/api_python`.  
";

%feature("docstring") vrna_pscore "
";

%feature("docstring") vrna_pscore_freq "
";

%feature("docstring") vrna_aln_slice "

Slice out a subalignment from a larger alignment.  

Parameters
----------
alignment : const char **
    The input alignment  
i : unsigned int
    The first column of the subalignment (1-based)  
j : unsigned int
    The last column of the subalignment (1-based)  

Returns
-------
char **  
    The subalignment between column :math:`i` and :math:`j`  

See Also
--------
RNA.aln_free()  

Note
----
The user is responsible to free the memory occupied by the returned subalignment  
";

%feature("docstring") vrna_aln_free "

Free memory occupied by a set of aligned sequences.  

Parameters
----------
alignment : char **
    The input alignment  
";

%feature("docstring") vrna_aln_uppercase "

Create a copy of an alignment with only uppercase letters in the sequences.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  

Returns
-------
char **  
    A copy of the input alignment where lowercase sequence letters are replaced by uppercase letters  

See Also
--------
RNA.aln_copy  
";

%feature("docstring") vrna_aln_toRNA "

Create a copy of an alignment where DNA alphabet is replaced by RNA alphabet.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  

Returns
-------
char **  
    A copy of the input alignment where DNA alphabet is replaced by RNA alphabet (T -> U)  

See Also
--------
RNA.aln_copy  
";

%feature("docstring") vrna_aln_copy "

Make a copy of a multiple sequence alignment.  

This function allows one to create a copy of a multiple sequence alignment. The `options` parameter
additionally allows for sequence manipulation, such as converting DNA to RNA alphabet, and
conversion to uppercase letters.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  
options : unsigned int
    Option flags indicating whether the aligned sequences should be converted  

Returns
-------
char **  
    A (manipulated) copy of the input alignment  

See Also
--------
RNA.aln_copy(), RNA.ALN_RNA, RNA.ALN_UPPERCASE, RNA.ALN_DEFAULT  
";

%feature("docstring") my_aln_conservation_struct "

Compute base pair conservation of a consensus structure.  

This function computes the base pair conservation (fraction of canonical base pairs) of a consensus
structure given a multiple sequence alignment. The base pair types that are considered canonical may
be specified using the RNA.md().pair array. Passing *NULL* as parameter `md` results in default
pairing rules, i.e. canonical Watson-Crick and GU Wobble pairs.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `aln_conservation_struct()` where the last
    parameter `md` may be omitted, indicating `md` = `NULL`. See, e.g.
    :py:func:`RNA.aln_conservation_struct()` in the :doc:`/api_python`.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  
structure : const char *
    The consensus structure in dot-bracket notation  
md : const RNA.md() *
    Model details that specify compatible base pairs (Maybe *NULL*)  

Returns
-------
float *  
    A 1-based vector of base pair conservations  
";

%feature("docstring") my_aln_conservation_col "

Compute nucleotide conservation in an alignment.  

This function computes the conservation of nucleotides in alignment columns. The simples measure is
Shannon Entropy and can be selected by passing the RNA.MEASURE_SHANNON_ENTROPY flag in the
`options` parameter.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `aln_conservation_col()` where the last two
    parameters may be omitted, indicating `md` = `NULL`, and `options` =
    RNA.MEASURE_SHANNON_ENTROPY, respectively. See e.g.  :py:func:`RNA.aln_conservation_col()` in
    the :doc:`/api_python`.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  
md :
    Model details that specify known nucleotides (Maybe *NULL*)  
options : unsigned int
    A flag indicating which measure of conservation should be applied  

Returns
-------
float *  
    A 1-based vector of column conservations  

See Also
--------
RNA.MEASURE_SHANNON_ENTROPY  

Note
----
Currently, only RNA.MEASURE_SHANNON_ENTROPY is supported as conservation measure.  
";

%feature("docstring") my_aln_consensus_sequence "

Compute the consensus sequence for a given multiple sequence alignment.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `aln_consensus_sequence()` where the last
    parameter may be omitted, indicating `md` = `NULL`. See e.g.
    :py:func:`RNA.aln_consensus_sequence()` in the :doc:`/api_python`.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  
md_p : const RNA.md() *
    Model details that specify known nucleotides (Maybe *NULL*)  

Returns
-------
char *  
    The consensus sequence of the alignment, i.e. the most frequent nucleotide for each alignment
    column  
";

%feature("docstring") my_aln_consensus_mis "

Compute the Most Informative Sequence (MIS) for a given multiple sequence alignment.  

The most informative sequence (MIS)   :cite:p:`freyhult:2005`  displays for each alignment column
the nucleotides with frequency greater than the background frequency, projected into IUPAC notation.
Columns where gaps are over-represented are in lower case.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `aln_consensus_mis()` where the last parameter
    may be omitted, indicating `md` = `NULL`. See e.g.  :py:func:`RNA.aln_consensus_mis()` in the
    :doc:`/api_python`.  

Parameters
----------
alignment : const char **
    The input sequence alignment (last entry must be *NULL* terminated)  
md_p : const RNA.md() *
    Model details that specify known nucleotides (Maybe *NULL*)  

Returns
-------
char *  
    The most informative sequence for the alignment  
";

%feature("docstring") ALN_DEFAULT "

Use default alignment settings.  
";

%feature("docstring") ALN_RNA "

Convert to RNA alphabet.  
";

%feature("docstring") ALN_DNA "

Convert to DNA alphabet.  
";

%feature("docstring") ALN_UPPERCASE "

Convert to uppercase nucleotide letters.  
";

%feature("docstring") ALN_LOWERCASE "

Convert to lowercase nucleotide letters.  
";

%feature("docstring") MEASURE_SHANNON_ENTROPY "

Flag indicating Shannon Entropy measure.  

Shannon Entropy is defined as :math:`H = - \\sum_{c} p_{c} \\cdot \\log_{2} p_{c}`  
";

// File: group__file__utils.xml

%feature("docstring") get_ribosum "

Retrieve a RiboSum Scoring Matrix for a given Alignment.  
";

%feature("docstring") readribosum "

Read a RiboSum or other user-defined Scoring Matrix and Store into global Memory.  
";

%feature("docstring") vrna_file_copy "

Inefficient cp.  
";

%feature("docstring") vrna_read_line "

Read a line of arbitrary length from a stream.  

Returns a pointer to the resulting string. The necessary memory is allocated and should be released
using *free()* when the string is no longer needed.  

Parameters
----------
fp : FILE *
    A file pointer to the stream where the function should read from  

Returns
-------
char *  
    A pointer to the resulting string  
";

%feature("docstring") vrna_mkdir_p "

Recursivly create a directory tree.  
";

%feature("docstring") vrna_basename "

Extract the filename from a file path.  
";

%feature("docstring") vrna_dirname "

Extract the directory part of a file path.  
";

%feature("docstring") my_filename_sanitize "

Sanitize a file name.  

Returns a new file name where all invalid characters are substituted by a replacement character. If
no replacement character is supplied, invalid characters are simply removed from the filename. File
names may also never exceed a length of 255 characters. Longer file names will undergo a 'smart'
truncation process, where the filenames suffix, i.e. everything after the last dot .', is attempted
to be kept intact. Hence, only the filename part before the suffix is reduced in such a way that the
total filename complies to the length restriction of 255 characters. If no suffix is present or the
suffix itself already exceeds the maximum length, the filename is simply truncated from the back of
the string.  

For now we consider the following characters invalid:  

*   backslash '\\'  
*   slash '/'  
*   question mark '?'  
*   percent sign ''  
*   asterisk '*'  
*   colon ':'  
*   pipe symbol '|'  
*   double quote '\"'  
*   triangular brackets '<' and '>'  

Furthermore, the (resulting) file name must not be a reserved file name, such as:  

*   '.'  
*   '..'  

Parameters
----------
name : const char *
    The input file name  
replacement : const char *
    The replacement character, or NULL  

Returns
-------
char *  
    The sanitized file name, or NULL  

Note
----
This function allocates a new block of memory for the sanitized string. It also may return (a) NULL
if the input is pointing to NULL, or (b) an empty string if the input only consists of invalid
characters which are simply removed!  
";

%feature("docstring") vrna_file_exists "

Check if a file already exists in the file system.  

Parameters
----------
filename : const char *
    The name of (path to) the file to check for existence  

Returns
-------
int  
    0 if it doesn't exists, 1 otherwise  
";

// File: group__file__formats.xml

%feature("docstring") vrna_file_helixlist "

Print a secondary structure as helix list.  

Parameters
----------
seq : const char *
    The RNA sequence  
db : const char *
    The structure in dot-bracket format  
energy : float
    Free energy of the structure in kcal/mol  
file : FILE *
    The file handle used to print to (print defaults to 'stdout' if(file == NULL) )  
";

%feature("docstring") vrna_file_connect "

Print a secondary structure as connect table.  

Connect table file format looks like this:  

    * 300  ENERGY = 7.0  example
    * 1 G       0    2   22    1
    * 2 G       1    3   21    2
    *  where the headerline is followed by 6 columns with:  

1.  Base number: index n  
2.  Base (A, C, G, T, U, X)  
3.  Index n-1 (0 if first nucleotide)  
4.  Index n+1 (0 if last nucleotide)  
5.  Number of the base to which n is paired. No pairing is indicated by 0 (zero).  
6.  Natural numbering.  

Parameters
----------
seq : const char *
    The RNA sequence  
db : const char *
    The structure in dot-bracket format  
energy : float
    The free energy of the structure  
identifier : const char *
    An optional identifier for the sequence  
file : FILE *
    The file handle used to print to (print defaults to 'stdout' if(file == NULL) )  
";

%feature("docstring") vrna_file_bpseq "

Print a secondary structure in bpseq format.  

Parameters
----------
seq : const char *
    The RNA sequence  
db : const char *
    The structure in dot-bracket format  
file : FILE *
    The file handle used to print to (print defaults to 'stdout' if(file == NULL) )  
";

%feature("docstring") vrna_file_json "

Print a secondary structure in jsonformat.  

Parameters
----------
seq : const char *
    The RNA sequence  
db : const char *
    The structure in dot-bracket format  
energy : double
    The free energy  
identifier : const char *
    An identifier for the sequence  
file : FILE *
    The file handle used to print to (print defaults to 'stdout' if(file == NULL) )  
";

%feature("docstring") my_file_fasta_read "

Get a (fasta) data set from a file or stdin.  

This function may be used to obtain complete datasets from a filehandle or stdin. A dataset is
always defined to contain at least a sequence. If data starts with a fasta header, i.e. a line like  

    >some header info  then RNA.file_fasta_read_record() will assume that the sequence that follows
the header may span over several lines. To disable this behavior and to assign a single line to the
argument 'sequence' one can pass RNA.INPUT_NO_SPAN in the 'options' argument. If no fasta header is
read in the beginning of a data block, a sequence must not span over multiple lines!  

Unless the options RNA.INPUT_NOSKIP_COMMENTS or RNA.INPUT_NOSKIP_BLANK_LINES are passed, a
sequence may be interrupted by lines starting with a comment character or empty lines.  
 A sequence is regarded as completely read if it was either assumed to not span over multiple lines,
a secondary structure or structure constraint follows the sequence on the next line, or a new header
marks the beginning of a new sequence...  

All lines following the sequence (this includes comments) that do not initiate a new dataset
according to the above definition are available through the line-array 'rest'. Here one can usually
find the structure constraint or other information belonging to the current dataset. Filling of
'rest' may be prevented by passing RNA.INPUT_NO_REST to the options argument.  

The main purpose of this function is to be able to easily parse blocks of data in the header of a
loop where all calculations for the appropriate data is done inside the loop. The loop may be then
left on certain return values, e.g.:  


In the example above, the while loop will be terminated when RNA.file_fasta_read_record() returns
either an error, EOF, or a user initiated quit request.  

As long as data is read from stdin (we are passing NULL as the file pointer), the id is printed if
it is available for the current block of data. The sequence will be printed in any case and if some
more lines belong to the current block of data each line will be printed as well.  

Parameters
----------
header : char **
    A pointer which will be set such that it points to the header of the record  
sequence : char **
    A pointer which will be set such that it points to the sequence of the record  
rest : char ***
    A pointer which will be set such that it points to an array of lines which also belong to the
    record  
file : FILE *
    A file handle to read from (if NULL, this function reads from stdin)  
options : unsigned int
    Some options which may be passed to alter the behavior of the function, use 0 for no options  

Returns
-------
unsigned int  
    A flag with information about what the function actually did read  

Note
----
This function will exit any program with an error message if no sequence could be read!  
 This function is NOT threadsafe! It uses a global variable to store information about the next data
block. Do not forget to free the memory occupied by header, sequence and rest!  
";

%feature("docstring") vrna_extract_record_rest_structure "

Extract a dot-bracket structure string from (multiline)character array.  

This function extracts a dot-bracket structure string from the 'rest' array as returned by
RNA.file_fasta_read_record() and returns it. All occurences of comments within the 'lines' array
will be skipped as long as they do not break the structure string. If no structure could be read,
this function returns NULL.  

**Precondition**
    The argument 'lines' has to be a 2-dimensional character array as obtained by
    RNA.file_fasta_read_record()  

Parameters
----------
lines : const char **
    The (multiline) character array to be parsed  
length : unsigned int
    The assumed length of the dot-bracket string (passing a value < 1 results in no length limit)  
option : unsigned int
    Some options which may be passed to alter the behavior of the function, use 0 for no options  

Returns
-------
char *  
    The dot-bracket string read from lines or NULL  

See Also
--------
RNA.file_fasta_read_record()  
";

%feature("docstring") my_file_SHAPE_read "

Read data from a given SHAPE reactivity input file.  

This function parses the informations from a given file and stores the result in the preallocated
string sequence and the double array values.  

Parameters
----------
file_name : const char *
    Path to the constraints file  
length : int
    Length of the sequence (file entries exceeding this limit will cause an error)  
default_value : double
    Value for missing indices  
sequence : char *
    Pointer to an array used for storing the sequence obtained from the SHAPE reactivity file  
values : double *
    Pointer to an array used for storing the values obtained from the SHAPE reactivity file  
";

%feature("docstring") my_file_connect_read_record "
";

%feature("docstring") my_file_RNAstrand_db_read_record "
";

%feature("docstring") vrna_extract_record_rest_constraint "

Extract a hard constraint encoded as pseudo dot-bracket string.  

.. deprecated:: 2.6.4
    Use RNA.extract_record_rest_structure() instead!  

**Precondition**
    The argument 'lines' has to be a 2-dimensional character array as obtained by
    RNA.file_fasta_read_record()  

Parameters
----------
cstruc : char **
    A pointer to a character array that is used as pseudo dot-bracket output  
lines : const char **
    A 2-dimensional character array with the extension lines from the FASTA input  
option : unsigned int
    The option flags that define the behavior and recognition pattern of this function  

See Also
--------
RNA.file_fasta_read_record(), RNA.CONSTRAINT_DB_PIPE, RNA.CONSTRAINT_DB_DOT,
RNA.CONSTRAINT_DB_XRNA.CONSTRAINT_DB_ANG_BRACK, RNA.CONSTRAINT_DB_RND_BRACK  
";

%feature("docstring") extract_record_rest_structure "
";

%feature("docstring") read_record "

Get a data record from stdin.  

.. deprecated:: 2.6.4
    This function is deprecated! Use RNA.file_fasta_read_record() as a replacment.  
";

%feature("docstring") get_multi_input_line "
";

%feature("docstring") OPTION_MULTILINE "

Tell a function that an input is assumed to span several lines.  

If used as input-option a function might also be returning this state telling that it has read data
from multiple lines.  

See Also
--------
RNA.extract_record_rest_structure(), RNA.file_fasta_read_record()  
";

%feature("docstring") CONSTRAINT_MULTILINE "

parse multiline constraint  

.. deprecated:: 2.6.4
    see RNA.extract_record_rest_structure()  
";

%feature("docstring") INPUT_VERBOSE "
";

// File: group__file__formats__msa.xml

%feature("docstring") my_file_msa_read "

Read a multiple sequence alignment from file.  

This function reads the (first) multiple sequence alignment from an input file. The read alignment
is split into the sequence id/name part and the actual sequence information and stored in memory as
arrays of ids/names and sequences. If the alignment file format allows for additional information,
such as an ID of the entire alignment or consensus structure information, this data is retrieved as
well and made available. The `options` parameter allows to specify the set of alignment file formats
that should be used to retrieve the data. If 0 is passed as option, the list of alignment file
formats defaults to RNA.FILE_FORMAT_MSA_DEFAULT.  

Currently, the list of parsable multiple sequence alignment file formats consists of:  

*   msa-formats-clustal  
*   msa-formats-stockholm  
*   msa-formats-fasta  
*   msa-formats-maf  

**SWIG Wrapper Notes**
    In the target scripting language, only the first and last argument, `filename` and `options`,
    are passed to the corresponding function. The other arguments, which serve as output in the
    C-library, are available as additional return values. This function exists as an overloaded
    version where the `options` parameter may be omitted! In that case, the `options` parameter
    defaults to RNA.FILE_FORMAT_MSA_STOCKHOLM. See, e.g.   :py:func:`RNA.file_msa_read()` in the
    :doc:`/api_python`  and   :ref:`examples/python:parsing alignments`  in the Python examples.  

Parameters
----------
filename : const char *
    The name of input file that contains the alignment  
names : char ***
    An address to the pointer where sequence identifiers should be written to  
aln : char ***
    An address to the pointer where aligned sequences should be written to  
id : char **
    An address to the pointer where the alignment ID should be written to (Maybe NULL)  
structure : char **
    An address to the pointer where consensus structure information should be written to (Maybe
    NULL)  
options : unsigned int
    Options to manipulate the behavior of this function  

Returns
-------
int  
    The number of sequences in the alignment, or -1 if no alignment record could be found  

See Also
--------
RNA.file_msa_read_record(), RNA.FILE_FORMAT_MSA_CLUSTAL, RNA.FILE_FORMAT_MSA_STOCKHOLM,
RNA.FILE_FORMAT_MSA_FASTA, RNA.FILE_FORMAT_MSA_MAF, RNA.FILE_FORMAT_MSA_DEFAULT,
RNA.FILE_FORMAT_MSA_NOCHECK  

Note
----
After successfully reading an alignment, this function performs a validation of the data that
includes uniqueness of the sequence identifiers, and equal sequence lengths. This check can be
deactivated by passing RNA.FILE_FORMAT_MSA_NOCHECK in the `options` parameter.  
 It is the users responsibility to free any memory occupied by the output arguments `names`, `aln`,
`id`, and `structure` after calling this function. The function automatically sets the latter two
arguments to `NULL` in case no corresponding data could be retrieved from the input alignment.  
";

%feature("docstring") my_file_msa_read_record "

Read a multiple sequence alignment from file handle.  

Similar to RNA.file_msa_read(), this function reads a multiple sequence alignment from an input
file handle. Since using a file handle, this function is not limited to the first alignment record,
but allows for looping over all alignments within the input.  

The read alignment is split into the sequence id/name part and the actual sequence information and
stored in memory as arrays of ids/names and sequences. If the alignment file format allows for
additional information, such as an ID of the entire alignment or consensus structure information,
this data is retrieved as well and made available. The `options` parameter allows to specify the
alignment file format used to retrieve the data. A single format must be specified here, see
RNA.file_msa_detect_format() for helping to determine the correct MSA file format.  

Currently, the list of parsable multiple sequence alignment file formats consists of:  

*   msa-formats-clustal  
*   msa-formats-stockholm  
*   msa-formats-fasta  
*   msa-formats-maf  

**SWIG Wrapper Notes**
    In the target scripting language, only the first and last argument, `fp` and `options`, are
    passed to the corresponding function. The other arguments, which serve as output in the
    C-library, are available as additional return values. This function exists as an overloaded
    version where the `options` parameter may be omitted! In that case, the `options` parameter
    defaults to RNA.FILE_FORMAT_MSA_STOCKHOLM. See, e.g.   :py:func:`RNA.file_msa_read_record()` in
    the :doc:`/api_python`  and   :ref:`examples/python:parsing alignments`  in the Python examples.  

Parameters
----------
fp : FILE *
    The file pointer the data will be retrieved from  
names : char ***
    An address to the pointer where sequence identifiers should be written to  
aln : char ***
    An address to the pointer where aligned sequences should be written to  
id : char **
    An address to the pointer where the alignment ID should be written to (Maybe NULL)  
structure : char **
    An address to the pointer where consensus structure information should be written to (Maybe
    NULL)  
options : unsigned int
    Options to manipulate the behavior of this function  

Returns
-------
int  
    The number of sequences in the alignment, or -1 if no alignment record could be found  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_detect_format(), RNA.FILE_FORMAT_MSA_CLUSTAL,
RNA.FILE_FORMAT_MSA_STOCKHOLM, RNA.FILE_FORMAT_MSA_FASTA, RNA.FILE_FORMAT_MSA_MAF,
RNA.FILE_FORMAT_MSA_DEFAULT, RNA.FILE_FORMAT_MSA_NOCHECK  

Note
----
After successfully reading an alignment, this function performs a validation of the data that
includes uniqueness of the sequence identifiers, and equal sequence lengths. This check can be
deactivated by passing RNA.FILE_FORMAT_MSA_NOCHECK in the `options` parameter.  
 It is the users responsibility to free any memory occupied by the output arguments `names`, `aln`,
`id`, and `structure` after calling this function. The function automatically sets the latter two
arguments to `NULL` in case no corresponding data could be retrieved from the input alignment.  
";

%feature("docstring") my_file_msa_detect_format "

Detect the format of a multiple sequence alignment file.  

This function attempts to determine the format of a file that supposedly contains a multiple
sequence alignment (MSA). This is useful in cases where a MSA file contains more than a single
record and therefore RNA.file_msa_read() can not be applied, since it only retrieves the first.
Here, one can try to guess the correct file format using this function and then loop over the file,
record by record using one of the low-level record retrieval functions for the corresponding MSA
file format.  

**SWIG Wrapper Notes**
    This function exists as an overloaded version where the `options` parameter may be omitted! In
    that case, the `options` parameter defaults to RNA.FILE_FORMAT_MSA_DEFAULT. See, e.g.
    :py:func:`RNA.file_msa_detect_format()` in the :doc:`/api_python` .  

Parameters
----------
filename : const char *
    The name of input file that contains the alignment  
options : unsigned int
    Options to manipulate the behavior of this function  

Returns
-------
unsigned int  
    The MSA file format, or RNA.FILE_FORMAT_MSA_UNKNOWN  

See Also
--------
RNA.file_msa_read(), RNA.file_stockholm_read_record(), RNA.file_clustal_read_record(),
RNA.file_fasta_read_record()  

Note
----
This function parses the entire first record within the specified file. As a result, it returns
RNA.FILE_FORMAT_MSA_UNKNOWN not only if it can't detect the file's format, but also in cases where
the file doesn't contain sequences!  
";

%feature("docstring") my_file_msa_write "

Write multiple sequence alignment file.  

**SWIG Wrapper Notes**
    In the target scripting language, this function exists as a set of overloaded versions, where
    the last four parameters may be omitted. If the `options` parameter is missing the options
    default to (RNA.FILE_FORMAT_MSA_STOCKHOLM | RNA.FILE_FORMAT_MSA_APPEND). See, e.g.
    :py:func:`RNA.file_msa_write()` in the :doc:`/api_python` .  

Parameters
----------
filename : const char *
    The output filename  
names : const char **
    The array of sequence names / identifies  
aln : const char **
    The array of aligned sequences  
id : const char *
    An optional ID for the alignment  
structure : const char *
    An optional consensus structure  
source : const char *
    A string describing the source of the alignment  
options : unsigned int
    Options to manipulate the behavior of this function  

Returns
-------
int  
    Non-null upon successfully writing the alignment to file  

See Also
--------
RNA.FILE_FORMAT_MSA_STOCKHOLM, RNA.FILE_FORMAT_MSA_APPEND, RNA.FILE_FORMAT_MSA_MIS  

Note
----
Currently, we only support msa-formats-stockholm output  
";

%feature("docstring") FILE_FORMAT_MSA_CLUSTAL "

Option flag indicating ClustalW formatted files.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record(), RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_STOCKHOLM "

Option flag indicating Stockholm 1.0 formatted files.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record(), RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_FASTA "

Option flag indicating FASTA (Pearson) formatted files.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record(), RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_MAF "

Option flag indicating MAF formatted files.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record(), RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_MIS "

Option flag indicating most informative sequence (MIS) output.  

The default reference sequence output for an alignment is simply a consensus sequence. This flag
allows to write the most informative equence (MIS) instead.  

See Also
--------
RNA.file_msa_write()  
";

%feature("docstring") FILE_FORMAT_MSA_DEFAULT "

Option flag indicating the set of default file formats.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record(), RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_NOCHECK "

Option flag to disable validation of the alignment.  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record()  
";

%feature("docstring") FILE_FORMAT_MSA_UNKNOWN "

Return flag of RNA.file_msa_detect_format() to indicate unknown or malformatted alignment.  

See Also
--------
RNA.file_msa_detect_format()  
";

%feature("docstring") FILE_FORMAT_MSA_APPEND "

Option flag indicating to append data to a multiple sequence alignment file rather than overwriting
it.  

See Also
--------
RNA.file_msa_write()  
";

%feature("docstring") FILE_FORMAT_MSA_QUIET "

Option flag to suppress unnecessary spam messages on `stderr`  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record()  
";

%feature("docstring") FILE_FORMAT_MSA_SILENT "

Option flag to completely silence any warnings on `stderr`  

See Also
--------
RNA.file_msa_read(), RNA.file_msa_read_record()  
";

// File: group__command__files.xml

%feature("docstring") my_file_commands_read "

Extract a list of commands from a command file.  

Read a list of commands specified in the input file and return them as list of abstract commands  

**SWIG Wrapper Notes**
    This function is available as global function `file_commands_read()`. See, e.g.
    :py:func:`RNA.file_commands_read()` in the :doc:`/api_python` .  

Parameters
----------
filename : const char *
    The filename  
options : unsigned int
    Options to limit the type of commands read from the file  

Returns
-------
RNA.cmd()  
    A list of abstract commands  

See Also
--------
RNA.fold_compound.commands_apply(), RNA.file_commands_apply(), RNA.commands_free()  
";

%feature("docstring") vrna_file_commands_apply "

Apply a list of commands from a command file.  

This function is a shortcut to directly parse a commands file and apply all successfully parsed
commands to a RNA.fold_compound() data structure. It is the same as:  

**SWIG Wrapper Notes**
    This function is attached as method `file_commands_apply()` to objects of type `fold_compound`.
    See, e.g.   :py:meth:`RNA.fold_compound.file_commands_apply()` in the :doc:`/api_python` .  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() the command list will be applied to  
filename : const char *
    The filename  
options : unsigned int
    Options to limit the type of commands read from the file  

Returns
-------
int  
    The number of commands successfully applied  
";

%feature("docstring") vrna_fold_compound_t::commands_apply "

Apply a list of commands to a RNA.fold_compound().  

**SWIG Wrapper Notes**
    This function is attached as method `commands_apply()` to objects of type `fold_compound`. See,
    e.g.   :py:meth:`RNA.fold_compound.commands_apply()` in the :doc:`/api_python` .  

Parameters
----------
commands : RNA.cmd()
    The commands to apply  
options : unsigned int
    Options to limit the type of commands read from the file  

Returns
-------
int  
    The number of commands successfully applied  
";

%feature("docstring") vrna_commands_free "

Free memory occupied by a list of commands.  

Release memory occupied by a list of commands  

Parameters
----------
commands : RNA.cmd()
    A pointer to a list of commands  
";

%feature("docstring") CMD_PARSE_HC "

Command parse/apply flag indicating hard constraints.  

See Also
--------
RNA.cmd(), RNA.file_commands_read(), RNA.file_commands_apply(), RNA.fold_compound.commands_apply()  
";

%feature("docstring") CMD_PARSE_SC "

Command parse/apply flag indicating soft constraints.  

See Also
--------
RNA.cmd(), RNA.file_commands_read(), RNA.file_commands_apply(), RNA.fold_compound.commands_apply()  
";

%feature("docstring") CMD_PARSE_UD "

Command parse/apply flag indicating unstructured domains.  

See Also
--------
RNA.cmd(), RNA.file_commands_read(), RNA.file_commands_apply(), RNA.fold_compound.commands_apply()  
";

%feature("docstring") CMD_PARSE_SD "

Command parse/apply flag indicating structured domains.  

See Also
--------
RNA.cmd(), RNA.file_commands_read(), RNA.file_commands_apply(), RNA.fold_compound.commands_apply()  
";

%feature("docstring") CMD_PARSE_DEFAULTS "

Command parse/apply flag indicating default set of commands.  

See Also
--------
RNA.cmd(), RNA.file_commands_read(), RNA.file_commands_apply(), RNA.fold_compound.commands_apply()  
";

%feature("docstring") CMD_PARSE_SILENT "
";

// File: group__plotting__utils.xml

%feature("docstring") vrna_file_PS_rnaplot "

Produce a secondary structure graph in PostScript and write it to 'filename'.  

Note that this function has changed from previous versions and now expects the structure to be
plotted in dot-bracket notation as an argument. It does not make use of the global base_pair array
anymore.  

Parameters
----------
seq : const char *
    The RNA sequence  
structure : const char *
    The secondary structure in dot-bracket notation  
file : const char *
    The filename of the postscript output  
md_p : RNA.md() *
    Model parameters used to generate a commandline option string in the output (Maybe NULL)  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") vrna_file_PS_rnaplot_a "

Produce a secondary structure graph in PostScript including additional annotation macros and write
it to 'filename'.  

Same as RNA.file_PS_rnaplot() but adds extra PostScript macros for various annotations (see
generated PS code). The 'pre' and 'post' variables contain PostScript code that is verbatim copied
in the resulting PS file just before and after the structure plot. If both arguments ('pre' and
'post') are NULL, no additional macros will be printed into the PostScript.  

Parameters
----------
seq : const char *
    The RNA sequence  
structure : const char *
    The secondary structure in dot-bracket notation  
file : const char *
    The filename of the postscript output  
pre : const char *
    PostScript code to appear before the secondary structure plot  
post : const char *
    PostScript code to appear after the secondary structure plot  
md_p : RNA.md() *
    Model parameters used to generate a commandline option string in the output (Maybe NULL)  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") vrna_file_PS_rnaplot_layout "
";

%feature("docstring") PS_rna_plot_snoop_a "
";

%feature("docstring") gmlRNA "

Produce a secondary structure graph in Graph Meta Language (gml) and write it to a file.  

If 'option' is an uppercase letter the RNA sequence is used to label nodes, if 'option' equals *'X'*
or *'x'* the resulting file will coordinates for an initial layout of the graph.  

Parameters
----------
string : char *
    The RNA sequence  
structure : char *
    The secondary structure in dot-bracket notation  
ssfile : char *
    The filename of the gml output  
option : char
    The option flag  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") ssv_rna_plot "

Produce a secondary structure graph in SStructView format.  

Write coord file for SStructView  

Parameters
----------
string : char *
    The RNA sequence  
structure : char *
    The secondary structure in dot-bracket notation  
ssfile : char *
    The filename of the ssv output  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") svg_rna_plot "

Produce a secondary structure plot in SVG format and write it to a file.  

Parameters
----------
string : char *
    The RNA sequence  
structure : char *
    The secondary structure in dot-bracket notation  
ssfile : char *
    The filename of the svg output  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") xrna_plot "

Produce a secondary structure plot for further editing in XRNA.  

Parameters
----------
string : char *
    The RNA sequence  
structure : char *
    The secondary structure in dot-bracket notation  
ssfile : char *
    The filename of the xrna output  

Returns
-------
int  
    1 on success, 0 otherwise  
";

%feature("docstring") PS_rna_plot "

Produce a secondary structure graph in PostScript and write it to 'filename'.  

.. deprecated:: 2.6.4
    Use RNA.file_PS_rnaplot() instead!  
";

%feature("docstring") PS_rna_plot_a "

Produce a secondary structure graph in PostScript including additional annotation macros and write
it to 'filename'.  

.. deprecated:: 2.6.4
    Use RNA.file_PS_rnaplot_a() instead!  
";

%feature("docstring") PS_rna_plot_a_gquad "

Produce a secondary structure graph in PostScript including additional annotation macros and write
it to 'filename' (detect and draw g-quadruplexes)  

.. deprecated:: 2.6.4
    Use RNA.file_PS_rnaplot_a() instead!  
";

// File: group__plot__layout__utils.xml

%feature("docstring") vrna_plot_layout "

Create a layout (coordinates, etc.) for a secondary structure plot.  

This function can be used to create a secondary structure nucleotide layout that is then further
processed by an actual plotting function. The layout algorithm can be specified using the
`plot_type` parameter, and the following algorithms are currently supported:  

*   RNA.PLOT_TYPE_SIMPLE  
*   RNA.PLOT_TYPE_NAVIEW  
*   RNA.PLOT_TYPE_CIRCULAR  
*   RNA.PLOT_TYPE_TURTLE  
*   RNA.PLOT_TYPE_PUZZLER  

Passing an unsupported selection leads to the default algorithm RNA.PLOT_TYPE_NAVIEW  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
plot_type : unsigned int
    The layout algorithm to be used  

Returns
-------
RNA.plot_layout() *  
    The layout data structure for the provided secondary structure  

See Also
--------
RNA.plot_layout_free(), RNA.plot_layout_simple(), RNA.plot_layout_naview(),
RNA.plot_layout_circular(), RNA.plot_layout_turtle(), RNA.plot_layout_puzzler(),
RNA.plot_coords(), RNA.file_PS_rnaplot_layout()  

Note
----
If only X-Y coordinates of the corresponding structure layout are required, consider using
RNA.plot_coords() instead!  
";

%feature("docstring") vrna_plot_layout_simple "

Create a layout (coordinates, etc.) for a *simple* secondary structure plot.  

This function basically is a wrapper to RNA.plot_layout() that passes the
`plot_type`RNA.PLOT_TYPE_SIMPLE.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
RNA.plot_layout() *  
    The layout data structure for the provided secondary structure  

See Also
--------
RNA.plot_layout_free(), RNA.plot_layout(), RNA.plot_layout_naview(), RNA.plot_layout_circular(),
RNA.plot_layout_turtle(), RNA.plot_layout_puzzler(), RNA.plot_coords_simple(),
RNA.file_PS_rnaplot_layout()  

Note
----
If only X-Y coordinates of the corresponding structure layout are required, consider using
RNA.plot_coords_simple() instead!  
";

%feature("docstring") vrna_plot_layout_circular "

Create a layout (coordinates, etc.) for a *circular* secondary structure plot.  

This function basically is a wrapper to RNA.plot_layout() that passes the
`plot_type`RNA.PLOT_TYPE_CIRCULAR.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
RNA.plot_layout() *  
    The layout data structure for the provided secondary structure  

See Also
--------
RNA.plot_layout_free(), RNA.plot_layout(), RNA.plot_layout_naview(), RNA.plot_layout_simple(),
RNA.plot_layout_turtle(), RNA.plot_layout_puzzler(), RNA.plot_coords_circular(),
RNA.file_PS_rnaplot_layout()  

Note
----
If only X-Y coordinates of the corresponding structure layout are required, consider using
RNA.plot_coords_circular() instead!  
";

%feature("docstring") vrna_plot_layout_turtle "

Create a layout (coordinates, etc.) for a secondary structure plot using the *Turtle Algorithm*
:cite:p:`wiegreffe:2018` .  

This function basically is a wrapper to RNA.plot_layout() that passes the
`plot_type`RNA.PLOT_TYPE_TURTLE.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
RNA.plot_layout() *  
    The layout data structure for the provided secondary structure  

See Also
--------
RNA.plot_layout_free(), RNA.plot_layout(), RNA.plot_layout_simple(), RNA.plot_layout_circular(),
RNA.plot_layout_naview(), RNA.plot_layout_puzzler(), RNA.plot_coords_turtle(),
RNA.file_PS_rnaplot_layout()  

Note
----
If only X-Y coordinates of the corresponding structure layout are required, consider using
RNA.plot_coords_turtle() instead!  
";

%feature("docstring") vrna_plot_layout_puzzler "

Create a layout (coordinates, etc.) for a secondary structure plot using the *RNApuzzler Algorithm*
:cite:p:`wiegreffe:2018` .  

This function basically is a wrapper to RNA.plot_layout() that passes the
`plot_type`RNA.PLOT_TYPE_PUZZLER.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
RNA.plot_layout() *  
    The layout data structure for the provided secondary structure  

See Also
--------
RNA.plot_layout_free(), RNA.plot_layout(), RNA.plot_layout_simple(), RNA.plot_layout_circular(),
RNA.plot_layout_naview(), RNA.plot_layout_turtle(), RNA.plot_coords_puzzler(),
RNA.file_PS_rnaplot_layout()  

Note
----
If only X-Y coordinates of the corresponding structure layout are required, consider using
RNA.plot_coords_puzzler() instead!  
";

%feature("docstring") vrna_plot_layout_free "

Free memory occupied by a figure layout data structure.  

Parameters
----------
layout : RNA.plot_layout() *
    The layout data structure to free  

See Also
--------
RNA.plot_layout(), RNA.plot_layout(), RNA.plot_layout_simple(), RNA.plot_layout_circular(),
RNA.plot_layout_naview(), RNA.plot_layout_turtle(), RNA.plot_layout_puzzler(),
RNA.file_PS_rnaplot_layout()  
";

%feature("docstring") get_xy_coordinates "

Compute nucleotide coordinates for secondary structure plot.  

This function takes a secondary structure and computes X-Y coordinates for each nucleotide that then
can be used to create a structure plot. The parameter `plot_type` is used to select the underlying
layout algorithm. Currently, the following selections are provided:  

*   RNA.PLOT_TYPE_SIMPLE  
*   RNA.PLOT_TYPE_NAVIEW  
*   RNA.PLOT_TYPE_CIRCULAR  
*   RNA.PLOT_TYPE_TURTLE  
*   RNA.PLOT_TYPE_PUZZLER  

Passing an unsupported selection leads to the default algorithm RNA.PLOT_TYPE_NAVIEW  

Here is a simple example how to use this function, assuming variable `structure` contains a valid
dot-bracket string:  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
plot_type : int
    The layout algorithm to be used  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords_pt(), RNA.plot_coords_simple(), RNA.plot_coords_naview()
RNA.plot_coords_circular(), RNA.plot_coords_turtle(), RNA.plot_coords_puzzler()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_pt "

Compute nucleotide coordinates for secondary structure plot.  

Same as RNA.plot_coords() but takes a pair table with the structure information as input.  

Parameters
----------
pt : const short *
    The pair table that holds the secondary structure  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
plot_type : int
    The layout algorithm to be used  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords(), RNA.plot_coords_simple_pt(), RNA.plot_coords_naview_pt()
RNA.plot_coords_circular_pt(), RNA.plot_coords_turtle_pt(), RNA.plot_coords_puzzler_pt()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_simple "

Compute nucleotide coordinates for secondary structure plot the *Simple way*  

This function basically is a wrapper to RNA.plot_coords() that passes the
`plot_type`RNA.PLOT_TYPE_SIMPLE.  

Here is a simple example how to use this function, assuming variable `structure` contains a valid
dot-bracket string:  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords(), RNA.plot_coords_simple_pt(), RNA.plot_coords_circular(),
RNA.plot_coords_naview(), RNA.plot_coords_turtle(), RNA.plot_coords_puzzler()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_simple_pt "

Compute nucleotide coordinates for secondary structure plot the *Simple way*  

Same as RNA.plot_coords_simple() but takes a pair table with the structure information as input.  

Parameters
----------
pt : const short *
    The pair table that holds the secondary structure  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords_pt(), RNA.plot_coords_simple(), RNA.plot_coords_circular_pt(),
RNA.plot_coords_naview_pt(), RNA.plot_coords_turtle_pt(), RNA.plot_coords_puzzler_pt()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_circular "

Compute coordinates of nucleotides mapped in equal distancies onto a unit circle.  

This function basically is a wrapper to RNA.plot_coords() that passes the
`plot_type`RNA.PLOT_TYPE_CIRCULAR.  

In order to draw nice arcs using quadratic bezier curves that connect base pairs one may calculate a
second tangential point :math:`P^t` in addition to the actual R2 coordinates. the simplest way to do
so may be to compute a radius scaling factor :math:`rs` in the interval :math:`[0,1]` that weights
the proportion of base pair span to the actual length of the sequence. This scaling factor can then
be used to calculate the coordinates for :math:`P^t`, i.e.  

.. math::

  P^{t}_{x}[i] = X[i] * rs  

and  

.. math::

  P^{t}_{y}[i] = Y[i] * rs
.  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords(), RNA.plot_coords_circular_pt(), RNA.plot_coords_simple(),
RNA.plot_coords_naview(), RNA.plot_coords_turtle(), RNA.plot_coords_puzzler()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_circular_pt "

Compute nucleotide coordinates for a *Circular Plot*  

Same as RNA.plot_coords_circular() but takes a pair table with the structure information as input.  

Parameters
----------
pt : const short *
    The pair table that holds the secondary structure  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords_pt(), RNA.plot_coords_circular(), RNA.plot_coords_simple_pt(),
RNA.plot_coords_naview_pt(), RNA.plot_coords_turtle_pt(), RNA.plot_coords_puzzler_pt()  

Note
----
On success, this function allocates memory for X and Y coordinates and assigns the pointers at
addressess `x` and `y` to the corresponding memory locations. It's the users responsibility to
cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_puzzler "

Compute nucleotide coordinates for secondary structure plot using the *RNApuzzler* algorithm
:cite:p:`wiegreffe:2018` .  

This function basically is a wrapper to RNA.plot_coords() that passes the
`plot_type`RNA.PLOT_TYPE_PUZZLER.  

Here is a simple example how to use this function, assuming variable `structure` contains a valid
dot-bracket string and using the default options (`options` = NULL):  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
arc_coords : double **
    The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL
    on failure)  
options : RNA.plot_options_puzzler() *
    The options for the RNApuzzler algorithm (or NULL)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords(), RNA.plot_coords_puzzler_pt(), RNA.plot_coords_circular(),
RNA.plot_coords_simple(), RNA.plot_coords_turtle(), RNA.plot_coords_naview(),
RNA.plot_options_puzzler()  

Note
----
On success, this function allocates memory for X, Y and arc coordinates and assigns the pointers at
addressess `x`, `y` and `arc_coords` to the corresponding memory locations. It's the users
responsibility to cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_puzzler_pt "

Compute nucleotide coordinates for secondary structure plot using the *RNApuzzler* algorithm
:cite:p:`wiegreffe:2018` .  

Same as RNA.plot_coords_puzzler() but takes a pair table with the structure information as input.  

Parameters
----------
pt :
    The pair table that holds the secondary structure  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
arc_coords : double **
    The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL
    on failure)  
options :
    The options for the RNApuzzler algorithm (or NULL)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords_pt(), RNA.plot_coords_puzzler(), RNA.plot_coords_circular_pt(),
RNA.plot_coords_simple_pt(), RNA.plot_coords_turtle_pt(), RNA.plot_coords_naview_pt()  

Note
----
On success, this function allocates memory for X, Y and arc coordinates and assigns the pointers at
addressess `x`, `y` and `arc_coords` to the corresponding memory locations. It's the users
responsibility to cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_options_puzzler "

Create an RNApuzzler options data structure.  

Returns
-------
RNA.plot_options_puzzler() *  
    An RNApuzzler options data structure with default settings  

See Also
--------
RNA.plot_options_puzzler_free(), RNA.plot_coords_puzzler(), RNA.plot_coords_puzzler_pt(),
RNA.plot_layout_puzzler()  
";

%feature("docstring") vrna_plot_options_puzzler_free "

Free memory occupied by an RNApuzzler options data structure.  

Parameters
----------
options : RNA.plot_options_puzzler() *
    A pointer to the options data structure to free  

See Also
--------
RNA.plot_options_puzzler(), RNA.plot_coords_puzzler(), RNA.plot_coords_puzzler_pt(),
RNA.plot_layout_puzzler()  
";

%feature("docstring") vrna_plot_coords_turtle "

Compute nucleotide coordinates for secondary structure plot using the *RNAturtle* algorithm
:cite:p:`wiegreffe:2018` .  

This function basically is a wrapper to RNA.plot_coords() that passes the
`plot_type`RNA.PLOT_TYPE_TURTLE.  

Here is a simple example how to use this function, assuming variable `structure` contains a valid
dot-bracket string:  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
arc_coords : double **
    The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL
    on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords(), RNA.plot_coords_turtle_pt(), RNA.plot_coords_circular(),
RNA.plot_coords_simple(), RNA.plot_coords_naview(), RNA.plot_coords_puzzler()  

Note
----
On success, this function allocates memory for X, Y and arc coordinates and assigns the pointers at
addressess `x`, `y` and `arc_coords` to the corresponding memory locations. It's the users
responsibility to cleanup this memory after usage!  
";

%feature("docstring") vrna_plot_coords_turtle_pt "

Compute nucleotide coordinates for secondary structure plot using the *RNAturtle* algorithm
:cite:p:`wiegreffe:2018` .  

Same as RNA.plot_coords_turtle() but takes a pair table with the structure information as input.  

Parameters
----------
pt :
    The pair table that holds the secondary structure  
x : float **
    The address of a pointer of X coordinates (pointer will point to memory, or NULL on failure)  
y : float **
    The address of a pointer of Y coordinates (pointer will point to memory, or NULL on failure)  
arc_coords : double **
    The address of a pointer that will hold arc coordinates (pointer will point to memory, or NULL
    on failure)  

Returns
-------
int  
    The length of the structure on success, 0 otherwise  

See Also
--------
RNA.plot_coords_pt(), RNA.plot_coords_turtle(), RNA.plot_coords_circular_pt(),
RNA.plot_coords_simple_pt(), RNA.plot_coords_puzzler_pt(), RNA.plot_coords_naview_pt()  

Note
----
On success, this function allocates memory for X, Y and arc coordinates and assigns the pointers at
addressess `x`, `y` and `arc_coords` to the corresponding memory locations. It's the users
responsibility to cleanup this memory after usage!  
";

%feature("docstring") PLOT_TYPE_SIMPLE "

Definition of Plot type *simple*  

This is the plot type definition for several RNA structure plotting functions telling them to use
**Simple** plotting algorithm  

See Also
--------
rna_plot_type, RNA.file_PS_rnaplot_a(), RNA.file_PS_rnaplot(), svg_rna_plot(), gmlRNA(),
ssv_rna_plot(), xrna_plot()  
";

%feature("docstring") PLOT_TYPE_NAVIEW "

Definition of Plot type *Naview*  

This is the plot type definition for several RNA structure plotting functions telling them to use
**Naview** plotting algorithm   :cite:p:`bruccoleri:1988` .  

See Also
--------
rna_plot_type, RNA.file_PS_rnaplot_a(), RNA.file_PS_rnaplot(), svg_rna_plot(), gmlRNA(),
ssv_rna_plot(), xrna_plot()  
";

%feature("docstring") PLOT_TYPE_CIRCULAR "

Definition of Plot type *Circular*  

This is the plot type definition for several RNA structure plotting functions telling them to
produce a **Circular plot**  

See Also
--------
rna_plot_type, RNA.file_PS_rnaplot_a(), RNA.file_PS_rnaplot(), svg_rna_plot(), gmlRNA(),
ssv_rna_plot(), xrna_plot()  
";

%feature("docstring") PLOT_TYPE_TURTLE "

Definition of Plot type *Turtle* :cite:p:`wiegreffe:2018` .  
";

%feature("docstring") PLOT_TYPE_PUZZLER "

Definition of Plot type *RNApuzzler* :cite:p:`wiegreffe:2018` .  
";

%feature("docstring") PLOT_TYPE_DEFAULT "
";

// File: group__annotation__utils.xml

%feature("docstring") vrna_annotate_covar_db "

Produce covariance annotation for an alignment given a secondary structure.  
";

%feature("docstring") vrna_annotate_covar_db_extended "
";

%feature("docstring") vrna_annotate_covar_pairs "

Produce covariance annotation for an alignment given a set of base pairs.  
";

// File: group__plot__probabilities.xml

%feature("docstring") plot_dp_EPS "

Produce an encapsulate PostScript (EPS) dot-plot from one or two lists of base pair probabilities.  

This function reads two RNA.ep() lists `upper` and `lower` (e.g. base pair probabilities and a
secondary structure) and produces an EPS \"dot plot\" with filename `'filename'` where data from
`upper` is placed in the upper-triangular and data from `lower` is placed in the lower triangular
part of the matrix.  
 For default output, provide the flag RNA.PLOT_PROBABILITIES_DEFAULT as `options` parameter.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `plot_dp_EPS()` where the last three
    parameters may be omitted. The default values for these parameters are `lower` = NULL, `auxdata`
    = NULL, `options` = RNA.PLOT_PROBABILITIES_DEFAULT. See, e.g.  :py:func:`RNA.plot_dp_EPS()` in
    the :doc:`/api_python`.  

Parameters
----------
filename : const char *
    A filename for the EPS output  
sequence : const char *
    The RNA sequence  
upper : RNA.ep() *
    The base pair probabilities for the upper triangular part  
lower : RNA.ep() *
    The base pair probabilities for the lower triangular part  
options : unsigned int
    Options indicating which of the input data should be included in the dot-plot  

Returns
-------
int  
    1 if EPS file was successfully written, 0 otherwise  

See Also
--------
RNA.plist(), RNA.fold_compound.plist_from_probs(), RNA.PLOT_PROBABILITIES_DEFAULT  
";

%feature("docstring") vrna_plot_dp_PS_list "

Produce a postscript dot-plot from two pair lists.  

This function reads two plist structures (e.g. base pair probabilities and a secondary structure) as
produced by RNA.fold_compound.plist_from_probs() and RNA.plist() and produces a postscript \"dot plot\" that is
written to 'filename'.  
Using base pair probabilities in the first and mfe structure in the second plist, the resulting
\"dot plot\" represents each base pairing probability by a square of corresponding area in a upper
triangle matrix. The lower part of the matrix contains the minimum free energy structure.  

Parameters
----------
seq : char *
    The RNA sequence  
filename : char *
    A filename for the postscript output  
pl : RNA.ep() *
    The base pair probability pairlist  
mf : RNA.ep() *
    The mfe secondary structure pairlist  
comment : char *
    A comment  

Returns
-------
int  
    1 if postscript was successfully written, 0 otherwise  

See Also
--------
RNA.fold_compound.plist_from_probs(), RNA.plist()  
";

%feature("docstring") PLOT_PROBABILITIES_BP "

Option flag for base pair probabilities in probability plot output functions.  
";

%feature("docstring") PLOT_PROBABILITIES_ACC "

Option flag for accessibilities in probability plot output functions.  
";

%feature("docstring") PLOT_PROBABILITIES_UD "

Option flag for unstructured domain probabilities in probability plot output functions.  
";

%feature("docstring") PLOT_PROBABILITIES_UD_LIN "

Option flag for unstructured domain probabilities (linear representation) in probability plot output
functions.  
";

%feature("docstring") PLOT_PROBABILITIES_SD "

Option flag for structured domain probabilities (such as G-quadruplexes) in probability plot output
functions.  
";

%feature("docstring") PLOT_PROBABILITIES_SC_MOTIF "

Option flag for soft-constraint motif probabilities in probability plot output functions.  
";

%feature("docstring") PLOT_PROBABILITIES_SC_UP "
";

%feature("docstring") PLOT_PROBABILITIES_SC_BP "
";

%feature("docstring") PLOT_PROBABILITIES_DEFAULT "

Default option flag for probability plot output functions.  

Default output includes actual base pair probabilties (RNA.PLOT_PROBABILITIES_BP), structured
domain probabilities such as G-quadruplexes (RNA.PLOT_PROBABILITIES_SD), probabilities obtained
from soft-constraint motif implementation (RNA.PLOT_PROBABILITIES_SC_MOTIF), and unstructured
domain probabilities (RNA.PLOT_PROBABILITIES_UD_LIN).  

See Also
--------
RNA.plot_dp_EPS()  
";

// File: group__alignment__plots.xml

%feature("docstring") vrna_file_PS_aln "

Create an annotated PostScript alignment plot.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `file_PS_aln()` with three additional
    parameters `start`, `end`, and `offset` before the `columns` argument. Thus, it resembles the
    `RNA.file_PS_aln_slice()` function. The last four arguments may be omitted, indicating the
    default of `start` = *0*, `end` = *0*, `offset` = *0*, and `columns` = *60*. See, e.g.
    :py:func:`RNA.file_PS_aln()` in the :doc:`/api_python`.  

Parameters
----------
filename : const char *
    The output file name  
seqs : const char **
    The aligned sequences  
names : const char **
    The names of the sequences  
structure : const char *
    The consensus structure in dot-bracket notation  
columns : unsigned int
    The number of columns before the alignment is wrapped as a new block (a value of 0 indicates no
    wrapping)  

See Also
--------
RNA.file_PS_aln_slice()  
";

%feature("docstring") file_PS_aln "

Create an annotated PostScript alignment plot.  

Similar to RNA.file_PS_aln() but allows the user to print a particular slice of the alignment by
specifying a `start` and `end` position. The additional `offset` parameter allows for adjusting the
alignment position ruler value.  

**SWIG Wrapper Notes**
    This function is available as overloaded function `file_PS_aln()` where the last four parameter
    may be omitted, indicating `start` = *0*, `end` = *0*, `offset` = *0*, and `columns` = *60*.
    See, e.g.  :py:func:`RNA.file_PS_aln()` in the :doc:`/api_python`.  

Parameters
----------
filename : const char *
    The output file name  
seqs : const char **
    The aligned sequences  
names : const char **
    The names of the sequences  
structure : const char *
    The consensus structure in dot-bracket notation  
start : unsigned int
    The start of the alignment slice (a value of 0 indicates the first position of the alignment,
    i.e. no slicing at 5' side)  
end : unsigned int
    The end of the alignment slice (a value of 0 indicates the last position of the alignment, i.e.
    no slicing at 3' side)  
offset : int
    The alignment coordinate offset for the position ruler.  
columns : unsigned int
    The number of columns before the alignment is wrapped as a new block (a value of 0 indicates no
    wrapping)  

See Also
--------
RNA.file_PS_aln_slice()  
";

// File: group__search__utils.xml

%feature("docstring") vrna_search_BMH_num "

Search for a string of elements in a larger string of elements using the Boyer-Moore-Horspool
algorithm.  

To speed-up subsequent searches with this function, the Bad Character Table should be precomputed
and passed as argument `badchars`.  

Parameters
----------
needle : const unsigned int *
    The pattern of object representations to search for  
needle_size : size()
    The size (length) of the pattern provided in `needle`  
haystack : const unsigned int *
    The string of objects the search will be performed on  
haystack_size : size()
    The size (length) of the `haystack` string  
start : size()
    The position within `haystack` where to start the search  
badchars : size() *
    A pre-computed Bad Character Table obtained from RNA.search_BM_BCT_num() (If NULL, a Bad
    Character Table will be generated automatically)  
cyclic : unsigned char
    Allow for cyclic matches if non-zero, stop search at end of haystack otherwise  

Returns
-------
const unsigned int *  
    A pointer to the first occurence of `needle` within `haystack` after position `start`  

See Also
--------
RNA.search_BM_BCT_num(), RNA.search_BMH()  
";

%feature("docstring") vrna_search_BMH "

Search for an ASCII pattern within a larger ASCII string using the Boyer-Moore-Horspool algorithm.  

To speed-up subsequent searches with this function, the Bad Character Table should be precomputed
and passed as argument `badchars`. Furthermore, both, the lengths of `needle` and the length of
`haystack` should be pre-computed and must be passed along with each call.  

Parameters
----------
needle : const char *
    The NULL-terminated ASCII pattern to search for  
needle_size : size()
    The size (length) of the pattern provided in `needle`  
haystack : const char *
    The NULL-terminated ASCII string of the search will be performed on  
haystack_size : size()
    The size (length) of the `haystack` string  
start : size()
    The position within `haystack` where to start the search  
badchars : size() *
    A pre-computed Bad Character Table obtained from RNA.search_BM_BCT() (If NULL, a Bad Character
    Table will be generated automatically)  
cyclic : unsigned char
    Allow for cyclic matches if non-zero, stop search at end of haystack otherwise  

Returns
-------
const char *  
    A pointer to the first occurence of `needle` within `haystack` after position `start`  

See Also
--------
RNA.search_BM_BCT(), RNA.search_BMH_num()  
";

%feature("docstring") vrna_search_BM_BCT_num "

Retrieve a Boyer-Moore Bad Character Table for a pattern of elements represented by natural numbers.  

Parameters
----------
pattern : const unsigned int *
    The pattern of element representations used in the subsequent search  
pattern_size : size()
    The size (length) of the pattern provided in `pattern`  
num_max : unsigned int
    The maximum number representation of an element, i.e. the size of the alphabet  

Returns
-------
size() *  
    A Bad Character Table for use in our Boyer-Moore search algorithm implementation(s)  

See Also
--------
RNA.search_BMH_num(), RNA.search_BM_BCT()  

Note
----
We store the maximum number representation of an element `num_max` at position `0`. So the actual
bad character table `T` starts at `T`[1] for an element represented by number `0`.  
";

%feature("docstring") vrna_search_BM_BCT "

Retrieve a Boyer-Moore Bad Character Table for a NULL-terminated pattern of ASCII characters.  

Parameters
----------
pattern : const char *
    The NULL-terminated pattern of ASCII characters used in the subsequent search  

Returns
-------
size() *  
    A Bad Character Table for use in our Boyer-Moore search algorithm implementation(s)  

See Also
--------
RNA.search_BMH(), RNA.search_BM_BCT_num()  

Note
----
We store the maximum number representation of an element, i.e. `127` at position `0`. So the actual
bad character table `T` starts at `T`[1] for an element represented by ASCII code `0`.  
";

// File: group__combinatorics__utils.xml

%feature("docstring") my_enumerate_necklaces "

Enumerate all necklaces with fixed content.  

This function implements *A fast algorithm to generate necklaces with fixed content* as published by
:cite:t:`sawada:2003` .  

The function receives a list of counts (the elements on the necklace) for each type of object within
a necklace. The list starts at index 0 and ends with an entry that has a count of 0. The algorithm
then enumerates all non-cyclic permutations of the content, returned as a list of necklaces. This
list, again, is zero-terminated, i.e. the last entry of the list is a `NULL` pointer.  

**SWIG Wrapper Notes**
    This function is available as global function `enumerate_necklaces()` which accepts lists input,
    an produces list of lists output. See, e.g.   :py:func:`RNA.enumerate_necklaces()` in the
    :doc:`/api_python` .  

Parameters
----------
type_counts : const unsigned int *
    A 0-terminated list of entity counts  

Returns
-------
unsigned int **  
    A list of all non-cyclic permutations of the entities  
";

%feature("docstring") vrna_rotational_symmetry_num "

Determine the order of rotational symmetry for a string of objects represented by natural numbers.  

The algorithm applies a fast search of the provided string within itself, assuming the end of the
string wraps around to connect with it's start. For example, a string of the form `011011` has
rotational symmetry of order `2`  

This is a simplified version of RNA.rotational_symmetry_pos_num() that may be useful if one is only
interested in the degree of rotational symmetry but not the actual set of rotational symmetric
strings.  

**SWIG Wrapper Notes**
    This function is available as global function `rotational_symmetry()`. See
    RNA.rotational_symmetry_pos() for details. Note, that in the target language the length of the
    list `string` is always known a-priori, so the parameter `string_length` must be omitted. See,
    e.g.   :py:func:`RNA.rotational_symmetry()` in the :doc:`/api_python` .  

Parameters
----------
string : const unsigned int *
    The string of elements encoded as natural numbers  
string_length : size()
    The length of the string  

Returns
-------
unsigned int  
    The order of rotational symmetry  

See Also
--------
RNA.rotational_symmetry_pos_num(), RNA.rotationa_symmetry()  
";

%feature("docstring") vrna_rotational_symmetry_pos_num "

Determine the order of rotational symmetry for a string of objects represented by natural numbers.  

The algorithm applies a fast search of the provided string within itself, assuming the end of the
string wraps around to connect with it's start. For example, a string of the form `011011` has
rotational symmetry of order `2`  

If the argument `positions` is not `NULL`, the function stores an array of string start positions
for rotational shifts that map the string back onto itself. This array has length of order of
rotational symmetry, i.e. the number returned by this function. The first element `positions`[0]
always contains a shift value of `0` representing the trivial rotation.  

**SWIG Wrapper Notes**
    This function is available as global function `rotational_symmetry()`. See
    RNA.rotational_symmetry_pos() for details. Note, that in the target language the length of the
    list `string` is always known a-priori, so the parameter `string_length` must be omitted. See,
    e.g.   :py:func:`RNA.rotational_symmetry()` in the :doc:`/api_python` .  

Parameters
----------
string : const unsigned int *
    The string of elements encoded as natural numbers  
string_length : size()
    The length of the string  
positions : unsigned int **
    A pointer to an (undefined) list of alternative string start positions that lead to an identity
    mapping (may be NULL)  

Returns
-------
unsigned int  
    The order of rotational symmetry  

See Also
--------
RNA.rotational_symmetry_num(), RNA.rotational_symmetry(), RNA.rotational_symmetry_pos()  

Note
----
Do not forget to release the memory occupied by `positions` after a successful execution of this
function.  
";

%feature("docstring") vrna_rotational_symmetry "

Determine the order of rotational symmetry for a NULL-terminated string of ASCII characters.  

The algorithm applies a fast search of the provided string within itself, assuming the end of the
string wraps around to connect with it's start. For example, a string of the form `AABAAB` has
rotational symmetry of order `2`  

This is a simplified version of RNA.rotational_symmetry_pos() that may be useful if one is only
interested in the degree of rotational symmetry but not the actual set of rotational symmetric
strings.  

**SWIG Wrapper Notes**
    This function is available as global function `rotational_symmetry()`. See
    RNA.rotational_symmetry_pos() for details. See, e.g.   :py:func:`RNA.rotational_symmetry()` in
    the :doc:`/api_python` .  

Parameters
----------
string : const char *
    A NULL-terminated string of characters  

Returns
-------
unsigned int  
    The order of rotational symmetry  

See Also
--------
RNA.rotational_symmetry_pos(), RNA.rotationa_symmetry_num()  
";

%feature("docstring") my_rotational_symmetry "

Determine the order of rotational symmetry for a NULL-terminated string of ASCII characters.  

The algorithm applies a fast search of the provided string within itself, assuming the end of the
string wraps around to connect with it's start. For example, a string of the form `AABAAB` has
rotational symmetry of order `2`  

If the argument `positions` is not `NULL`, the function stores an array of string start positions
for rotational shifts that map the string back onto itself. This array has length of order of
rotational symmetry, i.e. the number returned by this function. The first element `positions`[0]
always contains a shift value of `0` representing the trivial rotation.  

**SWIG Wrapper Notes**
    This function is available as overloaded global function `rotational_symmetry()`. It merges the
    functionalities of RNA.rotational_symmetry(), RNA.rotational_symmetry_pos(),
    RNA.rotational_symmetry_num(), and RNA.rotational_symmetry_pos_num(). In contrast to our
    C-implementation, this function doesn't return the order of rotational symmetry as a single
    value, but returns a list of cyclic permutation shifts that result in a rotationally symmetric
    string. The length of the list then determines the order of rotational symmetry. See, e.g.
    :py:func:`RNA.rotational_symmetry()` in the :doc:`/api_python` .  

Parameters
----------
string : const char *
    A NULL-terminated string of characters  
positions : unsigned int **
    A pointer to an (undefined) list of alternative string start positions that lead to an identity
    mapping (may be NULL)  

Returns
-------
unsigned int  
    The order of rotational symmetry  

See Also
--------
RNA.rotational_symmetry(), RNA.rotational_symmetry_num(), RNA.rotational_symmetry_num_pos()  

Note
----
Do not forget to release the memory occupied by `positions` after a successful execution of this
function.  
";

%feature("docstring") vrna_rotational_symmetry_db "

Determine the order of rotational symmetry for a dot-bracket structure.  

Given a (permutation of multiple) RNA strand(s) and a particular secondary structure in dot-bracket
notation, compute the degree of rotational symmetry. In case there is only a single linear RNA
strand, the structure always has degree 1, as there are no rotational symmetries due to the
direction of the nucleic acid sequence and the fixed positions of 5' and 3' ends. However, for
circular RNAs, rotational symmetries might arise if the sequence consists of a concatenation of
:math:`k` identical subsequences.  

This is a simplified version of RNA.fold_compound.rotational_symmetry_db() that may be useful if one is only
interested in the degree of rotational symmetry but not the actual set of rotational symmetric
strings.  

**SWIG Wrapper Notes**
    This function is attached as method `rotational_symmetry_db()` to objects of type
    `fold_compound` (i.e. RNA.fold_compound()). See RNA.fold_compound.rotational_symmetry_db() for
    details.
    See, e.g.   :py:meth:`RNA.fold_compound.rotational_symmetry_db()` in the :doc:`/api_python` .  

Parameters
----------
fc : RNA.fold_compound() *
    A fold_compound data structure containing the nucleic acid sequence(s), their order, and model
    settings  
structure : const char *
    The dot-bracket structure the degree of rotational symmetry is checked for  

Returns
-------
unsigned int  
    The degree of rotational symmetry of the `structure` (0 in case of any errors)  

See Also
--------
RNA.fold_compound.rotational_symmetry_db(), RNA.rotational_symmetry(), RNA.rotational_symmetry_num()  
";

%feature("docstring") vrna_fold_compound_t::rotational_symmetry_db "

Determine the order of rotational symmetry for a dot-bracket structure.  

Given a (permutation of multiple) RNA strand(s) and a particular secondary structure in dot-bracket
notation, compute the degree of rotational symmetry. In case there is only a single linear RNA
strand, the structure always has degree 1, as there are no rotational symmetries due to the
direction of the nucleic acid sequence and the fixed positions of 5' and 3' ends. However, for
circular RNAs, rotational symmetries might arise if the sequence consists of a concatenation of
:math:`k` identical subsequences.  

If the argument `positions` is not `NULL`, the function stores an array of string start positions
for rotational shifts that map the string back onto itself. This array has length of order of
rotational symmetry, i.e. the number returned by this function. The first element `positions`[0]
always contains a shift value of `0` representing the trivial rotation.  

**SWIG Wrapper Notes**
    This function is attached as method `rotational_symmetry_db()` to objects of type
    `fold_compound` (i.e. RNA.fold_compound()). Thus, the first argument must be omitted. In
    contrast to our C-implementation, this function doesn't simply return the order of rotational
    symmetry of the secondary structure, but returns the list `position` of cyclic permutation
    shifts that result in a rotationally symmetric structure. The length of the list then determines
    the order of rotational symmetry. See, e.g.
    :py:meth:`RNA.fold_compound.rotational_symmetry_db()` in the :doc:`/api_python` .  

Parameters
----------
structure : const char *
    The dot-bracket structure the degree of rotational symmetry is checked for  
positions : unsigned int **
    A pointer to an (undefined) list of alternative string start positions that lead to an identity
    mapping (may be NULL)  

Returns
-------
unsigned int  
    The degree of rotational symmetry of the `structure` (0 in case of any errors)  

See Also
--------
RNA.rotational_symmetry_db(), RNA.rotational_symmetry_pos(), RNA.rotational_symmetry_pos_num()  

Note
----
Do not forget to release the memory occupied by `positions` after a successful execution of this
function.  
";

%feature("docstring") vrna_n_multichoose_k "

Obtain a list of k-combinations with repetition (n multichoose k)  

This function compiles a list of k-combinations, or k-multicombination, i.e. a list of multisubsets
of size k from a set of integer values from 0 to n - 1. For that purpose, we enumerate n + k - 1
choose k and decrease each index position i by i to obtain n multichoose k.  

Parameters
----------
n : size()
    Maximum number to choose from (interval of integers from 0 to `n` - 1)  
k : size()
    Number of elements to choose, i.e. size of each multisubset  

Returns
-------
unsigned int **  
    A list of lists of elements of combinations (last entry is terminated by **NULL**  
";

%feature("docstring") boustrophedon "

Generate a sequence of Boustrophedon distributed numbers.  

This function generates a sequence of positive natural numbers within the interval :math:`[start,
end]` in a Boustrophedon fashion. That is, the numbers :math:`start, \\ldots, end` in the resulting
list are alternating between left and right ends of the interval while progressing to the inside,
i.e. the list consists of a sequence of natural numbers of the form:  

.. math::

  start, end, start + 1, end - 1, start + 2, end - 2, \\ldots  

The resulting list is 1-based and contains the length of the sequence of numbers at it's 0-th
position.  

Upon failure, the function returns **NULL**  

**SWIG Wrapper Notes**
    This function is available as overloaded global function `boustrophedon()`. See, e.g.
    :py:func:`RNA.boustrophedon()` in the :doc:`/api_python` .  

Parameters
----------
start : size()
    The first number of the list (left side of the interval)  
end : size()
    The last number of the list (right side of the interval)  

Returns
-------
unsigned int *  
    A list of alternating numbers from the interval :math:`[start, end]` (or **NULL** on error)  

See Also
--------
RNA.boustrophedon_pos()  
";

%feature("docstring") vrna_boustrophedon_pos "

Obtain the i-th element in a Boustrophedon distributed interval of natural numbers.  

**SWIG Wrapper Notes**
    This function is available as overloaded global function `boustrophedon()`. Omitting the `pos`
    argument yields the entire sequence from `start` to `end`. See, e.g.
    :py:func:`RNA.boustrophedon()` in the :doc:`/api_python` .  

Parameters
----------
start : size()
    The first number of the list (left side of the interval)  
end : size()
    The last number of the list (right side of the interval)  
pos : size()
    The index of the number within the Boustrophedon distributed sequence (1-based)  

Returns
-------
unsigned int  
    The `pos-th` element in the Boustrophedon distributed sequence of natural numbers of the
    interval  

See Also
--------
RNA.boustrophedon()  
";

// File: group__data__structures.xml

%feature("docstring") vrna_C11_features "

Dummy symbol to check whether the library was build using C11/C++11 features.  

By default, several data structures of our new v3.0 API use C11/C++11 features, such as unnamed
unions, unnamed structs. However, these features can be deactivated at compile time to allow
building the library and executables with compilers that do not support these features.  

Now, the problem arises that once our static library is compiled and a third-party application is
supposed to link against it, it needs to know, at compile time, how to correctly address particular
data structures. This is usually implicitely taken care of through the API exposed in our header
files. Unfortunately, we had some preprocessor directives in our header files that changed the API
depending on the capabilities of the compiler the third-party application is build with. This in
turn prohibited the use of an RNAlib compiled without C11/C++11 support in a program that
compiles/links with enabled C11/C++11 support and vice-versa.  

Therefore, we introduce this dummy symbol which can be used to check, whether the static library was
build with C11/C++11 features.  

since: v2.2.9  

Note
----
If the symbol is present, the library was build with enabled C11/C++11 features support and no
action is required. However, if the symbol is missing in RNAlib >= 2.2.9, programs that link to
RNAlib must define a pre-processor identifier *RNA.DISABLE_C11_FEATURES* before including any
ViennaRNA Package header file, for instance by adding a *CPPFLAG*  
";

// File: group__message__utils.xml

%feature("docstring") vrna_message_error "

Print an error message and die.  

This function is a wrapper to *fprintf(stderr, ...)* that puts a capital **ERROR:** in front of the
message and then exits the calling program.  

Parameters
----------
format : const char *
    The error message to be printed  
... :
    Optional arguments for the formatted message string  

See Also
--------
RNA.message_verror(), RNA.message_warning(), RNA.message_info()  
";

%feature("docstring") vrna_message_verror "

Print an error message and die.  

This function is a wrapper to *vfprintf(stderr, ...)* that puts a capital **ERROR:** in front of the
message and then exits the calling program.  

Parameters
----------
format : const char *
    The error message to be printed  
args : va_list
    The argument list for the formatted message string  

See Also
--------
RNA.message_error(), RNA.message_warning(), RNA.message_info()  
";

%feature("docstring") vrna_message_warning "

Print a warning message.  

This function is a wrapper to *fprintf(stderr, ...)* that puts a capital **WARNING:** in front of
the message.  

Parameters
----------
format : const char *
    The warning message to be printed  
... :
    Optional arguments for the formatted message string  

See Also
--------
RNA.message_vwarning(), RNA.message_error(), RNA.message_info()  
";

%feature("docstring") vrna_message_vwarning "

Print a warning message.  

This function is a wrapper to *fprintf(stderr, ...)* that puts a capital **WARNING:** in front of
the message.  

Parameters
----------
format : const char *
    The warning message to be printed  
args : va_list
    The argument list for the formatted message string  

See Also
--------
RNA.message_vwarning(), RNA.message_error(), RNA.message_info()  
";

%feature("docstring") vrna_message_info "

Print an info message.  

This function is a wrapper to *fprintf(...)*.  

Parameters
----------
fp : FILE *
    The file pointer where the message is printed to  
format : const char *
    The warning message to be printed  
... :
    Optional arguments for the formatted message string  

See Also
--------
RNA.message_vinfo(), RNA.message_error(), RNA.message_warning()  
";

%feature("docstring") vrna_message_vinfo "

Print an info message.  

This function is a wrapper to *fprintf(...)*.  

Parameters
----------
fp : FILE *
    The file pointer where the message is printed to  
format : const char *
    The info message to be printed  
args : va_list
    The argument list for the formatted message string  

See Also
--------
RNA.message_vinfo(), RNA.message_error(), RNA.message_warning()  
";

%feature("docstring") vrna_message_input_seq_simple "

Print a line to *stdout* that asks for an input sequence.  

There will also be a ruler (scale line) printed that helps orientation of the sequence positions  
";

%feature("docstring") vrna_message_input_seq "

Print a line with a user defined string and a ruler to stdout.  

(usually this is used to ask for user input) There will also be a ruler (scale line) printed that
helps orientation of the sequence positions  

Parameters
----------
s : const char *
    A user defined string that will be printed to stdout  
";

%feature("docstring") vrna_message_input_msa "
";

// File: group__units.xml

%feature("docstring") vrna_convert_energy "

Convert between energy / work units.  

Parameters
----------
energy : double
    Input energy value  
from : RNA.unit_energy
    Input unit  
to : RNA.unit_energy
    Output unit  

Returns
-------
double  
    Energy value in Output unit  

See Also
--------
RNA.unit_energy  
";

%feature("docstring") vrna_convert_temperature "

Convert between temperature units.  

Parameters
----------
temp : double
    Input temperature value  
from : RNA.unit_temperature
    Input unit  
to : RNA.unit_temperature
    Output unit  

Returns
-------
double  
    Temperature value in Output unit  

See Also
--------
RNA.unit_temperature  
";

%feature("docstring") vrna_convert_kcal_to_dcal "

Convert floating point energy value into integer representation.  

This function converts a floating point value in kcal/mol into its corresponding deka-cal/mol
integer representation as used throughout RNAlib.  

Parameters
----------
energy : double
    The energy value in kcal/mol  

Returns
-------
int  
    The energy value in deka-cal/mol  

See Also
--------
RNA.convert_dcal_to_kcal()  
";

%feature("docstring") vrna_convert_dcal_to_kcal "

Convert an integer representation of free energy in deka-cal/mol to kcal/mol.  

This function converts a free energy value given as integer in deka-cal/mol into the corresponding
floating point number in kcal/mol  

Parameters
----------
energy : int
    The energy in deka-cal/mol  

Returns
-------
double  
    The energy in kcal/mol  

See Also
--------
RNA.convert_kcal_to_dcal()  
";

// File: group__fold__compound.xml

%feature("docstring") vrna_fold_compound "

Retrieve a RNA.fold_compound() data structure for single sequences and hybridizing sequences.  

This function provides an easy interface to obtain a prefilled RNA.fold_compound() by passing a
single sequence, or two contatenated sequences as input. For the latter, sequences need to be
seperated by an '&' character like this:  

    char *sequence = \"GGGG&CCCC\";  

The optional parameter `md_p` can be used to specify the model details for successive computations
based on the content of the generated RNA.fold_compound(). Passing NULL will instruct the function
to use default model details. The third parameter `options` may be used to specify dynamic
programming (DP) matrix requirements.  

#### Options  

*   RNA.OPTION_DEFAULT -  Option flag to specify default settings/requirements.  
*   RNA.OPTION_MFE -  Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices
    and corresponding set of energy parameters.  
*   RNA.OPTION_PF -  Option flag to specify requirement of Partition Function (PF) DP matrices and
    corresponding set of Boltzmann factors.  
*   RNA.OPTION_WINDOW -  Option flag to specify requirement of DP matrices for local folding
    approaches.  

The above options may be OR-ed together.  

If you just need the folding compound serving as a container for your data, you can simply pass
RNA.OPTION_DEFAULT to the `option` parameter. This creates a RNA.fold_compound() without DP
matrices, thus saving memory. Subsequent calls of any structure prediction function will then take
care of allocating the memory required for the DP matrices. If you only intend to evaluate
structures instead of actually predicting them, you may use the RNA.OPTION_EVAL_ONLY macro. This
will seriously speedup the creation of the RNA.fold_compound().  

Parameters
----------
sequence : const char *
    A single sequence, or two concatenated sequences seperated by an '&' character  
md_p : const RNA.md() *
    An optional set of model details  
options : unsigned int
    The options for DP matrices memory allocation  

Returns
-------
RNA.fold_compound() *  
    A prefilled RNA.fold_compound() ready to be used for computations (may be `NULL` on error)  

See Also
--------
RNA.fold_compound_free(), RNA.fold_compound_comparative(), RNA.md()  

Note
----
The sequence string must be uppercase, and should contain only RNA (resp. DNA) alphabet depending on
what energy parameter set is used  
";

%feature("docstring") vrna_fold_compound_comparative "

Retrieve a RNA.fold_compound() data structure for sequence alignments.  

This function provides an easy interface to obtain a prefilled RNA.fold_compound() by passing an
alignment of sequences.  

The optional parameter `md_p` can be used to specify the model details for successive computations
based on the content of the generated RNA.fold_compound(). Passing NULL will instruct the function
to use default model details. The third parameter `options` may be used to specify dynamic
programming (DP) matrix requirements.  

#### Options  

*   RNA.OPTION_DEFAULT -  Option flag to specify default settings/requirements.  
*   RNA.OPTION_MFE -  Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices
    and corresponding set of energy parameters.  
*   RNA.OPTION_PF -  Option flag to specify requirement of Partition Function (PF) DP matrices and
    corresponding set of Boltzmann factors.  
*   RNA.OPTION_WINDOW -  Option flag to specify requirement of DP matrices for local folding
    approaches.  

The above options may be OR-ed together.  

If you just need the folding compound serving as a container for your data, you can simply pass
RNA.OPTION_DEFAULT to the `option` parameter. This creates a RNA.fold_compound() without DP
matrices, thus saving memory. Subsequent calls of any structure prediction function will then take
care of allocating the memory required for the DP matrices. If you only intend to evaluate
structures instead of actually predicting them, you may use the RNA.OPTION_EVAL_ONLY macro. This
will seriously speedup the creation of the RNA.fold_compound().  

Parameters
----------
sequences : const char **
    A sequence alignment including 'gap' characters  
md_p : RNA.md() *
    An optional set of model details  
options : unsigned int
    The options for DP matrices memory allocation  

Returns
-------
RNA.fold_compound() *  
    A prefilled RNA.fold_compound() ready to be used for computations (may be `NULL` on error)  

See Also
--------
RNA.fold_compound_free(), RNA.fold_compound(), RNA.md(), RNA.OPTION_MFE, RNA.OPTION_PF,
RNA.OPTION_EVAL_ONLY, read_clustal()  

Note
----
The sequence strings must be uppercase, and should contain only RNA (resp. DNA) alphabet including
gap characters depending on what energy parameter set is used.  
";

%feature("docstring") vrna_fold_compound_comparative2 "
";

%feature("docstring") vrna_fold_compound_TwoD "
";

%feature("docstring") vrna_fold_compound_prepare "
";

%feature("docstring") vrna_fold_compound_free "

Free memory occupied by a RNA.fold_compound().  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() that is to be erased from memory  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.mx_mfe_free(), RNA.mx_pf_free()  
";

%feature("docstring") vrna_fold_compound_t::add_auxdata "

Add auxiliary data to the RNA.fold_compound().  

This function allows one to bind arbitrary data to a RNA.fold_compound() which may later on be used
by one of the callback functions, e.g. RNA.recursion_status(). To allow for proper cleanup of the
memory occupied by this auxiliary data, the user may also provide a pointer to a cleanup function
that free's the corresponding memory. This function will be called automatically when the
RNA.fold_compound() is free'd with RNA.fold_compound_free().  

Parameters
----------
data : void *
    A pointer to an arbitrary data structure  
f : RNA.auxdata_free
    A pointer to function that free's memory occupied by the arbitrary data (May be NULL)  

See Also
--------
RNA.auxdata_free()  

Note
----
Before attaching the arbitrary data pointer, this function will call the RNA.auxdata_free() on
any pre-existing data that is already attached.  
";

%feature("docstring") vrna_fold_compound_t::add_callback "

Add a recursion status callback to the RNA.fold_compound().  

Binding a recursion status callback function to a RNA.fold_compound() allows one to perform
arbitrary operations just before, or after an actual recursive computations, e.g. MFE prediction, is
performed by the RNAlib. The callback function will be provided with a pointer to its
RNA.fold_compound(), and a status message. Hence, it has complete access to all variables that
incluence the recursive computations.  

Parameters
----------
f : RNA.recursion_status
    The pointer to the recursion status callback function  

See Also
--------
RNA.recursion_status(), RNA.fold_compound(), RNA.STATUS_MFE_PRE, RNA.STATUS_MFE_POST,
RNA.STATUS_PF_PRE, RNA.STATUS_PF_POST  
";

%feature("docstring") STATUS_MFE_PRE "

Status message indicating that MFE computations are about to begin.  

See Also
--------
RNA.fold_compound().stat_cb, RNA.recursion_status(), RNA.fold_compound.mfe(), RNA.fold(),
RNA.circfold(),
RNA.alifold(), RNA.circalifold(), RNA.cofold()  
";

%feature("docstring") STATUS_MFE_POST "

Status message indicating that MFE computations are finished.  

See Also
--------
RNA.fold_compound().stat_cb, RNA.recursion_status(), RNA.fold_compound.mfe(), RNA.fold(),
RNA.circfold(),
RNA.alifold(), RNA.circalifold(), RNA.cofold()  
";

%feature("docstring") STATUS_PF_PRE "

Status message indicating that Partition function computations are about to begin.  

See Also
--------
RNA.fold_compound().stat_cb, RNA.recursion_status(), RNA.fold_compound.pf()  
";

%feature("docstring") STATUS_PF_POST "

Status message indicating that Partition function computations are finished.  

See Also
--------
RNA.fold_compound().stat_cb, RNA.recursion_status(), RNA.fold_compound.pf()  
";

%feature("docstring") OPTION_DEFAULT "

Option flag to specify default settings/requirements.  
";

%feature("docstring") OPTION_MFE "

Option flag to specify requirement of Minimum Free Energy (MFE) DP matrices and corresponding set of
energy parameters.  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.OPTION_EVAL_ONLY  
";

%feature("docstring") OPTION_PF "

Option flag to specify requirement of Partition Function (PF) DP matrices and corresponding set of
Boltzmann factors.  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.OPTION_EVAL_ONLY  
";

%feature("docstring") OPTION_HYBRID "

Option flag to specify requirement of dimer DP matrices.  
";

%feature("docstring") OPTION_EVAL_ONLY "

Option flag to specify that neither MFE, nor PF DP matrices are required.  

Use this flag in conjuntion with RNA.OPTION_MFE, and RNA.OPTION_PF to save memory for a
RNA.fold_compound() obtained from RNA.fold_compound(), or RNA.fold_compound_comparative() in
cases where only energy evaluation but no structure prediction is required.  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.fold_compound.eval_structure()  
";

%feature("docstring") OPTION_WINDOW "

Option flag to specify requirement of DP matrices for local folding approaches.  
";

%feature("docstring") OPTION_F5 "
";

%feature("docstring") OPTION_F3 "
";

%feature("docstring") OPTION_WINDOW_F5 "
";

%feature("docstring") OPTION_WINDOW_F3 "
";

// File: group__dp__matrices.xml

%feature("docstring") vrna_mx_add "

Add Dynamic Programming (DP) matrices (allocate memory)  

This function adds DP matrices of a specific type to the provided RNA.fold_compound(), such that
successive DP recursion can be applied. The function caller has to specify which type of DP matrix
is requested, see RNA.mx_type, and what kind of recursive algorithm will be applied later on,
using the parameters type, and options, respectively. For the latter, Minimum free energy (MFE), and
Partition function (PF) computations are distinguished. A third option that may be passed is
RNA.OPTION_HYBRID, indicating that auxiliary DP arrays are required for RNA-RNA interaction
prediction.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() that holds pointers to the DP matrices  
type : RNA.mx_type
    The type of DP matrices requested  
options : unsigned int
    Option flags that specify the kind of DP matrices, such as MFE or PF arrays, and auxiliary
    requirements  

Returns
-------
int  
    1 if DP matrices were properly allocated and attached, 0 otherwise  

See Also
--------
RNA.mx_mfe_add(), RNA.mx_pf_add(), RNA.fold_compound(), RNA.fold_compound_comparative(),
RNA.fold_compound_free(), RNA.mx_pf_free(), RNA.mx_mfe_free(), RNA.mx_type, RNA.OPTION_MFE,
RNA.OPTION_PF, RNA.OPTION_HYBRID, RNA.OPTION_EVAL_ONLY  

Note
----
Usually, there is no need to call this function, since the constructors of RNA.fold_compound() are
handling all the DP matrix memory allocation.  
";

%feature("docstring") vrna_mx_mfe_add "
";

%feature("docstring") vrna_mx_pf_add "
";

%feature("docstring") vrna_mx_prepare "
";

%feature("docstring") vrna_mx_mfe_free "

Free memory occupied by the Minimum Free Energy (MFE) Dynamic Programming (DP) matrices.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() storing the MFE DP matrices that are to be erased from memory  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.fold_compound_free(), RNA.mx_pf_free()  
";

%feature("docstring") vrna_mx_pf_free "

Free memory occupied by the Partition Function (PF) Dynamic Programming (DP) matrices.  

Parameters
----------
fc : RNA.fold_compound() *
    The RNA.fold_compound() storing the PF DP matrices that are to be erased from memory  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.fold_compound_free(),
RNA.mx_mfe_free()  
";

// File: group__hash__table__utils.xml

%feature("docstring") vrna_ht_init "

Get an initialized hash table.  

This function returns a ready-to-use hash table with pre-allocated memory for a particular number of
entries.  

Note
----

If all function pointers are `NULL`, this function initializes the hash table with *default
functions*, i.e.  

*   RNA.ht_db_comp() for the `compare_function`,  
*   RNA.ht_db_hash_func() for the `hash_function`, and  
*   RNA.ht_db_free_entry() for the `free_hash_entry`  

arguments.  

Parameters
----------
b : unsigned int
    Number of bits for the hash table. This determines the size ( :math:`2^b -1`).  
compare_function : RNA.ht_cmp
    A function pointer to compare any two entries in the hash table (may be `NULL`)  
hash_function : RNA.ht_hashfunc
    A function pointer to retrieve the hash value of any entry (may be `NULL`)  
free_hash_entry : RNA.ht_free
    A function pointer to free the memory occupied by any entry (may be `NULL`)  

Returns
-------
RNA.hash_table()  
    An initialized, empty hash table, or `NULL` on any error  

Warnings
--------
If `hash_bits` is larger than 27 you have to compile it with the flag gcc -mcmodel=large.  
";

%feature("docstring") vrna_ht_size "

Get the size of the hash table.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  

Returns
-------
unsigned long  
    The size of the hash table, i.e. the maximum number of entries  
";

%feature("docstring") vrna_ht_collisions "

Get the number of collisions in the hash table.  

Parameters
----------
ht : struct RNA.hash_table() *
    The hash table  

Returns
-------
unsigned long  
    The number of collisions in the hash table  
";

%feature("docstring") vrna_ht_get "

Get an element from the hash table.  

This function takes an object `x` and performs a look-up whether the object is stored within the
hash table `ht`. If the object is already stored in `ht`, the function simply returns the entry,
otherwise it returns `NULL`.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  
x : void *
    The hash entry to look-up  

Returns
-------
void *  
    The entry `x` if it is stored in `ht`, `NULL` otherwise  

See Also
--------
RNA.ht_insert(), RNA.hash_delete(), RNA.ht_init()  
";

%feature("docstring") vrna_ht_insert "

Insert an object into a hash table.  

Writes the pointer to your hash entry into the table.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  
x : void *
    The hash entry  

Returns
-------
int  
    0 on success, 1 if the value is already in the hash table, -1 on error.  

Warnings
--------
In case of collisions, this function simply increments the hash key until a free entry in the hash
table is found.  

See Also
--------
RNA.ht_init(), RNA.hash_delete(), RNA.ht_clear()  
";

%feature("docstring") vrna_ht_remove "

Remove an object from the hash table.  

Deletes the pointer to your hash entry from the table.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  
x : void *
    The hash entry  

Note
----
This function doesn't free any memory occupied by the hash entry.  
";

%feature("docstring") vrna_ht_clear "

Clear the hash table.  

This function removes all entries from the hash table and automatically free's the memory occupied
by each entry using the bound RNA.ht_free() function.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  

See Also
--------
RNA.ht_free(), RNA.ht_init()  
";

%feature("docstring") vrna_ht_free "

Free all memory occupied by the hash table.  

This function removes all entries from the hash table by calling the RNA.ht_free() function for
each entry. Finally, the memory occupied by the hash table itself is free'd as well.  

Parameters
----------
ht : RNA.hash_table()
    The hash table  
";

%feature("docstring") vrna_ht_db_comp "

Default hash table entry comparison.  

This is the default comparison function for hash table entries. It assumes the both entries `x` and
`y` are of type RNA.ht_entry_db() and compares the `structure` attribute of both entries  

Parameters
----------
x : void *
    A hash table entry of type RNA.ht_entry_db()  
y : void *
    A hash table entry of type RNA.ht_entry_db()  

Returns
-------
int  
    -1 if x is smaller, +1 if x is larger than y. 0 if both are equal.  

See Also
--------
RNA.ht_entry_db(), RNA.ht_init(), RNA.ht_db_hash_func(), RNA.ht_db_free_entry()  
";

%feature("docstring") vrna_ht_db_hash_func "

Default hash function.  

This is the default hash function for hash table insertion/lookup. It assumes that entries are of
type RNA.ht_entry_db() and uses the Bob Jenkins 1996 mix function to create a hash key from the
`structure` attribute of the hash entry.  

Parameters
----------
x : void *
    A hash table entry to compute the key for  
hashtable_size : unsigned long
    The size of the hash table  

Returns
-------
unsigned int  
    The hash key for entry `x`  

See Also
--------
RNA.ht_entry_db(), RNA.ht_init(), RNA.ht_db_comp(), RNA.ht_db_free_entry()  
";

%feature("docstring") vrna_ht_db_free_entry "

Default function to free memory occupied by a hash entry.  

This function assumes that hash entries are of type RNA.ht_entry_db() and free's the memory
occupied by that entry.  

Parameters
----------
hash_entry : void *
    The hash entry to remove from memory  

Returns
-------
int  
    0 on success  

See Also
--------
RNA.ht_entry_db(), RNA.ht_init(), RNA.ht_db_comp(), RNA.ht_db_hash_func()  
";

// File: group__heap__utils.xml

%feature("docstring") vrna_heap_init "

Initialize a heap data structure.  

This function initializes a heap data structure. The implementation is based on a *min-heap*, i.e.
the minimal element is located at the root of the heap. However, by reversing the logic of the
compare function, one can easily transform this into a *max-heap* implementation.  

Beside the regular operations on a heap data structure, we implement removal and update of arbitrary
elements within the heap. For that purpose, however, one requires a reverse-index lookup system
that, (i) for a given element stores the current position in the heap, and (ii) allows for fast
lookup of an elements current position within the heap. The corresponding getter- and setter-
functions may be provided through the arguments `get_entry_pos` and `set_entry_pos`, respectively.  

Sometimes, it is difficult to simply compare two data structures without any context. Therefore, the
compare function is provided with a user-defined data pointer that can hold any context required.  

Parameters
----------
n : size()
    The initial size of the heap, i.e. the number of elements to store  
cmp : RNA.heap_cmp
    The address of a compare function that will be used to fullfill the partial order requirement  
get_entry_pos : RNA.heap_get_pos
    The address of a function that retrieves the position of an element within the heap (or NULL)  
set_entry_pos : RNA.heap_set_pos
    The address of a function that stores the position of an element within the heap (or NULL)  
data : void *
    An arbitrary data pointer passed through to the compare function `cmp`, and the set/get
    functions `get_entry_pos` / `set_entry_pos`  

Returns
-------
RNA.heap()  
    An initialized heap data structure, or NULL on error  

Warnings
--------
If any of the arguments `get_entry_pos` or `set_entry_pos` is NULL, the operations
RNA.heap_update() and RNA.heap_remove() won't work.  

See Also
--------
RNA.heap_free(), RNA.heap_insert(), RNA.heap_pop(), RNA.heap_top(), RNA.heap_remove(),
RNA.heap_update(), RNA.heap(), RNA.heap_cmp, RNA.heap_get_pos, RNA.heap_set_pos  
";

%feature("docstring") vrna_heap_free "

Free memory occupied by a heap data structure.  

Parameters
----------
h : RNA.heap()
    The heap that should be free'd  

See Also
--------
RNA.heap_init()  
";

%feature("docstring") vrna_heap_size "

Get the size of a heap data structure, i.e. the number of stored elements.  

Parameters
----------
h : struct RNA.heap() *
    The heap data structure  

Returns
-------
size()  
    The number of elements currently stored in the heap, or 0 upon any error  
";

%feature("docstring") vrna_heap_insert "

Insert an element into the heap.  

Parameters
----------
h : RNA.heap()
    The heap data structure  
v : void *
    A pointer to the object that is about to be inserted into the heap  

See Also
--------
RNA.heap_init(), RNA.heap_pop(), RNA.heap_top(), RNA.heap_free(), RNA.heap_remove(),
RNA.heap_update()  
";

%feature("docstring") vrna_heap_pop "

Pop (remove and return) the object at the root of the heap.  

This function removes the root from the heap and returns it to the caller.  

Parameters
----------
h : RNA.heap()
    The heap data structure  

Returns
-------
void *  
    The object at the root of the heap, i.e. the minimal element (or NULL if (a) the heap is empty
    or (b) any error occurred)  

See Also
--------
RNA.heap_init(), RNA.heap_top(), RNA.heap_insert(), RNA.heap_free()RNA.heap_remove(),
RNA.heap_update()  
";

%feature("docstring") vrna_heap_top "

Get the object at the root of the heap.  

Parameters
----------
h : RNA.heap()
    The heap data structure  

Returns
-------
const void *  
    The object at the root of the heap, i.e. the minimal element (or NULL if (a) the heap is empty
    or (b) any error occurred)  

See Also
--------
RNA.heap_init(), RNA.heap_pop(), RNA.heap_insert(), RNA.heap_free()RNA.heap_remove(),
RNA.heap_update()  
";

%feature("docstring") vrna_heap_remove "

Remove an arbitrary element within the heap.  

Parameters
----------
h : RNA.heap()
    The heap data structure  
v : const void *
    The object to remove from the heap  

Returns
-------
void *  
    The object that was removed from the heap (or NULL if (a) it wasn't found or (b) any error
    occurred)  

Warnings
--------
This function won't work if the heap was not properly initialized with callback functions for fast
reverse-index mapping!  

See Also
--------
RNA.heap_init(), RNA.heap_get_pos, RNA.heap_set_pos, RNA.heap_pop(), RNA.heap_free()  
";

%feature("docstring") vrna_heap_update "

Update an arbitrary element within the heap.  

Parameters
----------
h : RNA.heap()
    The heap data structure  
v : void *
    The object to update  

Returns
-------
void *  
    The 'previous' object within the heap that now got replaced by `v` (or NULL if (a) it wasn't
    found or (b) any error occurred)  

Warnings
--------
This function won't work if the heap was not properly initialized with callback functions for fast
reverse-index mapping!  

See Also
--------
RNA.heap_init(), RNA.heap_get_pos, RNA.heap_set_pos_fRNA.heap_pop(), RNA.heap_remove(),
RNA.heap_free()  

Note
----
If the object that is to be updated is not currently stored in the heap, it will be inserted. In
this case, the function returns NULL.  
";

// File: group__strings.xml

%feature("docstring") vrna_string_make "
";

%feature("docstring") vrna_string_free "
";

%feature("docstring") vrna_string_append "
";

%feature("docstring") vrna_string_append_cstring "
";

%feature("docstring") STRING_HEADER "
";

// File: group__array__utils.xml

%feature("docstring") vrna__array_set_capacity "

Explicitely set the capacity of an array.  

Note
----
Do not use this function. Rather resort to the RNA.array_set_capacity macro  
";

%feature("docstring") array "

Define an array.  
";

%feature("docstring") array_make "

Make an array `Name` of type `Type`.  
";

%feature("docstring") ARRAY_GROW_FORMULA "

The default growth formula for array.  
";

%feature("docstring") ARRAY_HEADER "

Retrieve a pointer to the header of an array `input`.  
";

%feature("docstring") array_size "

Get the number of elements of an array `input`.  
";

%feature("docstring") array_capacity "

Get the size of an array `input`, i.e. its actual capacity.  
";

%feature("docstring") array_set_capacity "

Explicitely set the capacity of an array `a`.  
";

%feature("docstring") array_init_size "

Initialize an array `a` with a particular pre-allocated size `init_size`.  
";

%feature("docstring") array_init "

Initialize an array `a`.  
";

%feature("docstring") array_free "

Release memory of an array `a`.  
";

%feature("docstring") array_append "

Safely append an item to an array `a`.  
";

%feature("docstring") array_grow "

Grow an array `a` to provide a minimum capacity `min_capacity`.  
";

// File: group__buffer__utils.xml

%feature("docstring") vrna_cstr "

Create a dynamic char * stream data structure.  

Parameters
----------
size : size()
    The initial size of the buffer in characters  
output : FILE *
    An optional output file stream handle that is used to write the collected data to (defaults to
    *stdout* if *NULL*)  

See Also
--------
RNA.cstr_free(), RNA.cstr_close(), RNA.cstr_fflush(), RNA.cstr_discard(), RNA.cstr_printf()  
";

%feature("docstring") vrna_cstr_discard "

Discard the current content of the dynamic char * stream data structure.  

Parameters
----------
buf : struct RNA.cstr() *
    The dynamic char * stream data structure to free  

See Also
--------
RNA.cstr_free(), RNA.cstr_close(), RNA.cstr_fflush(), RNA.cstr_printf()  
";

%feature("docstring") vrna_cstr_free "

Free the memory occupied by a dynamic char * stream data structure.  

This function first flushes any remaining character data within the stream and then free's the
memory occupied by the data structure.  

Parameters
----------
buf : RNA.cstr()
    The dynamic char * stream data structure to free  

See Also
--------
RNA.cstr_close(), RNA.cstr_fflush(), RNA.cstr()  
";

%feature("docstring") vrna_cstr_close "

Free the memory occupied by a dynamic char * stream and close the output stream.  

This function first flushes any remaining character data within the stream then closes the attached
output file stream (if any), and finally free's the memory occupied by the data structure.  

Parameters
----------
buf : RNA.cstr()
    The dynamic char * stream data structure to free  

See Also
--------
RNA.cstr_free(), RNA.cstr_fflush(), RNA.cstr()  
";

%feature("docstring") vrna_cstr_fflush "

Flush the dynamic char * output stream.  

This function flushes the collected char * stream, either by writing to the attached file handle, or
simply by writing to *stdout* if no file handle has been attached upon construction using
RNA.cstr().  

**Postcondition**
    The stream buffer is empty after execution of this function  

Parameters
----------
buf : struct RNA.cstr() *
    The dynamic char * stream data structure to flush  

See Also
--------
RNA.cstr(), RNA.cstr_close(), RNA.cstr_free()  
";

%feature("docstring") vrna_cstr_string "
";

%feature("docstring") vrna_cstr_vprintf "
";

%feature("docstring") vrna_cstr_printf "
";

%feature("docstring") vrna_cstr_message_info "
";

%feature("docstring") vrna_cstr_message_vinfo "
";

%feature("docstring") vrna_cstr_message_warning "
";

%feature("docstring") vrna_cstr_message_vwarning "
";

%feature("docstring") vrna_cstr_print_fasta_header "
";

%feature("docstring") vrna_cstr_printf_structure "
";

%feature("docstring") vrna_cstr_vprintf_structure "
";

%feature("docstring") vrna_cstr_printf_comment "
";

%feature("docstring") vrna_cstr_vprintf_comment "
";

%feature("docstring") vrna_cstr_printf_thead "
";

%feature("docstring") vrna_cstr_vprintf_thead "
";

%feature("docstring") vrna_cstr_printf_tbody "
";

%feature("docstring") vrna_cstr_vprintf_tbody "
";

%feature("docstring") vrna_cstr_print_eval_sd_corr "
";

%feature("docstring") vrna_cstr_print_eval_ext_loop "
";

%feature("docstring") vrna_cstr_print_eval_hp_loop "
";

%feature("docstring") vrna_cstr_print_eval_hp_loop_revert "
";

%feature("docstring") vrna_cstr_print_eval_int_loop "
";

%feature("docstring") vrna_cstr_print_eval_int_loop_revert "
";

%feature("docstring") vrna_cstr_print_eval_mb_loop "
";

%feature("docstring") vrna_cstr_print_eval_mb_loop_revert "
";

%feature("docstring") vrna_cstr_print_eval_gquad "
";

%feature("docstring") vrna_ostream_init "

Get an initialized ordered output stream.  

Parameters
----------
output : RNA.stream_output
    A callback function that processes and releases data in the stream  
auxdata : void *
    A pointer to auxiliary data passed as first argument to the `output` callback  

Returns
-------
RNA.ostream()  
    An initialized ordered output stream  

See Also
--------
RNA.ostream_free(), RNA.ostream_request(), RNA.ostream_provide()  
";

%feature("docstring") vrna_ostream_free "

Free an initialized ordered output stream.  

Parameters
----------
dat : RNA.ostream()
    The output stream for which occupied memory should be free'd  

See Also
--------
RNA.ostream_init()  
";

%feature("docstring") vrna_ostream_threadsafe "
";

%feature("docstring") vrna_ostream_request "

Request index in ordered output stream.  

This function must be called prior to RNA.ostream_provide() to indicate that data associted with a
certain index number is expected to be inserted into the stream in the future.  

Parameters
----------
dat : RNA.ostream()
    The output stream for which the index is requested  
num : unsigned int
    The index to request data for  

See Also
--------
RNA.ostream_init(), RNA.ostream_provide(), RNA.ostream_free()  
";

%feature("docstring") vrna_ostream_provide "

Provide output stream data for a particular index.  

**Precondition**
    The index data is provided for must have been requested using RNA.ostream_request() beforehand.  

Parameters
----------
dat : RNA.ostream()
    The output stream for which data is provided  
i : unsigned int
    The index of the provided data  
data : void *
    The data provided  

See Also
--------
RNA.ostream_request()  
";

// File: group__mfe__global__deprecated.xml

%feature("docstring") alifold "

Compute MFE and according consensus structure of an alignment of sequences.  

This function predicts the consensus structure for the aligned 'sequences' and returns the minimum
free energy; the mfe structure in bracket notation is returned in 'structure'.  

Sufficient space must be allocated for 'structure' before calling alifold().  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! Use RNA.alifold(), or RNA.fold_compound.mfe() instead!  

Parameters
----------
strings : const char **
    A pointer to a NULL terminated array of character arrays  
structure : char *
    A pointer to a character array that may contain a constraining consensus structure (will be
    overwritten by a consensus structure that exhibits the MFE)  

Returns
-------
float  
    The free energy score in kcal/mol  

See Also
--------
RNA.alifold(), RNA.fold_compound.mfe()  
";

%feature("docstring") circalifold "

Compute MFE and according structure of an alignment of sequences assuming the sequences are circular
instead of linear.  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! Use RNA.alicircfold(), and RNA.fold_compound.mfe()
instead!  

Parameters
----------
strings : const char **
    A pointer to a NULL terminated array of character arrays  
structure : char *
    A pointer to a character array that may contain a constraining consensus structure (will be
    overwritten by a consensus structure that exhibits the MFE)  

Returns
-------
float  
    The free energy score in kcal/mol  

See Also
--------
RNA.alicircfold(), RNA.alifold(), RNA.fold_compound.mfe()  
";

%feature("docstring") free_alifold_arrays "

Free the memory occupied by MFE alifold functions.  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! It only affects memory being free'd that was allocated by
    an old API function before. Release of memory occupied by the newly introduced
    RNA.fold_compound() is handled by RNA.fold_compound_free()  

See Also
--------
RNA.fold_compound_free()  
";

%feature("docstring") cofold "

Compute the minimum free energy of two interacting RNA molecules.  

The code is analog to the fold() function. If cut_point ==-1 results should be the same as with
fold().  

.. deprecated:: 2.6.4
    use RNA.fold_compound.mfe_dimer() instead  

Parameters
----------
sequence : const char *
    The two sequences concatenated  
structure : char *
    Will hold the barcket dot structure of the dimer molecule  

Returns
-------
float  
    minimum free energy of the structure  
";

%feature("docstring") cofold_par "

Compute the minimum free energy of two interacting RNA molecules.  

.. deprecated:: 2.6.4
    use RNA.fold_compound.mfe_dimer() instead  
";

%feature("docstring") free_co_arrays "

Free memory occupied by cofold()  

.. deprecated:: 2.6.4
    This function will only free memory allocated by a prior call of cofold() or cofold_par(). See
    RNA.fold_compound.mfe_dimer() for how to use the new API  

See Also
--------
RNA.fc_destroy(), RNA.fold_compound.mfe_dimer()  

Note
----
folding matrices now reside in the fold compound, and should be free'd there  
";

%feature("docstring") update_cofold_params "

Recalculate parameters.  

.. deprecated:: 2.6.4
    See RNA.fold_compound.params_subst() for an alternative using the new API  
";

%feature("docstring") update_cofold_params_par "

Recalculate parameters.  

.. deprecated:: 2.6.4
    See RNA.fold_compound.params_subst() for an alternative using the new API  
";

%feature("docstring") export_cofold_arrays_gq "

Export the arrays of partition function cofold (with gquadruplex support)  

Export the cofold arrays for use e.g. in the concentration Computations or suboptimal secondary
structure backtracking  

.. deprecated:: 2.6.4
    folding matrices now reside within the fold compound. Thus, this function will only work in
    conjunction with a prior call to cofold() or cofold_par()  

Parameters
----------
f5_p : int **
    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]  
c_p : int **
    A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that
    i pairs with j  
fML_p : int **
    A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any
    multiloop segment with at least one stem  
fM1_p : int **
    A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for
    multiloop segment with exactly one stem  
fc_p : int **
    A pointer to the 'fc' array, i.e. array ...  
ggg_p : int **
    A pointer to the 'ggg' array, i.e. array containing best free energy of a gquadruplex delimited
    by [i,j]  
indx_p : int **
    A pointer to the indexing array used for accessing the energy matrices  
ptype_p : char **
    A pointer to the ptype array containing the base pair types for each possibility (i,j)  

See Also
--------
RNA.fold_compound.mfe_dimer() for the new API  
";

%feature("docstring") export_cofold_arrays "

Export the arrays of partition function cofold.  

Export the cofold arrays for use e.g. in the concentration Computations or suboptimal secondary
structure backtracking  

.. deprecated:: 2.6.4
    folding matrices now reside within the RNA.fold_compound(). Thus, this function will only work
    in conjunction with a prior call to the deprecated functions cofold() or cofold_par()  

Parameters
----------
f5_p : int **
    A pointer to the 'f5' array, i.e. array conatining best free energy in interval [1,j]  
c_p : int **
    A pointer to the 'c' array, i.e. array containing best free energy in interval [i,j] given that
    i pairs with j  
fML_p : int **
    A pointer to the 'M' array, i.e. array containing best free energy in interval [i,j] for any
    multiloop segment with at least one stem  
fM1_p : int **
    A pointer to the 'M1' array, i.e. array containing best free energy in interval [i,j] for
    multiloop segment with exactly one stem  
fc_p : int **
    A pointer to the 'fc' array, i.e. array ...  
indx_p : int **
    A pointer to the indexing array used for accessing the energy matrices  
ptype_p : char **
    A pointer to the ptype array containing the base pair types for each possibility (i,j)  

See Also
--------
RNA.fold_compound.mfe_dimer() for the new API  
";

%feature("docstring") initialize_cofold "

allocate arrays for folding  

.. deprecated:: 2.6.4  
";

%feature("docstring") fold_par "

Compute minimum free energy and an appropriate secondary structure of an RNA sequence.  

The first parameter given, the RNA sequence, must be *uppercase* and should only contain an alphabet
:math:`\\Sigma` that is understood by the RNAlib  
(e.g. :math:`\\Sigma = \\{A,U,C,G\\}`)  
 The second parameter, *structure*, must always point to an allocated block of memory with a size of
at least :math:`\\mathrm{strlen}(\\mathrm{sequence})+1`  

If the third parameter is NULL, global model detail settings are assumed for the folding recursions.
Otherwise, the provided parameters are used.  

The fourth parameter indicates whether a secondary structure constraint in enhanced dot-bracket
notation is passed through the structure parameter or not. If so, the characters \" | x < > \" are
recognized to mark bases that are paired, unpaired, paired upstream, or downstream, respectively.
Matching brackets \" ( ) \" denote base pairs, dots \".\" are used for unconstrained bases.  

To indicate that the RNA sequence is circular and thus has to be post-processed, set the last
parameter to non-zero  

After a successful call of fold_par(), a backtracked secondary structure (in dot-bracket notation)
that exhibits the minimum of free energy will be written to the memory *structure* is pointing to.
The function returns the minimum of free energy for any fold of the sequence given.  

.. deprecated:: 2.6.4
    use RNA.fold_compound.mfe() instead!  

Note
----
OpenMP: Passing NULL to the 'parameters' argument involves access to several global model detail
variables and thus is not to be considered threadsafe  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  
parameters : RNA.param() *
    A data structure containing the pre-scaled energy contributions and the model details. (NULL may
    be passed, see OpenMP notes above)  
is_constrained : int
    Switch to indicate that a structure constraint is passed via the structure argument (0==off)  
is_circular : int
    Switch to (de-)activate post-processing steps in case RNA sequence is circular (0==off)  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
RNA.fold_compound.mfe(), fold(), circfold(), RNA.md(), set_energy_model(), get_scaled_parameters()  
";

%feature("docstring") fold "

Compute minimum free energy and an appropriate secondary structure of an RNA sequence.  

This function essentially does the same thing as fold_par(). However, it takes its model details,
i.e. temperature, dangles, tetra_loop, noGU, no_closingGU, fold_constrained, noLonelyPairs from the
current global settings within the library  

.. deprecated:: 2.6.4
    use RNA.fold(), or RNA.fold_compound.mfe() instead!  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
fold_par(), circfold()  
";

%feature("docstring") circfold "

Compute minimum free energy and an appropriate secondary structure of a circular RNA sequence.  

This function essentially does the same thing as fold_par(). However, it takes its model details,
i.e. temperature, dangles, tetra_loop, noGU, no_closingGU, fold_constrained, noLonelyPairs from the
current global settings within the library  

.. deprecated:: 2.6.4
    Use RNA.circfold(), or RNA.fold_compound.mfe() instead!  

Parameters
----------
sequence : const char *
    RNA sequence  
structure : char *
    A pointer to the character array where the secondary structure in dot-bracket notation will be
    written to  

Returns
-------
float  
    the minimum free energy (MFE) in kcal/mol  

See Also
--------
fold_par(), circfold()  
";

%feature("docstring") free_arrays "

Free arrays for mfe folding.  

.. deprecated:: 2.6.4
    See RNA.fold(), RNA.circfold(), or RNA.fold_compound.mfe() and RNA.fold_compound() for the usage
of the
    new API!  
";

%feature("docstring") update_fold_params "

Recalculate energy parameters.  

.. deprecated:: 2.6.4
    For non-default model settings use the new API with RNA.fold_compound.params_subst() and
RNA.fold_compound.mfe() instead!  
";

%feature("docstring") update_fold_params_par "

Recalculate energy parameters.  

.. deprecated:: 2.6.4
    For non-default model settings use the new API with RNA.fold_compound.params_subst() and
RNA.fold_compound.mfe() instead!  
";

%feature("docstring") export_fold_arrays "

.. deprecated:: 2.6.4
    See RNA.fold_compound.mfe() and RNA.fold_compound() for the usage of the new API!  
";

%feature("docstring") export_fold_arrays_par "

.. deprecated:: 2.6.4
    See RNA.fold_compound.mfe() and RNA.fold_compound() for the usage of the new API!  
";

%feature("docstring") export_circfold_arrays "

.. deprecated:: 2.6.4
    See RNA.fold_compound.mfe() and RNA.fold_compound() for the usage of the new API!  
";

%feature("docstring") export_circfold_arrays_par "

.. deprecated:: 2.6.4
    See RNA.fold_compound.mfe() and RNA.fold_compound() for the usage of the new API!  
";

%feature("docstring") LoopEnergy "

.. deprecated:: 2.6.4
    {This function is deprecated and will be removed soon. Use E_IntLoop() instead!}  
";

%feature("docstring") HairpinE "

.. deprecated:: 2.6.4
    {This function is deprecated and will be removed soon. Use E_Hairpin() instead!}  
";

%feature("docstring") initialize_fold "

Allocate arrays for folding  

.. deprecated:: 2.6.4
    See RNA.fold_compound.mfe() and RNA.fold_compound() for the usage of the new API!  
";

%feature("docstring") backtrack_fold_from_pair "
";

// File: group__mfe__window__deprecated.xml

%feature("docstring") Lfold "

The local analog to fold().  

Computes the minimum free energy structure including only base pairs with a span smaller than
'maxdist'  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.mfe_window() instead!  
";

%feature("docstring") Lfoldz "

.. deprecated:: 2.6.4
    Use RNA.fold_compound.mfe_window_zscore() instead!  
";

%feature("docstring") aliLfold "
";

%feature("docstring") aliLfold_cb "
";

// File: group__part__func__global__deprecated.xml

%feature("docstring") alipf_fold_par "

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf() instead  

Parameters
----------
sequences : const char **
structure : char *
pl : RNA.ep() **
parameters : RNA.exp_param() *
calculate_bppm : int
is_constrained : int
is_circular : int

Returns
-------
float  
";

%feature("docstring") alipf_fold "

The partition function version of alifold() works in analogy to pf_fold(). Pair probabilities and
information about sequence covariations are returned via the 'pi' variable as a list of RNA.pinfo()
structs. The list is terminated by the first entry with pi.i = 0.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf() instead  

Parameters
----------
sequences : const char **
structure : char *
pl : RNA.ep() **

Returns
-------
float  
";

%feature("docstring") alipf_circ_fold "

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf() instead  

Parameters
----------
sequences : const char **
structure : char *
pl : RNA.ep() **

Returns
-------
float  
";

%feature("docstring") export_ali_bppm "

Get a pointer to the base pair probability array.  

Accessing the base pair probabilities for a pair (i,j) is achieved by  

    FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j];  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! The new RNA.fold_compound() allows direct access to the
    folding matrices, including the pair probabilities! The pair probability array returned here
    reflects the one of the latest call to RNA.fold_compound.pf(), or any of the old API calls for
consensus
    structure partition function folding.  

Returns
-------
FLT_OR_DBL *  
    A pointer to the base pair probability array  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), and RNA.fold_compound.pf()  
";

%feature("docstring") free_alipf_arrays "

Free the memory occupied by folding matrices allocated by alipf_fold, alipf_circ_fold, etc.  

.. deprecated:: 2.6.4
    Usage of this function is discouraged! This function only free's memory allocated by old API
    function calls. Memory allocated by any of the new API calls (starting with RNA.) will be not
    affected!  

See Also
--------
RNA.fold_compound(), RNA.RNA.fold_compound_free()  
";

%feature("docstring") alipbacktrack "

Sample a consensus secondary structure from the Boltzmann ensemble according its probability.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pbacktrack() instead!  

Parameters
----------
prob : double *
    to be described (berni)  

Returns
-------
char *  
    A sampled consensus secondary structure in dot-bracket notation  
";

%feature("docstring") get_alipf_arrays "

Get pointers to (almost) all relavant arrays used in alifold's partition function computation.  

.. deprecated:: 2.6.4
    It is discouraged to use this function! The new RNA.fold_compound() allows direct access to all
    necessary consensus structure prediction related variables!  

See Also
--------
RNA.fold_compound(), RNA.fold_compound_comparative(), RNA.fold_compound.pf(), pf_alifold(),
alipf_circ_fold()  

Note
----
To obtain meaningful pointers, call alipf_fold first!  

Parameters
----------
S_p : short ***
    A pointer to the 'S' array (integer representation of nucleotides)  
S5_p : short ***
    A pointer to the 'S5' array  
S3_p : short ***
    A pointer to the 'S3' array  
a2s_p : unsigned short ***
    A pointer to the alignment-column to sequence position mapping array  
Ss_p : char ***
    A pointer to the 'Ss' array  
qb_p : FLT_OR_DBL **
    A pointer to the QB matrix  
qm_p : FLT_OR_DBL **
    A pointer to the QM matrix  
q1k_p : FLT_OR_DBL **
    A pointer to the 5' slice of the Q matrix ( :math:`q1k(k) = Q(1, k)`)  
qln_p : FLT_OR_DBL **
    A pointer to the 3' slice of the Q matrix ( :math:`qln(l) = Q(l, n)`)  
pscore : short **
    A pointer to the start of a pscore list  

Returns
-------
int  
    Non Zero if everything went fine, 0 otherwise  
";

%feature("docstring") pf_fold_par "

Compute the partition function :math:`Q` for a given RNA sequence.  

If *structure* is not a NULL pointer on input, it contains on return a string consisting of the
letters \" . , | { } ( ) \" denoting bases that are essentially unpaired, weakly paired, strongly
paired without preference, weakly upstream (downstream) paired, or strongly up- (down-)stream paired
bases, respectively. If fold_constrained is not 0, the *structure* string is interpreted on input as
a list of constraints for the folding. The character \"x\" marks bases that must be unpaired,
matching brackets \" ( ) \" denote base pairs, all other characters are ignored. Any pairs
conflicting with the constraint will be forbidden. This is usually sufficient to ensure the
constraints are honored. If the parameter calculate_bppm is set to 0 base pairing probabilities will
not be computed (saving CPU time), otherwise after calculations took place pr will contain the
probability that bases *i* and *j* pair.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf() instead  

**Postcondition**
    After successful run the hidden folding matrices are filled with the appropriate Boltzmann
    factors. Depending on whether the global variable do_backtrack was set the base pair
    probabilities are already computed and may be accessed for further usage via the export_bppm()
    function. A call of free_pf_arrays() will free all memory allocated by this function. Successive
    calls will first free previously allocated memory before starting the computation.  

Parameters
----------
sequence : const char *
    The RNA sequence input  
structure : char *
    A pointer to a char array where a base pair probability information can be stored in a pseudo-
    dot-bracket notation (may be NULL, too)  
parameters : RNA.exp_param() *
    Data structure containing the precalculated Boltzmann factors  
calculate_bppm : int
    Switch to Base pair probability calculations on/off (0==off)  
is_constrained : int
    Switch to indicate that a structure contraint is passed via the structure argument (0==off)  
is_circular : int
    Switch to (de-)activate postprocessing steps in case RNA sequence is circular (0==off)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
RNA.fold_compound.pf(), bppm_to_structure(), export_bppm(), RNA.exp_params(), free_pf_arrays()  

Note
----
The global array pr is deprecated and the user who wants the calculated base pair probabilities for
further computations is advised to use the function export_bppm()  
";

%feature("docstring") pf_fold "

Compute the partition function :math:`Q` of an RNA sequence.  

If *structure* is not a NULL pointer on input, it contains on return a string consisting of the
letters \" . , | { } ( ) \" denoting bases that are essentially unpaired, weakly paired, strongly
paired without preference, weakly upstream (downstream) paired, or strongly up- (down-)stream paired
bases, respectively. If fold_constrained is not 0, the *structure* string is interpreted on input as
a list of constraints for the folding. The character \"x\" marks bases that must be unpaired,
matching brackets \" ( ) \" denote base pairs, all other characters are ignored. Any pairs
conflicting with the constraint will be forbidden. This is usually sufficient to ensure the
constraints are honored. If do_backtrack has been set to 0 base pairing probabilities will not be
computed (saving CPU time), otherwise pr will contain the probability that bases *i* and *j* pair.  

**Precondition**
    This function takes its model details from the global variables provided in *RNAlib*  

**Postcondition**
    After successful run the hidden folding matrices are filled with the appropriate Boltzmann
    factors. Depending on whether the global variable do_backtrack was set the base pair
    probabilities are already computed and may be accessed for further usage via the export_bppm()
    function. A call of free_pf_arrays() will free all memory allocated by this function. Successive
    calls will first free previously allocated memory before starting the computation.  

Parameters
----------
sequence : const char *
    The RNA sequence input  
structure : char *
    A pointer to a char array where a base pair probability information can be stored in a pseudo-
    dot-bracket notation (may be NULL, too)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  

See Also
--------
pf_fold_par(), pf_circ_fold(), bppm_to_structure(), export_bppm()  

Note
----
The global array pr is deprecated and the user who wants the calculated base pair probabilities for
further computations is advised to use the function export_bppm().  **OpenMP:** This function is not
entirely threadsafe. While the recursions are working on their own copies of data the model details
for the recursions are determined from the global settings just before entering the recursions.
Consider using pf_fold_par() for a really threadsafe implementation.  
";

%feature("docstring") pf_circ_fold "

Compute the partition function of a circular RNA sequence.  

**Precondition**
    This function takes its model details from the global variables provided in *RNAlib*  

**Postcondition**
    After successful run the hidden folding matrices are filled with the appropriate Boltzmann
    factors. Depending on whether the global variable do_backtrack was set the base pair
    probabilities are already computed and may be accessed for further usage via the export_bppm()
    function. A call of free_pf_arrays() will free all memory allocated by this function. Successive
    calls will first free previously allocated memory before starting the computation.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf() instead!  

See Also
--------
RNA.fold_compound.pf()  

Note
----
The global array pr is deprecated and the user who wants the calculated base pair probabilities for
further computations is advised to use the function export_bppm().  **OpenMP:** This function is not
entirely threadsafe. While the recursions are working on their own copies of data the model details
for the recursions are determined from the global settings just before entering the recursions.
Consider using pf_fold_par() for a really threadsafe implementation.  

Parameters
----------
sequence : const char *
    The RNA sequence input  
structure : char *
    A pointer to a char array where a base pair probability information can be stored in a pseudo-
    dot-bracket notation (may be NULL, too)  

Returns
-------
float  
    The ensemble free energy :math:`G = -RT \\cdot \\log(Q)` in kcal/mol  
";

%feature("docstring") free_pf_arrays "

Free arrays for the partition function recursions.  

Call this function if you want to free all allocated memory associated with the partition function
forward recursion.  

.. deprecated:: 2.6.4
    See RNA.fold_compound() and its related functions for how to free memory occupied by the
    dynamic programming matrices  

Note
----
Successive calls of pf_fold(), pf_circ_fold() already check if they should free any memory from a
previous run.  **OpenMP notice:**  
 This function should be called before leaving a thread in order to avoid leaking memory  

**Postcondition**
    All memory allocated by pf_fold_par(), pf_fold() or pf_circ_fold() will be free'd  

See Also
--------
pf_fold_par(), pf_fold(), pf_circ_fold()  
";

%feature("docstring") update_pf_params "

Recalculate energy parameters.  

Call this function to recalculate the pair matrix and energy parameters after a change in folding
parameters like temperature  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.exp_params_subst() instead  
";

%feature("docstring") update_pf_params_par "

Recalculate energy parameters.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.exp_params_subst() instead  
";

%feature("docstring") export_bppm "

Get a pointer to the base pair probability array.  

Accessing the base pair probabilities for a pair (i,j) is achieved by  

**Precondition**
    Call pf_fold_par(), pf_fold() or pf_circ_fold() first to fill the base pair probability array  

Returns
-------
FLT_OR_DBL *  
    A pointer to the base pair probability array  

See Also
--------
pf_fold(), pf_circ_fold(), RNA.idx_row_wise()  
";

%feature("docstring") get_pf_arrays "

Get the pointers to (almost) all relavant computation arrays used in partition function computation.  

**Precondition**
    In order to assign meaningful pointers, you have to call pf_fold_par() or pf_fold() first!  

Parameters
----------
S_p : short **
    A pointer to the 'S' array (integer representation of nucleotides)  
S1_p : short **
    A pointer to the 'S1' array (2nd integer representation of nucleotides)  
ptype_p : char **
    A pointer to the pair type matrix  
qb_p : FLT_OR_DBL **
    A pointer to the QB matrix  
qm_p : FLT_OR_DBL **
    A pointer to the QM matrix  
q1k_p : FLT_OR_DBL **
    A pointer to the 5' slice of the Q matrix ( :math:`q1k(k) = Q(1, k)`)  
qln_p : FLT_OR_DBL **
    A pointer to the 3' slice of the Q matrix ( :math:`qln(l) = Q(l, n)`)  

Returns
-------
int  
    Non Zero if everything went fine, 0 otherwise  

See Also
--------
pf_fold_par(), pf_fold(), pf_circ_fold()  
";

%feature("docstring") get_subseq_F "

Get the free energy of a subsequence from the q[] array.  
";

%feature("docstring") mean_bp_distance "

Get the mean base pair distance of the last partition function computation.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.mean_bp_distance() or RNA.mean_bp_distance_pr() instead!  

Parameters
----------
length : int

Returns
-------
double  
    mean base pair distance in thermodynamic ensemble  

See Also
--------
RNA.fold_compound.mean_bp_distance(), RNA.mean_bp_distance_pr()  
";

%feature("docstring") mean_bp_distance_pr "

Get the mean base pair distance in the thermodynamic ensemble.  

This is a threadsafe implementation of mean_bp_dist() !  

:math:`<d> = \\sum_{a,b} p_{a} p_{b} d(S_{a},S_{b})`  
this can be computed from the pair probs :math:`p_{i}j` as  :math:`<d> = \\sum_{ij}
p_{ij}(1-p_{ij})`  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.mean_bp_distance() or RNA.mean_bp_distance_pr() instead!  

Parameters
----------
length : int
    The length of the sequence  
pr : FLT_OR_DBL *
    The matrix containing the base pair probabilities  

Returns
-------
double  
    The mean pair distance of the structure ensemble  
";

%feature("docstring") stackProb "

Get the probability of stacks.  

.. deprecated:: 2.6.4
    Use RNA.stack_prob() instead!  
";

%feature("docstring") init_pf_fold "

Allocate space for pf_fold()  

.. deprecated:: 2.6.4
    This function is obsolete and will be removed soon!  
";

%feature("docstring") co_pf_fold "

Calculate partition function and base pair probabilities.  

This is the cofold partition function folding. The second molecule starts at the cut_point
nucleotide.  

.. deprecated:: 2.6.4

Note
----
OpenMP: Since this function relies on the global parameters do_backtrack, dangles, temperature and
pf_scale it is not threadsafe according to concurrent changes in these variables! Use
co_pf_fold_par() instead to circumvent this issue.  

Parameters
----------
sequence : char *
    Concatenated RNA sequences  
structure : char *
    Will hold the structure or constraints  

Returns
-------
RNA.dimer_pf()  
    RNA.dimer_pf() structure containing a set of energies needed for concentration computations.  
";

%feature("docstring") co_pf_fold_par "

Calculate partition function and base pair probabilities.  

This is the cofold partition function folding. The second molecule starts at the cut_point
nucleotide.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pf_dimer() instead!  

Parameters
----------
sequence : char *
    Concatenated RNA sequences  
structure : char *
    Pointer to the structure constraint  
parameters : RNA.exp_param() *
    Data structure containing the precalculated Boltzmann factors  
calculate_bppm : int
    Switch to turn Base pair probability calculations on/off (0==off)  
is_constrained : int
    Switch to indicate that a structure contraint is passed via the structure argument (0==off)  

Returns
-------
RNA.dimer_pf()  
    RNA.dimer_pf() structure containing a set of energies needed for concentration computations.  

See Also
--------
get_boltzmann_factors(), co_pf_fold()  
";

%feature("docstring") compute_probabilities "

Compute Boltzmann probabilities of dimerization without homodimers.  

Given the pair probabilities and free energies (in the null model) for a dimer AB and the two
constituent monomers A and B, compute the conditional pair probabilities given that a dimer AB
actually forms. Null model pair probabilities are given as a list as produced by
assign_plist_from_pr(), the dimer probabilities 'prAB' are modified in place.  

.. deprecated:: 2.6.4  

Parameters
----------
FAB : double
    free energy of dimer AB  
FEA : double
    free energy of monomer A  
FEB : double
    free energy of monomer B  
prAB : RNA.ep() *
    pair probabilities for dimer  
prA : RNA.ep() *
    pair probabilities monomer  
prB : RNA.ep() *
    pair probabilities monomer  
Alength : int
    Length of molecule A  
";

%feature("docstring") init_co_pf_fold "

DO NOT USE THIS FUNCTION ANYMORE  

.. deprecated:: 2.6.4  
";

%feature("docstring") export_co_bppm "

Get a pointer to the base pair probability array.  

Accessing the base pair probabilities for a pair (i,j) is achieved by  

    FLT_OR_DBL *pr = export_bppm(); pr_ij = pr[iindx[i]-j];  

.. deprecated:: 2.6.4
    This function is deprecated and will be removed soon! The base pair probability array is
    available through the RNA.fold_compound() data structure, and its associated RNA.mx_pf()
    member.  

Returns
-------
FLT_OR_DBL *  
    A pointer to the base pair probability array  

See Also
--------
RNA.idx_row_wise()  
";

%feature("docstring") free_co_pf_arrays "

Free the memory occupied by co_pf_fold()  

.. deprecated:: 2.6.4
    This function will be removed for the new API soon! See RNA.fold_compound.pf_dimer(),
RNA.fold_compound(),
    and RNA.fold_compound_free() for an alternative  
";

%feature("docstring") update_co_pf_params "

Recalculate energy parameters.  

This function recalculates all energy parameters given the current model settings.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.exp_params_subst() instead!  

Parameters
----------
length : int
    Length of the current RNA sequence  
";

%feature("docstring") update_co_pf_params_par "

Recalculate energy parameters.  

This function recalculates all energy parameters given the current model settings. It's second
argument can either be NULL or a data structure containing the precomputed Boltzmann factors. In the
first scenario, the necessary data structure will be created automatically according to the current
global model settings, i.e. this mode might not be threadsafe. However, if the provided data
structure is not NULL, threadsafety for the model parameters dangles, pf_scale and temperature is
regained, since their values are taken from this data structure during subsequent calculations.  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.exp_params_subst() instead!  

Parameters
----------
length : int
    Length of the current RNA sequence  
parameters : RNA.exp_param() *
    data structure containing the precomputed Boltzmann factors  
";

%feature("docstring") assign_plist_from_db "

Create a RNA.ep() from a dot-bracket string.  

The dot-bracket string is parsed and for each base pair an entry in the plist is created. The
probability of each pair in the list is set by a function parameter.  

The end of the plist is marked by sequence positions i as well as j equal to 0. This condition
should be used to stop looping over its entries  

.. deprecated:: 2.6.4
    Use RNA.plist() instead  

Parameters
----------
pl : RNA.ep() **
    A pointer to the RNA.ep() that is to be created  
struc : const char *
    The secondary structure in dot-bracket notation  
pr : float
    The probability for each base pair  
";

%feature("docstring") assign_plist_from_pr "

Create a RNA.ep() from a probability matrix.  

The probability matrix given is parsed and all pair probabilities above the given threshold are used
to create an entry in the plist  

The end of the plist is marked by sequence positions i as well as j equal to 0. This condition
should be used to stop looping over its entries  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.plist_from_probs() instead!  

Note
----
This function is threadsafe  

Parameters
----------
pl : RNA.ep() **
    A pointer to the RNA.ep() that is to be created  
probs : FLT_OR_DBL *
    The probability matrix used for creating the plist  
length : int
    The length of the RNA sequence  
cutoff : double
    The cutoff value  
";

// File: group__part__func__window__deprecated.xml

%feature("docstring") update_pf_paramsLP "

Parameters
----------
length : int  
";

%feature("docstring") update_pf_paramsLP_par "
";

%feature("docstring") pfl_fold "

Compute partition functions for locally stable secondary structures.  

pfl_fold computes partition functions for every window of size 'winSize' possible in a RNA molecule,
allowing only pairs with a span smaller than 'pairSize'. It returns the mean pair probabilities
averaged over all windows containing the pair in 'pl'. 'winSize' should always be >= 'pairSize'.
Note that in contrast to Lfold(), bases outside of the window do not influence the structure at all.
Only probabilities higher than 'cutoffb' are kept.  

If 'pU' is supplied (i.e is not the NULL pointer), pfl_fold() will also compute the mean probability
that regions of length 'u' and smaller are unpaired. The parameter 'u' is supplied in 'pup[0][0]'.
On return the 'pup' array will contain these probabilities, with the entry on 'pup[x][y]' containing
the mean probability that x and the y-1 preceding bases are unpaired. The 'pU' array needs to be
large enough to hold n+1 float* entries, where n is the sequence length.  

If an array dpp2 is supplied, the probability of base pair (i,j) given that there already exists a
base pair (i+1,j-1) is also computed and saved in this array. If pUfp is given (i.e. not NULL), pU
is not saved but put out imediately. If spup is given (i.e. is not NULL), the pair probabilities in
pl are not saved but put out imediately.  

Parameters
----------
sequence : char *
    RNA sequence  
winSize : int
    size of the window  
pairSize : int
    maximum size of base pair  
cutoffb : float
    cutoffb for base pairs  
pU : double **
    array holding all unpaired probabilities  
dpp2 : RNA.ep() **
    array of dependent pair probabilities  
pUfp : FILE *
    file pointer for pU  
spup : FILE *
    file pointer for pair probabilities  

Returns
-------
RNA.ep() *  
    list of pair probabilities  
";

%feature("docstring") pfl_fold_par "

Compute partition functions for locally stable secondary structures.  
";

%feature("docstring") putoutpU_prob "

Writes the unpaired probabilities (pU) or opening energies into a file.  

Can write either the unpaired probabilities (accessibilities) pU or the opening energies -log(pU)kT
into a file  

Parameters
----------
pU : double **
    pair probabilities  
length : int
    length of RNA sequence  
ulength : int
    maximum length of unpaired stretch  
fp : FILE *
    file pointer of destination file  
energies : int
    switch to put out as opening energies  
";

%feature("docstring") putoutpU_prob_bin "

Writes the unpaired probabilities (pU) or opening energies into a binary file.  

Can write either the unpaired probabilities (accessibilities) pU or the opening energies -log(pU)kT
into a file  

Parameters
----------
pU : double **
    pair probabilities  
length : int
    length of RNA sequence  
ulength : int
    maximum length of unpaired stretch  
fp : FILE *
    file pointer of destination file  
energies : int
    switch to put out as opening energies  
";

// File: group__subopt__stochbt__deprecated.xml

%feature("docstring") pbacktrack "

Sample a secondary structure from the Boltzmann ensemble according its probability.  

**Precondition**
    st_back has to be set to 1 before calling pf_fold() or pf_fold_par()  pf_fold_par() or pf_fold()
    have to be called first to fill the partition function matrices  

Parameters
----------
sequence : char *
    The RNA sequence  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  
";

%feature("docstring") pbacktrack5 "

Sample a sub-structure from the Boltzmann ensemble according its probability.  
";

%feature("docstring") pbacktrack_circ "

Sample a secondary structure of a circular RNA from the Boltzmann ensemble according its
probability.  

This function does the same as pbacktrack() but assumes the RNA molecule to be circular  

**Precondition**
    st_back has to be set to 1 before calling pf_fold() or pf_fold_par()  pf_fold_par() or
    pf_circ_fold() have to be called first to fill the partition function matrices  

.. deprecated:: 2.6.4
    Use RNA.fold_compound.pbacktrack() instead.  

Parameters
----------
sequence : char *
    The RNA sequence  

Returns
-------
char *  
    A sampled secondary structure in dot-bracket notation  
";

// File: group__aln__utils__deprecated.xml

%feature("docstring") read_clustal "
";

%feature("docstring") consensus "
";

%feature("docstring") consens_mis "
";

%feature("docstring") get_ungapped_sequence "
";

%feature("docstring") get_mpi "

Get the mean pairwise identity in steps from ?to?(ident)  

.. deprecated:: 2.6.4
    Use RNA.aln_mpi() as a replacement  

Parameters
----------
Alseq : char *
n_seq : int
    The number of sequences in the alignment  
length : int
    The length of the alignment  
mini : int *

Returns
-------
int  
    The mean pairwise identity  
";

%feature("docstring") encode_ali_sequence "

Get arrays with encoded sequence of the alignment.  

this function assumes that in S, S5, s3, ss and as enough space is already allocated (size must be
at least sequence length+2)  

Parameters
----------
sequence : const char *
    The gapped sequence from the alignment  
S : short *
    pointer to an array that holds encoded sequence  
s5 : short *
    pointer to an array that holds the next base 5' of alignment position i  
s3 : short *
    pointer to an array that holds the next base 3' of alignment position i  
ss : char *
as : unsigned short *
circ : int
    assume the molecules to be circular instead of linear (circ=0)  
";

%feature("docstring") alloc_sequence_arrays "

Allocate memory for sequence array used to deal with aligned sequences.  

Note that these arrays will also be initialized according to the sequence alignment given  

Parameters
----------
sequences : const char **
    The aligned sequences  
S : short ***
    A pointer to the array of encoded sequences  
S5 : short ***
    A pointer to the array that contains the next 5' nucleotide of a sequence position  
S3 : short ***
    A pointer to the array that contains the next 3' nucleotide of a sequence position  
a2s : unsigned short ***
    A pointer to the array that contains the alignment to sequence position mapping  
Ss : char ***
    A pointer to the array that contains the ungapped sequence  
circ : int
    assume the molecules to be circular instead of linear (circ=0)  

See Also
--------
free_sequence_arrays()  
";

%feature("docstring") free_sequence_arrays "

Free the memory of the sequence arrays used to deal with aligned sequences.  

This function frees the memory previously allocated with alloc_sequence_arrays()  

Parameters
----------
n_seq : unsigned int
    The number of aligned sequences  
S : short ***
    A pointer to the array of encoded sequences  
S5 : short ***
    A pointer to the array that contains the next 5' nucleotide of a sequence position  
S3 : short ***
    A pointer to the array that contains the next 3' nucleotide of a sequence position  
a2s : unsigned short ***
    A pointer to the array that contains the alignment to sequence position mapping  
Ss : char ***
    A pointer to the array that contains the ungapped sequence  

See Also
--------
alloc_sequence_arrays()  
";

// File: group__struct__utils__deprecated.xml

%feature("docstring") b2HIT "

Converts the full structure from bracket notation to the HIT notation including root.  

.. deprecated:: 2.6.4
    See RNA.db_to_tree_string() and RNA.STRUCTURE_TREE_HIT for a replacement  

Parameters
----------
structure : const char *

Returns
-------
char *  
";

%feature("docstring") b2C "

Converts the full structure from bracket notation to the a coarse grained notation using the 'H' 'B'
'I' 'M' and 'R' identifiers.  

.. deprecated:: 2.6.4
    See RNA.db_to_tree_string() and RNA.STRUCTURE_TREE_SHAPIRO_SHORT for a replacement  

Parameters
----------
structure : const char *

Returns
-------
char *  
";

%feature("docstring") b2Shapiro "

Converts the full structure from bracket notation to the *weighted* coarse grained notation using
the 'H' 'B' 'I' 'M' 'S' 'E' and 'R' identifiers.  

.. deprecated:: 2.6.4
    See RNA.db_to_tree_string() and RNA.STRUCTURE_TREE_SHAPIRO_WEIGHT for a replacement  

Parameters
----------
structure : const char *

Returns
-------
char *  
";

%feature("docstring") add_root "

Adds a root to an un-rooted tree in any except bracket notation.  

Parameters
----------
structure : const char *

Returns
-------
char *  
";

%feature("docstring") expand_Shapiro "

Inserts missing 'S' identifiers in unweighted coarse grained structures as obtained from b2C().  

Parameters
----------
coarse : const char *

Returns
-------
char *  
";

%feature("docstring") expand_Full "

Convert the full structure from bracket notation to the expanded notation including root.  

Parameters
----------
structure : const char *

Returns
-------
char *  
";

%feature("docstring") unexpand_Full "

Restores the bracket notation from an expanded full or HIT tree, that is any tree using only
identifiers 'U' 'P' and 'R'.  

Parameters
----------
ffull : const char *

Returns
-------
char *  
";

%feature("docstring") unweight "

Strip weights from any weighted tree.  

Parameters
----------
wcoarse : const char *

Returns
-------
char *  
";

%feature("docstring") unexpand_aligned_F "

Converts two aligned structures in expanded notation.  

Takes two aligned structures as produced by tree_edit_distance() function back to bracket notation
with '_' as the gap character. The result overwrites the input.  

Parameters
----------
align : char *  
";

%feature("docstring") parse_structure "

Collects a statistic of structure elements of the full structure in bracket notation.  

The function writes to the following global variables: loop_size, loop_degree, helix_size, loops,
pairs, unpaired  

Parameters
----------
structure : const char *  
";

%feature("docstring") pack_structure "

Pack secondary secondary structure, 5:1 compression using base 3 encoding.  

Returns a binary string encoding of the secondary structure using a 5:1 compression scheme. The
string is NULL terminated and can therefore be used with standard string functions such as strcmp().
Useful for programs that need to keep many structures in memory.  

.. deprecated:: 2.6.4
    Use RNA.db_pack() as a replacement  

Parameters
----------
struc : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
char *  
    The binary encoded structure  
";

%feature("docstring") unpack_structure "

Unpack secondary structure previously packed with pack_structure()  

Translate a compressed binary string produced by pack_structure() back into the familiar dot-bracket
notation.  

.. deprecated:: 2.6.4
    Use RNA.db_unpack() as a replacement  

Parameters
----------
packed : const char *
    The binary encoded packed secondary structure  

Returns
-------
char *  
    The unpacked secondary structure in dot-bracket notation  
";

%feature("docstring") make_pair_table "

Create a pair table of a secondary structure.  

Returns a newly allocated table, such that table[i]=j if (i.j) pair or 0 if i is unpaired, table[0]
contains the length of the structure.  

.. deprecated:: 2.6.4
    Use RNA.ptable() instead  

Parameters
----------
structure : const char *
    The secondary structure in dot-bracket notation  

Returns
-------
short *  
    A pointer to the created pair_table  
";

%feature("docstring") copy_pair_table "

Get an exact copy of a pair table.  

.. deprecated:: 2.6.4
    Use RNA.ptable_copy() instead  

Parameters
----------
pt : const short *
    The pair table to be copied  

Returns
-------
short *  
    A pointer to the copy of 'pt'  
";

%feature("docstring") alimake_pair_table "

Pair table for snoop align  

.. deprecated:: 2.6.4
    Use RNA.pt_ali_get() instead!  
";

%feature("docstring") make_pair_table_snoop "

returns a newly allocated table, such that: table[i]=j if (i.j) pair or 0 if i is unpaired, table[0]
contains the length of the structure. The special pseudoknotted H/ACA-mRNA structure is taken into
account.  

.. deprecated:: 2.6.4
    Use RNA.pt_snoop_get() instead!  
";

%feature("docstring") bp_distance "

Compute the \"base pair\" distance between two secondary structures s1 and s2.  

The sequences should have the same length. dist = number of base pairs in one structure but not in
the other same as edit distance with open-pair close-pair as move-set  

.. deprecated:: 2.6.4
    Use RNA.bp_distance instead  

Parameters
----------
str1 : const char *
    First structure in dot-bracket notation  
str2 : const char *
    Second structure in dot-bracket notation  

Returns
-------
int  
    The base pair distance between str1 and str2  
";

%feature("docstring") make_referenceBP_array "

Make a reference base pair count matrix.  

Get an upper triangular matrix containing the number of basepairs of a reference structure for each
interval [i,j] with i<j. Access it via iindx!!!  

.. deprecated:: 2.6.4
    Use RNA.refBPcnt_matrix() instead  
";

%feature("docstring") compute_BPdifferences "

Make a reference base pair distance matrix.  

Get an upper triangular matrix containing the base pair distance of two reference structures for
each interval [i,j] with i<j. Access it via iindx!!!  

.. deprecated:: 2.6.4
    Use RNA.refBPdist_matrix() instead  
";

%feature("docstring") parenthesis_structure "

Create a dot-backet/parenthesis structure from backtracking stack.  

.. deprecated:: 2.6.4
    use RNA.parenthesis_structure() instead  

Note
----
This function is threadsafe  
";

%feature("docstring") parenthesis_zuker "

Create a dot-backet/parenthesis structure from backtracking stack obtained by zuker suboptimal
calculation in cofold.c.  

.. deprecated:: 2.6.4
    use RNA.parenthesis_zuker instead  

Note
----
This function is threadsafe  
";

%feature("docstring") bppm_to_structure "

Create a dot-bracket like structure string from base pair probability matrix.  

.. deprecated:: 2.6.4
    Use RNA.db_from_probs() instead!  
";

%feature("docstring") bppm_symbol "

Get a pseudo dot bracket notation for a given probability information.  

.. deprecated:: 2.6.4
    Use RNA.bpp_symbol() instead!  
";

%feature("docstring") STRUC "
";

// File: group__plotting__utils__deprecated.xml

%feature("docstring") PS_color_aln "

Produce PostScript sequence alignment color-annotated by consensus structure.  

.. deprecated:: 2.6.4
    Use RNA.file_PS_aln() instead!  
";

%feature("docstring") aliPS_color_aln "

PS_color_aln for duplexes.  

.. deprecated:: 2.6.4
    Use RNA.file_PS_aln() instead!  
";

%feature("docstring") simple_xy_coordinates "

Calculate nucleotide coordinates for secondary structure plot the *Simple way*  

.. deprecated:: 2.6.4
    Consider switching to RNA.plot_coords_simple_pt() instead!  

See Also
--------
make_pair_table(), rna_plot_type, simple_circplot_coordinates(), naview_xy_coordinates(),
RNA.file_PS_rnaplot_a(), RNA.file_PS_rnaplot, svg_rna_plot()  

Parameters
----------
pair_table : short *
    The pair table of the secondary structure  
X : float *
    a pointer to an array with enough allocated space to hold the x coordinates  
Y : float *
    a pointer to an array with enough allocated space to hold the y coordinates  

Returns
-------
int  
    length of sequence on success, 0 otherwise  
";

%feature("docstring") simple_circplot_coordinates "

Calculate nucleotide coordinates for *Circular Plot*  

This function calculates the coordinates of nucleotides mapped in equal distancies onto a unit
circle.  

.. deprecated:: 2.6.4
    Consider switching to RNA.plot_coords_circular_pt() instead!  

See Also
--------
make_pair_table(), rna_plot_type, simple_xy_coordinates(), naview_xy_coordinates(),
RNA.file_PS_rnaplot_a(), RNA.file_PS_rnaplot, svg_rna_plot()  

Note
----
In order to draw nice arcs using quadratic bezier curves that connect base pairs one may calculate a
second tangential point :math:`P^t` in addition to the actual R2 coordinates. the simplest way to do
so may be to compute a radius scaling factor :math:`rs` in the interval :math:`[0,1]` that weights
the proportion of base pair span to the actual length of the sequence. This scaling factor can then
be used to calculate the coordinates for :math:`P^t`, i.e. :math:`P^{t}_{x}[i] = X[i] * rs` and
:math:`P^{t}_{y}[i] = Y[i] * rs`.  

Parameters
----------
pair_table : short *
    The pair table of the secondary structure  
x : float *
    a pointer to an array with enough allocated space to hold the x coordinates  
y : float *
    a pointer to an array with enough allocated space to hold the y coordinates  

Returns
-------
int  
    length of sequence on success, 0 otherwise  
";

%feature("docstring") PS_color_dot_plot "
";

%feature("docstring") PS_color_dot_plot_turn "
";

%feature("docstring") PS_dot_plot_turn "
";

%feature("docstring") PS_dot_plot_list "

Produce a postscript dot-plot from two pair lists.  

This function reads two plist structures (e.g. base pair probabilities and a secondary structure) as
produced by assign_plist_from_pr() and assign_plist_from_db() and produces a postscript \"dot plot\"
that is written to 'filename'.  
Using base pair probabilities in the first and mfe structure in the second plist, the resulting
\"dot plot\" represents each base pairing probability by a square of corresponding area in a upper
triangle matrix. The lower part of the matrix contains the minimum free energy structure.  

Parameters
----------
seq : char *
    The RNA sequence  
filename : char *
    A filename for the postscript output  
pl : RNA.ep() *
    The base pair probability pairlist  
mf : RNA.ep() *
    The mfe secondary structure pairlist  
comment : char *
    A comment  

Returns
-------
int  
    1 if postscript was successfully written, 0 otherwise  

See Also
--------
assign_plist_from_pr(), assign_plist_from_db()  
";

%feature("docstring") PS_dot_plot "

Produce postscript dot-plot.  

Wrapper to PS_dot_plot_list  

Reads base pair probabilities produced by pf_fold() from the global array pr and the pair list
base_pair produced by fold() and produces a postscript \"dot plot\" that is written to 'filename'.
The \"dot plot\" represents each base pairing probability by a square of corresponding area in a
upper triangle matrix. The lower part of the matrix contains the minimum free energy  

.. deprecated:: 2.6.4
    This function is deprecated and will be removed soon! Use PS_dot_plot_list() instead!  

Note
----
DO NOT USE THIS FUNCTION ANYMORE SINCE IT IS NOT THREADSAFE  
";

// File: group__paths__deprecated.xml

%feature("docstring") find_saddle "

Find energy of a saddle point between 2 structures (search only direct path)  

.. deprecated:: 2.6.4
    Use RNA.path_findpath_saddle() instead!  

Parameters
----------
seq : const char *
    RNA sequence  
s1 : const char *
    A pointer to the character array where the first secondary structure in dot-bracket notation
    will be written to  
s2 : const char *
    A pointer to the character array where the second secondary structure in dot-bracket notation
    will be written to  
width : int
    integer how many strutures are being kept during the search  

Returns
-------
int  
    the saddle energy in 10cal/mol  
";

%feature("docstring") free_path "

Free memory allocated by get_path() function.  

.. deprecated:: 2.6.4
    Use RNA.path_free() instead!  

Parameters
----------
path : RNA.path() *
    pointer to memory to be freed  
";

%feature("docstring") get_path "

Find refolding path between 2 structures (search only direct path)  

.. deprecated:: 2.6.4
    Use RNA.path_findpath() instead!  

Parameters
----------
seq : const char *
    RNA sequence  
s1 : const char *
    A pointer to the character array where the first secondary structure in dot-bracket notation
    will be written to  
s2 : const char *
    A pointer to the character array where the second secondary structure in dot-bracket notation
    will be written to  
width : int
    integer how many strutures are being kept during the search  

Returns
-------
RNA.path() *  
    direct refolding path between two structures  
";

// File: distance_measures.xml

// File: plots.xml

// File: deprecated.xml

// File: callbacks.xml

// File: bug.xml

// File: wrappers.xml

// File: dir_e78bc76dd94a8744fef38d57611bd23c.xml

// File: dir_8a4ca9b4a73205797ab566a10a6011af.xml

// File: dir_e68e8157741866f444e17edd764ebbae.xml

// File: dir_04f2ecc425faf0d475a3caf484e551f3.xml

// File: dir_ff6b5900125bb0123025c1cb24bdc726.xml

// File: dir_45ac42044d8b814ecf694abd17d2de12.xml

// File: dir_3d74ce004720553e7871f65d12e07ff9.xml

// File: dir_54da8e06181df970e10d68c1aa7bbb33.xml

// File: dir_8b3d6293463b6781f7bf8859fc72e182.xml

// File: dir_1ca6f3a46323f79e1934000abf5d4efb.xml

// File: dir_1d2f085f11683a3b831dd25384489489.xml

// File: dir_511ee2b14d6879b7ada858ff02a7ab24.xml

// File: dir_ece6de16e10c3ab1bafd819fa99af76e.xml

// File: dir_1ddf03e95848b6488d60454f7a385ded.xml

// File: dir_94a42affd43e91f8a4e9ff1a4c7599c5.xml

// File: indexpage.xml


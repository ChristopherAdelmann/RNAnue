/*#######################################*/
/* Interface for static strings          */
/*#######################################*/

#include <string>

#ifdef SWIGPYTHON


%typemap(varout) const unsigned char parameter_set_dna_mathews1999[] {
  std::string str( parameter_set_dna_mathews1999, parameter_set_dna_mathews1999 + sizeof (parameter_set_dna_mathews1999) / sizeof (parameter_set_dna_mathews1999[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_dna_mathews2004[] {
  std::string str( parameter_set_dna_mathews2004, parameter_set_dna_mathews2004 + sizeof (parameter_set_dna_mathews2004) / sizeof (parameter_set_dna_mathews2004[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_andronescu2007[] {
  std::string str( parameter_set_rna_andronescu2007, parameter_set_rna_andronescu2007 + sizeof (parameter_set_rna_andronescu2007) / sizeof (parameter_set_rna_andronescu2007[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_langdon2018[] {
  std::string str( parameter_set_rna_langdon2018, parameter_set_rna_langdon2018 + sizeof (parameter_set_rna_langdon2018) / sizeof (parameter_set_rna_langdon2018[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_misc_special_hairpins[] {
  std::string str( parameter_set_rna_misc_special_hairpins, parameter_set_rna_misc_special_hairpins + sizeof (parameter_set_rna_misc_special_hairpins) / sizeof (parameter_set_rna_misc_special_hairpins[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_turner1999[] {
  std::string str( parameter_set_rna_turner1999, parameter_set_rna_turner1999 + sizeof (parameter_set_rna_turner1999) / sizeof (parameter_set_rna_turner1999[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_turner2004[] {
  std::string str( parameter_set_rna_turner2004, parameter_set_rna_turner2004 + sizeof (parameter_set_rna_turner2004) / sizeof (parameter_set_rna_turner2004[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_7DA_parameters[] {
  std::string str( parameter_set_rna_mod_7DA_parameters, parameter_set_rna_mod_7DA_parameters + sizeof (parameter_set_rna_mod_7DA_parameters) / sizeof (parameter_set_rna_mod_7DA_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_inosine_parameters[] {
  std::string str( parameter_set_rna_mod_inosine_parameters, parameter_set_rna_mod_inosine_parameters + sizeof (parameter_set_rna_mod_inosine_parameters) / sizeof (parameter_set_rna_mod_inosine_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_m6A_parameters[] {
  std::string str( parameter_set_rna_mod_m6A_parameters, parameter_set_rna_mod_m6A_parameters + sizeof (parameter_set_rna_mod_m6A_parameters) / sizeof (parameter_set_rna_mod_m6A_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_pseudouridine_parameters[] {
  std::string str( parameter_set_rna_mod_pseudouridine_parameters, parameter_set_rna_mod_pseudouridine_parameters + sizeof (parameter_set_rna_mod_pseudouridine_parameters) / sizeof (parameter_set_rna_mod_pseudouridine_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_purine_parameters[] {
  std::string str( parameter_set_rna_mod_purine_parameters, parameter_set_rna_mod_purine_parameters + sizeof (parameter_set_rna_mod_purine_parameters) / sizeof (parameter_set_rna_mod_purine_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_dihydrouridine_parameters[] {
  std::string str( parameter_set_rna_mod_dihydrouridine_parameters, parameter_set_rna_mod_dihydrouridine_parameters + sizeof (parameter_set_rna_mod_dihydrouridine_parameters) / sizeof (parameter_set_rna_mod_dihydrouridine_parameters[0]) );
  $result = PyUnicode_FromString(str.c_str());
}

#endif

#ifdef SWIGPERL5


%typemap(varout) const unsigned char parameter_set_dna_mathews1999[] {
  std::string str( parameter_set_dna_mathews1999, parameter_set_dna_mathews1999 + sizeof (parameter_set_dna_mathews1999) / sizeof (parameter_set_dna_mathews1999[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_dna_mathews2004[] {
  std::string str( parameter_set_dna_mathews2004, parameter_set_dna_mathews2004 + sizeof (parameter_set_dna_mathews2004) / sizeof (parameter_set_dna_mathews2004[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_andronescu2007[] {
  std::string str( parameter_set_rna_andronescu2007, parameter_set_rna_andronescu2007 + sizeof (parameter_set_rna_andronescu2007) / sizeof (parameter_set_rna_andronescu2007[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_langdon2018[] {
  std::string str( parameter_set_rna_langdon2018, parameter_set_rna_langdon2018 + sizeof (parameter_set_rna_langdon2018) / sizeof (parameter_set_rna_langdon2018[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_misc_special_hairpins[] {
  std::string str( parameter_set_rna_misc_special_hairpins, parameter_set_rna_misc_special_hairpins + sizeof (parameter_set_rna_misc_special_hairpins) / sizeof (parameter_set_rna_misc_special_hairpins[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_turner1999[] {
  std::string str( parameter_set_rna_turner1999, parameter_set_rna_turner1999 + sizeof (parameter_set_rna_turner1999) / sizeof (parameter_set_rna_turner1999[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_turner2004[] {
  std::string str( parameter_set_rna_turner2004, parameter_set_rna_turner2004 + sizeof (parameter_set_rna_turner2004) / sizeof (parameter_set_rna_turner2004[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_7DA_parameters[] {
  std::string str( parameter_set_rna_mod_7DA_parameters, parameter_set_rna_mod_7DA_parameters + sizeof (parameter_set_rna_mod_7DA_parameters) / sizeof (parameter_set_rna_mod_7DA_parameters[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_inosine_parameters[] {
  std::string str( parameter_set_rna_mod_inosine_parameters, parameter_set_rna_mod_inosine_parameters + sizeof (parameter_set_rna_mod_inosine_parameters) / sizeof (parameter_set_rna_mod_inosine_parameters[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_m6A_parameters[] {
  std::string str( parameter_set_rna_mod_m6A_parameters, parameter_set_rna_mod_m6A_parameters + sizeof (parameter_set_rna_mod_m6A_parameters) / sizeof (parameter_set_rna_mod_m6A_parameters[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_pseudouridine_parameters[] {
  std::string str( parameter_set_rna_mod_pseudouridine_parameters, parameter_set_rna_mod_pseudouridine_parameters + sizeof (parameter_set_rna_mod_pseudouridine_parameters) / sizeof (parameter_set_rna_mod_pseudouridine_parameters[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_purine_parameters[] {
  std::string str( parameter_set_rna_mod_purine_parameters, parameter_set_rna_mod_purine_parameters + sizeof (parameter_set_rna_mod_purine_parameters) / sizeof (parameter_set_rna_mod_purine_parameters[0]) );
  sv_setpv($result, str.c_str());
}

%typemap(varout) const unsigned char parameter_set_rna_mod_dihydrouridine_parameters[] {
  std::string str( parameter_set_rna_mod_dihydrouridine_parameters, parameter_set_rna_mod_dihydrouridine_parameters + sizeof (parameter_set_rna_mod_dihydrouridine_parameters) / sizeof (parameter_set_rna_mod_dihydrouridine_parameters[0]) );
  sv_setpv($result, str.c_str());
}

#endif

%include  <ViennaRNA/static/energy_parameter_sets.h>

### MAKEFILE:  TO MANAGE SATLAS_ODF AND ITS ASSOCIATED PROGRAMS
### 2017 JUN 15: REMOVED -assume byterecl FROM CLFAGS AND LFLAGS
### 2017 MAR 10: CHANGED EXTENSIONS FROM .f TO .for AND .x TO .exe

# MACRO DEFINITIONS

# THE INTEL FORTRAN COMPILER
FORT = ifort

# THESE COMPILER OPTIONS ARE FOR INTEL FORTRAN COMPILER RUNNING ON A
# DELL WORKSTATION

CFLAGS = -c -check bounds -check uninit -warn declarations -warn unused

# SOME POSSIBLE COMPILER OPTIONS:
# -assume byterecl = FORCE UNFORMATTED I/O TO USE byte AS UNIT
# -c = COMPILE ONLY, DO NOT LINK
# -check = -check all CHECKS EVERYTHING
# -check bounds = RUN-TIME CHECKS FOR SUBSCRIPT/STRING OUT OF BOUNDS
# -check uninit = CHECK FOR UNINITIALIZED VARIABLES
# -fast = MAXIMIZE SPEED
# -parallel = AUTO-PARALLELIZER GENERATES MULTITHREAD CODE FOR LOOPS
# -stand  CURRENTLY THE COMPILER DEFAULTS TO 2008 STANDARD IF NO STANDARD
#    IS EXPLICITLY GIVEN
# -warn all
# -warn declarations = WARN OF VARIABLES THAT ARE NOT DECLARED
# -warn unused = WARN OF DECLARED VARIABLES THAT ARE NEVER USED

# LINKER FLAGS
LFLAGS = -check bounds -check uninit

CONSTANTS = module_astro_parameters.o    \
            module_physical_constants.o

DIMS = module_code_dimensions.o

SYN_DIMS = module_synthe_dimensions.o

TYPES = module_var_types.o

SA_MODULES = module_abross_vars.o           \
             module_abtot_vars.o            \
             module_abundances.o            \
             module_atmosphere_parameters.o \
             module_conv_vars.o             \
             module_depart_vars.o           \
             module_edensity_vars.o         \
             module_elements_vars.o         \
             module_flux_vars.o             \
             module_freq_set.o              \
             module_freq_vars.o             \
             module_gravity.o               \
             module_if_vars.o               \
             module_intensity_vars.o        \
             module_isotope_vars.o          \
             module_iter_vars.o             \
             module_junk_vars.o             \
             module_molecular_vars.o        \
             module_odeint_vars.o           \
             module_odf_vars.o              \
             module_opacity.o               \
             module_opacity_switches.o      \
             module_potion_vars.o           \
             module_put_vars.o              \
             module_pzero_vars.o            \
             module_radius_vars.o           \
             module_rad_pressure.o          \
             module_rad_vars.o              \
             module_rhodr_var.o             \
             module_ross_tables.o           \
             module_state_vars.o            \
             module_tabex.o                 \
             module_tau_std.o               \
             module_temp_vars.o             \
             module_total_opacity.o         \
             module_total_pressure.o        \
             module_tsmooth.o               \
             module_turbpr_vars.o           \
             module_wave_vars.o             \
             module_xnf_vars.o              \
             module_xnmol.o
 
SA_OBJECTS = satlas_odf_main.o  \
             satlas_odf_conv.o  \
             satlas_odf_josh.o  \
             satlas_odf_kapp.o  \
             satlas_odf_pops.o  \
             satlas_odf_read.o  \
             satlas_odf_tcorr.o \
             satlas_odf_utils.o

SPECSYN_MODULES = module_continuum_edges.o  \
                  module_synth_lindat.o     \
                  module_synthe_nlines.o

SPECSYN_OBJECTS = specsyn_main.o    \
                  satlas_odf_josh.o \
                  satlas_odf_read.o \
                  satlas_odf_utils.o

SYNBEG_MODULES = module_list_vars.o     \
                 module_logtab.o        \
                 module_potion_vars.o   \
                 module_synbeg_syndat.o \
                 module_synth_lindat.o

SYNTH_MODULES = module_synth_ionpots.o   \
                module_synth_xnfh_vars.o \
                module_synth_xnfp_vars.o \
                module_synth_xnmol_vars.o

SYNTHE_MODULES = module_continuum_edges.o       \
                 module_synth_lindat.o          \
                 module_synth_xnfh_vars.o       \
                 module_synth_xnmol_vars.o      \
                 module_synthe_buffer.o         \
                 module_synthe_nlines.o         \
                 module_synthe_xnfdop.o

SYNTHE_OBJECTS = synthe_main.o     \
                 satlas_odf_utils.o

X_OBJECTS = xnfpelsyn_main.o  \
            synth_kapp.o      \
            synth_pops.o      \
            satlas_odf_read.o \
            satlas_odf_utils.o

# target definitions

# ------------------------------ BROADEN -------------------------------

broaden : $(TYPES) $(CONSTANTS) $(DIMS) module_intensity_vars.o \
         broaden.o
	$(FORT) $(LFLAGS) -o broaden.exe $(TYPES) $(CONSTANTS) \
                                         $(DIMS) \
                                         module_intensity_vars.o broaden.o

# ----------------------------- LIST_LINES -----------------------------

list_lines : list_lines.o
	$(FORT) $(LFLAGS) -o list_lines.exe list_lines.o

list_lines.o : list_lines.for
	$(FORT) $(CFLAGS) list_lines.for

# ------------------------------ PREDICT -------------------------------

predict : $(TYPES) predict_starting_model.o
	$(FORT) $(LFLAGS) -o predict_starting_model.exe $(TYPES) \
                             predict_starting_model.o

# ------------------------------ REFORM_F ------------------------------

reform_f : $(TYPES) $(CONSTANTS) $(DIMS) reform_f.o
	$(FORT) $(LFLAGS) -o reform_f.exe $(TYPES) $(CONSTANTS) $(DIMS) \
                             reform_f.o

# ------------------------------ REFORM_I ------------------------------

reform_i : $(TYPES) $(CONSTANTS) $(DIMS) reform_i.o
	$(FORT) $(LFLAGS) -o reform_i.exe $(TYPES) $(CONSTANTS) $(DIMS) \
                             reform_i.o

# ----------------------------- REFORM_ROT -----------------------------

reform_rot : $(TYPES) $(CONSTANTS) $(DIMS) reform_rot.o
	$(FORT) $(LFLAGS) -o reform_rot.exe $(TYPES) $(CONSTANTS) $(DIMS) \
                             reform_rot.o

# ------------------------------- ROTATE -------------------------------

rotate : $(TYPES) $(DIMS) $(SYN_DIMS) module_intensity_vars.o \
                                      module_synth_lindat.o rotate.o
	$(FORT) $(LFLAGS) -o rotate.exe $(TYPES) $(DIMS) $(SYN_DIMS) \
                                        module_intensity_vars.o \
                                        module_synth_lindat.o rotate.o

# ----------------------------- SATLAS_ODF -----------------------------

satlas_odf: $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES) \
            $(SA_OBJECTS)
	$(FORT) $(LFLAGS) -o satlas_odf.exe  $(TYPES) $(CONSTANTS) \
                                             $(DIMS) $(SA_MODULES) \
                                           $(SA_OBJECTS)

# ------------------------------- SPECSYN ------------------------------

specsyn : $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES) \
          $(SPECSYN_MODULES) $(SPECSYN_OBJECTS)
	$(FORT) $(LFLAGS) -o specsyn.exe $(TYPES) $(CONSTANTS) \
                                         $(DIMS) $(SA_MODULES) \
                                         $(SPECSYN_MODULES) \
                                         $(SPECSYN_OBJECTS)

# ------------------------------- SYNBEG -------------------------------

synbeg : $(TYPES) $(CONSTANTS) $(SYNBEG_MODULES) synbeg.o
	$(FORT) $(LFLAGS) -o synbeg.exe  $(TYPES) $(CONSTANTS) \
                                         $(SYNBEG_MODULES) synbeg.o

synbeg.o : $(TYPES) $(CONSTANTS) $(SYNBEG_MODULES) synbeg.for
	$(FORT) $(CFLAGS) synbeg.for

# ------------------------------- SYNTHE -------------------------------

synthe : $(TYPES) $(CONSTANTS) $(DIMS) $(SYN_DIMS) \
         $(SA_MODULES) $(SYNTHE_MODULES) $(SYNTHE_OBJECTS)
	$(FORT) $(LFLAGS) -o synthe.exe $(TYPES) $(CONSTANTS) \
                                        $(DIMS) $(SYN_DIMS) \
                                        $(SA_MODULES) $(SYNTHE_MODULES) \
                                        $(SYNTHE_OBJECTS)

# ------------------------------- VFIELD ------------------------------

vfield : $(TYPES) $(CONSTANTS) $(SYN_DIMS) module_vfield_wts.o \
         vfield.o
	$(FORT) $(LFLAGS) -o vfield.exe \
                $(TYPES) $(CONSTANTS) $(SYN_DIMS) \
                module_vfield_wts.o vfield.o

vfield.o : $(TYPES) $(CONSTANTS) $(SYN_DIMS) module_vfield_wts.o \
           vfield.for
	$(FORT) $(CFLAGS) vfield.for

# ------------------------------ XNFPELSYN -----------------------------

xnfpelsyn : $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES) $(SYNTH_MODULES) \
            module_continuum_edges.o module_synthe_xnfdop.o $(X_OBJECTS)
	$(FORT) $(LFLAGS) -o xnfpelsyn.exe $(TYPES) $(CONSTANTS) \
                                           $(DIMS) $(SA_MODULES) \
                                           $(SYNTH_MODULES) \
                                           module_continuum_edges.o \
                                           module_synthe_xnfdop.o \
                                           $(X_OBJECTS)

#
# ########################### DEPENDENCIES #############################
#

# ----------------------------- BROADEN ---------------------------------

broaden.o : $(TYPES) $(CONSTANTS) $(DIMS) module_intensity_vars.o \
            broaden.for
	$(FORT) $(CFLAGS) broaden.for

# ----------------------------- PREDICT --------------------------------

predict_starting_model.o : $(TYPES) predict_starting_model.for
	$(FORT) $(CFLAGS) predict_starting_model.for

# ----------------------------- REFORM_F -------------------------------

reform_f.o : $(TYPES) $(CONSTANTS) $(DIMS) reform_f.for
	$(FORT) $(CFLAGS) reform_f.for

# ----------------------------- REFORM_I -------------------------------

reform_i.o : $(TYPES) $(CONSTANTS) $(DIMS) reform_i.for
	$(FORT) $(CFLAGS) reform_i.for

# ----------------------------- REFORM_ROT -----------------------------

reform_rot.o : $(TYPES) $(CONSTANTS) $(DIMS) reform_rot.for
	$(FORT) $(CFLAGS) reform_rot.for

# ------------------------------- ROTATE -------------------------------

rotate.o : $(TYPES) $(CONSTANTS) $(DIMS) $(SYN_DIMS) \
           module_intensity_vars.o rotate.for
	$(FORT) $(CFLAGS) rotate.for

# ------------------------ SATLAS_ODF COMPONENTS -----------------------

satlas_odf_main.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_main.for
	$(FORT) $(CFLAGS) satlas_odf_main.for

satlas_odf_conv.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_conv.for
	$(FORT) $(CFLAGS) satlas_odf_conv.for

satlas_odf_josh.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_josh.for
	$(FORT) $(CFLAGS) satlas_odf_josh.for

satlas_odf_kapp.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_kapp.for
	$(FORT) $(CFLAGS) satlas_odf_kapp.for

satlas_odf_pops.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_pops.for
	$(FORT) $(CFLAGS) satlas_odf_pops.for

satlas_odf_read.o :   $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_read.for
	$(FORT) $(CFLAGS) satlas_odf_read.for

satlas_odf_tcorr.o :  $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_tcorr.for
	$(FORT) $(CFLAGS) satlas_odf_tcorr.for

satlas_odf_utils.o :  $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES)\
                          satlas_odf_utils.for
	$(FORT) $(CFLAGS) satlas_odf_utils.for

# ------------------------------ SPECSYN -------------------------------

specsyn_main.o :$(TYPES) $(CONSTANTS) $(DIMS) $(SPECSYN_MODULES) \
                          specsyn_main.for
	$(FORT) $(CFLAGS) specsyn_main.for

# ------------------------------ SYNTHE --------------------------------

synthe_main.o : $(TYPES) $(CONSTANTS) $(DIMS) $(SYN_DIMS) \
                $(SA_MODULES) $(SYNTHE_MODULES) synthe_main.for
	$(FORT) $(CFLAGS) synthe_main.for

# ------------------------ XNFPELSYN COMPONENTS ------------------------

xnfpelsyn_main.o : $(TYPES) $(CONSTANTS) $(DIMS) $(SYNTH_MODULES)\
                   module_continuum_edges.o xnfpelsyn_main.for
	$(FORT) $(CFLAGS)                   xnfpelsyn_main.for

synth_kapp.o : $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES) \
               $(SYNTH_MODULES) module_continuum_edges.o synth_kapp.for
	$(FORT) $(CFLAGS)                                synth_kapp.for

synth_pops.o : $(TYPES) $(CONSTANTS) $(DIMS) $(SA_MODULES) \
               $(SYNTH_MODULES) module_continuum_edges.o synth_pops.for
	$(FORT) $(CFLAGS)                                synth_pops.for

#
# ------------------------------ MODULES -------------------------------
#

# ----------------------------- CONSTANTS ------------------------------

module_astro_parameters.o :   $(TYPES) module_physical_constants.o \
                                       module_astro_parameters.for
	$(FORT) $(CFLAGS)              module_astro_parameters.for

module_physical_constants.o : $(TYPES) module_physical_constants.for
	$(FORT) $(CFLAGS)              module_physical_constants.for

# ----------------------------- DIMENSIONS -----------------------------

module_code_dimensions.o :    $(TYPES) module_code_dimensions.for
	$(FORT) $(CFLAGS)              module_code_dimensions.for

# ------------------------------ TYPES ---------------------------------

module_var_types.o :                   module_var_types.for
	$(FORT) $(CFLAGS)              module_var_types.for

# --------------------------- ROTATE MODULES ---------------------------

module_intensity_wts.o :      $(TYPES) module_intensity_wts.for
	$(FORT) $(CFLAGS)              module_intensity_wts.for

# --------------------------  SATLAS MODULES ---------------------------

module_abross_vars.o :        $(TYPES) $(DIMS) module_abross_vars.for
	$(FORT) $(CFLAGS)                      module_abross_vars.for

module_abtot_vars.o :         $(TYPES) $(DIMS) module_abtot_vars.for
	$(FORT) $(CFLAGS)                      module_abtot_vars.for

module_abundances.o :         $(TYPES) $(DIMS) module_abundances.for
	$(FORT) $(CFLAGS)                      module_abundances.for 

module_atmosphere_parameters.o : $(TYPES) $(DIMS) module_atmosphere_parameters.for
	$(FORT) $(CFLAGS)                         module_atmosphere_parameters.for 

module_conv_vars.o :          $(TYPES) $(DIMS) module_conv_vars.for
	$(FORT) $(CFLAGS)                      module_conv_vars.for

module_depart_vars.o :        $(TYPES) $(DIMS) module_depart_vars.for
	$(FORT) $(CFLAGS)                      module_depart_vars.for

module_edensity_vars.o :      $(TYPES) $(DIMS) module_edensity_vars.for
	$(FORT) $(CFLAGS)                      module_edensity_vars.for

module_elements_vars.o :      $(TYPES) $(DIMS) module_elements_vars.for
	$(FORT) $(CFLAGS)                      module_elements_vars.for

module_flux_vars.o :          $(TYPES) $(DIMS) module_flux_vars.for
	$(FORT) $(CFLAGS)                      module_flux_vars.for

module_freq_set.o :           $(TYPES) $(DIMS) module_freq_set.for
	$(FORT) $(CFLAGS)                      module_freq_set.for 

module_freq_vars.o :          $(TYPES) $(DIMS) module_freq_vars.for
	$(FORT) $(CFLAGS)                      module_freq_vars.for

module_gravity.o :            $(TYPES) $(DIMS) module_gravity.for
	$(FORT) $(CFLAGS)                      module_gravity.for

module_if_vars.o :            $(TYPES) $(DIMS) module_if_vars.for
	$(FORT) $(CFLAGS)                      module_if_vars.for

module_intensity_vars.o :     $(TYPES) $(DIMS) module_intensity_vars.for
	$(FORT) $(CFLAGS)                      module_intensity_vars.for

module_isotope_vars.o :       $(TYPES) $(DIMS) module_isotope_vars.for
	$(FORT) $(CFLAGS)                      module_isotope_vars.for

module_iter_vars.o :          $(TYPES) $(DIMS) module_iter_vars.for
	$(FORT) $(CFLAGS)                      module_iter_vars.for

module_junk_vars.o :          $(TYPES) $(DIMS) module_junk_vars.for
	$(FORT) $(CFLAGS)                      module_junk_vars.for

module_molecular_vars.o :     $(TYPES) $(DIMS) module_molecular_vars.for
	$(FORT) $(CFLAGS)                      module_molecular_vars.for

module_odeint_vars.o :        $(TYPES) $(DIMS) module_odeint_vars.for
	$(FORT) $(CFLAGS)                      module_odeint_vars.for

module_odf_vars.o :           $(TYPES) $(CONSTANTS) $(DIMS) module_odf_vars.for
	$(FORT) $(CFLAGS)                                   module_odf_vars.for

module_opacity.o :            $(TYPES) $(DIMS) module_opacity.for
	$(FORT) $(CFLAGS)                      module_opacity.for

module_opacity_switches.o :            $(DIMS) module_opacity_switches.for
	$(FORT) $(CFLAGS)                      module_opacity_switches.for

module_potion_vars.o :        $(TYPES)         module_potion_vars.for
	$(FORT) $(CFLAGS)                      module_potion_vars.for

module_put_vars.o :           $(TYPES) $(DIMS) module_put_vars.for
	$(FORT) $(CFLAGS)                      module_put_vars.for

module_pzero_vars.o :         $(TYPES) $(DIMS) module_pzero_vars.for
	$(FORT) $(CFLAGS)                      module_pzero_vars.for

module_radius_vars.o :        $(TYPES) $(DIMS) module_radius_vars.for
	$(FORT) $(CFLAGS)                      module_radius_vars.for

module_rad_pressure.o :       $(TYPES) $(DIMS) module_rad_pressure.for
	$(FORT) $(CFLAGS)                      module_rad_pressure.for

module_rad_vars.o :           $(TYPES) $(DIMS) module_rad_vars.for
	$(FORT) $(CFLAGS)                      module_rad_vars.for

module_rhodr_var.o :          $(TYPES) $(DIMS) module_rhodr_var.for
	$(FORT) $(CFLAGS)                      module_rhodr_var.for

module_ross_tables.o :        $(TYPES) $(DIMS) module_ross_tables.for
	$(FORT) $(CFLAGS)                      module_ross_tables.for

module_state_vars.o :         $(TYPES) $(DIMS) module_state_vars.for
	$(FORT) $(CFLAGS)                      module_state_vars.for

module_tabex.o :              $(TYPES)         module_tabex.for
	$(FORT) $(CFLAGS)                      module_tabex.for

module_tau_std.o :            $(TYPES) $(DIMS) module_tau_std.for
	$(FORT) $(CFLAGS)                      module_tau_std.for

module_temp_vars.o :          $(TYPES) $(DIMS) module_temp_vars.for
	$(FORT) $(CFLAGS)                      module_temp_vars.for

module_total_opacity.o :      $(TYPES) $(DIMS) module_total_opacity.for
	$(FORT) $(CFLAGS)                      module_total_opacity.for

module_total_pressure.o :     $(TYPES) $(DIMS) module_total_pressure.for
	$(FORT) $(CFLAGS)                      module_total_pressure.for

module_tsmooth.o :            $(TYPES) $(DIMS) module_tsmooth.for
	$(FORT) $(CFLAGS)                      module_tsmooth.for

module_turbpr_vars.o :        $(TYPES) $(DIMS) module_turbpr_vars.for
	$(FORT) $(CFLAGS)                      module_turbpr_vars.for

module_wave_vars.o :          $(TYPES) $(DIMS) module_wave_vars.for
	$(FORT) $(CFLAGS)                      module_wave_vars.for

module_xnf_vars.o :           $(TYPES) $(DIMS) module_xnf_vars.for
	$(FORT) $(CFLAGS)                      module_xnf_vars.for

module_xnmol.o :              $(TYPES) $(DIMS) module_xnmol.for
	$(FORT) $(CFLAGS)                      module_xnmol.for

# --------------------------- SYNBEG MODULES ---------------------------

module_list_vars.o :          $(TYPES) module_list_vars.for
	$(FORT) $(CFLAGS)              module_list_vars.for

module_logtab.o :             $(TYPES) module_logtab.for
	$(FORT) $(CFLAGS)              module_logtab.for

module_synbeg_syndat.o :      $(TYPES) module_synbeg_syndat.for
	$(FORT) $(CFLAGS)              module_synbeg_syndat.for

module_synth_lindat.o :       $(TYPES) module_synth_lindat.for
	$(FORT) $(CFLAGS)              module_synth_lindat.for

# --------------------------- SYNTH MODULES ----------------------------

module_synth_ionpots.o :      $(TYPES) module_synth_ionpots.for
	$(FORT) $(CFLAGS)              module_synth_ionpots.for

module_synth_xnf_vars.o :     $(TYPES) $(DIMS) module_synth_xnf_vars.for
	$(FORT) $(CFLAGS)                      module_synth_xnf_vars.for

module_synth_xnfh_vars.o :    $(TYPES) $(DIMS) module_synth_xnfh_vars.for
	$(FORT) $(CFLAGS)                      module_synth_xnfh_vars.for

module_synth_xnfp_vars.o :    $(TYPES) $(DIMS) module_synth_xnfp_vars.for
	$(FORT) $(CFLAGS)                      module_synth_xnfp_vars.for

module_synth_xnmol_vars.o :   $(TYPES) $(DIMS) module_synth_xnmol_vars.for
	$(FORT) $(CFLAGS)                      module_synth_xnmol_vars.for

# ---------------------------- SYN DIMENSION ---------------------------

module_synthe_dimensions.o :  $(TYPES) module_synthe_dimensions.for
	$(FORT) $(CFLAGS)              module_synthe_dimensions.for

# ----------------------------- SYN MODULE -----------------------------

module_continuum_edges.o :    $(TYPES)             module_continuum_edges.for
	$(FORT) $(CFLAGS)                          module_continuum_edges.for

module_synthe_buffer.o :      $(TYPES) $(SYN_DIMS) module_synthe_buffer.for
	$(FORT) $(CFLAGS)                          module_synthe_buffer.for

module_synthe_nlines.o :      $(TYPES) $(SYN_DIMS) module_synthe_nlines.for
	$(FORT) $(CFLAGS)                          module_synthe_nlines.for

module_synthe_xnfdop.o :      $(TYPES)             module_synthe_xnfdop.for
	$(FORT) $(CFLAGS)                          module_synthe_xnfdop.for

# ---------------------------- VFIELD MODULE ---------------------------

module_vfield_wts.o :         $(TYPES) module_vfield_wts.for
	$(FORT) $(CFLAGS)              module_vfield_wts.for


### MAKEFILE:  TO MANAGE SATLAS_OS
### 2017 MAR 10: CHANGED EXTENSIONS FROM .f TO .for AND .x TO .exe

# MACRO DEFINITIONS

# THE INTEL FORTRAN COMPILER
FORT = ifort

# THESE COMPILER OPTIONS ARE FOR:
# INTEL FORTRAN COMPILER RUNNING ON A DELL WORKSTATION

CFLAGS = -assume byterecl -c -check bounds -check uninit -warn declarations -warn unused

# SOME POSSIBLE COMPILER OPTIONS
# -assume byterecl = FORCE UNFORMATTED I/O TO USE byte AS UNIT
# -c = COMPILE ONLY, DO NOT LINK
# -check = -check all CHECKS EVERYTHING
# -check bounds = RUN-TIME CHECKS FOR SUBSCRIPT/STRING OUT OF BOUNDS
# -check uninit = CHECK FOR UNINITIALIZED VARIABLES
# -fast = MAXIMIZE SPEED
# -parallel = AUTO-PARALLELIZER GENERATES MULTITHREAD CODE FOR LOOPS
# -stand   CURRENTLY THE COMPILER DEFAULTS TO 2008 STANDARD IF NO STANDARD
#    IS EXPLICITLY GIVEN
# -warn all
# -warn declarations = WARN OF VARIABLES THAT ARE NOT DECLARED
# -warn unused = WARN OF DECLARED VARIABLES THAT ARE NEVER USED

# LINKER FLAGS
LFLAGS = -assume byterecl -check bounds -check uninit

CONSTANTS = module_astro_parameters.o    \
            module_physical_constants.o

DIMENSIONS = module_code_dimensions.o

TYPES = module_var_types.o

SA_MODULES = module_abross_vars.o           \
             module_abtot_vars.o            \
             module_abundances.o            \
             module_atmosphere_parameters.o \
             module_continuum_tables.o      \
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
             module_logtab.o                \
             module_molecular_vars.o        \
             module_odeint_vars.o           \
             module_opacity.o               \
             module_opacity_switches.o      \
             module_potion_vars.o           \
             module_put_vars.o              \
             module_pzero_vars.o            \
             module_radius_vars.o           \
             module_rad_pressure.o          \
             module_rad_vars.o              \
             module_rhodr_var.o             \
             module_state_vars.o            \
             module_tabex.o                 \
             module_tau_std.o               \
             module_temp_vars.o             \
             module_total_opacity.o         \
             module_total_pressure.o        \
             module_tsmooth.o               \
             module_turbpr_vars.o           \
             module_txnxn_vars.o            \
             module_xlines_vars.o           \
             module_xnf_dop.o               \
             module_xnf_vars.o              \
             module_xnmol.o
 
SA_OBJECTS = satlas_os_main.o  \
             satlas_os_conv.o  \
             satlas_os_josh.o  \
             satlas_os_kapp.o  \
             satlas_os_pops.o  \
             satlas_os_read.o  \
             satlas_os_tcorr.o \
             satlas_os_utils.o

# target definitions

# ----------------------------- SATLAS_OS ------------------------------

satlas_os: $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                                               $(SA_OBJECTS)
	$(FORT) $(LFLAGS) -o satlas_os.exe  $(TYPES) $(DIMENSIONS) \
                                            $(CONSTANTS) $(SA_MODULES) \
                                            $(SA_OBJECTS)

#
# ########################### DEPENDENCIES #############################
#

# ----------------------- SATLAS_OS COMPONENTS -------------------------

satlas_os_main.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_main.for
	$(FORT) $(CFLAGS) satlas_os_main.for

satlas_os_conv.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_conv.for
	$(FORT) $(CFLAGS) satlas_os_conv.for

satlas_os_josh.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_josh.for
	$(FORT) $(CFLAGS) satlas_os_josh.for

satlas_os_kapp.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_kapp.for
	$(FORT) $(CFLAGS) satlas_os_kapp.for

satlas_os_pops.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_pops.for
	$(FORT) $(CFLAGS) satlas_os_pops.for

satlas_os_read.o :  $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_read.for
	$(FORT) $(CFLAGS) satlas_os_read.for

satlas_os_tcorr.o : $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_tcorr.for
	$(FORT) $(CFLAGS) satlas_os_tcorr.for

satlas_os_utils.o : $(TYPES) $(DIMENSIONS) $(CONSTANTS) $(SA_MODULES) \
                          satlas_os_utils.for
	$(FORT) $(CFLAGS) satlas_os_utils.for

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

# ------------------------------- SATLAS -------------------------------

module_abross_vars.o :        $(TYPES) $(DIMENSIONS) module_abross_vars.for
	$(FORT) $(CFLAGS)                            module_abross_vars.for

module_abtot_vars.o :         $(TYPES) $(DIMENSIONS) module_abtot_vars.for
	$(FORT) $(CFLAGS)                            module_abtot_vars.for

module_abundances.o :         $(TYPES) $(DIMENSIONS) module_abundances.for
	$(FORT) $(CFLAGS)                            module_abundances.for 

module_atmosphere_parameters.o : $(TYPES) $(DIMENSIONS) module_atmosphere_parameters.for
	$(FORT) $(CFLAGS)                               module_atmosphere_parameters.for 

module_continuum_tables.o :   $(TYPES) $(DIMENSIONS) module_continuum_tables.for
	$(FORT) $(CFLAGS)                            module_continuum_tables.for

module_conv_vars.o :          $(TYPES) $(DIMENSIONS) module_conv_vars.for
	$(FORT) $(CFLAGS)                            module_conv_vars.for

module_depart_vars.o :        $(TYPES) $(DIMENSIONS) module_depart_vars.for
	$(FORT) $(CFLAGS)                            module_depart_vars.for

module_edensity_vars.o :      $(TYPES) $(DIMENSIONS) module_edensity_vars.for
	$(FORT) $(CFLAGS)                            module_edensity_vars.for

module_elements_vars.o :      $(TYPES) $(DIMENSIONS) module_elements_vars.for
	$(FORT) $(CFLAGS)                            module_elements_vars.for

module_flux_vars.o :          $(TYPES) $(DIMENSIONS) module_flux_vars.for
	$(FORT) $(CFLAGS)                            module_flux_vars.for

module_freq_set.o :           $(TYPES) $(DIMENSIONS) module_freq_set.for
	$(FORT) $(CFLAGS)                            module_freq_set.for 

module_freq_vars.o :          $(TYPES) $(DIMENSIONS) module_freq_vars.for
	$(FORT) $(CFLAGS)                            module_freq_vars.for

module_gravity.o :            $(TYPES) $(DIMENSIONS) module_gravity.for
	$(FORT) $(CFLAGS)                            module_gravity.for

module_if_vars.o :            $(TYPES) $(DIMENSIONS) module_if_vars.for
	$(FORT) $(CFLAGS)                            module_if_vars.for

module_intensity_vars.o :     $(TYPES) $(DIMENSIONS) module_intensity_vars.for
	$(FORT) $(CFLAGS)                            module_intensity_vars.for

module_isotope_vars.o :       $(TYPES) $(DIMENSIONS) module_isotope_vars.for
	$(FORT) $(CFLAGS)                            module_isotope_vars.for

module_iter_vars.o :          $(TYPES) $(DIMENSIONS) module_iter_vars.for
	$(FORT) $(CFLAGS)                            module_iter_vars.for

module_junk_vars.o :          $(TYPES) $(DIMENSIONS) module_junk_vars.for
	$(FORT) $(CFLAGS)                            module_junk_vars.for

module_logtab.o :             $(TYPES)               module_logtab.for
	$(FORT) $(CFLAGS)                            module_logtab.for

module_molecular_vars.o :     $(TYPES) $(DIMENSIONS) module_molecular_vars.for
	$(FORT) $(CFLAGS)                            module_molecular_vars.for

module_odeint_vars.o :        $(TYPES) $(DIMENSIONS) module_odeint_vars.for
	$(FORT) $(CFLAGS)                            module_odeint_vars.for

module_opacity.o :            $(TYPES) $(DIMENSIONS) module_opacity.for
	$(FORT) $(CFLAGS)                            module_opacity.for

module_opacity_switches.o :            $(DIMENSIONS) module_opacity_switches.for
	$(FORT) $(CFLAGS)                            module_opacity_switches.for

module_potion_vars.o :        $(TYPES)               module_potion_vars.for
	$(FORT) $(CFLAGS)                            module_potion_vars.for

module_put_vars.o :           $(TYPES) $(DIMENSIONS) module_put_vars.for
	$(FORT) $(CFLAGS)                            module_put_vars.for

module_pzero_vars.o :         $(TYPES) $(DIMENSIONS) module_pzero_vars.for
	$(FORT) $(CFLAGS)                            module_pzero_vars.for

module_radius_vars.o :        $(TYPES) $(DIMENSIONS) module_radius_vars.for
	$(FORT) $(CFLAGS)                            module_radius_vars.for

module_rad_pressure.o :       $(TYPES) $(DIMENSIONS) module_rad_pressure.for
	$(FORT) $(CFLAGS)                            module_rad_pressure.for

module_rad_vars.o :           $(TYPES) $(DIMENSIONS) module_rad_vars.for
	$(FORT) $(CFLAGS)                            module_rad_vars.for

module_rhodr_var.o :          $(TYPES) $(DIMENSIONS) module_rhodr_var.for
	$(FORT) $(CFLAGS)                            module_rhodr_var.for

module_state_vars.o :         $(TYPES) $(DIMENSIONS) module_state_vars.for
	$(FORT) $(CFLAGS)                            module_state_vars.for

module_tabex.o :              $(TYPES) $(DIMENSIONS) module_tabex.for
	$(FORT) $(CFLAGS)                            module_tabex.for

module_tau_std.o :            $(TYPES) $(DIMENSIONS) module_tau_std.for
	$(FORT) $(CFLAGS)                            module_tau_std.for

module_temp_vars.o :          $(TYPES) $(DIMENSIONS) module_temp_vars.for
	$(FORT) $(CFLAGS)                            module_temp_vars.for 

module_total_opacity.o :      $(TYPES) $(DIMENSIONS) module_total_opacity.for
	$(FORT) $(CFLAGS)                            module_total_opacity.for

module_total_pressure.o :     $(TYPES) $(DIMENSIONS) module_total_pressure.for
	$(FORT) $(CFLAGS)                            module_total_pressure.for

module_tsmooth.o :            $(TYPES) $(DIMENSIONS) module_tsmooth.for
	$(FORT) $(CFLAGS)                            module_tsmooth.for

module_turbpr_vars.o :        $(TYPES) $(DIMENSIONS) module_turbpr_vars.for
	$(FORT) $(CFLAGS)                            module_turbpr_vars.for

module_txnxn_vars.o :         $(TYPES) $(DIMENSIONS) module_txnxn_vars.for
	$(FORT) $(CFLAGS)                            module_txnxn_vars.for

module_xlines_vars.o :        $(TYPES) $(DIMENSIONS) module_xlines_vars.for
	$(FORT) $(CFLAGS)                            module_xlines_vars.for

module_xnf_dop.o :            $(TYPES) $(DIMENSIONS) module_xnf_dop.for
	$(FORT) $(CFLAGS)                            module_xnf_dop.for

module_xnf_vars.o :           $(TYPES) $(DIMENSIONS) module_xnf_vars.for
	$(FORT) $(CFLAGS)                            module_xnf_vars.for

module_xnmol.o :              $(TYPES) $(DIMENSIONS) module_xnmol.for
	$(FORT) $(CFLAGS)                            module_xnmol.for


SHELL = /bin/sh

include conf/configure.fv3

ifeq ($(strip $(FMS_DIR)),)
  FMS_DIR=$(realpath ../FMS/FMS_INSTALL)
endif

exist=$(wildcard $(FMS_DIR))
ifeq ($(strip $(exist)),)
  $(error ERROR: FMS_DIR variable is unset and FMS_INSTALL is not in ../FMS/FMS_INSTALL )
endif

ifeq ($(NAM_phys),Y)
  PHYSP  = nam
else
  PHYSP  = gfs
endif

FFLAGS   += -I$(FMS_DIR) -I$(PHYSP)physics -Iipd -Icpl -Iio -Iatmos_cubed_sphere

# Flag to CCPP build for 32bit dynamics
ifeq ($(32BIT),Y)
  DYN32 = Y
else
  DYN32 = N
endif

FV3_EXE  = fv3.exe
FV3CAP_LIB  = libfv3cap.a

all: libs
	$(MAKE) $(FV3_EXE) $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)

nems: libs
	$(MAKE) $(FV3CAP_LIB) $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) esmf_make_fragment FMS_DIR=$(FMS_DIR)

# Remove stochastic_physics from CCPP build
ifneq (,$(findstring CCPP,$(CPPDEFS)))
libs:
	$(MAKE) -C cpl                 $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) -C $(PHYSP)physics     $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  DYN32=$(DYN32) # force gfs physics to 64bit, flag to CCPP build for 32bit dynamics
	$(MAKE) -C ipd                 $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  # force gfs physics to 64bit
	$(MAKE) -C io                  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) -C atmos_cubed_sphere  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	#$(MAKE) -C stochastic_physics  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  # force gfs physics to 64bit

$(FV3_EXE): atmos_model.o coupler_main.o atmos_cubed_sphere/libfv3core.a io/libfv3io.a ipd/libipd.a $(PHYSP)physics/lib$(PHYSP)phys.a cpl/libfv3cpl.a
	$(LD) -o $@ $^ $(NCEPLIBS) $(LDFLAGS)

else
libs:
	$(MAKE) -C cpl                 $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) -C $(PHYSP)physics     $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  # force gfs physics to 64bit
	$(MAKE) -C ipd                 $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  # force gfs physics to 64bit
	$(MAKE) -C io                  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) -C atmos_cubed_sphere  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR)
	$(MAKE) -C stochastic_physics  $(MAKEFLAGS) FMS_DIR=$(FMS_DIR) 32BIT=N  # force gfs physics to 64bit

$(FV3_EXE): atmos_model.o coupler_main.o atmos_cubed_sphere/libfv3core.a io/libfv3io.a ipd/libipd.a $(PHYSP)physics/lib$(PHYSP)phys.a stochastic_physics/libstochastic_physics.a cpl/libfv3cpl.a
	$(LD) -o $@ $^ $(NCEPLIBS) $(LDFLAGS)

endif

$(FV3CAP_LIB): atmos_model.o module_fv3_config.o module_fcst_grid_comp.o time_utils.o fv3_cap.o
	ar rv $(FV3CAP_LIB) $?

atmos_model.o : atmos_model.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c atmos_model.F90

module_fv3_config.o: module_fv3_config.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c module_fv3_config.F90
module_fcst_grid_comp.o: module_fcst_grid_comp.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c module_fcst_grid_comp.F90
time_utils.o: time_utils.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c time_utils.F90
fv3_cap.o: fv3_cap.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) $(ESMF_INC) -c fv3_cap.F90

DEPEND_FILES = time_utils.F90 module_fv3_config.F90 atmos_model.F90 module_fcst_grid_comp.F90 fv3_cap.F90 coupler_main.F90

ifneq (,$(findstring SION,$(CPPDEFS)))
  SIONLIB_LINK_FLAGS = $(SIONLIB_LIB)
else
  SIONLIB_LINK_FLAGS =
endif

ifneq (,$(findstring CCPP,$(CPPDEFS)))
  ifneq (,$(findstring STATIC,$(CPPDEFS)))
    CCPP_LINK_FLAGS = -L$(PATH_CCPP)/lib -lccppphys -lccpp $(NCEPLIBS) -lxml2
  else
    CCPP_LINK_FLAGS = -L$(PATH_CCPP)/lib -lccpp
  endif
endif


# Remove stochastic_physics from CCPP build
ifneq (,$(findstring CCPP,$(CPPDEFS)))
esmf_make_fragment:
	@rm -rf nems_dir; mkdir nems_dir
	@cp $(FV3CAP_LIB) atmos_cubed_sphere/libfv3core.a io/libfv3io.a ipd/libipd.a $(PHYSP)physics/lib$(PHYSP)phys.a cpl/libfv3cpl.a nems_dir
	@cp fv3gfs_cap_mod.mod nems_dir
	@echo "# ESMF self-describing build dependency makefile fragment" > fv3.mk
	@echo "# src location $(PWD)" >> fv3.mk
	@echo  >> fv3.mk
	@echo "ESMF_DEP_FRONT     = fv3gfs_cap_mod"  >> fv3.mk
	# additional include files needed for PGI
	#@echo "ESMF_DEP_INCPATH   = $(PWD)/nems_dir" >> fv3.mk
	@echo "ESMF_DEP_INCPATH   = $(PWD) $(addprefix $(PWD)/, nems_dir atmos_cubed_sphere io gfsphysics cpl ipd)" >> fv3.mk
	@echo "ESMF_DEP_CMPL_OBJS ="                 >> fv3.mk
	@echo "ESMF_DEP_LINK_OBJS = $(addprefix $(PWD)/nems_dir/, libfv3cap.a libfv3core.a libfv3io.a libipd.a lib$(PHYSP)phys.a libfv3cpl.a) $(CCPP_LINK_FLAGS) $(SIONLIB_LINK_FLAGS)"  >> fv3.mk
	@echo "ESMF_DEP_SHRD_PATH ="                 >> fv3.mk
	@echo "ESMF_DEP_SHRD_LIBS ="                 >> fv3.mk
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" fv3.mk
	@echo
else
esmf_make_fragment:
	@rm -rf nems_dir; mkdir nems_dir
	@cp $(FV3CAP_LIB) atmos_cubed_sphere/libfv3core.a io/libfv3io.a ipd/libipd.a $(PHYSP)physics/lib$(PHYSP)phys.a cpl/libfv3cpl.a stochastic_physics/libstochastic_physics.a nems_dir
	@cp fv3gfs_cap_mod.mod nems_dir
	@echo "# ESMF self-describing build dependency makefile fragment" > fv3.mk
	@echo "# src location $(PWD)" >> fv3.mk
	@echo  >> fv3.mk
	@echo "ESMF_DEP_FRONT     = fv3gfs_cap_mod"  >> fv3.mk
	# additional include files needed for PGI
	#@echo "ESMF_DEP_INCPATH   = $(PWD)/nems_dir" >> fv3.mk
	@echo "ESMF_DEP_INCPATH   = $(PWD) $(addprefix $(PWD)/, nems_dir atmos_cubed_sphere io gfsphysics cpl ipd)" >> fv3.mk
	@echo "ESMF_DEP_CMPL_OBJS ="                 >> fv3.mk
	@echo "ESMF_DEP_LINK_OBJS = $(addprefix $(PWD)/nems_dir/, libfv3cap.a libfv3core.a libfv3io.a libipd.a lib$(PHYSP)phys.a libfv3cpl.a libstochastic_physics.a) $(SIONLIB_LINK_FLAGS)"  >> fv3.mk
	@echo "ESMF_DEP_SHRD_PATH ="                 >> fv3.mk
	@echo "ESMF_DEP_SHRD_LIBS ="                 >> fv3.mk
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" fv3.mk
	@echo
endif

# fv3 library installation defaults (for NEMS):
DESTDIR  := $(PWD)
INSTDIR  := FV3_INSTALL

nemsinstall: nems
	@mkdir -p $(DESTDIR)/$(INSTDIR)
	@cp nems_dir/* $(DESTDIR)/$(INSTDIR)
	@sed -e 's;'$(PWD)/nems_dir';'$(DESTDIR)/$(INSTDIR)';g' fv3.mk > $(DESTDIR)/$(INSTDIR)/fv3.mk
	@echo Installation into \"$(DESTDIR)/$(INSTDIR)\" complete!
	@echo

.PHONY: clean cleanall
clean:
	@echo "Cleaning ... "
	@echo
	(cd $(PHYSP)physics     && make clean)
	(cd ipd                 && make clean)
	(cd stochastic_physics  && make clean)
	(cd io                  && make clean)
	(cd atmos_cubed_sphere  && make clean)
	(cd cpl                 && make clean)
	$(RM) -f $(FV3_EXE) $(FV3CAP_LIB) *.o *.mod *.i90 *.lst depend

cleanall: clean
	$(RM) -rf nems_dir fv3.mk $(INSTDIR)
	$(RM) -f conf/modules.fv3
	$(RM) -f conf/configure.fv3

MKDEPENDS = ./mkDepends.pl
include conf/make.rules

# do not include 'depend' file if the target contains string 'clean'
ifneq (clean,$(findstring clean,$(MAKECMDGOALS)))
    -include depend
endif


# Source lists, compilation rules and targets for FastEarth3D.
# Included by config/Makefile after the flag sets are assembled.

# --- Library object list (in module-dependency order) ------------------------
obj_fastearth = \
	$(objdir)/fe_precision.o \
	$(objdir)/fe_constants.o \
	$(objdir)/fe_params.o \
	$(objdir)/fe_control.o \
	$(objdir)/fe_sht.o \
	$(objdir)/fe_tensor_sh.o \
	$(objdir)/fe_field.o \
	$(objdir)/fe_earth_structure.o \
	$(objdir)/fe_radial_integrals.o \
	$(objdir)/fe_band.o \
	$(objdir)/fe_radial_fe.o \
	$(objdir)/fe_viscoelastic.o \
	$(objdir)/fe_modal.o \
	$(objdir)/fe_response.o \
	$(objdir)/fe_sle.o \
	$(objdir)/fe_timestep.o \
	$(objdir)/fe_rotation.o \
	$(objdir)/fe_remap.o \
	$(objdir)/fe_vilma.o \
	$(objdir)/fe_coupling.o \
	$(objdir)/fe_io.o \
	$(objdir)/fe_drive.o \
	$(objdir)/fastearth3d.o

# --- Inter-module dependencies (so `make -j` stays correct) ------------------
$(objdir)/fe_constants.o:        $(objdir)/fe_precision.o
$(objdir)/fe_params.o:           $(objdir)/fe_precision.o $(objdir)/fe_constants.o
$(objdir)/fe_control.o:          $(objdir)/fe_precision.o $(objdir)/fe_constants.o $(objdir)/fe_params.o
$(objdir)/fe_sht.o:              $(objdir)/fe_precision.o
$(objdir)/fe_tensor_sh.o:        $(objdir)/fe_precision.o $(objdir)/fe_sht.o \
                                 $(objdir)/fe_constants.o
$(objdir)/fe_field.o:            $(objdir)/fe_precision.o $(objdir)/fe_sht.o
$(objdir)/fe_earth_structure.o:  $(objdir)/fe_precision.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_params.o $(objdir)/fe_sht.o
$(objdir)/fe_radial_integrals.o: $(objdir)/fe_precision.o
$(objdir)/fe_band.o:             $(objdir)/fe_precision.o
$(objdir)/fe_radial_fe.o:        $(objdir)/fe_constants.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_radial_integrals.o $(objdir)/fe_band.o
$(objdir)/fe_viscoelastic.o:     $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o
$(objdir)/fe_modal.o:            $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_viscoelastic.o $(objdir)/fe_precision.o
$(objdir)/fe_response.o:         $(objdir)/fe_radial_fe.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_tensor_sh.o \
                                 $(objdir)/fe_constants.o $(objdir)/fe_viscoelastic.o \
                                 $(objdir)/fe_modal.o
$(objdir)/fe_sle.o:              $(objdir)/fe_sht.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_response.o
$(objdir)/fe_timestep.o:         $(objdir)/fe_response.o $(objdir)/fe_sle.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_viscoelastic.o \
                                 $(objdir)/fe_precision.o
$(objdir)/fe_rotation.o:         $(objdir)/fe_sht.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_earth_structure.o $(objdir)/fe_radial_fe.o \
                                 $(objdir)/fe_viscoelastic.o
# Optional VILMA backend. Compiled ALWAYS; with vilma=0 (the default) -DVILMA is
# absent and this is a pure-Fortran stub that references no VILMA symbol.
$(objdir)/fe_vilma.o:            $(objdir)/fe_precision.o $(objdir)/fe_constants.o \
                                 $(objdir)/fe_params.o $(objdir)/fe_sht.o
$(objdir)/fe_coupling.o:         $(objdir)/fe_vilma.o \
                                 $(objdir)/fe_response.o $(objdir)/fe_sle.o \
                                 $(objdir)/fe_rotation.o $(objdir)/fe_earth_structure.o \
                                 $(objdir)/fe_sht.o $(objdir)/fe_params.o $(objdir)/fe_remap.o \
                                 $(objdir)/fe_timestep.o $(objdir)/fe_viscoelastic.o
$(objdir)/fe_remap.o:            $(objdir)/fe_precision.o $(objdir)/fe_sht.o
$(objdir)/fe_io.o:               $(objdir)/fe_coupling.o $(objdir)/fe_response.o \
                                 $(objdir)/fe_viscoelastic.o $(objdir)/fe_sht.o \
                                 $(objdir)/fe_constants.o
$(objdir)/fe_drive.o:            $(objdir)/fe_params.o $(objdir)/fe_control.o \
                                 $(objdir)/fe_sht.o \
                                 $(objdir)/fe_coupling.o $(objdir)/fe_io.o \
                                 $(objdir)/fe_remap.o \
                                 $(objdir)/fe_constants.o $(objdir)/fe_precision.o
# The umbrella module re-exports every component, so it compiles last.
$(objdir)/fastearth3d.o:         $(objdir)/fe_coupling.o $(objdir)/fe_io.o \
                                 $(objdir)/fe_drive.o $(objdir)/fe_params.o \
                                 $(objdir)/fe_control.o

# --- Pattern rule ------------------------------------------------------------
# Every object depends on FESMUTILS_LIB (libfesmutils.a): fesm-utils modules
# (coords/ncio/nml) are used across FastEarth3D, and a rebuilt fesm-utils with a
# changed module interface must force a recompile against the new .mod rather than
# silently relinking a stale object (ABI mismatch -> segfault, the i_geo=3 coords
# crash). Blanket dep keeps this correct without tracking which objects use which
# fesm-utils module. FESMUTILS_LIB (config/common.mk) is required by every build.
$(objdir)/%.o: $(srcdir)/%.f90 $(FESMUTILS_LIB) | $(objdir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) -c $< -o $@

$(objdir):
	mkdir -p $(objdir)

$(bindir):
	mkdir -p $(bindir)

# --- Library -----------------------------------------------------------------
fastearth-static: $(obj_fastearth)
	ar rcs $(objdir)/libfastearth.a $(obj_fastearth)
	@echo ""
	@echo "    $(objdir)/libfastearth.a is ready."
	@echo ""

# --- Standalone driver -------------------------------------------------------
fastearth: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/fastearth_main.f90 \
		-o $(bindir)/fastearth.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/fastearth.x is ready."

# --- Offline lon-lat -> Gauss remapper ---------------------------------------
fastearth_remap: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/fastearth_remap_main.f90 \
		-o $(bindir)/fastearth_remap.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/fastearth_remap.x is ready."

# --- Offline reference (bed + ice) generation onto the Gauss grid ------------
fastearth_mkref: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(srcdir)/fastearth_mkref_main.f90 \
		-o $(bindir)/fastearth_mkref.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/fastearth_mkref.x is ready."

# --- Tests -------------------------------------------------------------------
test_params: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_params.f90 \
		-o $(bindir)/test_params.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_params.x is ready."

test_drive: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_drive.f90 \
		-o $(bindir)/test_drive.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_drive.x is ready."

test_band: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_band.f90 \
		-o $(bindir)/test_band.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_band.x is ready."

test_sht: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sht.f90 \
		-o $(bindir)/test_sht.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sht.x is ready."

test_earth: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_earth.f90 \
		-o $(bindir)/test_earth.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_earth.x is ready."

test_mesh: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_mesh.f90 \
		-o $(bindir)/test_mesh.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_mesh.x is ready."

test_integrals: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_integrals.f90 \
		-o $(bindir)/test_integrals.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_integrals.x is ready."

test_assembly: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_assembly.f90 \
		-o $(bindir)/test_assembly.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_assembly.x is ready."

test_love: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_love.f90 \
		-o $(bindir)/test_love.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_love.x is ready."

# Reference dumper for the Julia port (writes NetCDF oracle data; not in TESTS).
dump_reference: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/dump_reference.f90 \
		-o $(bindir)/dump_reference.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/dump_reference.x is ready."

test_relax: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_relax.f90 \
		-o $(bindir)/test_relax.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_relax.x is ready."

test_tidal: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_tidal.f90 \
		-o $(bindir)/test_tidal.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_tidal.x is ready."

test_rotation: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotation.f90 \
		-o $(bindir)/test_rotation.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_rotation.x is ready."

test_rotation_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotation_sle.f90 \
		-o $(bindir)/test_rotation_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_rotation_sle.x is ready."

test_etd1: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_etd1.f90 \
		-o $(bindir)/test_etd1.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_etd1.x is ready."

# Coupling-order characterization (§3c): measures the convergence order of the
# strain<->memory coupling. Standalone like test_etd1 -- a dt-sweep diagnostic,
# intentionally NOT in TESTS / `make check`. Build + run directly.
test_couple_order: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_couple_order.f90 \
		-o $(bindir)/test_couple_order.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_couple_order.x is ready."

# SLE<->memory coupling-order characterization (§3c 3b): drives a fast-evolving load
# through the full sea-level driver and measures the order restored by co-converging
# the ocean load σ and the end-of-step memory τ. Standalone dt-sweep diagnostic, like
# test_couple_order -- intentionally NOT in TESTS / `make check`. Build + run directly.
test_sle_couple_order: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_couple_order.f90 \
		-o $(bindir)/test_sle_couple_order.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_couple_order.x is ready."

# Adaptive-Δt controller (§3c): field step-doubling estimate order + the adaptive
# stepper converging to a fine reference with far fewer steps. Standalone diagnostic,
# NOT in `make check`.
test_timestep: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_timestep.f90 \
		-o $(bindir)/test_timestep.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_timestep.x is ready."

test_response: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response.f90 \
		-o $(bindir)/test_response.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_response.x is ready."

test_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle.f90 \
		-o $(bindir)/test_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle.x is ready."

test_flotation: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation.f90 \
		-o $(bindir)/test_flotation.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation.x is ready."

test_ve_response: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_ve_response.f90 \
		-o $(bindir)/test_ve_response.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_ve_response.x is ready."

test_tensor_sh: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_tensor_sh.f90 \
		-o $(bindir)/test_tensor_sh.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_tensor_sh.x is ready."

test_response_3d: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_response_3d.f90 \
		-o $(bindir)/test_response_3d.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_response_3d.x is ready."

test_visc_load: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_visc_load.f90 \
		-o $(bindir)/test_visc_load.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_visc_load.x is ready."

test_rotinv: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_rotinv.f90 \
		-o $(bindir)/test_rotinv.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_rotinv.x is ready."

test_benchmark_lvz: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_lvz.f90 \
		-o $(bindir)/test_benchmark_lvz.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_lvz.x is ready."

test_sle_ve: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_ve.f90 \
		-o $(bindir)/test_sle_ve.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_ve.x is ready."

test_benchmark_love: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_love.f90 \
		-o $(bindir)/test_benchmark_love.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_love.x is ready."

test_coupling: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_coupling.f90 \
		-o $(bindir)/test_coupling.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_coupling.x is ready."

test_restart: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_restart.f90 \
		-o $(bindir)/test_restart.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_restart.x is ready."

test_couple_remap: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_couple_remap.f90 \
		-o $(bindir)/test_couple_remap.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_couple_remap.x is ready."

test_spinup: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_spinup.f90 \
		-o $(bindir)/test_spinup.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_spinup.x is ready."

test_benchmark_disc: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_disc.f90 \
		-o $(bindir)/test_benchmark_disc.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_disc.x is ready."

test_benchmark_martinec: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_martinec.f90 \
		-o $(bindir)/test_benchmark_martinec.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_martinec.x is ready."

test_field: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_field.f90 \
		-o $(bindir)/test_field.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_field.x is ready."

test_flotation_load: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_flotation_load.f90 \
		-o $(bindir)/test_flotation_load.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_flotation_load.x is ready."

test_marine_reference: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_marine_reference.f90 \
		-o $(bindir)/test_marine_reference.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_marine_reference.x is ready."

# Standalone SLE benchmark (Martinec-2018 case E2): ~750 steps at lmax=128, runs
# in minutes -- intentionally NOT in TESTS / `make check`. Build with `make
# openmp=1 test_benchmark_sle` and run $(bindir)/test_benchmark_sle.x directly.
test_benchmark_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_benchmark_sle.f90 \
		-o $(bindir)/test_benchmark_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_benchmark_sle.x is ready."

test_sle_eustatic: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_eustatic.f90 \
		-o $(bindir)/test_sle_eustatic.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_eustatic.x is ready."

test_sle_subgrid: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_sle_subgrid.f90 \
		-o $(bindir)/test_sle_subgrid.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_sle_subgrid.x is ready."

test_remap: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_remap.f90 \
		-o $(bindir)/test_remap.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_remap.x is ready."

test_modal: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_modal.f90 \
		-o $(bindir)/test_modal.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_modal.x is ready."

test_modal_resp: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_modal_resp.f90 \
		-o $(bindir)/test_modal_resp.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_modal_resp.x is ready."

test_modal_visc3d: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/test_modal_visc3d.f90 \
		-o $(bindir)/test_modal_visc3d.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/test_modal_visc3d.x is ready."

# Diagnostic (not in `check`): sweep the modal Krylov block size n_krylov and
# report convergence of n_modes=all toward RESP_VE over the full degree spectrum.
diag_modal_pblock: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_modal_pblock.f90 \
		-o $(bindir)/diag_modal_pblock.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/diag_modal_pblock.x is ready."

diag_modal_ramp: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_modal_ramp.f90 \
		-o $(bindir)/diag_modal_ramp.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/diag_modal_ramp.x is ready."

diag_modal_sle: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_modal_sle.f90 \
		-o $(bindir)/diag_modal_sle.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/diag_modal_sle.x is ready."

diag_modal_latvisc: fastearth-static | $(bindir)
	$(FC) $(DFLAGS) $(CPPFLAGS) $(FFLAGS) $(testdir)/diag_modal_latvisc.f90 \
		-o $(bindir)/diag_modal_latvisc.x $(objdir)/libfastearth.a $(LFLAGS)
	@echo "    $(bindir)/diag_modal_latvisc.x is ready."

TESTS = test_params test_drive test_band test_sht test_earth test_mesh test_integrals test_assembly test_love test_relax test_tidal test_rotation test_rotation_sle test_response test_sle test_flotation test_flotation_load test_marine_reference test_etd1 test_ve_response test_tensor_sh test_response_3d test_sle_ve test_benchmark_love test_coupling test_couple_remap test_spinup test_restart test_benchmark_disc test_benchmark_martinec test_field test_sle_subgrid test_visc_load test_rotinv test_remap test_modal test_modal_resp test_modal_visc3d

check: $(TESTS)
	@echo ""
	@echo "=== Running FastEarth3D test suite ==="
	@for t in $(TESTS); do \
		echo "--- $$t ---"; \
		$(bindir)/$$t.x || exit 1; \
	done
	@echo ""
	@echo "=== All tests passed ==="

# Slow quantitative benchmarks: full-resolution, multi-minute runs kept OUT of
# `make check`. Build with `make openmp=1 check-slow` (the per-step solve is ~100x
# faster threaded). The SLE benchmark sweeps all four Martinec migrating-coast cases;
# test_rotinv re-runs the off-pole rotational-invariance check at full resolution
# (lmax 128, vs lmax 16 in `make check`).
SLOW = test_benchmark_sle test_rotinv test_benchmark_lvz
check-slow: $(SLOW)
	@echo ""
	@echo "=== Running FastEarth3D slow benchmarks ==="
	@for c in C2 D3 E2 F1; do \
		echo "--- test_benchmark_sle $$c ---"; \
		$(bindir)/test_benchmark_sle.x $$c || exit 1; \
	done
	@echo "--- test_rotinv (lmax 128) ---"; $(bindir)/test_rotinv.x 128 || exit 1
	@echo ""
	@echo "=== All slow benchmarks passed ==="

# --- Housekeeping ------------------------------------------------------------
.PHONY: usage check check-slow clean showconfig fastearth fastearth-static

usage:
	@echo ""
	@echo "    * FastEarth3D build *"
	@echo ""
	@echo " make fastearth-static : build libfastearth.a"
	@echo " make fastearth        : build the standalone driver (bin/fastearth.x)"
	@echo " make check            : build + run the test suite"
	@echo " make openmp=1 check-slow : build + run the slow full-res benchmarks"
	@echo " make test_sht         : build the SHT round-trip test"
	@echo " make clean            : remove objects and binaries"
	@echo " make showconfig       : show the active build configuration"
	@echo ""
	@echo "   switches:  debug=0|1|2   openmp=0|1   vilma=0|1"
	@echo ""
	@echo "   vilma=1 additionally links the optional VILMA backend (solver=\"vilma\"):"
	@echo "     make fastearth vilma=1 VILMAROOT=/path/to/vilma"
	@echo "   It is OFF by default and is not a dependency; see doc/vilma-backend.md."
	@echo ""

showconfig:
	@echo "----------------------"
	@echo "FastEarth3D build configuration"
	@echo "----------------------"
	@echo "compiler  : $(FC)"
	@echo "host      : $(shell hostname)"
	@echo "openmp    : $(openmp)"
	@echo "debug     : $(debug)"
	@echo "vilma     : $(vilma)   (VILMAROOT=$(VILMAROOT))"
	@echo "FFLAGS    : $(FFLAGS)"
	@echo "LFLAGS    : $(LFLAGS)"

clean:
	rm -f $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a
	rm -f $(bindir)/*.x
	rm -rf $(bindir)/*.dSYM

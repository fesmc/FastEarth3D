module fastearth3d
   !! FastEarth3D umbrella module — single `use fastearth3d` entry point that
   !! re-exports the public API of every component (including the parameter record
   !! fe_param_class / fe_par_load and the standalone driver fastearth_run). Host
   !! models and tests should depend on this rather than the individual fe_* modules.
   use fe_precision,       only: wp, sp, dp
   use fe_constants
   use fe_params,          only: fe_param_class, fe_par_load, fe_par_print
   use fe_control,         only: fe_ctl_class, fe_ctl_load, fe_ctl_print, DEFAULTS_FILE
   use fe_sht,             only: sht_grid
   use fe_earth_structure, only: earth_model, build_earth, build_M3L70V01
   use fe_radial_fe,       only: radial_operator
   use fe_viscoelastic,    only: ve_degree
   use fe_sle,             only: sle_solver
   use fe_rotation,        only: rotation_state
   use fe_coupling,        only: solid_earth
   use fe_io,              only: fe_restart_write, fe_restart_read, fe_write_step, &
                                 fe_io_set_table
   use fe_drive,           only: fastearth_run
   implicit none
   public

#ifndef VERSION
#define VERSION "unknown"
#endif
   character(len=*), parameter :: fastearth_version = VERSION

end module fastearth3d

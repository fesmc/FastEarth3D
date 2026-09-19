program fastearth
   !! Standalone FastEarth3D driver.
   !!
   !!   ./bin/fastearth.x [run-config.nml]      (default: fastearth.nml)
   !!
   !! The run config supplies the &ctl run-control group and any &fe3d physics
   !! overrides; its &fe3d is overlaid on the complete physics defaults in
   !! input/fastearth3d_defaults.nml (loaded automatically — see DEFAULTS_FILE).
   use fe_drive,   only: fastearth_run
   use fe_control, only: DEFAULTS_FILE
   implicit none
   character(len=512) :: cfg

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "fastearth.nml"
   end if

   call fastearth_run(trim(cfg), DEFAULTS_FILE)
end program fastearth

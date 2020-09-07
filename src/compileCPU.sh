ifort -fast -fpp -DDIPSLIP -ofd3d_pt_CPU_TSN dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversion.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90
#ifort -fast -fpp -ofd3d_pt_CPU_SS_TSN dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversion.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90

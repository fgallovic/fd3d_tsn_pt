source load_nv
#ifort -fast -fpp -ofd3d_pt_CPU_TSN_SS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90
nvfortran -Mpreprocess -Mbackslash -fast -ofd3d_pt_CPU_TSN_SS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90

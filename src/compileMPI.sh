source load_nv2
mpif90 -gpu:ccall -acc:gpu -fast -Mpreprocess -Mbackslash -DMPI -DGPUMPI -DDIPSLIP -ofd3d_pt_MPI_TSN_DS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90
mpif90 -gpu:ccall -acc:gpu -fast -Mpreprocess -Mbackslash -DMPI -DGPUMPI -ofd3d_pt_MPI_TSN_SS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90

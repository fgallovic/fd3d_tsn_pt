source load_nv
#nvfortran -gpu:ccall -acc:gpu -fast -Mpreprocess -DDIPSLIP -acc -Mbackslash -ofd3d_pt_GPU_TSN_DS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90
nvfortran -gpu:ccall,cc89 -acc:gpu -fast -Mpreprocess -Mbackslash -ofd3d_pt_GPU_TSN_SS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionSW.f90 filters.for mod_pt.f90 qdyn.f90 randomdynmod.f90 PGAmisf.f90

#FVW:
#nvfortran -lfftw3 -gpu:ccall -acc:gpu -fast -Mpreprocess -DFVW -r8 -Mbackslash -ofd3d_pt_GPU_FVW_SS dynamicsolver.f90 fd3d_deriv.f90 fd3d_init.f90 fd3d_theo.f90 waveforms.f90 inversionRS.f90 filters.for mod_pt.f90 qdyn.f90 PGAmisf.f90 mFFT_FFTW.f90

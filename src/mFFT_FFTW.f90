! Forward (negative exponent) and backward (positive exponent) FFT by FFTW3
!
! to compile and run...
! (GCC 1 thread)    gfortran -O3 mFFT_FFTW.f90 -lfftw3; time ./a.out
! (GCC & OpenMP)    gfortran -O3 -fopenmp mFFT_FFTW.f90 -lfftw3_omp -lfftw3; time OMP_NUM_THREADS=4 ./a.out
! (GCC & MKL)       gfortran -O3 -fopenmp mFFT_FFTW.f90 -L/opt/intel/oneapi/mkl/2023.0.0/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
!                   time LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2023.0.0/lib/intel64 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 ./a.out
! (Intel 1 thread)  ifort -O3 -heap-arrays mFFT_FFTW.f90 -lfftw3; time ./a.out
! (Intel OMP FFTW)  ifort -O3 -qopenmp -heap-arrays mFFT_FFTW.f90 -lfftw3_omp -lfftw3; time OMP_NUM_THREADS=4 ./a.out
! (Intel OMP MKL)   ifort -O3 -qopenmp -heap-arrays mFFT_FFTW.f90 -qmkl
!                   time OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 ./a.out
! (Nvidia 1 thread) nvfortran -O3 mFFT_FFTW.f90 -lfftw3; time ./a.out
! (Nvidia OMP FFTW) nvfortran -O3 -mp mFFT_FFTW.f90 -lfftw3_omp -lfftw3; time OMP_NUM_THREADS=4 ./a.out
! (Nvidia OMP MKL)  nvfortran -O3 -mp mFFT_FFTW.f90 -L/opt/intel/oneapi/mkl/2023.0.0/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
!                   time LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2023.0.0/lib/intel64 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 ./a.out
!
! a source code for running a single thread
! use mFFT_FFTW
! call init_fourn(arr,nn,ndim)
! call fourn(arr,nn,ndim,-1); call fourn(arr,nn,ndim,1)
! call destroy_fourn

MODULE FFTW3
USE,INTRINSIC :: iso_c_binding
INCLUDE '/usr/include/fftw3.f03'
END MODULE FFTW3

! ----------------------------------------------------------------------

MODULE mFFT_FFTW
IMPLICIT NONE
INTEGER(8),SAVE,PRIVATE :: planN,planP
!$OMP THREADPRIVATE (planN,planP)

CONTAINS

SUBROUTINE init_fourn(data,nn,ndim)
use FFTW3
INTEGER ndim,nn(ndim)
COMPLEX(8) data(product(nn))
select case (ndim)
case (1); call dfftw_plan_dft_1d(planN,nn(1),data,data,FFTW_FORWARD,FFTW_ESTIMATE)
case (2); call dfftw_plan_dft_2d(planN,nn(1),nn(2),data,data,FFTW_FORWARD,FFTW_ESTIMATE)
case (3); call dfftw_plan_dft_3d(planN,nn(1),nn(2),nn(3),data,data,FFTW_FORWARD,FFTW_ESTIMATE)
case default; call dfftw_plan_dft(planN,nn,data,data,FFTW_FORWARD,FFTW_ESTIMATE)
end select
select case (ndim)
case (1); call dfftw_plan_dft_1d(planP,nn(1),data,data,FFTW_BACKWARD,FFTW_ESTIMATE)
case (2); call dfftw_plan_dft_2d(planP,nn(1),nn(2),data,data,FFTW_BACKWARD,FFTW_ESTIMATE)
case (3); call dfftw_plan_dft_3d(planP,nn(1),nn(2),nn(3),data,data,FFTW_BACKWARD,FFTW_ESTIMATE)
case default; call dfftw_plan_dft(planP,nn,data,data,FFTW_BACKWARD,FFTW_ESTIMATE)
end select
END SUBROUTINE

SUBROUTINE fourn(data,nn,ndim,isign)
use FFTW3
INTEGER isign,ndim,nn(ndim)
COMPLEX(8) data(product(nn))
! select case (ndim)
! case (1); call dfftw_plan_dft_1d(plan,nn(1),data,data,isign,FFTW_ESTIMATE)
! case (2); call dfftw_plan_dft_2d(plan,nn(1),nn(2),data,data,isign,FFTW_ESTIMATE)
! case (3); call dfftw_plan_dft_3d(plan,nn(1),nn(2),nn(3),data,data,isign,FFTW_ESTIMATE)
! case default; call dfftw_plan_dft(plan,nn,data,data,isign,FFTW_ESTIMATE)
! end select
if (isign==-1) then
call dfftw_execute_dft(planN,data,data)
else
call dfftw_execute_dft(planP,data,data)
endif
! call dfftw_destroy_plan(plan)
END SUBROUTINE

SUBROUTINE destroy_fourn
use FFTW3
call dfftw_destroy_plan(planN)
call dfftw_destroy_plan(planP)
END SUBROUTINE

SUBROUTINE init_FFTW
END SUBROUTINE

SUBROUTINE cleanup_FFTW
END SUBROUTINE

END MODULE mFFT_FFTW

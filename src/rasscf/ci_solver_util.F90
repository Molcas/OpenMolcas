!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Giovanni Li Manni                                *
!               2020, Oskar Weser                                      *
!***********************************************************************

#include "macros.fh"

module CI_solver_util

use Para_Info, only: MyRank
use linalg_mod, only: verify_
use rasscf_global, only: nAc, nAcPar, nAcpr2, nroots
use general_data, only: JobIPH
#ifdef _MOLCAS_MPP_
use mpi
use Para_Info, only: Is_Real_Par
use Definitions, only: MPIInt
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
private

public :: wait_and_read, RDM_to_runfile, rdm_from_runfile, cleanMat, triangular_number, inv_triang_number, write_RDM

#ifdef _MOLCAS_MPP_
#include "global.fh"
integer(kind=MPIInt) :: error
integer(kind=MPIInt), parameter :: ROOT = 0_MPIInt
#include "mpi_interfaces.fh"
#endif

interface
  function isfreeunit(iseed)
    import :: iwp
    integer(kind=iwp) :: isfreeunit
    integer(kind=iwp), intent(in) :: iseed
  end function
end interface

#include "intent.fh"

contains

subroutine wait_and_read(filename,energy)

  character(len=*), intent(in) :: filename
  real(kind=wp), intent(out) :: energy(nroots)
  logical(kind=iwp) :: newcycle_found
  integer(kind=iwp) :: i, LuNewC

  newcycle_found = .false.
  do while (.not. newcycle_found)
    call sleepf(1)
    if (myrank == 0) call f_Inquire(trim(filename),newcycle_found)
#   ifdef _MOLCAS_MPP_
    if (is_real_par()) call MPI_Bcast(newcycle_found,1_MPIInt,MPI_LOGICAL,ROOT,MPI_COMM_WORLD,error)
#   endif
  end do
  if (myrank == 0) then
    write(u6,*) 'NEWCYCLE file found. Proceeding with SuperCI'
    LuNewC = isFreeUnit(12)
    call molcas_open(LuNewC,'NEWCYCLE')
    read(LuNewC,*) (energy(i),i=1,nroots)
    close(LuNewC,status='delete')
    write(u6,*) 'I read the following energies:',energy
  end if
# ifdef _MOLCAS_MPP_
  if (is_real_par()) call MPI_Bcast(energy,1_MPIInt,MPI_REAL8,ROOT,MPI_COMM_WORLD,error)
# endif

end subroutine wait_and_read

!>  @brief
!>  Put RDMs into runfile.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @param[out] DMAT Average 1 body density matrix
!>  @param[out] D1S_MO Average spin 1-dens matrix
!>  @param[out] PSMAT Average symm. 2-dens matrix
!>  @param[out] PAMAT Average antisymm. 2-dens matrix
!>  @param[in,out] jDisk
subroutine RDM_to_runfile(DMAT,D1S_MO,PSMAT,PAMAT,jDisk)

  ! _IN_ is not a semantic IN, since DDAFILE is both a read and
  ! write routine. Redefinition to suppress compiler warning.
  real(kind=wp), intent(_IN_) :: DMAT(nAcpar), D1S_MO(nAcPar), PSMAT(nAcpr2), PAMAT(nAcpr2)
  integer(kind=iwp), intent(inout), optional :: jDisk

  ! Put it on the RUNFILE
  call Put_dArray('D1mo',DMAT,NACPAR)
  call Put_dArray('P2mo',PSMAT,NACPR2)
  ! Save density matrices on disk
  ! DDAFILE calls BDAFile, iOpt option code
  ! 1 = synchronous write
  ! 2 = synchronous read
  ! BUF = array carrying data
  ! lBUF = length of array carrying data
  ! jDisk = memory address (automatically incremented upon
  !         repeated DDAFILE call)
  call DDafile(JOBIPH,1,DMAT,NACPAR,jDisk)
  call DDafile(JOBIPH,1,D1S_MO,NACPAR,jDisk)
  call DDafile(JOBIPH,1,PSMAT,NACPR2,jDisk)
  call DDafile(JOBIPH,1,PAMAT,NACPR2,jDisk)

end subroutine RDM_to_runfile

subroutine rdm_from_runfile(dmat,d1s_mo,psmat,pamat,jdisk)

  real(kind=wp), intent(out) :: dmat(nacpar), d1s_mo(nacpar), psmat(nacpr2), pamat(nacpr2)
  integer(kind=iwp), intent(inout), optional :: jdisk

  call ddafile(jobiph,2,dmat,nacpar,jdisk)
  call ddafile(jobiph,2,d1s_mo,nacpar,jdisk)
  call ddafile(jobiph,2,psmat,nacpr2,jdisk)
  call ddafile(jobiph,2,pamat,nacpr2,jdisk)

end subroutine rdm_from_runfile

subroutine CleanMat(MAT)
  !************ by G. Li Manni Stuttgart April 2016 *************
  !
  ! MAT: One-body density matrix in MO basis as passed by QMC calculation.

  ! It could well be an average matrix in SA calculation.
  !
  ! It has following shape:
  !        11
  !        12 22
  !        ** ** 33
  !        ** ** ** 44
  !        ** ** ** 45 55
  !        ** ** ** 46 56 66
  !        ** ** ** 47 57 67 77
  !        ** ** ** ** ** ** ** 88
  !        ** ** ** ** ** ** ** 89  99
  !        ** ** ** ** ** ** ** 810 910 1010
  !        """""""""""""""""""""""""""""""""""
  ! mimicking a system with (2 0 0 1 4 3 0 0)  active orbitals (blocked by Irreps)

  ! DMAT will be destroyed and replaced with a positive semi-definite one.
  ! N-representability will be preserved.

  use Constants, only: Zero, One, Two

  real(kind=wp), intent(inout) :: MAT(NacPar)
  integer(kind=iwp) :: i, j
  real(kind=wp) :: trace
  logical(kind=iwp) :: cleanup_required
  real(kind=wp), allocatable :: EVC(:), MAT_copy(:), Tmp(:), Tmp2(:)

  if (nacpar < 1) then
    write(u6,*) 'matrix size < 1.'
    Go To 10
  end if

  call mma_allocate(MAT_copy,NacPar)
  MAT_copy(:) = MAT(:)

  ! Allocate memory for eigenvectors and new DMAT
  call mma_allocate(EVC,NAC**2)
  ! Initialize eigenvectors
  ! set eigenvector array to identity for this version of JACOB
  call unitmat(EVC,NAC)

  ! Step 1: Diagonalize MAT. Eigenvalues are stored in diagonal of MAT
  trace = Zero
  do i=1,nac
    trace = trace+mat(i*(i+1)/2)
  end do
  call JACOB(MAT_copy,EVC,NAC,NAC)

# ifdef _DEBUGPRINT_
  write(u6,*) 'eigenvalues: '
  do i=1,nac
    write(u6,*) MAT_copy(I*(I+1)/2)
  end do
  write(u6,*) 'eigenvectors: '
  do i=1,nac
    write(u6,*) (EVC(i*NAC+j),j=0,NAC)
  end do
# endif
  ! Set to zero negative eigenvalue and to TWO values larger than 2.0.
  cleanup_required = .false.
  do j=1,nac
    if (MAT_copy(j*(j+1)/2) > Two) then
      MAT_copy(j*(j+1)/2) = Two
      cleanup_required = .true.
    end if
    if (MAT_copy(j*(j+1)/2) < 1.0e-12_wp) then
      MAT_copy(j*(j+1)/2) = Zero
      cleanup_required = .true.
    end if
  end do

  if (cleanup_required) then
    trace = Zero
    do i=1,nac
      trace = trace+MAT_copy(I*(I+1)/2)
    end do
    write(u6,*) 'trace after removing negative eigenvalues =',trace
    ! Combine pieced to form the output MAT
    ! blas routine for square*triangular operation
    call mma_allocate(Tmp,nac**2)
    call mma_allocate(Tmp2,nac**2)
    call dCopy_(nac**2,[Zero],0,Tmp,1)
    call dCopy_(nac**2,[Zero],0,Tmp2,1)
    do i=1,nac
      do j=1,nac
        Tmp(j+(i-1)*nac) = EVC(j+(i-1)*NAC)*MAT_copy(I*(I+1)/2)
      end do
    end do
    call DGEMM_('N','T',nac,nac,nac,One,Tmp,nac,EVC,nac,Zero,Tmp2,nac)
    ! Copy back to MAT
    do i=1,nac
      do j=1,i
        MAT(j+(i-1)*i/2) = Tmp2(j+(i-1)*nac)
      end do
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'trace after recombination:'
    trace = Zero
    do i=1,nac
      trace = trace+MAT(i*(i+1)/2)
    end do
#   endif
    call mma_deallocate(tmp)
    call mma_deallocate(tmp2)
  end if
  call mma_deallocate(MAT_copy)
  call mma_deallocate(EVC)
!***************** Exit ****************
10 continue

  return

end subroutine cleanMat

elemental function triangular_number(n)

  integer(kind=iwp) :: triangular_number
  integer(kind=iwp), intent(in) :: n

  triangular_number = n*(n+1)/2

end function

!> This is the inverse function of triangular_number
elemental function inv_triang_number(n) result(res)

  use Constants, only: Half, Quart

  integer(kind=iwp) :: res
  integer(kind=iwp), intent(in) :: n

  res = nint(sqrt(Quart+real(2*n,kind=wp))-Half)

end function

subroutine write_RDM(RDM,i_unit)

  real(kind=wp), intent(in) :: RDM(:)
  integer(kind=iwp), intent(in) :: i_unit
  integer(kind=iwp) :: curr_line, i, io_err, j, n_lines

  if (myrank == 0) then
    n_lines = inv_triang_number(size(RDM))

    i = 1
    do curr_line=1,n_lines
      do j=i,i+curr_line-1
        write(i_unit,'(ES25.15)',advance='no',iostat=io_err) RDM(j)
        call verify_(io_err == 0,'Error on writing RDM.')
      end do
      write(i_unit,*)
      i = i+curr_line
    end do
  end if

end subroutine write_RDM

end module CI_solver_util

************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2019, Giovanni Li Manni                                *
*               2020, Oskar Weser                                      *
************************************************************************
      module CI_solver_util
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_data, only: lRoots, nRoots, iAdr15,
     &                       iRoot, Weight, nAc, nAcPar, nAcpr2
      use general_data, only: JobIPH
      implicit none
      private
      public :: wait_and_read, abort_, assert_, RDM_to_runfile, dp,
     &      cleanMat
      integer, parameter :: dp = kind(1.0d0)
#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
      integer*4 :: error
      integer*4, parameter :: one4=1, root4=0
#endif

      interface
        integer function isfreeunit(iseed)
          integer, intent(in) :: iseed
        end function
      end interface

      contains


      subroutine wait_and_read(filename, energy)
        character(*), intent(in) :: filename
        real*8, intent(out) :: energy
        logical :: newcycle_found
        integer :: LuNewC
        newcycle_found = .false.
        do while(.not. newcycle_found)
          call sleep(1)
          if (myrank == 0) call f_Inquire(trim(filename),newcycle_found)
#ifdef _MOLCAS_MPP_
          if (is_real_par()) then
            call MPI_Bcast(newcycle_found, one4, MPI_LOGICAL,
     &                     root4, MPI_COMM_WORLD, error)
          end if
#endif
        end do
        if (myrank == 0) then
          write(6, *) 'NEWCYCLE file found. Proceding with SuperCI'
          LuNewC = isFreeUnit(12)
          call molcas_open(LuNewC, 'NEWCYCLE')
            read(LuNewC,*) energy
          close(LuNewC, status='delete')
          write(6, *) 'I read the following energy:', energy
        end if
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          call MPI_Bcast(energy, one4, MPI_REAL8,
     &                   root4, MPI_COMM_WORLD, error)
        end if
#endif
      end subroutine wait_and_read

      subroutine abort_(message)
        character(*), intent(in) :: message
        call WarningMessage(2, message)
        call QTrace()
        call Abend()
      end subroutine

      subroutine assert_(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) call abort_(message)
      end subroutine


!>  @brief
!>    State Average RDMs and put into runfile.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] DSPN Average spin 1-dens matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
      subroutine RDM_to_runfile(DMAT, D1S_MO, PSMAT, PAMAT)
        real*8, intent(in) :: DMAT(nAcpar), D1S_MO(nAcPar),
     &                        PSMAT(nAcpr2), PAMAT(nAcpr2)
        integer :: iDisk, jDisk

! Put it on the RUNFILE
        call Put_D1MO(DMAT,NACPAR)
        call Put_P2MO(PSMAT,NACPR2)
! Save density matrices on disk
        iDisk = IADR15(4)
        jDisk = IADR15(3)
        call DDafile(JOBIPH, 1, DMAT, NACPAR, jDisk)
        call DDafile(JOBIPH, 1, D1S_MO, NACPAR, jDisk)
        call DDafile(JOBIPH, 1, PSMAT, NACPR2, jDisk)
        call DDafile(JOBIPH, 1, PAMAT, NACPR2, jDisk)
      end subroutine RDM_to_runfile


      Subroutine CleanMat(MAT)
************* by G. Li Manni Stuttgart April 2016 *************
*
* MAT: One-body density matrix in MO basis as passed by QMC calculation.

* It could well be an average matrix in SA calculation.
*
* It has following shape:
*        11
*        12 22
*        ** ** 33
*        ** ** ** 44
*        ** ** ** 45 55
*        ** ** ** 46 56 66
*        ** ** ** 47 57 67 77
*        ** ** ** ** ** ** ** 88
*        ** ** ** ** ** ** ** 89  99
*        ** ** ** ** ** ** ** 810 910 1010
*        """""""""""""""""""""""""""""""""""
* mimicking a system with (2 0 0 1 4 3 0 0)  actice orbitals (blocked by Irreps)

*           DMAT will be destroyed and replaced with a positive semi-definite one.
*           N-representability will be preserved.

      real*8, intent(inout) :: MAT(NacPar)
      real*8, allocatable :: EVC(:), Tmp(:), Tmp2(:), MAT_copy(:)
      integer :: rc, i, j
      real*8 :: trace
      character(12), parameter :: routine = 'CleanMat'
      logical :: cleanup_required

      Call qEnter(routine)

      rc = 0
      If (nacpar .lt. 1) then
        rc= -1
        write(6,*) 'matrix size < 1.'
        Go To 10
      end if

      call mma_allocate(MAT_copy, NacPar)
      MAT_copy(:) = MAT(:)

* Allocate memory for eigenvectors and new DMAT
      call mma_allocate(EVC, NAC**2)
* Initialize eigenvectors
      Call dCopy_(NAC**2, [0.0d0], 0, EVC, 1)
* set eigenvector array to identity for this version of JACOB
      Call dCopy_(NAC, [1.0d0], 0, EVC, NAC + 1)

* Step 1: Diagonalize MAT. Eigenvalues are stored in diagonal of MAT
      trace = 0.0d0
      do i = 1, nac
         trace = trace + mat(i * (i + 1) / 2)
      end do
      CALL JACOB(MAT_copy, EVC, NAC, NAC)

#ifdef _DEBUG_
      write(6,*) 'eigenvalues: '
      do i=1,nac
         write(6,*) MAT_copy(I*(I+1)/2)
      end do
      write(6,*) 'eigenvectors: '
      do i=1, nac
        write(6,*) (EVC(i * NAC + j), j = 0, NAC)
      end do
#endif
* Set to zero negative eigenvalue and to TWO values larger than 2.0d0.
      cleanup_required = .false.
      do j = 1, nac
        if (MAT_copy(j * (j + 1) / 2) > 2.0d0) then
          MAT_copy(j * (j + 1) / 2) = 2.0d0
          cleanup_required = .true.
        end if
        if (MAT_copy(j * (j + 1) / 2) < 1.0d-12) then
          MAT_copy(j * (j + 1) / 2) = 0.0d0
          cleanup_required = .true.
        end if
      end do

      if (cleanup_required) then
        trace = 0.0d0
        do i = 1, nac
          trace = trace + MAT_copy(I * (I + 1) / 2)
        end do
        write(6,*) 'trace after removing negative eigenvalues =', trace
* Combine pieced to form the output MAT
* blas routine for square*triangular operation
        call mma_allocate(Tmp, nac**2)
        call mma_allocate(Tmp2, nac**2)
        Call dCopy_(nac**2, [0.0d0], 0, Tmp, 1)
        Call dCopy_(nac**2, [0.0d0], 0, Tmp2, 1)
        do i = 1, nac
          do j = 1, nac
            Tmp(j + (i - 1) * nac) =
     &          EVC(j + (i - 1) * NAC) * MAT_copy(I * (I + 1) / 2)
          end do
        end do
        Call DGEMM_('N','T',nac,nac,nac,
     &              1.0d0, Tmp, nac, EVC, nac,
     &              0.0d0, Tmp2, nac)
* Copy back to MAT
        do i = 1, nac
          do j = 1, i
            MAT(j + (i - 1) * i / 2) = Tmp2(j + (i - 1) * nac)
          end do
        end do
#ifdef _DEBUG_
        write(6,*) 'trace after recombination:'
        trace = 0.0d0
        do i = 1, nac
           trace = trace + MAT(i * (i + 1) / 2)
        end do
#endif
        call mma_deallocate(tmp)
        call mma_deallocate(tmp2)
      end if
      call mma_deallocate(MAT_copy)
      call mma_deallocate(EVC)
****************** Exit ****************
10    Continue
      Call qExit(routine)
      return
      end subroutine cleanMat



      end module CI_solver_util

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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

module Vec_Backup

use caspt2_global, only: LUSBT
use caspt2_module, only: nASup, nInDep, nISup, nSym
use EQSOLV, only: IDBMAT
use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
use GA_Wrapper, only: DBL_MB
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
private

public :: BackupVecs

contains

!-----------------------------------------------------------------------

subroutine BackupVecs(iMode,iVec,VecKeep)

  integer(kind=iwp), intent(in) :: iMode, iVec
  real(kind=wp), allocatable, intent(inout) :: VecKeep(:)
  integer(kind=iwp) :: iCase, iHi, iLo, iPass, iPos, iSym, jD, jHi, jLo, lg_V, mV, NAS, NIN, NIS, nKeep
  real(kind=wp), allocatable :: BD(:)
  real(kind=wp), parameter :: ThrB = 1.0e+3_wp

  ! The transformation of a vector, SR -> contrav C -> SR, is not numerically exact,
  ! and the difference matters for the functions with a large eigenvalue of B.
  ! iMode = 1 stores these problematic vectors of iVec in VecKeep, and iMode = 2 restores them.

  ! iPass = 0, 1, 2: count, store, and restore
  nKeep = 0
  do iPass=2*(iMode-1),iMode
    if (iPass == 1) call mma_allocate(VecKeep,nKeep,Label='VecKeep')
    iPos = 0
    do iCase=1,11
      do iSym=1,nSym
        NIN = nInDep(iSym,iCase)
        NIS = nISup(iSym,iCase)
        NAS = nASup(iSym,iCase)
        if (NIN*NIS == 0) cycle
        call mma_allocate(BD,NAS,Label='BD')
        jD = IDBMAT(iSym,iCase)
        call DDAFILE(LUSBT,2,BD,NAS,jD)
        if (any(BD(1:NIN) > ThrB)) then
          if (iPass == 0) then
            call RHS_DISTRIBUTION(NIN,NIS,iLo,iHi,jLo,jHi)
            nKeep = nKeep+count(BD(1:NIN) > ThrB)*(jHi-jLo+1)
          else
            call RHS_ALLO(NIN,NIS,lg_V)
            call RHS_READ(NIN,NIS,lg_V,iCase,iSym,iVec)
            call RHS_ACCESS(NIN,NIS,lg_V,iLo,iHi,jLo,jHi,mV)
#           ifdef _MOLCAS_MPP_
            if (Is_Real_Par()) then
              call CopyVecs(jHi-jLo+1,DBL_MB(mV))
            else
#           endif
              call CopyVecs(jHi-jLo+1,GA_Arrays(mV)%A)
#           ifdef _MOLCAS_MPP_
            end if
#           endif
            if (iPass == 2) then
              call RHS_RELEASE_UPDATE(lg_V,iLo,iHi,jLo,jHi)
              call RHS_SAVE(NIN,NIS,lg_V,iCase,iSym,iVec)
            else
              call RHS_RELEASE(lg_V,iLo,iHi,jLo,jHi)
            end if
            call RHS_FREE(lg_V)
          end if
        end if
        call mma_deallocate(BD)
      end do
    end do
  end do
  if (iMode == 2) call mma_deallocate(VecKeep)

  contains

  subroutine CopyVecs(NCOL,W)

    integer(kind=iwp), intent(in) :: NCOL
    real(kind=wp), intent(inout) :: W(NIN,NCOL)
    integer(kind=iwp) :: i

    ! Copy the rows of W with BD > ThrB to VecKeep (iPass = 1) or back (iPass = 2)

    do i=1,NIN
      if (BD(i) <= ThrB) cycle
      if (iPass == 1) then
        VecKeep(iPos+1:iPos+NCOL) = W(i,:)
      else
        W(i,:) = VecKeep(iPos+1:iPos+NCOL)
      end if
      iPos = iPos+NCOL
    end do

  end subroutine CopyVecs

end subroutine BackupVecs

end module Vec_Backup

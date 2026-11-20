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

use caspt2_global, only: do_grad, LUSBT, nStpGrd
use caspt2_module, only: NASUP, NINDEP, NISUP, nSym
use EQSOLV, only: IDBMAT, IDTMAT
use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
use GA_Wrapper, only: DBL_MB, GA_Destroy, GA_NodeId
use Para_Info, only: Is_Real_Par
use Constants, only: Zero
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
  integer(kind=iwp) :: iCase, iHi, iLo, iOff, iPass, iPos, iSym, jHi, jLo, lg_V, mV, NIN, NIS, nKeep
  logical(kind=iwp), allocatable :: Keep(:)

  ! The transformation of a vector, SR --> contrav C --> SR, is not numerically exact.
  ! iMode = 1 stores the rows of iVec that are sensitive to it (see SelectRows) in VecKeep, and iMode = 2 restores them.

  if (.not. (do_grad .or. (nStpGrd == 2))) return

  call mma_allocate(Keep,sum(NINDEP(1:nSym,1:11)),Label='Keep')
  call SelectRows(sum(NINDEP(1:nSym,1:11)),Keep)

  ! iPass = 0, 1, 2: count, store, and restore
  nKeep = 0
  do iPass=2*(iMode-1),iMode
    if (iPass == 1) call mma_allocate(VecKeep,nKeep,Label='VecKeep')
    iPos = 0
    iOff = 0
    do iCase=1,11
      do iSym=1,nSym
        NIN = NINDEP(iSym,iCase)
        NIS = NISUP(iSym,iCase)
        iOff = iOff+NIN
        if (.not. any(Keep(iOff-NIN+1:iOff))) cycle
        if (iPass == 0) then
          call RHS_DISTRIBUTION(NIN,NIS,iLo,iHi,jLo,jHi)
          nKeep = nKeep+count(Keep(iOff-NIN+1:iOff))*(jHi-jLo+1)
        else
          call RHS_ALLO(NIN,NIS,lg_V)
          call RHS_READ(NIN,NIS,lg_V,iCase,iSym,iVec)
          call RHS_ACCESS(NIN,NIS,lg_V,iLo,iHi,jLo,jHi,mV)
#         ifdef _MOLCAS_MPP_
          if (Is_Real_Par()) then
            call CopyVecs(jHi-jLo+1,DBL_MB(mV))
          else
#         endif
            call CopyVecs(jHi-jLo+1,GA_Arrays(mV)%A)
#         ifdef _MOLCAS_MPP_
          end if
#         endif
          if (iPass == 2) then
            call RHS_RELEASE_UPDATE(lg_V,iLo,iHi,jLo,jHi)
            call RHS_SAVE(NIN,NIS,lg_V,iCase,iSym,iVec)
          else
            call RHS_RELEASE(lg_V,iLo,iHi,jLo,jHi)
          end if
          call RHS_FREE(lg_V)
        end if
      end do
    end do
  end do
  call mma_deallocate(Keep)
  if (iMode == 2) call mma_deallocate(VecKeep)

  contains

  subroutine CopyVecs(NCOL,W)

    integer(kind=iwp), intent(in) :: NCOL
    real(kind=wp), intent(inout) :: W(NIN,NCOL)
    integer(kind=iwp) :: i

    ! Copy the selected rows of W to VecKeep (iPass = 1) or back (iPass = 2)

    do i=1,NIN
      if (.not. Keep(iOff-NIN+i)) cycle
      if (iPass == 1) then
        VecKeep(iPos+1:iPos+NCOL) = W(i,:)
      else
        W(i,:) = VecKeep(iPos+1:iPos+NCOL)
      end if
      iPos = iPos+NCOL
    end do

  end subroutine CopyVecs

end subroutine BackupVecs

!-----------------------------------------------------------------------

subroutine SelectRows(nDimKeep,Keep)

  integer(kind=iwp), intent(in) :: nDimKeep
  logical(kind=iwp), intent(out) :: Keep(nDimKeep)
  integer(kind=iwp) :: iCase, iOff, iSym, j, jD, NAS, NIN
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: i, iHi, iLo, jHi, jLo, LDT, lg_T, mT, myRank
  logical(kind=iwp) :: bStat
# endif
  real(kind=wp), allocatable :: BD(:), TN(:), TRANS(:,:)
  real(kind=wp), parameter :: ThrB = 1.0e+3_wp, ThrT = 1.0e+6_wp
  real(kind=wp), external :: DDot_

  ! Select the rows (functions in the SR basis) for which the difference in the transformation matters.
  ! Large eigenvalue of B (BD) amplifies the difference in the residual.
  ! Large norm of the column of the transformation matrix (TN) means a small eigenvalue of S,
  ! for which the transformation matrix is less accurately orthonormal.

  Keep(:) = .false.
  iOff = 0
  do iCase=1,11
    do iSym=1,nSym
      NIN = NINDEP(iSym,iCase)
      NAS = NASUP(iSym,iCase)
      iOff = iOff+NIN
      if (NIN*NISUP(iSym,iCase) == 0) cycle
      call mma_allocate(BD,NIN,Label='BD')
      call mma_allocate(TN,NIN,Label='TN')
      jD = IDBMAT(iSym,iCase)
      call DDAFILE(LUSBT,2,BD,NIN,jD)
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par() .and. ((iCase == 1) .or. (iCase == 4))) then
        TN(:) = Zero
        call GA_CREATE_STRIPED('H',NAS,NIN,'TMAT',lg_T)
        call PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
        myRank = GA_NodeId()
        call GA_Distribution(lg_T,myRank,iLo,iHi,jLo,jHi)
        if ((iLo > 0) .and. (jLo > 0)) then
          call GA_Access(lg_T,iLo,iHi,jLo,jHi,mT,LDT)
          do j=jLo,jHi
            do i=1,iHi-iLo+1
              TN(j) = TN(j)+DBL_MB(mT+i-1+LDT*(j-jLo))**2
            end do
          end do
          call GA_Release(lg_T,iLo,iHi,jLo,jHi)
        end if
        call GADGOP(TN,NIN,'+')
        bStat = GA_Destroy(lg_T)
#       include "macros.fh"
        unused_var(bStat)
      else
#     endif
        call mma_allocate(TRANS,NAS,NIN,Label='TRANS')
        jD = IDTMAT(iSym,iCase)
        call DDAFILE(LUSBT,2,TRANS,NAS*NIN,jD)
        do j=1,NIN
          TN(j) = DDot_(NAS,TRANS(:,j),1,TRANS(:,j),1)
        end do
        call mma_deallocate(TRANS)
#     ifdef _MOLCAS_MPP_
      end if
#     endif
      Keep(iOff-NIN+1:iOff) = (BD(:) > ThrB) .or. (TN(:) > ThrT)
      call mma_deallocate(BD)
      call mma_deallocate(TN)
    end do
  end do

end subroutine SelectRows

end module Vec_Backup

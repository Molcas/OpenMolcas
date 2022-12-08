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
! Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine GetRawPAOs(R,C,nBas,nOrb,nFro,nOrb2Loc,nSym,Normalize)
! Thomas Bondo Pedersen, December 2005.
! - revised january 2006 (Thomas Bondo Pedersen).
!
! Purpose: compute projected AOs spanning the primary space from the
!          formula
!
!          R = 1 - Do*S = D*S
!
!          where Do is the "AO density" matrix of the orthogonal
!          complement orbitals and D that of the orbital space to be
!          localised (primary space).
!          S is the AO overlap matrix (which is read from disk).
!
!          Which formula is used depends on the dimensions of the
!          two spaces (such that the most economical is computed).
!
!          If (Normalize): normalize each PAO (recommended).
!                          Note that if the norm of the PAO is
!                          smaller than 1.0e-6_wp, it will not be
!                          normalized (it is left unchanged).
!
!-------------------------------------------------------------
!-TODO/FIXME: it is in most cases faster to use Do*S=C*(C^T*S)
!-------------------------------------------------------------

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: R(*)
real(kind=wp), intent(in) :: C(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nOrb(nSym), nFro(nSym), nOrb2Loc(nSym)
logical(kind=iwp), intent(in) :: Normalize
integer(kind=iwp) :: i, iSym, lDo, lOff, mu, nB, nF, nO2L, nOrth, nRest
real(kind=wp) :: OvlpR
character(len=80) :: Txt
type(DSBA_Type) :: CB, Ovlp, RB
real(kind=wp), allocatable, target :: DoRt(:)
real(kind=wp), pointer :: DoR(:,:)
character(len=*), parameter :: SecNam = 'GetRawPAOs'
real(kind=wp), external :: ddot_

! Read the overlap matrix from disk.
! ----------------------------------

call Allocate_DT(Ovlp,nBas,nBas,nSym,label='Ovlp')
call GetOvlp_Localisation(Ovlp%A0,'Sqr',nBas,nSym)

call Allocate_DT(RB,nBas,nBas,nSym,label='RB',Ref=R)
call Allocate_DT(CB,nBas,nBas,nSym,label='RB',Ref=C)

! Compute R.
! ----------

lDo = nBas(1)**2
do iSym=2,nSym
  lDo = max(lDo,nBas(iSym)**2)
end do
call mma_allocate(DoRt,lDo,label='Do')

do iSym=1,nSym

  nB = nBas(iSym)
  DoR(1:nB,1:nB) => DoRt(1:nB**2)
  if (nB > 0) then

    nF = nFro(iSym)
    nO2L = nOrb2Loc(iSym)
    nRest = nOrb(iSym)-nF-nO2L
    nOrth = nF+nRest ! dim. of orthogonal complement

    if (nO2L < 1) then ! R = 0
      RB%SB(iSym)%A2(:,:) = Zero
    else if (nOrth < 0) then ! error
      call SysAbendMsg(SecNam,'Dimension of orthogonal complement space < 0',' ')
    else if (nOrth == 0) then ! R = 1
      call unitmat(RB%SB(iSym)%A2,nB)
    else if (nOrth < nO2L) then ! R = 1 - Do*S
      if (nRest > 0) then
        lOff = nF+nO2L+1
        call GetDens_Localisation(DoR,CB%SB(iSym)%A2(:,lOff:),nB,nRest)
      else
        DoR(:,:) = Zero
      end if
      if (nF > 0) then
        call GetDens_Localisation(RB%SB(iSym)%A2,CB%SB(iSym)%A2,nB,nF)
        DoR(:,:) = DoR(:,:)+RB%SB(iSym)%A2(:,:)
      end if
      call DGEMM_('N','N',nB,nB,nB,-One,DoR,nB,Ovlp%SB(iSym)%A2,nB,Zero,RB%SB(iSym)%A2,nB)
      do i=1,nB
        RB%SB(iSym)%A2(i,i) = RB%SB(iSym)%A2(i,i)+One
      end do
    else ! R = D*S
      lOff = nF+1
      call GetDens_Localisation(DoR,CB%SB(iSym)%A2(:,lOff:),nB,nO2L)
      call DGEMM_('N','N',nB,nB,nB,One,DoR,nB,Ovlp%SB(iSym)%A2,nB,Zero,RB%SB(iSym)%A2,nB)
    end if

  end if

end do

! If requested, normalize the PAOs.
! ---------------------------------

if (Normalize) then
  do iSym=1,nSym
    nB = nBas(iSym)
    DoR(1:nB,1:nB) => DoRt(1:nB**2)
    if (nB > 0) then
      call DGEMM_('N','N',nB,nB,nB,One,Ovlp%SB(iSym)%A2,nB,RB%SB(iSym)%A2,nB,Zero,DoR,nB)
      do mu=1,nB
        OvlpR = dDot_(nB,RB%SB(iSym)%A2(:,mu),1,DoR(:,mu),1)
        if (OvlpR > 1.0e-6_wp) then
          RB%SB(iSym)%A2(:,mu) = RB%SB(iSym)%A2(:,mu)/sqrt(OvlpR)
        else if (OvlpR < Zero) then
          write(Txt,'(A,1P,D15.5)') 'Overlap = ',OvlpR
          call SysAbendMsg(SecNam,'Negative raw PAO overlap!',Txt)
        end if
      end do
    end if
  end do
end if

nullify(DoR)
call mma_deallocate(DoRt)
call Deallocate_DT(Ovlp)
call Deallocate_DT(RB)
call Deallocate_DT(CB)

end subroutine GetRawPAOs

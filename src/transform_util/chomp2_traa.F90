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
! Copyright (C) 2004,2005, Giovanni Ghigo                              *
!***********************************************************************

subroutine ChoMP2_TraA(iSym,jSym,NumV,CMO,NCMO,lUCHFV,iStrtVec_AB,nFVec)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden                                    *
! Written:  October 2004                                               *
! Modified for Cholesky-MP2 May 2005                                   *
!***********************************************************************

use Cho_Tra, only: nBas, nFro, nIsh, nSsh, TCVX, TCVXist
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, jSym, NumV, NCMO, lUCHFV, iStrtVec_AB, nFVec
real(kind=wp), intent(in) :: CMO(NCMO)
integer(kind=iwp) :: i, iFBatch, iiVec, iStrt, iStrt0MO, iStrtVec_FAB, iVec, j, jStrt, jStrt0MO, jVec, Len_XAj, Len_XBi, Naj, Nbi, &
                     NFAB, NumFV
logical(kind=iwp) :: TCVC, TCVCt
real(kind=wp), allocatable :: XAj(:), XBi(:), FAB(:,:)

! Memory to allocate & Nr. of Cholesky vectors transformable
! A=Alpha(AO);  B=Beta(AO)
TCVC = .false.
TCVCt = .false.

NFAB = 0
Naj = 0  ! C
Nbi = 0  ! C"
Len_XAj = 0   ! C
Len_XBi = 0   ! C"
NFAB = nBas(iSym)*nBas(jSym)

! Allocate memory for Transformed Cholesky Vectors - TCVx
! TCV-C :
if (TCVXist(3,iSym,jSym)) then
  TCVC = .true.
  Len_XAj = nBas(iSym)*nIsh(jSym)
  Naj = nSsh(iSym)*nIsh(jSym)

  call mma_allocate(TCVX(3,iSym,jSym)%A,Naj,NumV,Label='TCVX')
end if
! TCV-Ct:
if (TCVXist(3,jSym,iSym)) then
  TCVCt = .true.
  Len_XBi = nBas(jSym)*nIsh(iSym)
  Nbi = nSsh(jSym)*nIsh(iSym)

  call mma_allocate(TCVX(3,jSym,iSym)%A,Nbi,NumV,Label='TCVX')
end if

iStrt = 1
do i=1,iSym-1
  iStrt = iStrt+nBas(i)*nBas(i)
end do

jStrt = 1
do j=1,jSym-1
  jStrt = jStrt+nBas(j)*nBas(j)
end do

! START LOOP iiVec
do iiVec=1,NumV,nFVec
  NumFV = max(nFVec,NumV-iiVec+1)
  iFBatch = (iiVec+nFVec-1)/nFVec

  ! Allocate memory & Load Full Cholesky Vectors - CHFV

  iStrtVec_FAB = iStrtVec_AB+nFVec*(iFBatch-1)

  call mma_allocate(FAB,nFAB,NumFV,Label='FAB')
  call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

  ! Start Loop jVec
  do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
    iVec = jVec-iiVec+1

    ! 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
    ! From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
    if (TCVC) then
      jStrt0MO = jStrt+nFro(jSym)*nBas(jSym)
      call mma_allocate(XAj,Len_XAj,Label='XAj')
      call ProdsA_2(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(jStrt0MO),nIsh(jSym),XAj)
    end if
    ! From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
    if (TCVCt) then
      iStrt0MO = iStrt+nFro(iSym)*nBas(iSym)
      call mma_allocate(XBi,Len_XBi,Label='XBi')
      call ProdsA_2t(FAB(:,iVec),nBas(iSym),nBas(jSym),CMO(iStrt0MO),nIsh(iSym),XBi)
    end if

    ! 2nd Half-Transformation  iAlpha(AO) -> p(MO)
    ! From XAj(Alpha,jMO) to aj(a,j)
    if (TCVC) then
      iStrt0MO = iStrt+(nFro(iSym)+nIsh(iSym))*nBas(iSym)
      call ProdsA_1(XAj,nBas(iSym),nIsh(jSym),CMO(iStrt0MO),nSsh(iSym),TCVX(3,iSym,jSym)%A(:,jVec))
    end if

    ! From XBi(Beta,jMO) to bi(b,i)
    if (TCVCt) then
      jStrt0MO = jStrt+(nFro(jSym)+nIsh(jSym))*nBas(jSym)
      call ProdsA_1(XBi,nBas(jSym),nIsh(iSym),CMO(jStrt0MO),nSsh(jSym),TCVX(3,jSym,iSym)%A(:,jVec))
    end if

    ! End of Transformations

    if (allocated(XAj)) call mma_deallocate(XAj)
    if (allocated(XBi)) call mma_deallocate(XBi)

  end do
  ! End Loop jVec

  call mma_deallocate(FAB)
end do
! END LOOP iiVec

return

end subroutine ChoMP2_TraA

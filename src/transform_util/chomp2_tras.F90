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

subroutine ChoMP2_TraS(iSym,jSym,NumV,CMO,NCMO,lUCHFV,iStrtVec_AB,nFVec)
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
integer(kind=iwp) :: i, iFBatch, iiVec, iStrt, iStrt0MO, iStrtVec_FAB, iVec, j, jStrt, jStrt0MO, jVec, Len_XAj, Naj, NFAB, NumFV
real(kind=wp), allocatable :: XAj(:), FAB(:,:)

! Memory to allocate & Nr. of Cholesky vectors transformable
! A=Alpha(AO);  B=Beta(AO)

if (.not. TCVXist(3,iSym,jSym)) return

Naj = 0
Len_XAj = 0
NFAB = nBas(iSym)*(nBas(jSym)+1)/2
Len_XAj = nBas(iSym)*nIsh(jSym)
Naj = nSsh(iSym)*nIsh(jSym)

! Allocate memory for Transformed Cholesky Vectors - TCVx
call mma_allocate(TCVX(3,iSym,jSym)%A,Naj,NumV,Label='TCVX')

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

  call mma_allocate(FAB,NFAB,NumFV,Label='FAB')
  call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

  ! Start Loop jVec
  do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
    iVec = jVec-iiVec+1

    ! 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
    jStrt0MO = jStrt+nFro(jSym)*nBas(jSym)

    ! From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
    call mma_allocate(XAj,Len_XAj,Label='XAj')
    call ProdsS_1(FAB(:,iVec),nBas(iSym),CMO(jStrt0MO),nIsh(jSym),XAj)

    ! 2nd Half-Transformation  iAlpha(AO) -> p(MO)
    iStrt0MO = iStrt+(nFro(iSym)+nIsh(iSym))*nBas(iSym)

    ! From XAj(Alpha,jMO) to aj(a,j)
    call ProdsA_1(XAj,nBas(iSym),nIsh(jSym),CMO(iStrt0MO),nSsh(iSym),TCVX(3,iSym,jSym)%A(:,jVec))

    ! End of Transformations

    call mma_deallocate(XAj)

  end do
  ! End Loop  jVec

  call mma_deallocate(FAB)
end do
! END LOOP iiVec

return

end subroutine ChoMP2_TraS

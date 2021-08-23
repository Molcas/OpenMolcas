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

subroutine ChoMP2_TwoEl(iBatch,numV,LUINTM,iAddrIAD2M,iSymI,iSymJ,iSymA,iSymB)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden                                    *
! Written:  October-November 2004                                      *
! Modified for Cholesky-MP2 May 2005                                   *
!***********************************************************************

use Cho_Tra, only: IAD2M, nOsh, nSsh, nSym, SubBlocks
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iBatch, numV, LUINTM, iSymI, iSymJ, iSymA, iSymB
integer(kind=iwp), intent(inout) :: iAddrIAD2M
integer(kind=iwp) :: iAddrIAD2Mij, iEndJ, iI, iIJAB, iJ, nA, nB, nN_AB, nN_Ex1, nN_Ex2, nN_IJ, nSymP
real(kind=wp), allocatable :: AddEx1(:), AddEx2(:), AddEx2t(:)

nSymP = (nSym**2+nSym)/2
call LenInt(iSymI,iSymJ,iSymA,iSymB,nN_IJ,nN_AB,nN_Ex1,nN_Ex2)

!**** START GENERATION of EXCHANGE-1 INTEGRALS *************************
if (nN_IJ*nN_Ex1 > 0) then
  iIJAB = ((iSymI**2-iSymI)/2+iSymJ-1)*nSymP+(iSymA**2-iSymA)/2+iSymB
  SubBlocks(3,3) = .true.
  if (iBatch == 1) then
    IAD2M(2,iIJAB) = iAddrIAD2M
  else
    iAddrIAD2M = IAD2M(2,iIJAB)
  end if
  ! Start Loop on i, j
  iAddrIAD2Mij = iAddrIAD2M
  do iI=1,nOsh(iSymI)
    if (iSymI == iSymJ) then
      iEndJ = iI
    else
      iEndJ = nOsh(iSymJ)
    end if
    do iJ=1,iEndJ
      call mma_allocate(AddEx1,nN_Ex1,Label='AddEx1')
      if (iBatch > 1) then
        ! Reload Int
        call dDaFile(LUINTM,2,AddEx1,nN_Ex1,iAddrIAD2Mij)
        iAddrIAD2Mij = iAddrIAD2Mij-nN_Ex1
      else
        AddEx1(:) = Zero
      end if
      call ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx1,nN_Ex1)
      call dDaFile(LUINTM,1,AddEx1,nN_Ex1,iAddrIAD2Mij)
      call mma_deallocate(AddEx1)
    end do
  end do
  ! End Loop on i, j
  iAddrIAD2M = iAddrIAD2Mij
end if
!**** END GENERATION of EXCHANGE-1 INTEGRALS ***************************

!**** START GENERATION of EXCHANGE-2 INTEGRALS *************************
if (nN_IJ*nN_Ex2 > 0) then
  iIJAB = ((iSymI**2-iSymI)/2+iSymJ-1)*nSymP+(iSymB**2-iSymB)/2+iSymA
  SubBlocks(3,3) = .true.
  if (iBatch == 1) then
    IAD2M(3,iIJAB) = iAddrIAD2M
  else
    iAddrIAD2M = IAD2M(3,iIJAB)
  end if
  ! Start Loop on i, j
  iAddrIAD2Mij = iAddrIAD2M
  do iI=1,nOsh(iSymI)
    if (iSymI == iSymJ) then
      iEndJ = iI
    else
      iEndJ = nOsh(iSymJ)
    end if
    do iJ=1,iEndJ
      nA = nSsh(iSymA)
      nB = nSsh(iSymB)
      call mma_allocate(AddEx2,nN_Ex2,Label='AddEx2')
      call mma_allocate(AddEx2t,nN_Ex2,Label='AddEx2t')
      if (iBatch > 1) then
        ! Reload Int
        call dDaFile(LUINTM,2,AddEx2,nN_Ex2,iAddrIAD2Mij)
        iAddrIAD2Mij = iAddrIAD2Mij-nN_Ex2
        call Trnsps(nA,nB,AddEx2,AddEx2t)
      else
        AddEx2t(:) = Zero
      end if
      call ChoMP2_GenE(iSymI,iSymJ,iSymA,iSymB,iI,iJ,numV,AddEx2t,nN_Ex2)
      call Trnsps(nB,nA,AddEx2t,AddEx2)
      call dDaFile(LUINTM,1,AddEx2,nN_Ex2,iAddrIAD2Mij)
      call mma_deallocate(AddEx2t)
      call mma_deallocate(AddEx2)
    end do
  end do
  ! End Loop on i, j
  iAddrIAD2M = iAddrIAD2Mij
end if
!**** END GENERATION of EXCHANGE-2 INTEGRALS ***************************

return

end subroutine ChoMP2_TwoEl

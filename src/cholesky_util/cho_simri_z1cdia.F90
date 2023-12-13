!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Cho_SimRI_Z1CDia(Diag,Thr,Indx)
!
! Purpose: Zero 1-center diagonals smaller then Thr.
!          Diag is the diagonal (1st reduced set).
!          On exit, Indx(i)=1 if diagonal i was zeroed, else
!          Indx(i)=0 (thus, Indx must have same dimension as Diag).

use Cholesky, only: iAtomShl, iiBstR, iiBstRSh, IPRINT, iSP2F, LuPri, nnBstR, nnBstRSh, nnShl
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: Diag(*)
real(kind=wp), intent(in) :: Thr
integer(kind=iwp), intent(_OUT_) :: Indx(*)
integer(kind=iwp) :: iAB, iAB1, iAB2, iShlA, iShlAB, iShlB, n
real(kind=wp) :: zmx
integer(kind=iwp), parameter :: Inf_SimRI = 0

Indx(1:nnBstR(1,1)) = 0

zmx = Zero
n = 0
do iShlAB=1,nnShl
  call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
  if (iAtomShl(iShlA) == iAtomShl(iShlB)) then
    iAB1 = iiBstR(1,1)+iiBstRSh(1,iShlAB,1)+1
    iAB2 = iAB1+nnBstRSh(1,iShlAB,1)-1
    do iAB=iAB1,iAB2
      if (Diag(iAB) < Thr) then
        zmx = max(zmx,Diag(iAB))
        n = n+1
        Indx(iAB) = 1
        Diag(iAB) = Zero
      end if
    end do
  end if
end do

if (iPrint > Inf_SimRI) then
  write(LuPri,'(/,A,I7,A,ES10.2,A)') 'Simulating RI:',n,' 1-center diagonals < ',Thr,' have been zeroed'
  if (n > 0) write(LuPri,'(A,ES15.7)') 'Largest zeroed diagonal: ',zmx
end if

end subroutine Cho_SimRI_Z1CDia

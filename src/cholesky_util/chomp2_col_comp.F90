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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_Col_Comp(Col,nDim,iCol,nCol,Vec,nVec,Buf,lBuf,Fac,irc)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: compute columns from a set of vectors.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim, nCol, iCol(nCol), nVec, lBuf
real(kind=wp), intent(inout) :: Col(nDim,nCol)
real(kind=wp), intent(in) :: Vec(nDim,nVec), Fac
real(kind=wp), intent(out) :: Buf(lBuf)
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: iBat, iC, iC1, NumBat, NumC, NumCol

irc = 0
if ((nDim < 1) .or. (nCol < 1)) return
if (nVec < 1) then
  if (Fac /= One) Col(:,:) = Fac*Col(:,:)
  return
end if

if ((nCol > 1) .and. (lBuf >= nVec)) then

  NumCol = min(lBuf/nVec,nCol)
  NumBat = (nCol-1)/NumCol+1

  do iBat=1,NumBat

    if (iBat == NumBat) then
      NumC = nCol-NumCol*(NumBat-1)
    else
      NumC = NumCol
    end if
    iC1 = NumCol*(iBat-1)+1

    if (lBuf < NumC*nVec) then
      irc = -1
      return
    end if

    call ChoMP2_Col_cp(Vec,nDim,nVec,Buf,NumC,iCol(iC1))
    call DGEMM_('N','T',nDim,NumC,nVec,One,Vec,nDim,Buf,NumC,Fac,Col(:,iC1),nDim)

  end do

else

  do iC=1,nCol
    call dGeMV_('N',nDim,nVec,One,Vec(:,:),nDim,Vec(iCol(iC),1),nDim,Fac,Col(:,iC),1)
  end do

end if

end subroutine ChoMP2_Col_Comp

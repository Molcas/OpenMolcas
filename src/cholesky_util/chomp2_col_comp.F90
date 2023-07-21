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

implicit real*8(a-h,o-z)
real*8 Col(nDim,nCol), Vec(nDim,nVec), Buf(lBuf)
integer iCol(nCol)

irc = 0
if ((nDim < 1) .or. (nCol < 1)) return
if (nVec < 1) then
  if (Fac /= 1.0d0) call dScal_(nDim*nCol,Fac,Col,1)
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
    call DGEMM_('N','T',nDim,NumC,nVec,1.0d0,Vec,nDim,Buf,NumC,Fac,Col(1,iC1),nDim)

  end do

else

  do iC=1,nCol
    call dGeMV_('N',nDim,nVec,1.0d0,Vec(1,1),nDim,Vec(iCol(iC),1),nDim,Fac,Col(1,iC),1)
  end do

end if

end subroutine ChoMP2_Col_Comp

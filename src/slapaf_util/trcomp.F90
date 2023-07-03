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

subroutine TRComp(TRVec,nTR,nX,Smmtrc)

implicit real*8(a-h,o-z)
real*8 TRVec(nTR,nX)
logical Smmtrc(nX)

if (nTR == 0) return
iDim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    iDim = iDim+1
    if (iDim /= iX) call dcopy_(nTR,TRVec(1,iX),1,TRVec(1,iDim),1)
  end if
end do

return

end subroutine TRComp

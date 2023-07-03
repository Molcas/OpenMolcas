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

subroutine TROrder(TRVec,nTR,nX)

implicit real*8(a-h,o-z)
real*8 TRVec(6*nX)

if (nTR == 6) return
do i=2,nX
  do iTR=1,nTR
    iFrom = (i-1)*6+iTR
    iTo = (i-1)*nTR+iTR
    TRVec(iTo) = TRVec(iFrom)
  end do
end do

return

end subroutine TROrder

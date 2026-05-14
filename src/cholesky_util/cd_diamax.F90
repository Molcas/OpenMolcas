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

subroutine CD_DiaMax(Diag,nDim,iPivot,iQual,nQual,DiaMin)
!
! Thomas Bondo Pedersen, October 2004.
!
! Purpose: find nQual largest elements in Diag and leave pointers to
!          them in iQual. Only elements larger than DiaMin are
!          returned (thus, nQual may be reduced here).

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(in) :: Diag(nDim), DiaMin
integer(kind=iwp), intent(inout) :: nQual
integer(kind=iwp), intent(out) :: iPivot(nDim), iQual(nQual)
integer(kind=iwp) :: i, iMax, iTmp, j

do i=1,nDim
  iPivot(i) = i
end do

do j=1,nQual
  do i=nDim,j+1,-1
    if (Diag(iPivot(i)) > Diag(iPivot(i-1))) then
      iTmp = iPivot(i-1)
      iPivot(i-1) = iPivot(i)
      iPivot(i) = iTmp
    end if
  end do
end do

iQual(:) = 0
iMax = nQual
i = 0
nQual = 0
do while (i < iMax)
  i = i+1
  if (Diag(iPivot(i)) >= DiaMin) then
    nQual = nQual+1
    iQual(nQual) = iPivot(i)
  else
    i = iMax+1
  end if
end do

end subroutine CD_DiaMax

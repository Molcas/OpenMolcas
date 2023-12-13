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
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine SortDiag(HH,EigVec,nVec,nDim)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nVec, nDim
real(kind=wp), intent(inout) :: HH(nTri_Elem(nDim)), EigVec(nDim,nVec)
integer(kind=iwp) :: i, iMax, ii, jj
integer(kind=iwp), external :: idAMax_

do i=1,nVec-1
  iMax = idAMax_(nVec-i+1,EigVec(i,i),nDim)
  if (iMax > 1) then
    iMax = iMax+i-1
    ii = nTri_Elem(i)
    jj = nTri_Elem(iMax)
    call dSwap_(1,HH(ii),1,HH(jj),1)
    call dSwap_(nDim,EigVec(:,i),1,EigVec(:,iMax),1)
  end if
end do

end subroutine SortDiag

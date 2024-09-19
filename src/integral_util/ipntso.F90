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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

function iPntSO(j1,j2,lOper,nbas)
!***********************************************************************
!                                                                      *
! Object: to compute the offset to an one-electron integral block.     *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iPntSO
integer(kind=iwp), intent(in) :: j1, j2, lOper, nbas(0:7)
integer(kind=iwp) :: iIrrep, ij, iSmLbl, jIrrep, jMax

iPntSO = 0
iSmLbl = lOper
do iIrrep=0,j1
  jMax = iIrrep
  if (iIrrep == j1) jMax = j2-1
  do jIrrep=0,jMax
    ij = Mul(iIrrep+1,jIrrep+1)-1
    if (.not. btest(iSmLbl,ij)) cycle
    if (iIrrep == jIrrep) then
      iPntSO = iPntSO+nTri_Elem(nBas(iIrrep))
    else
      iPntSO = iPntSO+nBas(iIrrep)*nBas(jIrrep)
    end if
  end do
end do

return

end function iPntSO

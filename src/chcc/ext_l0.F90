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

subroutine Ext_L0(V1,V2,no,dimij,dimc,nbs)
! this routine does:
! V2(i,j,m') <- V1(p,q,m')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: no, dimij, dimc, nbs
real(kind=wp) :: V1(nbs,nbs,dimc), V2(dimij,dimc)
integer(kind=iwp) :: i, ij, j, m

do m=1,dimc
  ij = 0
  do i=1,no
    do j=1,i
      ij = ij+1
      V2(ij,m) = V1(i,j,m)
    end do
  end do
end do

return

end subroutine Ext_L0

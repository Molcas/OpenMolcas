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

subroutine Ext_L1(V1,V2,no,dima,dimc,adda,nbs)
! this routine does:
! V2(i,a',m') <- V1(p,q,m')

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: no, dima, dimc, adda, nbs
real(kind=wp) :: V1(nbs,nbs,dimc), V2(no,dima,dimc)
integer(kind=iwp) :: a, i, m

do m=1,dimc
  do a=1,dima
    do i=1,no
      V2(i,a,m) = V1(i,adda+a,m)
    end do
  end do
end do

return

end subroutine Ext_L1

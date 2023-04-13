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

subroutine Ext_L2s(V1,V2,dima,dimab,dimc,adda,addb,nbs)
! this routine does:
! V2(a'b',m') <- V1(p,q,m') for aGrp=bGrp

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimab, dimc, adda, addb, nbs
real(kind=wp) :: V1(nbs,nbs,dimc), V2(dimab,dimc)
integer(kind=iwp) :: a, ab, b, m

do m=1,dimc
  ab = 0
  do a=1,dima
    do b=1,a
      ab = ab+1
      V2(ab,m) = V1(adda+a,addb+b,m)
    end do
  end do
end do

return

end subroutine Ext_L2s

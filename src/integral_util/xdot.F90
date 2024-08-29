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

function XDot(A,n,m,k,l)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: XDot
integer(kind=iwp), intent(in) :: n, m, k, l
real(kind=wp), intent(in) :: A(n*m+1,k,l)
integer(kind=iwp) :: ik, il
real(kind=wp), external :: DDot_

XDot = Zero
do ik=1,k
  do il=1,l
    XDot = XDot+DDot_(n*m,[One],0,A(1:n*m,ik,il),1)
  end do
end do

end function XDot

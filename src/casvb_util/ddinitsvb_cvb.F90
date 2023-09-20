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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine ddinitsvb_cvb(ifollow1,isaddle1,ip1)

use casvb_global, only: ifollow, ipdd, isaddledd, nroot

implicit real*8(a-h,o-z)

ifollow = ifollow1
isaddledd = isaddle1
nroot = max(1,isaddledd+1)
ipdd = ip1

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(dum)
  call Unused_integer(nfrdim1)
end if

end subroutine ddinitsvb_cvb

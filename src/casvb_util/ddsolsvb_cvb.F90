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

subroutine ddsolsvb_cvb(dum,rhsp,itdav,maxdav,nfrdim1,solp,solp_res,eig,eig_res)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: itdav, maxdav, nfrdim1
real(kind=wp) :: dum, rhsp(maxdav), solp(maxdav), solp_res(maxdav), eig, eig_res
real(kind=wp), external :: dnrm2_

call fmove_cvb(rhsp,solp,itdav)
eig = dnrm2_(itdav,solp,1)
call dscal_(itdav,One/eig,solp,1)
eig_res = eig
call fmove_cvb(solp,solp_res,itdav)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(dum)
  call Unused_integer(nfrdim1)
end if

end subroutine ddsolsvb_cvb

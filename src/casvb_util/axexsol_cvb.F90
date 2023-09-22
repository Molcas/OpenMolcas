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

subroutine axexsol_cvb(ap,dum,itdav,maxdav,nfrdim1,solp,solp_res,eig,eig_res)
! Diagonalize Hamiltonian in Davidson subspace:

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: itdav, maxdav, nfrdim1
real(kind=wp) :: ap(maxdav,maxdav), dum, solp(maxdav), solp_res(maxdav), eig, eig_res
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2
integer(kind=iwp), external :: mstackr_cvb

i1 = mstackr_cvb(itdav)
i2 = mstackr_cvb(itdav*itdav)

call axexsol2_cvb(ap,work(i1),work(i2),dum,itdav,maxdav,solp,solp_res,eig,eig_res)

call mfreer_cvb(i1)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(nfrdim1)

end subroutine axexsol_cvb

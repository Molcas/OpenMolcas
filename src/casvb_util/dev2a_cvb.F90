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

!***********************************************************************
!*                                                                     *
!*  DEV2A  := calculate two-electron Hessian                           *
!*                                                                     *
!***********************************************************************
subroutine dev2a_cvb(v1,v2,cfrom,hessorb,oaa2,aa1)
! Calculates V1 EijEkl CFROM and V2 EijEkl CFROM

use casvb_global, only: i1alf, i1bet, iafrm, iapr, iato, ibfrm, ibpr, ibto, ixapr, ixbpr, phato, phbto
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
real(kind=wp) :: v1(0:ndet), v2(0:ndet), cfrom(0:ndet), hessorb(nprorb,nprorb), oaa2, aa1
integer(kind=iwp) :: icfrom

icfrom = nint(cfrom(0))
n_2el = n_2el+2
if (iform_ci(icfrom) /= 0) then
  write(u6,*) ' Unsupported format in DEV2A :',iform_ci(icfrom)
  call abend_cvb()
end if

call dev2a_2_cvb(v1(1:),v2(1:),cfrom(1:),hessorb,oaa2,aa1,nprorb,i1alf,i1bet,iafrm,ibfrm,iato,ibto,phato,phbto,iapr,ixapr,ibpr, &
                 ixbpr,npvb,nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3))

return

end subroutine dev2a_cvb

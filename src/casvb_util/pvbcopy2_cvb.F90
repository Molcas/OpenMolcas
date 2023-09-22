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

subroutine pvbcopy2_cvb(cfrom,cto,iapr,ixapr,ret,ic)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cfrom(nda,ndb), cto(nda,ndb), ret
integer(kind=iwp) :: iapr(ndetvb), ixapr(nda+1), ic
integer(kind=iwp) :: ia, ib, idetvb, ixa

if (ic == 0) then
  call fzero(cto,nda*ndb)
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      cto(ia,ib) = cfrom(ia,ib)
    end do
  end do
else if (ic == 1) then
  ret = Zero
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      ret = ret+cto(ia,iapr(ixa))*cfrom(ia,iapr(ixa))
    end do
  end do
end if

return

end subroutine pvbcopy2_cvb

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

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension cfrom(nda,ndb), cto(nda,ndb)
dimension iapr(ndetvb), ixapr(nda+1)

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
  ret = zero
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      ret = ret+cto(ia,iapr(ixa))*cfrom(ia,iapr(ixa))
    end do
  end do
end if

return

end subroutine pvbcopy2_cvb

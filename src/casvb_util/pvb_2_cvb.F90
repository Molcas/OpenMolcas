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
!*  PVB    := Zero parts of CI vector not in VB wfn.                   *
!*                                                                     *
!***********************************************************************
subroutine pvb_2_cvb(cfrom,cto,csk,iapr,ixapr,mult)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: cfrom(nda,ndb), cto(nda,ndb), csk(ndetvb)
integer(kind=iwp) :: iapr(ndetvb), ixapr(nda+1), mult
integer(kind=iwp) :: ia, ib, idetvb, ixa

if (mult == -1) then
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      csk(idetvb) = cfrom(ia,ib)
    end do
  end do
else if (mult == 0) then
  call fzero(cto,nda*ndb)
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      cto(ia,ib) = cfrom(ia,ib)
      csk(idetvb) = cfrom(ia,ib)
    end do
  end do
else if (mult == 1) then
  csk(1) = Zero
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      csk(1) = csk(1)+cto(ia,iapr(ixa))*cfrom(ia,iapr(ixa))
    end do
  end do
else if (mult == 2) then
  call fzero(cto,nda*ndb)
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      cto(ia,ib) = csk(idetvb)
    end do
  end do
else if (mult == 3) then
  csk(1) = Zero
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ! CFROM is really CDETVB
      csk(1) = csk(1)+cto(ia,iapr(ixa))*cfrom(idetvb,1)
    end do
  end do
end if

return

end subroutine pvb_2_cvb

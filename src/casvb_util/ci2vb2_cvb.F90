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

subroutine ci2vb2_cvb(civec,cvbdet,iapr,ixapr,ret,ic)

use casvb_global, only: nda, ndb, ndetvb
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: civec(nda,ndb), cvbdet(ndetvb)
integer(kind=iwp), intent(in) :: iapr(ndetvb), ixapr(nda+1), ic
real(kind=wp), intent(_OUT_) :: ret
integer(kind=iwp) :: ia, ib, idetvb, ixa

if (ic == 0) then
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      cvbdet(idetvb) = civec(ia,ib)
    end do
  end do
else if (ic == 1) then
  civec(:,:) = Zero
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      civec(ia,ib) = cvbdet(idetvb)
    end do
  end do
else if (ic == 2) then
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      civec(ia,ib) = civec(ia,ib)+cvbdet(idetvb)
    end do
  end do
else if (ic == 3) then
  ret = Zero
  idetvb = 0
  do ia=1,nda
    do ixa=ixapr(ia),ixapr(ia+1)-1
      idetvb = idetvb+1
      ib = iapr(ixa)
      ret = ret+civec(ia,ib)*cvbdet(idetvb)
    end do
  end do
end if

return

end subroutine ci2vb2_cvb

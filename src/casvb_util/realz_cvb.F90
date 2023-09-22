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

subroutine realz_cvb(arr,nmax,nread,ifc)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nmax, nread, ifc
real(kind=wp) :: arr(nmax)
#include "WrkSpc.fh"
integer(kind=iwp), parameter :: nbuf = 100
integer(kind=iwp) :: i1, nleft, nread1
integer(kind=iwp), external :: mstackr_cvb

i1 = mstackr_cvb(nbuf)
nread = 0
do
  call fzero(work(i1),nbuf)
  nleft = nmax-nread
  call real_cvb(work(i1),min(nbuf,nleft),nread1,ifc)
  call fmove_cvb(work(i1),arr(1+nread),nread1)
  nread = nread+nread1
  if (nread1 <= 0) exit
end do
call mfreer_cvb(i1)

return

end subroutine realz_cvb

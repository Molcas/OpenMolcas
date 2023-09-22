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

subroutine fxdx_cvb(fx,fast,dx)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: fx, dx(*)
logical(kind=iwp) :: fast
#include "main_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iwf1, iwf2

dxmove = .false.
iwf1 = lw(12)
iwf2 = lw(13)
call upd_cvb(dx,work(iwf1),work(iwf2))
if (.not. memplenty) then
  call ciwr_cvb(work(lc(2)),61002.2_wp)
  call ciwr_cvb(work(lc(3)),61003.2_wp)
  call ciwr_cvb(work(lc(4)),61004.2_wp)
  call setcnt2_cvb(2,0)
  call setcnt2_cvb(3,0)
  call setcnt2_cvb(4,0)
end if
call setcnt2_cvb(6,0)
call setcnt2_cvb(7,0)
call setcnt2_cvb(8,0)
if (icrit == 1) then
  call fx_svb1_cvb(fx,fast,work(iwf1),work(iwf2),work(lc(1)),work(lc(6)),work(lc(7)),work(lc(8)),work(lw(4)),work(lw(5)), &
                   work(lw(6)),work(lw(9)))
else if (icrit == 2) then
  call fx_evb1_cvb(fx,fast,work(iwf1),work(iwf2),work(lc(1)),work(lc(6)),work(lc(7)),work(lc(8)),work(lw(4)),work(lw(5)), &
                   work(lw(6)),work(lw(9)))
end if
if (.not. memplenty) then
  call ciwr_cvb(work(lc(6)),61006.2_wp)
  call ciwr_cvb(work(lc(7)),61007.2_wp)
  call ciwr_cvb(work(lc(8)),61008.2_wp)
  call cird_cvb(work(lc(2)),61002.2_wp)
  call cird_cvb(work(lc(3)),61003.2_wp)
  call cird_cvb(work(lc(4)),61004.2_wp)
end if
! Figure out what we just calculated, and make it up2date:
if (dxmove) then
  if (icrit == 1) then
    call make_cvb('SVB')
  else if (icrit == 2) then
    call make_cvb('EVB')
  end if
else
  if (icrit == 1) then
    call make_cvb('SVBTRY')
  else if (icrit == 2) then
    call make_cvb('EVBTRY')
  end if
end if

return

end subroutine fxdx_cvb

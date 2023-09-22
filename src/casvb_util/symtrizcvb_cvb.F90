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

subroutine symtrizcvb_cvb(vecstr)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: vecstr(nvb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1, i2
real(kind=wp) :: dum(1)
integer(kind=iwp), external :: mstackr_cvb

if (iconstruc == 0) then
  return
else if (iconstruc == 1) then
  i1 = mstackr_cvb(ndetvb)
  i2 = mstackr_cvb(nvb)
  call symtrizcvb2_cvb(vecstr,iwork(ls(13)),iwork(ls(16)),work(i1),work(i2))
  call mfreer_cvb(i1)
  call symtrizcvb3_cvb(vecstr,iwork(ls(10)))
else if (iconstruc == 2) then
  call schmidtd_cvb(work(ls(15)),nconstr,vecstr,1,dum,nvb,0)
end if

return

end subroutine symtrizcvb_cvb

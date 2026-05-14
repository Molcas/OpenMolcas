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

subroutine o12sa_cvb( &
#                    define _CALLING_
#                    include "opta_interface.fh"
                    )

use casvb_global, only: civb1, civb2, civb3, civb4, civb5, civb6, civb7, civb8, cvb, cvbdet, have_solved_it, icnt_ci, nprorb, nv, &
                        nvb, nvguess, nvrestart, nvrhs, odx, orbs, strucopt
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "opta_interface.fh"
real(kind=wp), pointer :: vuse(:), vuse2(:)
integer(kind=iwp) :: icivb

nvrestart = 0
nvguess = 0
nvrhs = 0
have_solved_it = .false.

! Find CIVBS:
nullify(vuse)
vuse2 => civb3
if (nv >= 1) then
  icivb = nint(civb1(0))
  if (icnt_ci(icivb) == 4) vuse => civb1
end if
if (nv >= 2) then
  icivb = nint(civb2(0))
  if (icnt_ci(icivb) == 4) vuse => civb2
end if
if (nv >= 3) then
  icivb = nint(civb3(0))
  if (icnt_ci(icivb) == 4) vuse => civb3
end if
if (nv >= 4) then
  icivb = nint(civb4(0))
  if (icnt_ci(icivb) == 4) vuse => civb4
end if
if (nv >= 5) then
  icivb = nint(civb5(0))
  if (icnt_ci(icivb) == 4) vuse => civb5
end if
if (nv >= 6) then
  icivb = nint(civb6(0))
  if (icnt_ci(icivb) == 4) vuse => civb6
end if
if (nv >= 7) then
  icivb = nint(civb7(0))
  if (icnt_ci(icivb) == 4) vuse => civb7
end if
if (nv >= 8) then
  icivb = nint(civb8(0))
  if (icnt_ci(icivb) == 4) vuse => civb8
end if
if (associated(vuse)) then
  if (associated(vuse,civb3)) vuse2 => civb2
  call o12sa2_cvb(nparam,vuse2,vuse,cvbdet,cvb)
else
  if (strucopt) then
    call ddguess_cvb(cvb,nvb,nprorb)
  else
    call ddguess_cvb([One],1,0)
  end if
end if

call o12sa3_cvb(odx,cvb,orbs,civb1,civb2,civb3,cvbdet,nparam)

return

end subroutine o12sa_cvb

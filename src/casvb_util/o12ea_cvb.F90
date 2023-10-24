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

subroutine o12ea_cvb( &
#                    define _CALLING_
#                    include "opta_interface.fh"
                    )

use casvb_global, only: civb1, civb2, civb3, civb4, civb5, civb6, civb7, civb8, cvb, cvbdet, have_solved_it, nprorb, nv, nvb, &
                        nvguess, nvrestart, nvrhs, strucopt
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "opta_interface.fh"
real(kind=wp), pointer :: vuse2(:), vuse_h(:), vuse_s(:)
logical(kind=iwp), external :: tstcnt_cvb ! ... Content of CI vectors ...

nvrestart = 0
nvguess = 0
nvrhs = 0
have_solved_it = .false.

! Find CIVBS & CIVBH:
nullify(vuse_s)
nullify(vuse_h)
vuse2 => civb3
if (nv >= 1) then
  if (tstcnt_cvb(civb1,4)) vuse_s => civb1
  if (tstcnt_cvb(civb1,5)) vuse_h => civb1
end if
if (nv >= 2) then
  if (tstcnt_cvb(civb2,4)) vuse_s => civb2
  if (tstcnt_cvb(civb2,5)) vuse_h => civb2
end if
if (nv >= 3) then
  if (tstcnt_cvb(civb3,4)) vuse_s => civb3
  if (tstcnt_cvb(civb3,5)) vuse_h => civb3
end if
if (nv >= 4) then
  if (tstcnt_cvb(civb4,4)) vuse_s => civb4
  if (tstcnt_cvb(civb4,5)) vuse_h => civb4
end if
if (nv >= 5) then
  if (tstcnt_cvb(civb5,4)) vuse_s => civb5
  if (tstcnt_cvb(civb5,5)) vuse_h => civb5
end if
if (nv >= 6) then
  if (tstcnt_cvb(civb6,4)) vuse_s => civb6
  if (tstcnt_cvb(civb6,5)) vuse_h => civb6
end if
if (nv >= 7) then
  if (tstcnt_cvb(civb7,4)) vuse_s => civb7
  if (tstcnt_cvb(civb7,5)) vuse_h => civb7
end if
if (nv >= 8) then
  if (tstcnt_cvb(civb8,4)) vuse_s => civb8
  if (tstcnt_cvb(civb8,5)) vuse_h => civb8
end if
if (associated(vuse_h,civb3) .or. associated(vuse_s,civb3)) vuse2 => civb2
if (associated(vuse_h,civb2) .or. associated(vuse_s,civb2)) vuse2 => civb4
if (associated(vuse_s) .and. associated(vuse_h)) then
  call o12ea2_cvb(nparam,vuse2,vuse_s,vuse_h,cvbdet,cvb)
else if (strucopt) then
  call ddguess_cvb(cvb,nvb,nprorb)
else
  call ddguess_cvb([One],1,0)
end if
call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb3)

return

end subroutine o12ea_cvb

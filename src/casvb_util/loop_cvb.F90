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

subroutine loop_cvb(nel,nk,nkmin,nkmax,rc)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nel, nkmin(0:nel), nkmax(0:nel)
integer(kind=iwp), intent(inout) :: nk(0:nel)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: iel, ik

rc = 0
do iel=1,nel-1
  ik = nk(iel)
  ! Tests: (1) IK+1 OCCUPIED; (2) IK NOT OCCUPIED; (3) IK MINIMAL
  if (.not. ((nk(iel+1)-ik == 1) .or. (ik == nk(iel-1)) .or. (ik == nkmin(iel)))) then
    ! SITUATION IS :  IEL       \     <= NOT MINIMAL
    !                 IEL+1     |
    nk(1:iel-1) = min(nkmax(1:iel-1),ik-1)
    nk(iel) = ik-1
    rc = 1
    return
  end if
end do
! Maximize the loop on exit
nk(0:nel-1) = nkmax(0:nel-1)

return

end subroutine loop_cvb

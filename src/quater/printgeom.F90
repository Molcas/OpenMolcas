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
! Copyright (C) Yannick Carissan                                       *
!***********************************************************************
!  PrintGeom
!
!> @brief
!>   Print the geometry
!> @author Y. Carissan
!>
!> @details
!> Print the geometry.
!>
!> @param[in] iLU   Logic unit number
!> @param[in] Nat   number of atoms
!> @param[in] title title to be printed
!> @param[in] Geom  XYZ coordinates
!> @param[in] lbl   Atom labels
!***********************************************************************

#include "macros.fh"

subroutine PrintGeom(iLU,Nat,title,Geom,lbl)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iLU, nat
character(len=20), intent(in) :: title, lbl(nat)
real(kind=wp), intent(in) :: Geom(3,nat)
integer(kind=iwp) :: iat
unused_var(lbl)

write(u6,'(a8,i1)') '--- GEOM'
write(u6,'(i4)') nat
write(u6,*) title
do iat=1,nat
  !if (debug) then
  !  ! lbl is not defined. I'll comment out this statement
  !  ! until someone needs it and replace it with the
  !  ! the standard printout.
  !  write(u6,'(i3,5x,a10,8x,3f16.8)') iat,lbl(iat),Geom(:,iat)
  !else
  !  write(u6,   '(5x,a10,8x,3f16.8)')     lbl(iat),Geom(:,iat)
  !end if
  write(u6,'(i3,5x,3f16.8)') iat,Geom(:,iat)
end do

if (iLU /= -1) then
  write(iLU,'(a8,i1)') '--- GEOM'
  write(iLU,'(i4)') nat
  write(iLU,*) title
  do iat=1,nat
    !write(iLU,     '(a10,8x,3f16.8)')     lbl(iat),Geom(:,iat)
    write(iLU,'(i3,5x,3f16.8)') iat,Geom(:,iat)
  end do
end if

return

end subroutine PrintGeom

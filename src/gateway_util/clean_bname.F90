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
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Clean_BName
!
!> @brief
!>   Clean a basis function name for printing.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Clean up a basis function name for printing purposes, by removing
!> unneeded zeros and padding.
!>
!> @param[in] BName  Basis function name
!> @param[in] Offset Initial characters not to be touched
!>
!> @return Cleaned basis function name
!***********************************************************************

function Clean_BName(BName,Offset)

use Definitions, only: iwp

implicit none
#include "Molcas.fh"
character(len=LenIn8) :: Clean_BName
character(len=*), intent(in) :: BName
integer(kind=iwp), intent(in) :: Offset
integer(kind=iwp) :: i, Err
logical(kind=iwp) :: Cart
character(len=8) :: Clean

Clean = BName(Offset+1:)
! For spherical functions, the 3rd character is a letter
!   (counting s and p as spherical),
! for Cartesian functions, it is a number
read(Clean(3:3),'(I1)',iostat=Err) i
Cart = Err == 0
#include "macros.fh"
unused_var(i)

if (Cart) then
  ! For Cartesian functions, remove zeros if all indices are below 10,
  ! and shift everything one space to the right
  if ((Clean(2:2) == '0') .and. (Clean(4:4) == '0') .and. (Clean(6:6) == '0')) Clean(2:) = Clean(3:3)//Clean(5:5)//Clean(7:7)
  Clean(2:) = Clean(1:)
  Clean(1:1) = ' '
else
  ! For spherical functions, remove leading zero in m number,
  ! replace leading zero in shell number with space,
  ! replace "*0" tag with a single *
  if (Clean(1:1) == '0') Clean(1:1) = ' '
  if (Clean(1:2) == '*0') Clean(1:2) = ' *'
  if (Clean(4:4) == '0') Clean(4:) = Clean(5:)
end if

Clean_BName = BName(1:Offset)//Clean

return

end function Clean_BName

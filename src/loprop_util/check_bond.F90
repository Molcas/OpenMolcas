!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function Check_Bond(CXi,CXj,iANr,jANr,Factor)
! Returns true if the bond should be included, otherwise false

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: Check_Bond
real(kind=wp), intent(in) :: CXi(3), CXj(3), Factor
integer(kind=iwp), intent(in) :: iANr, jANr
real(kind=wp) :: Bond_Length, Bond_Max, Radius_i, Radius_j
real(kind=wp), external :: Bragg_Slater

Check_Bond = .true.
if (Factor < Zero) then
  Check_Bond = .true.
else
  Radius_i = Bragg_Slater(iANr)
  Radius_j = Bragg_Slater(jANr)
  Bond_Length = sqrt((CXi(1)-CXj(1))**2+(CXi(2)-CXj(2))**2+(CXi(3)-CXj(3))**2)
  Bond_Max = Factor*(Radius_i+Radius_j)
  if (Bond_Length > Bond_Max) Check_Bond = .false.
end if

return

end function Check_Bond

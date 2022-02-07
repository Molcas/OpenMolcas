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

subroutine ComputeFunc(nAtoms,nOrb2Loc,PA,Functional,Debug)
! Author: Y. Carissan

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nOrb2Loc
real(kind=wp), intent(in) :: PA(nOrb2Loc,nOrb2Loc,nAtoms)
real(kind=wp), intent(out) :: Functional
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: iAt, iMO_s

Functional = Zero
do iAt=1,nAtoms
  do iMO_s=1,nOrb2Loc
    Functional = Functional+PA(iMO_s,iMO_s,iAt)**2
  end do
end do

if (Debug) then
  write(u6,*) 'ComputeFunc: Functional: ',Functional
end if

return

end subroutine ComputeFunc

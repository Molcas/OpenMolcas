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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine GradPrep(nState,UEFF,VECROT)

use caspt2_global, only: iRoot1, iRoot2, jStLag
use caspt2_module, only: IFMSCOUP, JSTATE
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nState
real(kind=wp), intent(in) :: UEFF(nState,nState)
real(kind=wp), intent(inout) :: VECROT(nState)

! If i == j, UIi*dHij/dx*UJj
! If i \= j, (UIi*UJj+UJi*UIj)*dHij/dx*0.5

!! Construct the rotation vector
if (IFMSCOUP) then
  VECROT(:) = Half*(UEFF(:,iRoot1)*UEFF(jState,iRoot2)+UEFF(:,iRoot2)*UEFF(jState,iRoot1))
else
  !write(u6,*) 'jState in gradprep: ',jstate
  VECROT(1:nState) = Zero
  VECROT(jState) = One
end if
jStLag = jState

end subroutine GradPrep

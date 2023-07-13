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

subroutine NACInt(xyz,nCent,H12,Bf,lWrite_,Label,dBf,ldB,lIter)

use Slapaf_Info, only: NAC
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent, lIter
real(kind=wp), intent(in) :: xyz(3,nCent)
real(kind=wp), intent(out) :: H12, Bf(3,nCent)
logical(kind=iwp), intent(in) :: lWrite_, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBf(3,nCent,3,nCent)
integer(kind=iwp) :: iCent
real(kind=wp) :: Fact
integer(kind=iwp), external :: iDeg

H12 = Zero
if (lWrite_) write(u6,'(2A,F18.8,A,F18.8,A)') Label,' : H12               = ',H12,' hartree '

! Compute the WDC B-matrix

do iCent=1,nCent
  Fact = real(iDeg(xyz(1,iCent)),kind=wp)
  Bf(:,iCent) = Fact*NAC(:,iCent,lIter)
end do
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(u6,*) 'NACInt, lIter:',lIter
call RecPrt('NAC',' ',NAC(1,1,lIter),3,nCent)
call RecPrt('Bf',' ',Bf,3,nCent)
#endif

! Compute the cartesian derivative of the B-Matrix.

if (ldB) dBf(:,:,:,:) = Zero

return

end subroutine NACInt

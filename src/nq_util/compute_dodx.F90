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

subroutine Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx,Px)

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, iAtom, iCar
real(kind=wp), intent(in) :: ZA(nAtoms), RA(3,nAtoms), T(3), O(3,3), EVal(3), dTdRAi
logical(kind=iwp), intent(in) :: Rot_Corr
real(kind=wp), intent(out) :: dMdx(3,3), dOdx(3,3), Px(3,3)
real(kind=wp) :: Alpha, Beta, Gmma, OtMx(3,3), OtMxO(3,3)
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Rot_Corr) then
  dOdx(:,:) = Zero
  return
end if

! Differentiate the nuclear charge moment tensor, M.

call Compute_dMdx(ZA,RA,nAtoms,T,iAtom,iCar,dTdRAi,dMdx)

! Form O(t) dMdRAi O

call DGEMM_('T','N',3,3,3,One,O,3,dMdx,3,Zero,OtMx,3)
call DGEMM_('N','N',3,3,3,One,OtMx,3,O,3,Zero,OtMxO,3)

! Compute the off diagonal elements in Px

if (abs(OtMxO(2,3)+OtMxO(3,2)) < Thrs) then
  if (abs(EVal(2)-EVal(3)) < Thrs) then
    Alpha = One
  else
    Alpha = Zero
  end if
else
  if (abs(EVal(2)-EVal(3)) < Thrs) then
    Alpha = huge(Alpha)
  else
    Alpha = -(OtMxO(2,3)+OtMxO(3,2))/(Two*(EVal(2)-EVal(3)))
  end if
end if

if (abs(OtMxO(1,3)+OtMxO(3,1)) < Thrs) then
  if (abs(EVal(3)-EVal(1)) < Thrs) then
    Beta = One
  else
    Beta = Zero
  end if
else
  if (abs(EVal(3)-EVal(1)) < Thrs) then
    Beta = huge(Beta)
  else
    Beta = -(OtMxO(1,3)+OtMxO(3,1))/(Two*(EVal(3)-EVal(1)))
  end if
end if
if (abs(OtMxO(1,2)+OtMxO(2,1)) < Thrs) then
  if (abs(EVal(1)-EVal(2)) < Thrs) then
    Gmma = One
  else
    Gmma = Zero
  end if
else
  if (abs(EVal(1)-EVal(2)) < Thrs) then
    Gmma = huge(Gmma)
  else
    Gmma = -(OtMxO(1,2)+OtMxO(2,1))/(Two*(EVal(1)-EVal(2)))
  end if
end if

Px(:,:) = Zero
Px(1,2) = Gmma
Px(2,1) = -Gmma
Px(1,3) = -Beta
Px(3,1) = Beta
Px(2,3) = Alpha
Px(3,2) = -Alpha

! Finally evaluate dO/dRAi

call DGEMM_('N','N',3,3,3,One,O,3,Px,3,Zero,dOdx(1,1),3)
#ifdef _DEBUGPRINT_
!call RecPrt('M','(3G20.10)',M,3,3)
!call RecPrt('RotGrd: O','(3G20.10)',O,3,3)
!call RecPrt('RotGrd: EVal',' ',EVal,3,1)
!call RecPrt('RotGrd: dMdx','(3G20.10)',dMdx,3,3)
!call RecPrt('RotGrd: OtMxO','(3G20.10)',OtMxO,3,3)
!write(u6,*) 'A,B,G=',Alpha,Beta,Gmma
!call RecPrt('RotGrd: Px','(3G20.10)',Px,3,3)
!call RecPrt('dOdx','(3G20.10)',dOdx,3,3)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_dOdx

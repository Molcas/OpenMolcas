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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 ZA(nAtoms), RA(3,nAtoms), T(3), O(3,3), EVal(3), dOdx(3,3), dMdx(3,3), Px(3,3)
logical Rot_Corr
! Local arrays
real*8 OtMx(3,3), OtMxO(3,3)

!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Rot_Corr) then
  call FZero(dOdx,9)
  return
end if

! Differentiate the nuclear charge moment tensor, M.

call Compute_dMdx(ZA,RA,nAtoms,T,iAtom,iCar,dTdRAi,dMdx)

! Form O(t) dMdRAi O

call DGEMM_('T','N',3,3,3,1.0d0,O,3,dMdx,3,0.0d0,OtMx,3)
call DGEMM_('N','N',3,3,3,1.0d0,OtMx,3,O,3,0.0d0,OtMxO,3)

! Compute the off diagonal elements in Px

if (abs(OtMxO(2,3)+OtMxO(3,2)) < 1.0D-14) then
  if (abs(EVal(2)-EVal(3)) < 1.0D-14) then
    Alpha = One
  else
    Alpha = Zero
  end if
else
  if (abs(EVal(2)-EVal(3)) < 1.0D-14) then
    Alpha = 1.0d99
  else
    Alpha = -(OtMxO(2,3)+OtMxO(3,2))/(Two*(EVal(2)-EVal(3)))
  end if
end if

if (abs(OtMxO(1,3)+OtMxO(3,1)) < 1.0D-14) then
  if (abs(EVal(3)-EVal(1)) < 1.0D-14) then
    Beta = One
  else
    Beta = Zero
  end if
else
  if (abs(EVal(3)-EVal(1)) < 1.0D-14) then
    Beta = 1.0d99
  else
    Beta = -(OtMxO(1,3)+OtMxO(3,1))/(Two*(EVal(3)-EVal(1)))
  end if
end if
if (abs(OtMxO(1,2)+OtMxO(2,1)) < 1.0D-14) then
  if (abs(EVal(1)-EVal(2)) < 1.0D-14) then
    Gamma = One
  else
    Gamma = Zero
  end if
else
  if (abs(EVal(1)-EVal(2)) < 1.0D-14) then
    Gamma = 1.0d99
  else
    Gamma = -(OtMxO(1,2)+OtMxO(2,1))/(Two*(EVal(1)-EVal(2)))
  end if
end if

call FZero(Px,9)
Px(1,2) = Gamma
Px(2,1) = -Gamma
Px(1,3) = -Beta
Px(3,1) = Beta
Px(2,3) = Alpha
Px(3,2) = -Alpha

! Finally evaluate dO/dRAi

call DGEMM_('N','N',3,3,3,1.0d0,O,3,Px,3,0.0d0,dOdx(1,1),3)
#ifdef _DEBUGPRINT_
!call RecPrt('M','(3G20.10)',M,3,3)
!call RecPrt('RotGrd: O','(3G20.10)',O,3,3)
!call RecPrt('RotGrd: EVal',' ',EVal,3,1)
!call RecPrt('RotGrd: dMdx','(3G20.10)',dMdx,3,3)
!call RecPrt('RotGrd: OtMxO','(3G20.10)',OtMxO,3,3)
!write(6,*) 'A,B,G=',Alpha,Beta,Gamma
!call RecPrt('RotGrd: Px','(3G20.10)',Px,3,3)
!call RecPrt('dOdx','(3G20.10)',dOdx,3,3)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_dOdx

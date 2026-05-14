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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine CrSph1(Win,ijkla,Scrt,nScrt,Coeff3,kCar,kSph,Tr3,Coeff4,lCar,lSph,Tr4,Wout,mcd)
!***********************************************************************
!                                                                      *
! Object : to transform the two-electron integrals from cartesian      *
!          gaussians to spherical gaussians.                           *
!                                                                      *
!          Observe that most of the time Win and Wout will overlap.    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ijkla, nScrt, kCar, kSph, lCar, lSph, mcd
real(kind=wp), intent(in) :: Win(ijkla*kCar*lCar), Coeff3(kCar,kCar), Coeff4(lCar,lCar)
real(kind=wp), intent(out) :: Scrt(nScrt)
real(kind=wp), intent(inout) :: Wout(mcd*ijkla)
logical(kind=iwp), intent(in) :: Tr3, Tr4

if (Tr3 .and. Tr4) then
  !call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
  ! Starting with IJKL,a,cd transforming to D,IJKL,a,c
  !call xxDGeMul(Coeff4,lCar,'T',Win,ijkla*kCar,'T',Scrt,lSph,lSph,lCar,ijkla*kCar)
  call TTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)

  !call RecPrt(' In CrSph: (a0|cD) ',' ',Scrt,lSph*ijkla,kCar)
  !call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
  ! Transform D,IJKL,a,c to CD,IJKL,a
  !call xxDGeMul(Coeff3,kCar,'T',Scrt,lSph*ijkla,'T',Wout,kSph,kSph,kCar,lSph*ijkla)
  call TTMul(Coeff3,Scrt,Wout,kCar,kSph,lSph*ijkla)
else if (Tr4) then
  !call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
  ! Starting with IJKL,a,cd transforming to D,IJKL,a,c
  !call xxDGeMul(Coeff4,lCar,'T',Win,ijkla*kCar,'T',Scrt,lSph,lSph,lCar,ijkla*kCar)
  call TTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
  ! Transpose D,IJKL,a,c to cD,IJKL,a
  call DGeTMO(Scrt,lSph*ijkla,lSph*ijkla,kCar,Wout,kCar)
else if (Tr3) then
  ! Transpose IJKL,a,c,d to d,IJKL,a,c
  call DGeTMO(Win,ijkla*kCar,ijkla*kCar,lCar,Scrt,lCar)

  !call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
  ! Transform D,IJKL,a,c to CD,IJKL,a
  !call xxDGeMul(Coeff3,kCar,'T',Scrt,lSph*ijkla,'T',Wout,kSph,kSph,kCar,lSph*ijkla)
  call TTMul(Coeff3,Scrt,Wout,kCar,kSph,lSph*ijkla)
else
  ! Transpose IJKL,a,c,d to c,d,IJKL,a
  Scrt(1:ijkla*kCar*lCar) = Win(:)
  if (kCar*lCar /= 1) then
    call DGeTMO(Scrt,ijkla,ijkla,kCar*lCar,Wout,kCar*lCar)
  else
    Wout(1:ijkla*kCar*lCar) = Scrt(1:ijkla*kCar*lCar)
  end if
end if

!call RecPrt(' In CrSph1: (a0|CD)  ',' ',Wout,mcd,ijkla)

return

end subroutine CrSph1

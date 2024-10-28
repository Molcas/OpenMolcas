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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine SphCr1(Win,ijkla,Scrt,nScrt,Coeff3,kCar,kSph,Tr3,Pr3,Coeff4,lCar,lSph,Tr4,Pr4,Wout,mcd)
!***********************************************************************
!                                                                      *
! Object : to transform the two-electron integrals from cartesian      *
!          gaussians to real spherical harmonic gaussians.             *
!                                                                      *
!          Observe that most of the time Win and Wout will overlap.    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified to back projection to cartesian gaussians,      *
!             January '92.                                             *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ijkla, nScrt, kCar, kSph, lCar, lSph, mcd
real(kind=wp), intent(in) :: Win(ijkla*kSph*lSph), Coeff3(kCar,kCar), Coeff4(lCar,lCar)
real(kind=wp), intent(out) :: Scrt(nScrt)
real(kind=wp), intent(inout) :: Wout(mcd*ijkla)
logical(kind=iwp), intent(in) :: Tr3, Pr3, Tr4, Pr4

#include "macros.fh"
unused_var(Pr3)
unused_var(Pr4)

!call RecPrt(' In SphCr1: P(AB|CD) ',' ',Win,ijkla,kSph*lSph)
if (Tr3 .and. Tr4) then
  !call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
  ! Starting with IJKL,AB,CD transforming to d,IJKL,AB,C
  !call xxDGeMul(Coeff4,lCar,'N',Win,ijkla*kSph,'T',Scrt,lCar,lCar,lSph,ijkla*kSph)
  call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kSph)

  !call RecPrt(' In SphCr: P(AB|Cd) ',' ',Scrt,lCar*ijkla,kSph)
  !call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
  ! Transform d,IJKL,AB,C to cd,IJKL,AB
  !call xxDGeMul(Coeff3,kCar,'N',Scrt,lCar*ijkla,'T',Wout,kCar,kCar,kSph,lCar*ijkla)
  call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
else if (Tr4) then
  !call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
  ! Starting with IJKL,AB,cD transforming to d,IJKL,AB,c
  !call xxDGeMul(Coeff4,lCar,'N',Win,ijkla*kCar,'T',Scrt,lCar,lCar,lSph,ijkla*kCar)
  call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
  ! Transpose d,IJKL,AB,c to cd,IJKL,AB
  call DGeTMO(Scrt,lCar*ijkla,lCar*ijkla,kCar,Wout,kCar)
else if (Tr3) then
  ! Transpose IJKL,AB,C,d to d,IJKL,AB,C
  call DGeTMO(Win,ijkla*kSph,ijkla*kSph,lCar,Scrt,lCar)

  !call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
  ! Transform d,IJKL,AB,c to cd,IJKL,AB
  !call xxDGeMul(Coeff3,kCar,'N',Scrt,lCar*ijkla,'T',Wout,kCar,kCar,kSph,lCar*ijkla)
  call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
else
  ! Transpose IJKL,AB,cd to cd,IJKL,AB
  Scrt(1:ijkla*kCar*lCar) = Win(1:ijkla*kCar*lCar)
  if (kCar*lCar /= 1) then
    call DGeTMO(Scrt,ijkla,ijkla,kCar*lCar,Wout,kCar*lCar)
  else
    Wout(1:ijkla*kCar*lCar) = Scrt(1:ijkla*kCar*lCar)
  end if
end if

!call RecPrt(' In SphCr1: P(AB|cd)  ',' ',Wout,mcd,ijkla)

return

end subroutine SphCr1

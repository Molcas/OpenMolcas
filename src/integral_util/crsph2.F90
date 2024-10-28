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

subroutine CrSph2(Win,ijkl,ncd,Scrt,nScrt,Coeff1,iCar,iSph,Tr1,Coeff2,jCar,jSph,Tr2,Wout,mab)
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
integer(kind=iwp), intent(in) :: ijkl, ncd, nScrt, iCar, iSph, jCar, jSph, mab
real(kind=wp), intent(in) :: Win(ijkl*ncd*iCar*jCar), Coeff1(iCar,iCar), Coeff2(jCar,jCar)
real(kind=wp), intent(out) :: Scrt(nScrt)
real(kind=wp), intent(inout) :: Wout(ijkl*ncd*mab)
logical(kind=iwp), intent(in) :: Tr1, Tr2

if (Tr1 .and. Tr2) then
  !call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
  ! Starting with CD,IJKL,a,b transform to B,CD,IJKL,a
  !call xxDGeMul(Coeff2,jCar,'T',Win,ijkl*ncd*iCar,'T',Scrt,jSph,jSph,jCar,ijkl*ncd*iCar)
  call TTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iCar)

  !call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
  ! Transform B,CD,IJKL,a to AB,CD,IJKL
  !call xxDGeMul(Coeff1,iCar,'T',Scrt,jSph*ncd*ijkl,'T',Wout,iSph,iSph,iCar,jSph*ncd*ijkl)
  call TTMul(Coeff1,Scrt,Wout,iCar,iSph,jSph*ncd*ijkl)
  ! Transpose AB,CD,IJKL to IJKL,AB,CD
  Scrt(1:mab*ncd*ijkl) = Wout(:)
  call dGeTMO(Scrt,mab*ncd,mab*ncd,ijkl,Wout,ijkl)
else if (Tr2) then
  !call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
  ! Starting with CD,IJKL,a,b transform to B,CD,IJKL,a
  !call xxDGeMul(Coeff2,jCar,'T',Win,ncd*ijkl*iCar,'T',Scrt,jSph,jSph,jCar,ncd*ijkl*iCar)
  call TTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iCar)
  ! Transpose B,CD,IJKL,a to IJKL,aB,CD
  call DGeTMO(Scrt,jSph*ncd,jSph*ncd,ijkl*iCar,Wout,ijkl*iCar)
else if (Tr1) then
  ! Transpose CD,IJKL,a,b to b,CD,IJKL,a
  call DGeTMO(Win,ncd*ijkl*iCar,ncd*ijkl*iCar,jCar,Scrt,jCar)

  !call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
  ! Transform b,CD,IJKL,a to Ab,CD,IJKL
  !call xxDGeMul(Coeff1,iCar,'T',Scrt,jCar*ncd*ijkl,'T',Wout,iSph,iSph,iCar,jCar*ncd*ijkl)
  call TTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
  ! Transpose Ab,CD,IJKL to IJKL,Ab,CD
  Scrt(1:iSph*jCar*ncd*ijkl) = Wout(1:iSph*jCar*ncd*ijkl)
  call DGeTMO(Scrt,iSph*jCar*ncd,iSph*jCar*ncd,ijkl,Wout,ijkl)
else
  ! Transpose CD,IJKL,a,b to IJKL,ab,CD
  Scrt(1:ncd*ijkl*iCar*jCar) = Win(:)
  if (ncd /= 1) then
    call DGeTMO(Scrt,ncd,ncd,ijkl*iCar*jCar,Wout,ijkl*iCar*jCar)
  else
    Wout(1:ncd*ijkl*iCar*jCar) = Scrt(1:ncd*ijkl*iCar*jCar)
  end if
end if

!call RecPrt(' In CrSph2: (AB|CD)',' ',Wout,ijkl,ncd*mab)

return

end subroutine CrSph2

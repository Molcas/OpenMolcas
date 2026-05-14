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

subroutine SphCr2(Win,ijkl,ncd,Scrt,nScrt,Coeff1,iCar,iSph,Tr1,Pr1,Coeff2,jCar,jSph,Tr2,Pr2,Wout,mab)
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
integer(kind=iwp), intent(in) :: ijkl, ncd, nScrt, iCar, iSph, jCar, jSph, mab
real(kind=wp), intent(in) :: Win(ijkl*ncd*iSph*jSph), Coeff1(iCar,iCar), Coeff2(jCar,jCar)
real(kind=wp), intent(out) :: Scrt(nScrt)
real(kind=wp), intent(inout) :: Wout(ijkl*ncd*mab)
logical(kind=iwp), intent(in) :: Tr1, Pr1, Tr2, Pr2

#include "macros.fh"
unused_var(Pr1)
unused_var(Pr2)

!call RecPrt(' In SphCr2: P(AB|cd) ',' ',Win,ncd*ijkl,iSph*jSph)
if (Tr1 .and. Tr2) then
  !call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
  ! Starting with cd,IJKL,A,B transform to b,cd,IJKL,A
  !call xxDGeMul(Coeff2,jCar,'N',Win,ijkl*ncd*iSph,'T',Scrt,jCar,jCar,jSph,ijkl*ncd*iSph)
  call NTMul(Coeff2,Win,Scrt,jCar,jSph,ijkl*ncd*iSph)

  !call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
  ! Transform b,cd,IJKL,A to ab,cd,IJKL
  !call xxDGeMul(Coeff1,iCar,'N',Scrt,jCar*ncd*ijkl,'T',Wout,iCar,iCar,iSph,jCar*ncd*ijkl)
  call NTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
  ! Transpose ab,cd,IJKL to IJKL,ab,cd
  Scrt(1:mab*ncd*ijkl) = Wout(1:mab*ncd*ijkl)
  call DGeTMO(Scrt,mab*ncd,mab*ncd,ijkl,Wout,ijkl)
else if (Tr2) then
  !call RecPrt(' Right contraction',' ',Coeff2,jCar,jSph)
  ! Starting with cd,IJKL,a,B transform to b,cd,IJKL,a
  !call xxDGeMul(Coeff2,jCar,'N',Win,ncd*ijkl*iCar,'T',Scrt,jCar,jCar,jSph,ncd*ijkl*iCar)
  call NTMul(Coeff2,Win,Scrt,jCar,jSph,ncd*ijkl*iCar)
  ! Transpose b,cd,IJKL,a to IJKL,ab,cd
  call DGeTMO(Scrt,jCar*ncd,jCar*ncd,ijkl*iCar,Wout,ijkl*iCar)
else if (Tr1) then
  ! Transpose cd,IJKL,A,b to b,cd,IJKL,A
  call DGeTMO(Win,ncd*ijkl*iSph,ncd*ijkl*iSph,jCar,Scrt,jCar)

  !call RecPrt(' Left contraction',' ',Coeff1,iCar,iSph)
  ! Transform b,cd,IJKL,A to ab,cd,IJKL
  !call xxDGeMul(Coeff1,iCar,'N',Scrt,jCar*ncd*ijkl,'T',Wout,iCar,iCar,iSph,jCar*ncd*ijkl)
  call NTMul(Coeff1,Scrt,Wout,iCar,iSph,jCar*ncd*ijkl)
  ! Transpose ab,cd,IJKL to IJKL,ab,cd
  Scrt(1:iCar*jCar*ncd*ijkl) = Wout(1:iCar*jCar*ncd*ijkl)
  call DGeTMO(Scrt,iCar*jCar*ncd,iCar*jCar*ncd,ijkl,Wout,ijkl)
else
  ! Transpose cd,IJKL,ab to IJKL,ab,cd
  Scrt(1:ncd*ijkl*iCar*jCar) = Win(1:ncd*ijkl*iCar*jCar)
  if (ncd /= 1) then
    call DGeTMO(Scrt,ncd,ncd,ijkl*iCar*jCar,Wout,ijkl*iCar*jCar)
  else
    Wout(1:ncd*ijkl*iCar*jCar) = Scrt(1:ncd*ijkl*iCar*jCar)
  end if
end if

!call RecPrt(' In SphCr2: P(ab|cd)',' ',Wout,ijkl,ncd*mab)

return

end subroutine SphCr2

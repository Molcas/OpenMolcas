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

subroutine MSP(B,rGamma,Delta,nDim)

use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Constants, only: Two
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(inout) :: B(nDim,nDim)
real(kind=wp), intent(in) :: rGamma(nDim), Delta(nDim)
integer(kind=iwp) :: i
real(kind=wp) :: dd, gd, gg, phi
real(kind=wp), external :: DDot_

!                          T       T            ( T)
!                |(1-phi)/d g phi/d d|        | (g )
!                |     T       T    T         | ( T)
! B = B + (g  d )|phi/d d    -Phi*d g/(d d)**2| (d )

gd = DDot_(nDim,rGamma,1,Delta,1)
dd = DDot_(nDim,Delta,1,Delta,1)
gg = DDot_(nDim,rGamma,1,rGamma,1)
phi = (One-((gd**2)/(dd*gg)))

#ifdef _DEBUGPRINT_
call RecPrt(' MSP: Hessian',' ',B,nDim,nDim)
call RecPrt(' MSP: Delta',' ',Delta,nDim,1)
call RecPrt(' MSP: Gamma',' ',rGamma,nDim,1)
write(u6,*) 'MSP: Phi=',Phi
write(u6,*) 'gd,dd,gg=',gd,dd,gg
write(u6,*) 'MSP: a=',sqrt(Phi)
write(u6,*) 'MSP: E_msp=',(gd/dd)**2*((Two/(One-Phi*sqrt(Phi)))-One)
#endif

do i=1,nDim
  B(:,i) = B(:,i)+((One-phi)/gd)*rGamma(:)*rGamma(i)+phi*((rGamma(:)*Delta(i)+Delta(:)*rGamma(i))/dd-gd*Delta(:)*Delta(i)/dd**2)
end do

#ifdef _DEBUGPRINT_
call RecPrt(' MSP: Updated Hessian',' ',B,nDim,nDim)
#endif

end subroutine MSP

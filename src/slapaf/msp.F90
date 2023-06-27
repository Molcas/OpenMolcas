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

subroutine MSP(B,Gamma,Delta,nDim)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 B(nDim,nDim), gamma(nDim), Delta(nDim)

!                          T       T            ( T)
!                |(1-phi)/d g phi/d d|        | (g )
!                |     T       T    T         | ( T)
! B = B + (g  d )|phi/d d    -Phi*d g/(d d)**2| (d )

iRout = 212
iPrint = nPrint(iRout)

gd = DDot_(nDim,Gamma,1,Delta,1)
dd = DDot_(nDim,Delta,1,Delta,1)
gg = DDot_(nDim,Gamma,1,Gamma,1)
phi = (One-((gd**2)/(dd*gg)))
e_msp = (gd/dd)**2*((Two/(One-Phi*sqrt(Phi)))-One)
if (iPrint >= 99) then
  call RecPrt(' MSP: Hessian',' ',B,nDim,nDim)
  call RecPrt(' MSP: Delta',' ',Delta,nDim,1)
  call RecPrt(' MSP: Gamma',' ',Gamma,nDim,1)
  write(6,*) 'MSP: Phi=',Phi
  write(6,*) 'gd,dd,gg=',gd,dd,gg
  write(6,*) 'MSP: a=',sqrt(Phi)
  write(6,*) 'MSP: E_msp=',E_msp
end if
do i=1,nDim
  do j=1,nDim
    B(i,j) = B(i,j)+((One-phi)/gd)*gamma(i)*gamma(j)+phi*((gamma(i)*Delta(j)+Delta(i)*gamma(j))/dd-gd*Delta(i)*Delta(j)/dd**2)
  end do
end do

if (iPrint >= 99) call RecPrt(' MSP: Updated Hessian',' ',B,nDim,nDim)

return

end subroutine MSP

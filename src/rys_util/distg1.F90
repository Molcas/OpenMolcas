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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Distg1(Temp,mVec,Grad,nGrad,IfGrad,IndGrd,iStab,kOp)

!***********************************************************************
!                                                                      *
! Object: trace the gradient of the ERI's with the second order        *
!         density matrix                                               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iChBas

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
real*8 Grad(nGrad), Temp(9), PAOg1(12), Prmt(0:7)
logical IfGrad(3,4)
integer IndGrd(3,4), kOp(4), iStab(4)
data Prmt/1.d0,-1.d0,-1.d0,1.d0,-1.d0,1.d0,1.d0,-1.d0/
! Statement Function
xPrmt(i,j) = Prmt(iand(i,j))

#ifdef _DEBUGPRINT_
iRout = 239
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt('Accumulated gradient on entrance',' ',Grad,nGrad,1)
  write(6,*) ' kOp=',kOp
  write(6,*) ' iStab=',iStab
  call RecPrt('Distg1: Temp',' ',Temp,9,1)
end if
if (iPrint >= 49) write(6,*) IndGrd
#endif

! Distribute Temp in PAOg1

nVec = 0
#ifdef __INTEL_COMPILER
do kl=1,12
  iCar = (kl-1)/4+1
  iCent = kl-(iCar-1)*4
  ij = 3*(iCent-1)+iCar
  if (IfGrad(iCar,iCent)) then
    nVec = nVec+1
    PAOg1(ij) = Temp(nVec)
  else
    PAOg1(ij) = Zero
  end if
end do
#else

! Original code didn't work for Intel compiler with -O3
! options since it swaps the loops.

do iCar=1,3
  do iCent=1,4
    ij = 3*(iCent-1)+iCar
    if (IfGrad(iCar,iCent)) then
      nVec = nVec+1
      PAOg1(ij) = Temp(nVec)
    else
      PAOg1(ij) = Zero
    end if
  end do
end do
#endif

! a) Compute some of the contributions via the translational invariance
! b) Distribute contribution to the gradient.

do iCn=1,4
  do iCar=1,3
    ij = 3*(iCn-1)+iCar

    ! a)

    if (IndGrd(iCar,iCn) < 0) then
      do jCn=1,4
        if ((iCn /= jCn) .and. IfGrad(iCar,jCn)) then
          kl = 3*(jCn-1)+iCar
          PAOg1(ij) = PAOg1(ij)-PAOg1(kl)
        end if
      end do
    end if

    ! b)

    if (IndGrd(iCar,iCn) /= 0) then
      iGrad = abs(IndGrd(iCar,iCn))
      ! Parity due to integration direction
      ps = xPrmt(kOp(iCn),iChBas(1+iCar))
      Fact = ps*dble(iStab(iCn))/dble(nIrrep)
      Grad(iGrad) = Grad(iGrad)+Fact*PAOg1(ij)
    end if

  end do
end do

#ifdef _DEBUGPRINT_
if (iPrint >= 49) then
  call RecPrt('PAOg1',' ',PAOg1,12,1)
  call RecPrt('Accumulated gradient on exit',' ',Grad,nGrad,1)
end if
#endif

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer(mVec)
#endif

end subroutine Distg1

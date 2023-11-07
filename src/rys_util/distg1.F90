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

subroutine Distg1(Temp,Grad,nGrad,IfGrad,IndGrd,iStab,kOp)
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
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(in) :: Temp(9)
integer(kind=iwp), intent(in) :: nGrad, IndGrd(3,4), iStab(4), kOp(4)
real(kind=wp), intent(inout) :: Grad(nGrad)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp) :: iCar, iCent, iCn, iGrad, ij, jCn, kl, nVec
real(kind=wp) :: Fact, PAOg1(12), ps
real(kind=wp), parameter :: Prmt(0:7) = [One,-One,-One,One,-One,One,One,-One]

#ifdef _DEBUGPRINT_
call RecPrt('Accumulated gradient on entrance',' ',Grad,nGrad,1)
write(u6,*) ' kOp=',kOp
write(u6,*) ' iStab=',iStab
call RecPrt('Distg1: Temp',' ',Temp,9,1)
write(u6,*) IndGrd
#endif

! Distribute Temp in PAOg1

nVec = 0
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
      ps = Prmt(iand(kOp(iCn),iChBas(1+iCar)))
      Fact = ps*real(iStab(iCn),kind=wp)/real(nIrrep,kind=wp)
      Grad(iGrad) = Grad(iGrad)+Fact*PAOg1(ij)
    end if

  end do
end do

#ifdef _DEBUGPRINT_
call RecPrt('PAOg1',' ',PAOg1,12,1)
call RecPrt('Accumulated gradient on exit',' ',Grad,nGrad,1)
#endif

return

end subroutine Distg1

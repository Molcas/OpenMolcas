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

subroutine Distg1X(g1,PAO,nT,mPAO,mVec,Grad,nGrad,IfGrad,IndGrd,iStab,kOp)
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
integer(kind=iwp), intent(in) :: nT, mPAO, mVec, nGrad, IndGrd(3,4), iStab(4), kOp(4)
real(kind=wp), intent(in) :: g1(nT,mPAO,mVec), PAO(nT,mPAO)
real(kind=wp), intent(inout) :: Grad(nGrad)
logical(kind=iwp), intent(in) :: IfGrad(3,4)
integer(kind=iwp) :: iCar, iCent, iCn, iGrad, ij, jCn, kl, nVec
real(kind=wp) :: Fact, PAOg1(12), ps, Temp(9)
#ifdef _DEBUGPRINT_
character(len=80) :: Label
#endif
real(kind=wp), parameter :: Prmt(0:7) = [One,-One,-One,One,-One,One,One,-One]

#ifdef _DEBUGPRINT_
iRout = 239
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt('PAO',' ',PAO,nT,mPAO)
  do iVec=1,mVec
    write(Label,'(A,I2,A)') ' g1(',iVec,')'
    call RecPrt(Label,' ',g1(:,:,iVec),nT,mPAO)
  end do
  call RecPrt('Accumulated gradient on entrance',' ',Grad,nGrad,1)
end if
if (iPrint >= 49) write(u6,*) IndGrd
#endif

! Trace the integral derivatives with the second order density matrix.

call dGeMV_('T',nT*mPAO,mVec,One,g1,nT*mPAO,PAO,1,Zero,Temp,1)
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

! Compute some of the contributions via the translational invariance

do iCn=1,4
  do iCar=1,3
    if (IndGrd(iCar,iCn) < 0) then
      ij = 3*(iCn-1)+iCar
      do jCn=1,4
        if (iCn == jCn) cycle
        if (IfGrad(iCar,jCn)) then
          kl = 3*(jCn-1)+iCar
          PAOg1(ij) = PAOg1(ij)-PAOg1(kl)
        end if
      end do
    end if
  end do
end do
#ifdef _DEBUGPRINT_
if (iPrint >= 49) call RecPrt('PAOg1',' ',PAOg1,12,1)
#endif

! Distribute contribution to the gradient.

do iCn=1,4
  do iCar=1,3
    ij = 3*(iCn-1)+iCar
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
if (iPrint >= 49) then
  call RecPrt('Accumulated gradient on exit',' ',Grad,nGrad,1)
end if
#endif

return

end subroutine Distg1X

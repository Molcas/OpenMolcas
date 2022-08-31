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

subroutine SmAdNa(ArrIn,nb,ArrOut,nop,lOper,IndGrd,iuv,Indx,iDCar,rf,tr)

use Symmetry_Info, only: iChBas, iChTbl, nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nb, nOp(3), lOper, IndGrd(3,4,0:nIrrep-1), iuv(3), Indx(3,4), iDCar
real(kind=wp), intent(in) :: ArrIn(nb,*), rf
real(kind=wp), intent(inout) :: ArrOut(nb,*)
logical(kind=iwp), intent(in) :: tr(4)
integer(kind=iwp) :: i1, i2, iCn, iComp, iIrrep
real(kind=wp) :: Fact, ps
integer(kind=iwp), external :: iPrmt

!iRout = 200
!iPrint = nPrint(iRout)

! Accumulate contributions

iComp = 0
do iIrrep=0,nIrrep-1
  if (btest(lOper,iIrrep)) then
    iComp = iComp+1
    do iCn=1,3
      !if (Indx(idCar,iCn) /= 0) then
      if ((Indgrd(idCar,iCn,iIrrep) /= 0) .and. ((Indx(idcar,icn) > 0) .or. tr(icn))) then
        ! Accumulate contribution to the gradient
        i1 = 0
        i2 = 0
        if (iCn == 1) then
          ps = real(iPrmt(nOp(1),iChBas(1+idCar)),kind=wp)
          Fact = rf*real(iuv(1),kind=wp)/real(nIrrep,kind=wp)
          if (.not. tr(iCn)) then
            i1 = Indx(idCar,iCn)
          else
            if (Indx(idcar,2) > 0) i1 = Indx(idCar,2)
            if (Indx(idCar,3) > 0) i2 = Indx(idCar,3)
            Fact = -Fact
          end if
        else if (iCn == 2) then
          ps = real(iChTbl(iIrrep,nOp(2)),kind=wp)
          ps = ps*real(iPrmt(nOp(2),iChBas(1+idCar)),kind=wp)
          Fact = rf*ps*real(iuv(2),kind=wp)/real(nIrrep,kind=wp)
          if (.not. tr(iCn)) then
            i1 = Indx(idCar,iCn)
          else
            if (Indx(idcar,1) > 0) i1 = Indx(idCar,1)
            if (Indx(idCar,3) > 0) i2 = Indx(idCar,3)
            Fact = -Fact
          end if
        else
          ps = real(iChTbl(iIrrep,nOp(3)),kind=wp)
          ps = ps*real(iPrmt(nOp(3),iChBas(1+idCar)),kind=wp)
          Fact = rf*ps*real(iuv(3),kind=wp)/real(nIrrep,kind=wp)
          if (.not. tr(iCn)) then
            i1 = Indx(idCar,iCn)
          else
            if (Indx(idcar,1) > 0) i1 = Indx(idCar,1)
            if (Indx(idCar,2) > 0) i2 = Indx(idCar,2)
            Fact = -Fact
          end if
        end if
        if (i1 /= 0) ArrOut(:,iComp) = ArrOut(:,iComp)+Fact*ArrIn(:,i1)
        if (i2 /= 0) ArrOut(:,iComp) = ArrOut(:,iComp)+Fact*ArrIn(:,i2)
      end if
    end do
  end if
end do

return

end subroutine SmAdNa

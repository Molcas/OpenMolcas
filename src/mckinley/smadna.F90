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

subroutine SmAdNa(ArrIn,nb,ArrOut,nop,lOper,IndGrd,iuv,IfGrd,Index,iDCar,rf,IFG,tr)

use Symmetry_Info, only: nIrrep, iChTbl, iChBas

implicit real*8(A-H,O-Z)
#include "real.fh"
!#include "print.fh"
real*8 ArrIn(nb,*), ArrOut(nb,*)
integer lOper, IndGrd(3,4,0:nIrrep-1), iuv(3), index(3,4), nOp(3)
logical IfGrd(3,4), IFG(4), tr(4)

!iRout = 200
!iPrint = nPrint(iRout)

! Accumulate contributions

iComp = 0
do iIrrep=0,nIrrep-1
  if (iand(lOper,2**iIrrep) /= 0) then
    iComp = iComp+1
    do iCn=1,3
      !if (Index(idCar,iCn) /= 0) then
      if ((Indgrd(idCar,iCn,iIrrep) /= 0) .and. ((index(idcar,icn) > 0) .or. tr(icn))) then
        ! Accumulate contribution to the gradient
        i1 = 0
        i2 = 0
        if (iCn == 1) then
          ps = dble(iPrmt(nOp(1),iChBas(1+idCar)))
          Fact = rf*dble(iuv(1))/dble(nIrrep)
          if (.not. tr(iCn)) then
            i1 = index(idCar,iCn)
          else
            if (index(idcar,2) > 0) i1 = index(idCar,2)
            if (index(idCar,3) > 0) i2 = index(idCar,3)
            Fact = -Fact
          end if
        else if (iCn == 2) then
          ps = dble(iChTbl(iIrrep,nOp(2)))
          ps = ps*dble(iPrmt(nOp(2),iChBas(1+idCar)))
          Fact = rf*ps*dble(iuv(2))/dble(nIrrep)
          if (.not. tr(iCn)) then
            i1 = index(idCar,iCn)
          else
            if (index(idcar,1) > 0) i1 = index(idCar,1)
            if (index(idCar,3) > 0) i2 = index(idCar,3)
            Fact = -Fact
          end if
        else
          ps = dble(iChTbl(iIrrep,nOp(3)))
          ps = ps*dble(iPrmt(nOp(3),iChBas(1+idCar)))
          Fact = rf*ps*dble(iuv(3))/dble(nIrrep)
          if (.not. tr(iCn)) then
            i1 = index(idCar,iCn)
          else
            if (index(idcar,1) > 0) i1 = index(idCar,1)
            if (index(idCar,2) > 0) i2 = index(idCar,2)
            Fact = -Fact
          end if
        end if
        if (i1 /= 0) call DaXpY_(nb,Fact,ArrIn(1,i1),1,ArrOut(1,iComp),1)
        if (i2 /= 0) call DaXpY_(nb,Fact,ArrIn(1,i2),1,ArrOut(1,iComp),1)
      end if
    end do
  end if
end do

!call GetMem(' Exit SymAdO','LIST','REAL',iDum,iDum)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_logical_array(IfGrd)
  call Unused_logical_array(IFG)
end if

end subroutine SmAdNa

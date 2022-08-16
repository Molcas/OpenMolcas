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

subroutine SymAdO_mck(ArrIn,nB,ArrOut,nrOp,nop,lOper,IndGrd,iu,iv,ifgrd,idCar,trans)

use Symmetry_Info, only: nIrrep, iChTbl, iChBas

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 ArrIn(nB,2), ArrOut(nB,nrOp)
integer IndGrd(0:7), nop(2)
logical IfGrd(3,2), trans(2)

! Accumulate contributions

n = 0
do iIrrep=0,nIrrep-1
  if (Indgrd(iIrrep) /= 0) then
    n = n+1
    do iCn=1,2
      if ((Trans(iCn) .or. IfGrd(idCar,iCn)) .and. (IndGrd(iIrrep) /= 0)) then

        ! Accumulate contribution to the gradient

        if (iCn == 1) then
          ps = dble(iPrmt(nOp(1),iChBas(1+idCar)))
          Fact = dble(iu)/dble(nIrrep)
          if (trans(1)) then
            Fact = -Fact
          end if
        else
          ps = dble(iChTbl(iIrrep,nOp(2)))
          ps = ps*dble(iPrmt(nOp(2),iChBas(1+idCar)))
          Fact = ps*dble(iv)/dble(nIrrep)
          if (trans(2)) then
            Fact = -Fact
          end if
        end if
        call DaXpY_(nB,Fact,ArrIn(1,1),1,ArrOut(1,n),1)
      end if
    end do
  end if
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(lOper)

end subroutine SymAdO_mck

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

subroutine SymAdO_mck2(ArrIn,nB,ArrOut,nrOp,nop,IndGrd,ksym,iu,iv,ifgrd,idCar,trans)

use Symmetry_Info, only: iChBas, iChTbl, iOper, nIrrep
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nB, nrOp, nop(2), IndGrd(0:nIrrep-1), ksym, iu, iv, idCar
real(kind=wp), intent(in) :: ArrIn(nB,2)
real(kind=wp), intent(inout) :: ArrOut(nB,nrOp)
logical(kind=iwp), intent(in) :: IfGrd(3,2), trans(2)
integer(kind=iwp) :: iCn, iirrep, jIrrep, n
real(kind=wp) :: Fact, ps
integer(kind=iwp), external :: iPrmt, NrOpr

! Accumulate contributions

n = 0
do jIrrep=0,nIrrep-1
  iirrep = nropr(ieor(ioper(jirrep),ioper(ksym)))
  if (Indgrd(jIrrep) /= 0) then
    n = n+1
    do iCn=1,2
      if ((Trans(iCn) .or. IfGrd(idCar,iCn)) .and. (IndGrd(jIrrep) /= 0)) then

        ! Accumulate contribution to the gradient

        if (iCn == 1) then
          ps = real(iPrmt(nOp(1),iChBas(1+idCar)),kind=wp)
          Fact = real(iu,kind=wp)/real(nIrrep,kind=wp)
          if (trans(1)) then
            Fact = -Fact
          end if
        else
          ps = real(iChTbl(iIrrep,nOp(2)),kind=wp)
          ps = ps*real(iPrmt(nOp(2),iChBas(1+idCar)),kind=wp)
          Fact = ps*real(iv,kind=wp)/real(nIrrep,kind=wp)
          if (trans(2)) then
            Fact = -Fact
          end if
        end if
        ArrOut(:,n) = ArrOut(:,n)+Fact*ArrIn(:,iCn)
      end if
    end do
  end if
end do

return

end subroutine SymAdO_mck2

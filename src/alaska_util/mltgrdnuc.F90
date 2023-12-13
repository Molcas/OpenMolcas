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

subroutine MltGrdNuc(Grad,nGrad,nOrdOp)

use Basis_Info, only: dbsc, nCnttp
use finfld, only: force
use Index_Functions, only: C_Ind
use Disp, only: Dirct, IndDsp
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nGrad, nOrdOp
real(kind=wp), intent(inout) :: Grad(nGrad)
integer(kind=iwp) :: iCar, icomp, iIrrep, ixop, iyop, izop, kCnt, kCnttp, kdc, ndc, nDisp
real(kind=wp) :: C(3), Fact, ff, XGrad
logical(kind=iwp), external :: TF

iIrrep = 0
do ixop=0,nOrdOp
  do iyop=0,nOrdOp-ixop
    izop = nOrdOp-ixop-iyop
    icomp = C_Ind(nOrdOp,ixop,izop)
    ff = Force(icomp)
    if (ff == Zero) cycle
    kdc = 0
    do kCnttp=1,nCnttp
      if (kCnttp > 1) kdc = kdc+dbsc(kCnttp-1)%nCntr
      if (dbsc(kCnttp)%Charge == Zero) cycle
      do kCnt=1,dbsc(kCnttp)%nCntr
        C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
        ndc = kdc+kCnt
        Fact = -dbsc(kCnttp)%Charge*ff
        nDisp = IndDsp(ndc,iIrrep)
        do iCar=0,2
          iComp = 2**iCar
          if (TF(ndc,iIrrep,iComp) .and. (.not. dbsc(kCnttp)%pChrg)) then
            nDisp = nDisp+1
            if (Dirct(nDisp)) then
              XGrad = Zero
              if (iCar == 0) then
                if (ixop > 0) XGrad = Fact*real(ixop,kind=wp)*C(1)**(ixop-1)*C(2)**iyop*C(3)**izop
              else if (iCar == 1) then
                if (iyop > 0) XGrad = Fact*real(iyop,kind=wp)*C(1)**ixop*C(2)**(iyop-1)*C(3)**izop
              else
                if (izop > 0) XGrad = Fact*real(izop,kind=wp)*C(1)**ixop*C(2)**iyop*C(3)**(izop-1)
              end if
              Grad(nDisp) = Grad(nDisp)+XGrad
            end if
          end if
        end do
      end do
    end do
  end do
end do

return

end subroutine MltGrdNuc

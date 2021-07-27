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

use Basis_Info
use Center_Info

implicit real*8(A-H,O-Z)
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "disp.fh"
real*8 Grad(nGrad), C(3)
#include "finfld.fh"
logical, external :: TF
! Statement function
Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2+iz+1

iIrrep = 0
do ixop=0,nOrdOp
  do iyop=0,nOrdOp-ixop
    izop = nOrdOp-ixop-iyop
    icomp = Ind(nOrdOp,ixop,izop)
    ff = Force(icomp)
    if (ff == 0.d0) goto 801
    kdc = 0
    do kCnttp=1,nCnttp
      if (dbsc(kCnttp)%Charge == 0.d0) Go To 411
      do kCnt=1,dbsc(kCnttp)%nCntr
        C(1:3) = dbsc(kCnttp)%Coor(1:3,kCnt)
        ndc = kdc+kCnt
        Fact = -dbsc(kCnttp)%Charge*ff
        nDisp = IndDsp(ndc,iIrrep)
        do iCar=0,2
          iComp = 2**iCar
          if (TF(ndc,iIrrep,iComp) .and. (.not. dbsc(kCnttp)%pChrg)) then
            nDisp = nDisp+1
            if (Direct(nDisp)) then
              XGrad = 0.d0
              if (iCar == 0) then
                if (ixop > 0) XGrad = Fact*dble(ixop)*C(1)**(ixop-1)*C(2)**iyop*C(3)**izop
              else if (iCar == 1) then
                if (iyop > 0) XGrad = Fact*dble(iyop)*C(1)**ixop*C(2)**(iyop-1)*C(3)**izop
              else
                if (izop > 0) XGrad = Fact*dble(izop)*C(1)**ixop*C(2)**iyop*C(3)**(izop-1)
              end if
              Grad(nDisp) = Grad(nDisp)+XGrad
            end if
          end if
        end do
      end do
411   continue
      kdc = kdc+dbsc(kCnttp)%nCntr
    end do
801 continue
  end do
end do

return

end subroutine MltGrdNuc

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

subroutine DMinvCI_sa(ipSigma,rout,S)

use ipPage, only: ipin, W
use MCLR_Data, only: nConf1, ipCI
use MCLR_Data, only: ipDia
use input_mclr, only: nRoots, ERASSCF, nCSF, State_Sym
use Constants, only: Zero

implicit none
integer ipSigma
real*8 rout(nCSF(State_Sym),nRoots), S(nroots,nroots,nroots)
#include "rasdim.fh"
real*8 rcoeff(mxroot), alpha(mxRoot)
real*8, external :: DDot_
integer k, iR, jR

!                                  -1           -1
!                             (H -E) |0><0|(H -E) |Sigma>
!                -1             0            0
! |rNew> = (H - E) |Sigma> - -----------------------------
!            0                               -1
!                                    <0|(H -E) |0>
!                                         0

if (nconf1 > 1) then
  call ipin(ipdia)
  call ipin(ipsigma)
  k = 0
  do iR=1,nroots
    rout(:,iR) = W(ipSigma)%A(k+1:k+nCSF(State_Sym))/(W(ipdia)%A(1:nCSF(State_Sym))-ERASSCF(iR))
    k = k+nCSF(State_Sym)
  end do
  do iR=1,nroots

    !We = weight(iR)
    call ipin(ipCI)
    do jR=1,nroots
      rcoeff(jR) = ddot_(nconf1,rout(:,iR),1,W(ipCI)%A(1+(jR-1)*nCSF(State_Sym)),1)
    end do

    do jR=1,nroots
      alpha(jR) = sum(S(jR,:,iR)*rcoeff(1:nroots))
    end do

    do jR=1,nroots
      rout(:,iR) = rout(:,iR)- &
                   W(ipCI)%A((jR-1)*nCSF(State_Sym)+1:jR*nCSF(State_Sym))*alpha(jR)/(W(ipdia)%A(1:nCSF(State_Sym))-ERASSCF(iR))
    end do

  end do
else
  rout(:,:) = Zero
end if

end subroutine DMinvCI_sa

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
use MCLR_Data, only: ipCI, ipDia, nConf1
use input_mclr, only: ERASSCF, nCSF, nRoots, State_Sym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipSigma
real(kind=wp), intent(out) :: rout(nCSF(State_Sym),nRoots)
real(kind=wp), intent(in) :: S(nRoots,nRoots,nRoots)
integer(kind=iwp) :: iR, jR, k
real(kind=wp) :: alpha(nRoots), rcoeff(nRoots)
real(kind=wp), external :: DDot_

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
  do iR=1,nRoots
    rout(:,iR) = W(ipSigma)%A(k+1:k+nCSF(State_Sym))/(W(ipdia)%A(1:nCSF(State_Sym))-ERASSCF(iR))
    k = k+nCSF(State_Sym)
  end do
  do iR=1,nRoots

    !We = weight(iR)
    call ipin(ipCI)
    do jR=1,nRoots
      rcoeff(jR) = ddot_(nconf1,rout(:,iR),1,W(ipCI)%A(1+(jR-1)*nCSF(State_Sym)),1)
    end do

    do jR=1,nRoots
      alpha(jR) = sum(S(jR,:,iR)*rcoeff(:))
    end do

    do jR=1,nRoots
      rout(:,iR) = rout(:,iR)- &
                   W(ipCI)%A((jR-1)*nCSF(State_Sym)+1:jR*nCSF(State_Sym))*alpha(jR)/(W(ipdia)%A(1:nCSF(State_Sym))-ERASSCF(iR))
    end do
    !if (abs(Weight(iR)) > 1.0e-9_wp) rout(:,iR) = rout(:,iR)/Weight(iR)

  end do
else
  rout(:,:) = Zero
end if

end subroutine DMinvCI_sa

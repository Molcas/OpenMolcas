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

use CI_Info, only: CI_St
use ipPage, only: ipin, W
use MCLR_Data, only: FANCY_PRECONDITIONER, ipCI, ipDia, nConf1
use input_mclr, only: ERASSCF, nCSF, nRoots, Weight
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: ipSigma
real(kind=wp), intent(out) :: rout(nCSF(CI_St(iState)%STSym),nRoots)
real(kind=wp), intent(in) :: S(nRoots,nRoots,nRoots)
integer(kind=iwp) :: iR, jR, k
real(kind=wp) :: alpha(nRoots), rcoeff(nRoots), den(nCSF(CI_St(1)%STSym),nRoots)
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
  ! Denominators |H0-E_i|: the CI block of the SA Hessian, 2 w_i (H-E_i), is
  ! positive definite on the space orthogonal to the SA states, but H0-E_i
  ! can be negative for excited roots, which made the preconditioner
  ! indefinite. The same denominators are used in SA_PREC2. The simpler
  ! (NODEfault) preconditioner, which only projects against the own state,
  ! keeps the signed denominators of CIDIA_sa.
  do iR=1,nRoots
    den(:,iR) = W(ipdia)%A(1:nCSF(CI_St(iState)%STSym))-ERASSCF(iR)
    if (FANCY_PRECONDITIONER) den(:,iR) = max(abs(den(:,iR)),1.0e-16_wp)
  end do
  k = 0
  do iR=1,nRoots
    rout(:,iR) = W(ipSigma)%A(k+1:k+nCSF(CI_St(iState)%STSym))/den(:,iR)
    k = k+nCSF(CI_St(iState)%STSym)
  end do
  do iR=1,nRoots

    !We = weight(iR)
    call ipin(ipCI)
    do jR=1,nRoots
      rcoeff(jR) = ddot_(nconf1,rout(:,iR),1,W(ipCI)%A(1+(jR-1)*nCSF(CI_St(iState)%STSym)),1)
    end do

    do jR=1,nRoots
      alpha(jR) = sum(S(jR,:,iR)*rcoeff(:))
    end do

    do jR=1,nRoots
      rout(:,iR) = rout(:,iR)-W(ipCI)%A((jR-1)*nCSF(CI_St(iState)%STSym)+1:jR*nCSF(CI_St(iState)%STSym))*alpha(jR)/den(:,iR)
    end do
    ! the CI block of the Hessian is 2 w_i (H-E_i), cf. the factor 1/2 in DMinvCI
    if (FANCY_PRECONDITIONER .and. (Weight(iR) > 1.0e-12_wp)) rout(:,iR) = rout(:,iR)/(Two*Weight(iR))

  end do
else
  rout(:,:) = Zero
end if

end subroutine DMinvCI_sa

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
real*8 rout(*), S(nroots,nroots,nroots)
#include "rasdim.fh"
real*8 rcoeff(mxroot), alpha(mxRoot)
real*8, external :: DDot_
integer i, j, k, iR, jR
real*8 E

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
  do i=1,nroots
    E = ERASSCF(i)
    do j=1,ncsf(state_SYM)
      k = k+1
      rout(k) = W(ipSigma)%A(k)/(W(ipdia)%A(j)-E)
    end do
  end do
  do iR=1,nroots

    !We = weight(iR)
    E = ERASSCF(iR)
    call ipin(ipCI)
    do jR=1,nroots
      rcoeff(jR) = ddot_(nconf1,rout(1+(iR-1)*ncsf(State_Sym)),1,W(ipCI)%A(1+(jR-1)*ncsf(State_Sym)),1)
    end do

    do i=1,nroots
      alpha(i) = Zero
      do j=1,nroots
        alpha(i) = alpha(i)+S(i,j,iR)*rcoeff(j)
      end do
    end do

    do i=1,nroots
      do j=1,ncsf(State_Sym)
        rout(j+(iR-1)*ncsf(State_Sym)) = rout(j+(iR-1)*ncsf(State_Sym))- &
                                         W(ipCI)%A(j+(i-1)*ncsf(State_Sym))*alpha(i)/(W(ipdia)%A(j)-E)
      end do
    end do

  end do
else
  rout(1:nconf1*nroots) = Zero
end if

end subroutine DMinvCI_sa

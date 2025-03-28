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

subroutine DMinvCI(ipSigma,rout,rC_HE_C,idsym)

use Exp, only: NewPre
use ipPage, only: W
use MCLR_Data, only: ngp
use Constants, only: Zero, One, Half
use MCLR_Data, only: nConf1, ipCI
use MCLR_Data, only: ipDia

implicit none
integer ipSIgma, idSym
real*8 rout(*), rC_HE_C
real*8 rCoeff
real*8, external :: DDot_

!                                  -1           -1
!                             (H -E) |0><0|(H -E) |Sigma>
!                -1             0            0
! |rNew> = (H - E) |Sigma> - -----------------------------
!            0                               -1
!                                    <0|(H -E) |0>
!                                         0

if (nconf1 > 1) then

  call ipin(ipdia)
  call ipin(ipSigma)
  call exphinvv(W(ipdia)%Vec,W(ipsigma)%Vec,rout,Zero,One)
  call ipout(ipsigma)
  call opout(ipdia)

  ! OBS <0|(H-E)|Sigma>=0 if idsym=/=1

  if (NewPre .and. (idsym == 1)) then
    !                  -1
    ! rcoeff = <0|(H -E) |Sigma>
    !               0
    !         -------------------
    !                    -1
    !            <0|(H -E) |0>
    !                 0

    if (.not. ngp) then
      call ipin(ipCI)
      rcoeff = ddot_(nconf1,rout,1,W(ipCI)%Vec,1)/rC_HE_C

      !                          -1
      ! rout = rout-rocoeff*(H -E) |0>
      !                       0
      call ipin(ipdia)
      call exphinvv(W(ipdia)%Vec,W(ipci)%Vec,rOUT,One,-rcoeff)
      call opout(ipCI)
    else
      call NEGP(ipdia,ipSigma,rout)
    end if

  end if

  call DSCAL_(nconf1,Half,rout,1)

else

  call ipin(ipsigma)
  rout(1:nConf1) = W(ipSigma)%Vec(:)

end if

end subroutine DMinvCI

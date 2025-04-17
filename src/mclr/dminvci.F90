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

use ipPage, only: ipin, ipout, opout, W
use MCLR_Data, only: ipCI, ipDia, nConf1, NewPre, ngp
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ipSIgma, idSym
real(kind=wp), intent(inout) :: rout(nConf1)
real(kind=wp), intent(in) :: rC_HE_C
real(kind=wp) :: rCoeff
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
  call ipin(ipSigma)
  call exphinvv(W(ipdia)%A,W(ipsigma)%A,rout,Zero,One)
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
      rcoeff = ddot_(nconf1,rout,1,W(ipCI)%A,1)/rC_HE_C

      !                          -1
      ! rout = rout-rocoeff*(H -E) |0>
      !                       0
      call ipin(ipdia)
      call exphinvv(W(ipdia)%A,W(ipci)%A,rOUT,One,-rcoeff)
      call opout(ipCI)
    else
      call NEGP(ipdia,ipSigma,rout)
    end if

  end if

  rout(:) = Half*rout(:)

else

  call ipin(ipsigma)
  rout(:) = W(ipSigma)%A(:)

end if

end subroutine DMinvCI

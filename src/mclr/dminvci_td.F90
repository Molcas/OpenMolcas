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

subroutine DMinvCI_td(rin,rout,rome,idsym)

use ipPage, only: ipin, W
use MCLR_Data, only: ipCI, ipDia, nConf1
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rin(nConf1), rome
real(kind=wp), intent(out) :: rout(nConf1)
integer(kind=iwp), intent(in) :: idSym
real(kind=wp) :: r1, r2
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
  rout(:) = rin(:)/(W(ipdia)%A(1:nConf1)+rome)

  ! To ensure orthogonal response if response is in same symmetry as wavefunction

  if (idsym == 1) then
    call ipin(ipCI)
    r1 = ddot_(nConf1,W(ipCI)%A,1,rout,1)

    call ipin(ipDia)
    r2 = sum(W(ipCI)%A(1:nConf1)**2/(W(ipDia)%A(1:nConf1)+rome))

    rout(:) = rout(:)-r1/r2*W(ipCI)%A(1:nConf1)/(W(ipDia)%A(1:nConf1)+rome)
  end if

else

  rout(:) = rin(:)

end if

rout(:) = Half*rout(:)

end subroutine DMinvCI_td

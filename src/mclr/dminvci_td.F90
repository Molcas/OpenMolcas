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
use Constants, only: Zero, Half
use MCLR_Data, only: nConf1, ipCI
use MCLR_Data, only: ipDia

implicit none
integer idSym
real*8 rout(*), rin(*), rome
integer i
real*8 r1, r2
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
  do i=1,nconf1
    rout(i) = rin(i)/(W(ipdia)%A(i)+rome)
  end do

  ! To ensure orthogonal response if response is in same symmetry as wavefunction

  if (idsym == 1) then
    call ipin(ipCI)
    r1 = ddot_(nconf1,W(ipCI)%A,1,rout,1)

    r2 = Zero
    call ipin(ipDia)
    do i=1,nconf1
      r2 = r2+W(ipCI)%A(i)**2/(W(ipDia)%A(i)+rome)
    end do

    do i=1,nconf1
      rout(i) = rout(i)-r1/r2*W(ipCI)%A(i)/(W(ipDia)%A(i)+rome)
    end do
  end if

else

  rout(1:nConf1) = rin(1:nConf1)

end if

rout(1:nConf1) = Half*rout(1:nConf1)

end subroutine DMinvCI_td

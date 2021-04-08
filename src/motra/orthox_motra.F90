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

subroutine ORTHOX_MOTRA(S,C,NORB,NBAS)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NORB, NBAS
real(kind=wp), intent(inout) :: S(NORB,NORB), C(NBAS,NORB)
integer(kind=iwp) :: IBAS, IORB, JORB, KORB
real(kind=wp) :: A, F

do IORB=1,NORB
  F = One/sqrt(S(IORB,IORB))
  do IBAS=1,NBAS
    C(IBAS,IORB) = F*C(IBAS,IORB)
  end do
  do JORB=1,NORB
    S(IORB,JORB) = F*S(IORB,JORB)
    S(JORB,IORB) = F*S(JORB,IORB)
  end do
  do JORB=IORB+1,NORB
    A = S(IORB,JORB)
    do IBAS=1,NBAS
      C(IBAS,JORB) = C(IBAS,JORB)-A*C(IBAS,IORB)
    end do
    do KORB=1,NORB
      S(JORB,KORB) = S(JORB,KORB)-A*S(IORB,KORB)
    end do
    do KORB=1,NORB
      S(KORB,JORB) = S(KORB,JORB)-A*S(KORB,IORB)
    end do
  end do
end do

return

end subroutine ORTHOX_MOTRA

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

subroutine ASIND(IAS,ISYM,ICASE,IP,IQ,IR)

use SUPERINDEX, only: MAGEB, MAGTB
use caspt2_module, only: NAGEBES, NAGTBES, IEXTIS
use definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IAS, ISYM, ICASE
integer(kind=iwp), intent(out) :: IP, IQ, IR
integer(kind=iwp) :: IABABS, IAABS, IBABS

if (ICASE == 2) then
  IABABS = IAS+NAGTBES(ISYM)
  IAABS = MAGTB(1,IABABS)
  IBABS = MAGTB(2,IABABS)
else
  IABABS = IAS+NAGEBES(ISYM)
  IAABS = MAGEB(1,IABABS)
  IBABS = MAGEB(2,IABABS)
end if

IP = IEXTIS(IAABS)
IQ = IEXTIS(IBABS)
IR = 0

end subroutine ASIND

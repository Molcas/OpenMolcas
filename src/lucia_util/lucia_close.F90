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

subroutine Lucia_Close()

use lucia_data, only: LUC, LUDIA, LUHC, LUMOUT, LUSC1, LUSC2, LUSC3, LUSC34, LUSC35, LUSC36, LUSC37, LUSC38, LUSC39, LUSC40

implicit none

! Free memory allocated by Lucia

call FREESTR_GAS()
call DeAlloc_Lucia()

! Close any files opened by Lucia

call DAClos(LUDIA)
call DAClos(LUC)
call DAClos(LUHC)
call DAClos(LUSC1)
call DAClos(LUSC2)
call DAClos(LUSC3)
call DAClos(LUSC34)
call DAClos(LUSC35)
call DAClos(LUSC36)
call DAClos(LUSC37)
call DAClos(LUSC38)
call DAClos(LUSC39)
call DAClos(LUSC40)
call DAClos(LUMOUT)

end subroutine Lucia_Close

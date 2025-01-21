!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2023, Roland Lindh                                     *
!***********************************************************************

subroutine Set_Breit(n)

use Breit, only: nOrdOp, nComp, Do_BP_Integrals, PSO, PAO
use Definitions, only: iwp, u6

integer(kind=iwp), intent(in) :: n

nOrdOp = n
if (nOrdOp == 0) then
  nComp = 1
  Do_BP_Integrals=.False.
  PSO=>Null()
  PAO=>Null()
else if (nOrdOp == 1) then
  nComp = 6
  Do_BP_Integrals=.True.
else
  Write(u6,*) 'Set_Breit: Illegal n value, n:',n
  Call abend()
end if

end subroutine Set_Breit


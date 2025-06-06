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

module timers

use Constants, only: Zero
use Definitions, only: wp

implicit none
private

real(kind=wp) :: C_Dress = Zero, C_get_Cm = Zero, TimeAoMo = Zero, TimeCIOpt = Zero, TimeDavid = Zero, TimeDens = Zero, &
                 TimeFock = Zero, TimeHCSCE = Zero, TimeHDiag = Zero, TimeHSel = Zero, TimeInput = Zero, TimeOrb = Zero, &
                 TimeOutput = Zero, TimePage = Zero, TimeRelax = Zero, TimeSigma = Zero, TimeTotal = Zero, TimeTrans = Zero, &
                 TimeWfn = Zero, W_Dress = Zero, W_get_Cm = Zero

public :: C_Dress, C_get_Cm, TimeAoMo, TimeCIOpt, TimeDavid, TimeDens, TimeFock, TimeHCSCE, TimeHDiag, TimeHSel, TimeInput, &
          TimeOrb, TimeOutput, TimePage, TimeRelax, TimeSigma, TimeTotal, TimeTrans, TimeWfn, W_Dress, W_get_Cm

end module timers

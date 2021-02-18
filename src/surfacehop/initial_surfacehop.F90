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

subroutine initial_surfacehop()

use Tully_variables, only: tullyL, decoherence, DECO, Ethreshold, RandThreshold, tullySubVerb, fixedrandL, FixedRand, NSUBSTEPS
use Surfacehop_globals, only: lH5Restart
use Constants, only: Zero, One, auTofs
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: Found
real(kind=wp) :: DT

tullyL = .false.
decoherence = .false.
DECO = Zero
Ethreshold = huge(Ethreshold)
RandThreshold = Zero
tullySubVerb = .false.
fixedrandL = .false.
FixedRand = -One
lH5Restart = .false.

call qpg_dscalar('Timestep',Found)
if (Found) then
  call Get_dScalar('Timestep',DT)
  NSUBSTEPS = int(200*DT*auTofs)
else
  NSUBSTEPS = 0
end if

return

end subroutine initial_surfacehop

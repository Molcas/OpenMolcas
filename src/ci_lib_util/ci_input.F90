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

!> The CI space as read from the input or the wave-function file (spin,
!> active electrons, state symmetry, RAS and GAS description). Only the
!> input readers and the CI setups use it; after the setup the programs
!> read the state from CI_Lib_Info (CI_Lib_Space(iState)).
module CI_Input

use Molcas, only: MxSym, MxGas
use CI_Space_Def, only: CI_Space_Input
use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: ISPIN, NACTEL, NELEC3, NHOLE1, NLEV, NRS1(mxSym), NRS1T, NRS2(mxSym), NRS2T, NRS3(mxSym), NRS3T, STSYM
integer(kind=iwp) :: IGSOCCX(mxGAS,2) = 0, NGAS = 0, NGSSH(mxGAS,mxSym) = 0
logical(kind=iwp) :: iDoGas = .false.
! Explicit CI Hamiltonian of the CI start from SGUGA (keyword SGEX)
logical(kind=iwp) :: DoSGEx = .false.
! The HEXS/DEXS requests of the input (keywords HEXS, DEXS), from which
! CI_Space_Setup builds the removals of the CI space
type(CI_Space_Input) :: HEXS_DEXS_Inp

public :: DoSGEx, HEXS_DEXS_Inp, iDoGas, IGSOCCX, ISPIN, NACTEL, NELEC3, NGAS, NGSSH, NHOLE1, NLEV, NRS1, NRS1T, NRS2, NRS2T, &
  NRS3, NRS3T, STSYM

end module CI_Input

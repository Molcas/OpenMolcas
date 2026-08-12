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

module rasscf_files

! common logical unit numbers
!
! LUStartOrb : MO-coefficients and occupation numbers
!              (formatted ASCI file, input)
! JOBIPH     : MO-coefficients and occupation numbers etc.
!              (binary, output)
! JOBOLD     : MO-coefficients and occupation numbers etc.
!              (binary, input)
! LUONEL     : one-electron integrals in AO basis
!              (binary, input)
! LUINTA     : two-electron integrals in AO basis
!              (binary, input)
! LUINTM     : two-electron integrals in MO basis
!              (binary, temporary)
! LUQUNE     : orbital gradients
!              (binary, temporary)
! LUDAVID    : Intermediate results of the diagonalization
!              (binary, temporary)

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: ITERFILE, JOBIPH, JOBOLD, LUDAVID, LUINTA, LUINTM, LUONEL, LUQUNE, LUStartOrb
character(len=256) :: StartOrbFile
public :: ITERFILE, JOBIPH, JOBOLD, LUDAVID, LUINTA, LUINTM, LUONEL, LUQUNE, &
          LUStartOrb, StartOrbFile

end module rasscf_files

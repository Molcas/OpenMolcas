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

module general_data

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

use Molcas, only: MxSym
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: MAXALTER = 16
integer(kind=iwp) :: INVEC, ISPIN, ITERFILE, JOBIPH, JOBOLD, LUDAVID, LUINTA, LUINTM, LUONEL, LUQUNE, LUStartOrb, &
                     MALTER(MAXALTER,3), NACTEL, NALTER, NASH(mxSym), NBAS(mxSym), NCONF, NCRVEC, NDEL(mxSym), NDELT, NELEC3, &
                     NFRO(mxSym), NFROT, NHOLE1, NISH(mxSym), NORB(mxSym), NRS1(mxSym), NRS1T, NRS2(mxSym), NRS2T, NRS3(mxSym), &
                     NRS3T, NSEL, NSSH(mxSym), NSYM, NTOT, NTOT1, NTOT2, NTOTSP, STSYM
real(kind=wp) :: SXDAMP
character(len=256) :: StartOrbFile
logical(kind=iwp) :: Lowdin_ON
integer(kind=iwp), allocatable :: CleanMask(:)
real(kind=wp), allocatable :: CRPROJ(:), CRVEC(:)

public :: CleanMask, CRPROJ, CRVEC, INVEC, ISPIN, ITERFILE, JOBIPH, JOBOLD, Lowdin_ON, LUDAVID, LUINTA, LUINTM, LUONEL, LUQUNE, &
          LUStartOrb, MALTER, MAXALTER, NACTEL, NALTER, NASH, NBAS, NCONF, NCRVEC, NDEL, NDELT, NELEC3, NFRO, NFROT, NHOLE1, NISH, &
          NORB, NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T, NSEL, NSSH, NSYM, NTOT, NTOT1, NTOT2, NTOTSP, StartOrbFile, STSYM, SXDAMP

end module general_data

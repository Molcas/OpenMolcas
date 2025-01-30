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

subroutine MKNSM()
! PURPOSE: CREATE THE SYMMETRY INDEX VECTOR

use gugx, only: SGS
use gas_data, only: NGAS, NGSSH
use rasscf_global, only: NSM
use general_data, only: NSYM
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IGAS, ISYM, NLEV, NSTA

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do

if (SGS%nSym /= 0) then
  SGS%nLev = nLev
  call mma_allocate(SGS%ISM,nLev,Label='SGS%ISM')
  SGS%ISM(1:nLev) = NSM(1:nLev)
end if

end subroutine MKNSM


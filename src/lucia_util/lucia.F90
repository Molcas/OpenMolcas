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

subroutine LUCIA()

use stdalloc, only: mma_allocate
use GLBBAS, only: CI_VEC, SIGMA_VEC
use lucia_data, only: MXSOOB, XISPSM
use lucia_data, only: IPRCIX, IPRORB, IPRSTR
use lucia_data, only: NOINT, LCSBLK
use lucia_data, only: IREFSM, PSSIGN
use Constants, only: Zero, Two
use Definitions, only: u6

implicit none
! Parameters for dimensioning
#include "warnings.h"
!Scratch : A character line
integer LBLOCK

! No floating point underflow
!call XUFLOW()
! Assign diskunits
!if (ENVIRO(1:6) == 'RASSCF') then
call DISKUN2()
!else
!  call DISKUN()
!end if
! Header
! From shells to orbitals
call ORBINF(IPRORB)
! Number of string types
call STRTYP_GAS(IPRSTR)
! Divide orbital spaces into inactive/active/secondary
call GASSPC()
! Symmetry information
call SYMINF_LUCIA(IPRORB)
! Number of integrals
call INTDIM(IPRORB)
! Static memory, initialization and  allocation
!if (ENVIRO(1:6) /= 'RASSCF') then
!  KBASE = 1
!  KADD = MXPWRD
!  call MEMMAN(KBASE,KADD,'INI   ',IDUMMY,'DUMMY')
!end if
call ALLOC_LUCIA()
! Read in integrals
if (NOINT == 0) then
  call INTIM()
else
  write(u6,*) ' No integrals imported'
end if
! READ in MO-AO matrix
!if (NOMOFL == 0) call GET_CMOAO(MOAOIN)
! Internal string information
call STRINF_GAS(IPRSTR)
! Internal subspaces
call LCISPC(IPRCIX)

! Symmetry of reference
!if (PNTGRP > 1) call MLSM(IREFSM,IREFPA,IREFSM,'CI',1)
if (NOINT == 1) then
  write(u6,*) ' End of calculation without integrals'
  call QUIT(_RC_ALL_IS_WELL_)
end if

LBLOCK = MXSOOB
LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
!if (ENVIRO(1:6) == 'RASSCF') then
LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

call mma_allocate(CI_VEC,LBLOCK,Label='CI_VEC')
call mma_allocate(SIGMA_VEC,LBLOCK,Label='SIGMA_VEC')

end subroutine LUCIA

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

subroutine MV7(C,HC,LUC,LUHC)
! Outer routine for sigma vector generation
! GAS version !!!!
!
! Written in terms of RASG3/SBLOCK, May 1997

use stdalloc, only: mma_allocate, mma_deallocate
use strbas, only: NSTSO
! Definition of c and sigma
use CandS, only: ISSPC, ISSM
use lucia_data, only: MXNTTS, MXSOOB, ISMOST, XISPSM
use lucia_data, only: ICISTR, ENVIRO, ISIMSYM, LCSBLK
use lucia_data, only: IDC, IREFSM, PSSIGN
use lucia_data, only: I_AM_OUT, N_ELIMINATED_BATCHES
use lucia_data, only: NOCTYP
use csm_data, only: NSMST
use Constants, only: Zero
use Definitions, only: u6

implicit none
integer LUC, LUHC
real*8 C(*), HC(*)
integer, allocatable :: SIOIO(:)
integer, allocatable :: SVST(:)
! this is a the same structure as for local_arrays but it can not be
! used since lower level routines will use it.
integer, allocatable :: CBLTP(:), CLBT(:), CLEBT(:), CI1BT(:), CIBT(:)
integer IATP, IBTP, NOCTPA, NOCTPB, NTTS, LBLOCK, LLUC, LLUHC, NBATCH

if (ICISTR == 1) then
  write(u6,*) ' MV7 does not work for ICISTR = 1'
  write(u6,*) ' SWITCH to ICISTR = 2,3 or program'
  !stop ' MV7 does not work for ICISTR = 1'
  call SYSABENDMSG('lucia_util/mv7','Internal error','')
end if

IATP = 1
IBTP = 2

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Arrays giving allowed type combinations
call mma_allocate(SIOIO,NOCTPA*NOCTPB,Label='SIOIO')
call IAIBCM(ISSPC,SIOIO)
! Arrays for additional symmetry operation
if ((IDC == 3) .or. (IDC == 4)) then
  call mma_allocate(SVST,NSMST,Label='SVST')
  call SIGVST(SVST,NSMST)
else
  call mma_allocate(SVST,1,Label='SVST')
end if
! Arrays giving block type
call mma_allocate(CBLTP,NSMST,Label='CBLTP')
call ZBLTP(ISMOST(1,ISSM),NSMST,IDC,CBLTP,SVST)
call mma_deallocate(SVST)
! Arrays for partitioning of sigma
NTTS = MXNTTS
call mma_allocate(CLBT,NTTS,Label='CLBT')
call mma_allocate(CLEBT,NTTS,Label='CLEBT')
call mma_allocate(CI1BT,NTTS,Label='CI1BT')
call mma_allocate(CIBT,8*NTTS,Label='CIBT')
! Batches  of C vector
!if (ISIMSYM == 0) then
LBLOCK = MXSOOB
!else
!  LBLOCK = MXSOOB_AS
!end if
LBLOCK = max(LBLOCK,LCSBLK)
! JESPER : Should reduce I/O
if (ENVIRO(1:6) == 'RASSCF') then
  LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
  if (PSSIGN /= Zero) LBLOCK = int(2*XISPSM(IREFSM,1))
end if
call PART_CIV2(IDC,CBLTP,NSTSO(IATP)%I,NSTSO(IBTP)%I,NOCTPA,NOCTPB,NSMST,LBLOCK,SIOIO,ISMOST(1,ISSM),NBATCH,CLBT,CLEBT,CI1BT,CIBT, &
               0,ISIMSYM)
call mma_deallocate(SIOIO)
call mma_deallocate(CBLTP)

if (ICISTR == 1) then
  LLUC = 0
  LLUHC = 0
else
  LLUC = LUC
  LLUHC = LUHC
end if

call RASSG3(C,HC,NBATCH,CLBT,CI1BT,CIBT,LLUC,LLUHC,I_AM_OUT,N_ELIMINATED_BATCHES)
!write(u6,*) ' LSCMAX_MX = ',LSCMAX_MX
! Eliminate local memory
call mma_deallocate(CLBT)
call mma_deallocate(CLEBT)
call mma_deallocate(CI1BT)
call mma_deallocate(CIBT)

end subroutine MV7

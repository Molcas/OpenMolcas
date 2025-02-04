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
! Copyright (C) 1994,1995,1999, Jeppe Olsen                            *
!***********************************************************************

subroutine LCISPC(IPRNT)
! Number of dets and combinations
! per symmetry for each type of internal space
!
! Jeppe Olsen, Winter 1994/1995 ( woops !)
!              MXSOOB_AS added,MXSB removed May 1999
!
! GAS VERSION

use stdalloc, only: mma_allocate, mma_deallocate
use strbas
use lucia_data, only: NCMBSPC, IGSOCCX, NGAS
use lucia_data, only: NICISP, MXSB, MXSOOB, MXSOOB_AS, MXNTTS, ISMOST, LCOLIC, NBLKIC, XISPSM
use lucia_data, only: IDC
use lucia_data, only: IBSPGPFTP, ISPGPFTP
use lucia_data, only: NOCTYP
use lucia_data, only: MXPCSM, MXPNGAS
use csm_data, only: NSMST, NSMCI
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer IPRNT
integer, allocatable :: LBLTP(:), LIOIO(:), CVST(:)
integer IATP, IBTP, NOCTPA, NOCTPB, ICI, ISYM, NTTSBL, LCOL, ISM, MXS, MXSOO, MXSOO_AS, NCOMB
real*8 XNCOMB
#ifdef _DEBUGPRINT_
integer NTEST, II
#endif

! Number of spaces
NICISP = NCMBSPC
!write(u6,*) ' LCISPC : NICISP ',NICISP
! Type of alpha- and beta strings
IATP = 1
IBTP = 2

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Local memory
call mma_allocate(LBLTP,NSMST,Label='LBLTP')
call mma_allocate(CVST,NSMST,Label='CVST')
call mma_allocate(LIOIO,NOCTPA*NOCTPB,Label='LIOIO')
! Obtain array giving symmetry of sigma v reflection times string
! symmetry.
!if ((IDC == 3) .or. (IDC == 4)) call SIGVST(CVST,NSMST)

! Array defining symmetry combinations of internal strings
! Number of internal dets for each symmetry
call SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
! MXSB is not calculated anymore, set to 0
MXSB = 0

MXSOOB = 0
MXSOOB_AS = 0
do ICI=1,NICISP
  ! allowed combination of types
  call IAIBCM(ICI,LIOIO)

  do ISYM=1,NSMCI
    call ZBLTP(ISMOST(1,ISYM),NSMST,IDC,LBLTP,CVST)
    call NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),NGAS,ISYM,NSMST,NOCTPA,NOCTPB,NSTSO(IATP)%I,NSTSO(IBTP)%I, &
                ISPGPFTP(1,IBSPGPFTP(IATP)),ISPGPFTP(1,IBSPGPFTP(IBTP)),MXPNGAS,NCOMB,XNCOMB,MXS,MXSOO,LBLTP,NTTSBL,LCOL,LIOIO, &
                MXSOO_AS)

    XISPSM(ISYM,ICI) = XNCOMB
    MXSOOB = max(MXSOOB,MXSOO)
    MXSB = max(MXSB,MXS)
    MXSOOB_AS = max(MXSOO_AS,MXSOOB_AS)
    NBLKIC(ISYM,ICI) = NTTSBL
    LCOLIC(ISYM,ICI) = LCOL
  end do
  call mma_deallocate(LBLTP)
  call mma_deallocate(CVST)
  call mma_deallocate(LIOIO)
end do

#ifdef _DEBUGPRINT_
NTEST = 0
NTEST = max(NTEST,IPRNT)
if (NTEST >= 5) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' Number of internal combinations per symmetry'
  write(u6,*) ' ==========================================='

  do ICI=1,NCMBSPC
    write(u6,*) ' CI space ',ICI
    write(u6,'(1X, 4ES22.15)') (XISPSM(II,ICI),II=1,NSMCI)
    !call WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
  end do
  write(u6,*)
  write(u6,*) ' Largest Symmetry-type-type block ',MXSOOB
  write(u6,*) ' Largest type-type block (all symmetries) ',MXSOOB_AS
  write(u6,*)

  write(u6,*) ' Number of TTS subblocks per CI expansion'
  write(u6,*) ' ========================================'

  do ICI=1,NCMBSPC
    write(u6,*) ' Internal CI space ',ICI
    call IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
  end do
end if
#else
call Unused_Integer(IPRNT)
#endif
! Largest number of BLOCKS in a CI expansion
MXNTTS = 0
do ICI=1,NCMBSPC
  do ISM=1,NSMCI
    MXNTTS = max(MXNTTS,NBLKIC(ISM,ICI))
  end do
end do

#ifdef _DEBUGPRINT_
if (NTEST >= 5) then
  write(u6,*) ' Largest number of blocks in CI expansion',MXNTTS

  write(u6,*) ' Number of columns per CI expansion'
  write(u6,*) ' =================================='

  do ICI=1,NCMBSPC
    write(u6,*) ' Internal CI space ',ICI
    call IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
  end do
end if
#endif

end subroutine LCISPC

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

!#define _DEBUGPRINT_
subroutine LCISPC()
! Number of dets and combinations
! per symmetry for each type of internal space
!
! Jeppe Olsen, Winter 1994/1995 (woops!)
!              MXSOOB_AS added, MXSB removed May 1999
!
! GAS VERSION

use lucia_data, only: IDC, MXNTTS, MXSOOB, NCMBSPC, NIRREP, NOCTYP, NSTSO, XISPSM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: IATP, IBTP, ICI, ISYM, LCOL, MXSOO, MXSOO_AS, NCOMB, NOCTPA, NOCTPB, NTTSBL
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: II
#endif
real(kind=wp) :: XNCOMB
integer(kind=iwp), allocatable :: CVST(:), LBLTP(:), LIOIO(:), NBLKIC(:,:)

!write(u6,*) ' LCISPC : NCMBSPC ',NCMBSPC
! Type of alpha- and beta strings
IATP = 1
IBTP = 2

NOCTPA = NOCTYP(IATP)
NOCTPB = NOCTYP(IBTP)
! Local memory
call mma_allocate(LBLTP,NIRREP,Label='LBLTP')
call mma_allocate(CVST,NIRREP,Label='CVST')
call mma_allocate(LIOIO,NOCTPA*NOCTPB,Label='LIOIO')
call mma_allocate(NBLKIC,NIRREP,NCMBSPC,Label='NBLKIC')
! Obtain array giving symmetry of sigma v reflection times string
! symmetry.
!if ((IDC == 3) .or. (IDC == 4)) call SIGVST(CVST,NIRREP)

! Number of internal dets for each symmetry

MXSOOB = 0
do ICI=1,NCMBSPC
  ! allowed combination of types
  call IAIBCM(ICI,LIOIO)

  do ISYM=1,NIRREP
    call ZBLTP(ISYM,NIRREP,IDC,LBLTP,CVST)
    call NGASDT(ISYM,NIRREP,NOCTPA,NOCTPB,NSTSO(IATP)%A,NSTSO(IBTP)%A,NCOMB,XNCOMB,MXSOO,LBLTP,NTTSBL,LCOL,LIOIO,MXSOO_AS)

    XISPSM(ISYM,ICI) = XNCOMB
    MXSOOB = max(MXSOOB,MXSOO)
    !MXSOOB_AS = max(MXSOO_AS,MXSOOB_AS)
    NBLKIC(ISYM,ICI) = NTTSBL
    !LCOLIC(ISYM,ICI) = LCOL
  end do
end do
call mma_deallocate(LBLTP)
call mma_deallocate(CVST)
call mma_deallocate(LIOIO)

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*)
write(u6,*) ' Number of internal combinations per symmetry'
write(u6,*) ' ============================================'

do ICI=1,NCMBSPC
  write(u6,*) ' CI space ',ICI
  write(u6,'(1X, 4ES22.15)') (XISPSM(II,ICI),II=1,NIRREP)
  !call WRTMAT(XISPSM(1,ICI),1,NIRREP,1,NIRREP)
end do
write(u6,*)
write(u6,*) ' Largest Symmetry-type-type block ',MXSOOB
!write(u6,*) ' Largest type-type block (all symmetries) ',MXSOOB_AS
write(u6,*)

write(u6,*) ' Number of TTS subblocks per CI expansion'
write(u6,*) ' ========================================'

do ICI=1,NCMBSPC
  write(u6,*) ' Internal CI space ',ICI
  call IWRTMA(NBLKIC(:,ICI),1,NIRREP,1,NIRREP)
end do
#endif
! Largest number of BLOCKS in a CI expansion
MXNTTS = max(0,maxval(NBLKIC(:,:)))

#ifdef _DEBUGPRINT_
write(u6,*) ' Largest number of blocks in CI expansion',MXNTTS

write(u6,*) ' Number of columns per CI expansion'
write(u6,*) ' =================================='

!do ICI=1,NCMBSPC
!  write(u6,*) ' Internal CI space ',ICI
!  call IWRTMA(LCOLIC(:,ICI),1,NIRREP,1,NIRREP)
!end do
#endif

call mma_deallocate(NBLKIC)

end subroutine LCISPC

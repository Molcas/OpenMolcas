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

subroutine FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,NO_ZEROING)
! TRANSFER ARRAY FROM DISC FILE IFILE
!
! Version allowing zero and packed blocks
!
! If NO_ZEROING = 1, the elements of zero blocks
! are not set to zero, the routine just returns with
! IMZERO = 1

use lucia_data, only: IDISK
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NDIM, MBLOCK, IFILE, NO_ZEROING
real(kind=wp), intent(out) :: ARRAY(NDIM)
integer(kind=iwp), intent(out) :: IMZERO
integer(kind=iwp), intent(inout) :: I_AM_PACKED
integer(kind=iwp), parameter :: LPBLK = 50000
integer(kind=iwp) :: IBASE, IDUMMY(1), IELMNT, IPACK, IPAK(LPBLK), IREST, ISCR(2), ISTOP, LBATCH, LBATCHP, NBATCH, NBLOCK
real(kind=wp) :: XPAK(LPBLK)

IMZERO = 0

IPACK = 1
if (IPACK /= 0) then
  ! Read if ARRAY is zero
  !MMBLOCK = MBLOCK
  !if (MMBLOCK >= 2) MMBLOCK = 2
  !    IFRMDS(ISCR,2,MMBLOCK,IFILE)
  call IFRMDS(ISCR,2,2,IFILE)
  IMZERO = ISCR(1)
  I_AM_PACKED = ISCR(2)
  if (IMZERO == 1) then
    if (NO_ZEROING == 0) ARRAY(:) = Zero
    return
  end if
end if

if (I_AM_PACKED == 1) then
  ARRAY(:) = Zero
  ! Loop over packed records of dimension LPBLK
  NBATCH = 0
  !1000 continue
  ! The next LPBLK elements
  LBATCH = -2**30
  do
    NBATCH = NBATCH+1
    if (NBATCH /= 1) LBATCHP = LBATCH
    ! Read next batch
    call IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
    LBATCH = IDUMMY(1)
    if (LBATCH > 0) then
      call IDAFILE(IFILE,2,IPAK,LBATCH,IDISK(IFILE))
      call DDAFILE(IFILE,2,XPAK,LBATCH,IDISK(IFILE))
    end if
    call IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
    ISTOP = IDUMMY(1)
    do IELMNT=1,LBATCH
      if ((IPAK(IELMNT) <= 0) .or. (IPAK(IELMNT) > NDIM)) then
        write(u6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT
        write(u6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT)
        write(u6,*) ' LBATCH IFILE  = ',LBATCH,IFILE
        if (NBATCH == 1) then
          write(u6,*) ' NBATCH = 1'
        else
          write(u6,*) ' NBATCH, LBATCHP',NBATCH,LBATCHP
        end if
        write(u6,*) ' NDIM,IMZERO = ',NDIM,IMZERO
        !stop ' problem in FRMDSC'
        call SYSABENDMSG('lucia_util/frmdsc','Internal error','')
      end if
      ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
    end do
    if (ISTOP /= 0) exit
  end do
  ! End of loop over records of truncated elements
else if (I_AM_PACKED == 0) then
  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  IREST = NDIM
  IBASE = 0
  do
    if (IREST > NBLOCK) then
      call DDAFILE(IFILE,2,ARRAY(IBASE+1),NBLOCK,IDISK(IFILE))
      IBASE = IBASE+NBLOCK
      IREST = IREST-NBLOCK
    else
      call DDAFILE(IFILE,2,ARRAY(IBASE+1),IREST,IDISK(IFILE))
      IREST = 0
    end if
    call IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
    if (IREST <= 0) exit
  end do

end if

end subroutine FRMDSC2

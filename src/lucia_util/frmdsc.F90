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

subroutine FRMDSC(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
! TRANSFER ARRAY FROM DISC FILE IFILE
!
! Version allowing zero and packed blocks

use lucia_data, only: IDISK
use Constants, only: Zero
use Definitions, only: u6

implicit none
real*8 ARRAY(*)
integer ISCR(2)
integer, parameter :: LPBLK = 50000
integer IPAK(LPBLK)
real*8 XPAK(LPBLK)
integer IDUMMY(1)
integer IPACK, IFILE, IMZERO, I_AM_PACKED, NDIM, NBATCH, LBATCH, LBATCHP, ISTOP, IELMNT, NBLOCK, MBLOCK, IREST, IBASE

IPACK = 1
if (IPACK /= 0) then
  ! Read if ARRAY is zero
  !MMBLOCK = MBLOCK
  !if (MMBLOCK >= 2) MMBLOCK = 2
  !    IFRMDS(ISCR,2,MMBLOCK,IFILE)
  call IFRMDS(ISCR,2,2,IFILE)
  IMZERO = ISCR(1)
  I_AM_PACKED = ISCR(2)
  !if (I_AM_PACKED /= 0) write(u6,*) ' File is packed, file number = ',IFILE
  if (IMZERO == 1) then
    !write(u6,*) ' frmdsc, length of zero block',NDIM
    call SETVEC(ARRAY,Zero,NDIM)
    goto 1001
  end if
end if
!write(u6,*) ' IMZERO I_AM_PACKED',IMZERO,I_AM_PACKED

if (I_AM_PACKED == 1) then
  call SETVEC(ARRAY,Zero,NDIM)
  ! Loop over packed records of dimension LPBLK
  NBATCH = 0
  !1000 CONTINUE
  ! The next LPBLK elements
  LBATCH = -2**30
999 continue
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
  if (ISTOP == 0) goto 999
  ! End of loop over records of truncated elements
else if (I_AM_PACKED == 0) then
  NBLOCK = MBLOCK
  if (MBLOCK <= 0) NBLOCK = NDIM
  IREST = NDIM
  IBASE = 0
100 continue
  if (IREST > NBLOCK) then
    call DDAFILE(IFILE,2,ARRAY(IBASE+1),NBLOCK,IDISK(IFILE))
    IBASE = IBASE+NBLOCK
    IREST = IREST-NBLOCK
  else
    call DDAFILE(IFILE,2,ARRAY(IBASE+1),IREST,IDISK(IFILE))
    IREST = 0
  end if
  call IDAFILE(IFILE,2,IDUMMY,1,IDISK(IFILE))
  if (IREST > 0) goto 100
end if

1001 continue

end subroutine FRMDSC

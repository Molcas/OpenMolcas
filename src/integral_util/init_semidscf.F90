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
! Copyright (C) 1990,1991,1993,1996, Roland Lindh                      *
!               1990, IBM                                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine Init_SemiDSCF(FstItr,Thize,Cutint)

use dEAF, only: dEAFARead, dEAFAWrite, dEAFRead
use IOBUF, only: Buffer, Disk, Disk_1, Disk_2, iBuf, ID, IODone, ipos, iStatIO, lBuf, LuTmp, Mode_Read, Mode_Write, nBuf, OnDisk
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, RtoI

implicit none
logical(kind=iwp), intent(in) :: FstItr
real(kind=wp), intent(inout) :: Thize
real(kind=wp), intent(in) :: CutInt
integer(kind=iwp) :: lBufOld, nBufOld
real(kind=wp) :: control(4), CutIntOld, ThizeOld

!write(u6,*) 'Enter: Init_SemiDSCF'
!write(u6,*) 'Ondisk=',Ondisk
!write(u6,*) 'lBuf=',lBuf

! Initiate asynchronous double buffer I/O.

IODone = .false.
Disk = Zero
iBuf = 1
iPos = 1
if (FstItr) then
  iStatIO = Mode_Write
  !write(u6,*) 'write istatio=',istatio
  control(1) = real(lbuf,kind=wp)
  control(2) = real(nbuf,kind=wp)
  control(3) = thize
  control(4) = cutint
  !write(u6,*) 'control written:',control
  !write(u6,*) ' Initiate write @',Disk,'iBuf=',iBuf
  if (OnDisk) call dEAFAwrite(LuTmp,control,4*RtoI,Disk,id)
else
  iStatIO = Mode_Read
  !write(u6,*) 'read istatio=',istatio

  ! Initiate first read ahead of time.

  !write(u6,*) 'lBuf*RtoI=',lbuf*RtoI,' rtoi=',Rtoi
  if (OnDisk) then
    !write(u6,*) ' Initiate read @',Disk,'iBuf=',iBuf
    call dEAFread(LuTmp,control,4*RtoI,Disk)
    Disk_2 = Disk
    Disk_1 = Disk
    !write(u6,*) 'control read:',control
    lbufold = nint(control(1))
    nbufold = nint(control(2))
    thizeold = control(3)
    cutintold = control(4)
    if (lbufold < lbuf) then
      write(u6,*) 'Reducing the buffer size from ',lbuf,' to',lbufold
      lbuf = lbufold
    else if (lbufold > lbuf) then
      write(u6,*) 'Inconsistent buffer lengths. Old:',lbufold,'  current:',lbuf
      call Abend()
    end if
    if (nbuf /= nbufold) then
      write(u6,*) 'Inconsistent buffer number. Old:',nbufold,'  current:',nbuf
      call Abend()
    end if
    if (abs(thize-thizeold) > 1.0e-10_wp) then
      write(u6,*) 'Resetting thize from',thize,' to',thizeold
      thize = thizeold
    end if
    if (cutintold > cutint) then
      write(u6,*) 'Inconsistent Cutint. Old:',cutintold,'  current:',cutint
      call Abend()
    end if
    !write(u6,*) ' Initiate read @',Disk,'iBuf=',iBuf
    !if (OnDisk) write(u6,*) ' Initial EAFARead'
    call dEAFARead(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk,id)
  end if
end if

!write(u6,*) 'Exit: Init_SemiDSCF'

end subroutine Init_SemiDSCF

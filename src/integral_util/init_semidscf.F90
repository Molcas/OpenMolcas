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

use dEAF
use IOBUF, only: IODone, Disk, iBuf, ipos, iStatIO, Mode_Write, OnDisk, Mode_Read, Disk_1, Disk_2, lBuf, nBuf, Buffer, ID, LuTmp

implicit none
#include "SysDef.fh"
logical FstItr
real*8 Thize, CutInt
integer lBufOld, nBufOld
real*8 ThizeOld, CutIntOld
real*8 control(4)

!write(6,*) 'Enter: Init_SemiDSCF'
!write(6,*) 'Ondisk=',Ondisk
!write(6,*) 'lBuf=',lBuf

! Initiate asynchronous double buffer I/O.

IODone = .false.
Disk = 0.0d0
iBuf = 1
iPos = 1
if (FstItr) then
  iStatIO = Mode_Write
  !write(6,*) 'write istatio=',istatio
  control(1) = dble(lbuf)
  control(2) = dble(nbuf)
  control(3) = thize
  control(4) = cutint
  !write(6,*) 'control written:',control
  !write(6,*) ' Initiate write @',Disk,'iBuf=',iBuf
  if (OnDisk) call dEAFAwrite(LuTmp,control,4*RtoI,Disk,id)
else
  iStatIO = Mode_Read
  !write(6,*) 'read istatio=',istatio

  ! Initiate first read ahead of time.

  !write(6,*) 'lBuf*RtoI=',lbuf*RtoI,' rtoi=',Rtoi
  if (OnDisk) then
    !write(6,*) ' Initiate read @',Disk,'iBuf=',iBuf
    call dEAFread(LuTmp,control,4*RtoI,Disk)
    Disk_2 = Disk
    Disk_1 = Disk
    !write(6,*) 'control read:',control
    lbufold = nint(control(1))
    nbufold = nint(control(2))
    thizeold = control(3)
    cutintold = control(4)
    if (lbufold < lbuf) then
      write(6,*) 'Reducing the buffer size from ',lbuf,' to',lbufold
      lbuf = lbufold
    else if (lbufold > lbuf) then
      write(6,*) 'Inconsistent buffer lengths. Old:',lbufold,'  current:',lbuf
      call Abend()
    end if
    if (nbuf /= nbufold) then
      write(6,*) 'Inconsistent buffer number. Old:',nbufold,'  current:',nbuf
      call Abend()
    end if
    if (abs(thize-thizeold) > 1.d-10) then
      write(6,*) 'Resetting thize from',thize,' to',thizeold
      thize = thizeold
    end if
    if (cutintold > cutint) then
      write(6,*) 'Inconsistent Cutint. Old:',cutintold,'  current:',cutint
      call Abend()
    end if
    !write(6,*) ' Initiate read @',Disk,'iBuf=',iBuf
    !if (OnDisk) write(6,*) ' Initial EAFARead'
    call dEAFARead(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk,id)
  end if
end if

!write(6,*) 'Exit: Init_SemiDSCF'

end subroutine Init_SemiDSCF

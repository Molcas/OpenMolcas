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

subroutine dWBuf(Array,nArray)

use dEAF, only: dEAFAWrite
use IOBUF, only: InCore, iBuf, IODone, lBuf, iPos, Disk, OnDisk, DiskMx_Byte, Disk_1, Disk_2, Buffer, iD, LuTmp

implicit none
#include "SysDef.fh"
integer nArray
real*8 Array(nArray)
integer iArray, mArray, Left
real*8 Temp

!write(6,*) 'Enter WBuf: iPos @',iPos,' iBuf,lBuf=',iBuf,lBuf
if (InCore .and. (iBuf == 2)) then
  call WarningMessage(2,'Error in in-core semi-direct implementation')
  call Abend()
end if
IODone = .true.
iArray = 1
mArray = nArray
10 continue
Left = lBuf-iPos+1
if (mArray > Left) then
  call dCopy_(Left,Array(iArray),1,Buffer(iPos,iBuf),1)
  iArray = iArray+Left
  mArray = mArray-Left
  iPos = 1

  ! Wait for previous buffer to complete I/O.
  ! Disk=32 after writing the control!

  !write(6,*) 'In dwbuf'
  if ((Disk /= 32.0d0) .and. OnDisk) call EAFWait(LuTmp,id)

  ! Put current buffer on disk and change buffer.

  temp = Disk+dble(lBuf*RtoB)
  !write(6,*) 'temp=',temp
  if (temp <= DiskMx_Byte) then
    Disk_2 = Disk_1
    Disk_1 = Disk
    !write(6,*) 'WBuf write on disk @',Disk,'iBuf=',iBuf
    if (OnDisk) call dEAFAWrite(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk,id)
    if (iBuf == 1) then
      iBuf = 2
    else
      iBuf = 1
    end if
  else
    call WarningMessage(2,'WBuf: Disc is full!!')
    call Abend()
  end if
else
  !write(6,*) ' Add ',mArray,'elements to buffer',iPos,ibuf
  call dCopy_(mArray,Array(iArray),1,Buffer(iPos,iBuf),1)
  iPos = iPos+mArray
  mArray = 0
end if
if (mArray > 0) goto 10

!write(6,*) 'Exit WBuf: iPos @',iPos,'iBuf=',iBuf

return

end subroutine dWBuf

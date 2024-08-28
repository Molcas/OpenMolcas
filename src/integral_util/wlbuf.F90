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

subroutine WLBuf()

use dEAF, only: dEAFWrite
use IOBUF, only: iStatIO, Mode_Read, OnDisk, InCore, iBuf, iPos, Disk, lBuf, DiskMx_Byte, Disk_1, Disk_2, Buffer, ID, LuTmp
use Constants, only: Zero
use Definitions, only: wp, u6

implicit none
#include "SysDef.fh"
real*8 Temp

if (iStatIO == Mode_Read) then
  !write(u6,*) 'In WLbuf'
  if (OnDisk) call EAFWait(LuTmp,id)
  return
end if
!Disk_Save = Disk
!write(u6,*) 'Enter WLBuf: Disk,iPos,iBuf=',Disk,iPos,iBuf
if (InCore .and. (iBuf == 2)) then
  call WarningMessage(2,'Error in in-core semi-direct implementation')
  call Abend()
end if

! If any data in buffer write buffer to disk.

!write(u6,*) 'In WLbuf'
if (OnDisk) call EAFWait(LuTmp,id)
if (iPos /= 1) then
  temp = Disk+real(lBuf*RtoB,kind=wp)
  !write(u6,*) 'temp,DiskMx_Byte=',temp,DiskMx_Byte
  if (temp <= DiskMx_Byte) then
    Disk_2 = Disk_1
    Disk_1 = Disk
    !if (OnDisk) write(u6,*) 'Disk=',Disk,' lBuf*RtoI=',lBuf*RtoI
    !write(u6,*) 'WLBuf write on disk @',Disk,'iBuf=',iBuf
    if (OnDisk) call dEAFWrite(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk)
    ! Put a dummy record at the end
    temp = Disk+real(lBuf*RtoB,kind=wp)
    !write(u6,*) 'temp,DiskMx_Byte=',temp,DiskMx_Byte
    if ((temp <= DiskMx_Byte) .and. OnDisk) then
      !write(u6,*) 'WLBuf write on disk @',Disk,'iBuf=',iBuf
      call dCopy_(lBuf,[Zero],0,Buffer(1,iBuf),1)
      call dEAFWrite(LuTmp,Buffer(1,iBuf),lBuf*RtoI,Disk)
    end if
  else
    call WarningMessage(2,'WLBuf: Disc is full!')
    write(u6,*) 'temp           =',temp
    write(u6,*) 'DiskMx_Byte    =',DiskMx_Byte
    call FastIO('STATUS')
    call Abend()
  end if
end if
iPos = 1

!if (Disk_save /= Disk) then
!  write(u6,*) 'Enter WLBuf: Disk @:',Disk_Save
!  write(u6,*) 'Exit  WLBuf: Disk @:',Disk
!end if
!write(u6,*) 'Exit WLBuf: Disk,iPos,iBuf=',Disk,iPos,iBuf

end subroutine WLBuf

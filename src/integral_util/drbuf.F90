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

subroutine dRBuf(Array,nArray,Copy)

use dEAF, only: dEAFARead
use IOBUF, only: Buffer, Disk, Disk_1, Disk_2, DiskMX_Byte, iBuf, iD, InCore, iPos, lBuf, LuTmp, OnDisk
use Definitions, only: wp, iwp, RtoB, RtoI

implicit none
integer(kind=iwp), intent(in) :: nArray
real(kind=wp), intent(inout) :: Array(nArray)
logical(kind=iwp), intent(in) :: Copy
integer(kind=iwp) :: iArray, jBuf, Left, mArray
real(kind=wp) :: Temp

!write(u6,*) 'Enter RBuf: iPos @',iPos,'iBuf=,lBuf',iBuf,lBuf
if (InCore .and. (iBuf == 2)) then
  call WarningMessage(2,'Error in in-core semi-direct implementation')
  call Abend()
end if
iArray = 1
mArray = nArray
do while (mArray >= 1)
  if (iPos == 1) then

    ! Wait for pending buffer.

    !write(u6,*) 'In drbuf.'
    if (OnDisk) call EAFWait(LuTmp,id)

    ! Get the next buffer, make sure that request is not beyond
    ! the disc limitation.

    temp = Disk+real(lBuf*RtoB,kind=wp)
    !write(u6,*) 'temp=',temp
    if (temp <= DiskMx_Byte) then
      if (iBuf == 1) then
        jBuf = 2
      else
        jBuf = 1
      end if
      Disk_2 = Disk_1
      Disk_1 = Disk
      !write(u6,*) 'RBuf aread on disk @',Disk,'jBuf=',jBuf
      if (OnDisk) call dEAFARead(LuTmp,Buffer(1,jBuf),lBuf*RtoI,Disk,id)
    end if
  end if
  Left = lBuf-iPos+1
  if (mArray > Left) then
    if (Copy) Array(iArray:iArray+Left-1) = Buffer(iPos:,iBuf)
    iArray = iArray+Left
    mArray = mArray-Left
    iPos = 1
    !write(u6,*) 'LuTmp,Disk=',LuTmp,Disk
    ! Swap buffer
    if (iBuf == 1) then
      iBuf = 2
    else
      iBuf = 1
    end if
  else
    !write(u6,*) ' Copy ',mArray,'elements from buffer',iPos
    if (Copy) Array(iArray:iArray+mArray-1) = Buffer(iPos:iPos+mArray-1,iBuf)
    iPos = iPos+mArray
    mArray = 0
  end if
end do

!write(u6,*) 'Exit RBuf: iPos @',iPos,'iBuf=',iBuf
return

end subroutine dRBuf

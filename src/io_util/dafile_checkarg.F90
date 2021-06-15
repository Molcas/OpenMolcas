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
! Copyright (C) 2012, Victor P. Vysotskiy                              *
!***********************************************************************

subroutine DaFile_checkarg(Lu,iOpt,lBuf,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Check arguments to iDafile                                       *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
!     iOpt    : integer, input                                         *
!               valid values: 0,1,2,5,6,7,8,10,99                      *
!     lBuf    : integer, input                                         *
!               length of the buffer Buf:   0<lBuf<2**29(2**60)        *
!     iDisk   : integer, input/output                                  *
!               disk address  >=0                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
!                                                                      *
!***********************************************************************

use Fast_IO, only: isOpen, MxFile
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf, iDisk
character(len=*), parameter :: TheName = 'DaFile_checkarg'

! 2012
! VpV: a lot of checking is here.
! Check arguments
if ((Lu <= 0) .or. (Lu > MxFile)) call SysFileMsg(TheName,'MSG: unit',Lu,' ')

if (isOpen(Lu) == 0) call SysFileMsg(TheName,'MSG: not opened',Lu,' ')

if (lBuf < 0) then
  write(u6,*) 'Invalid buffer size ',lBuf
  call Error()
end if

if (iDisk < 0) then
  write(u6,*) 'Invalid disk address ',iDisk
  call Error()
end if

if ((iOpt < 0) .or. ((iOpt > 10) .and. (iOpt /= 99))) then
  write(u6,*) 'Invalid action code ',iOpt
  call Error()
end if

if ((iOpt == 3) .or. (iOpt == 4) .or. (iOpt == 9)) then
  write(u6,*) 'DaFile: GSlist option is not in operation!'
  call Error()
end if

return

contains

subroutine Error()
  write(u6,*) 'I/O error in ',TheName
  write(u6,*) 'Unit = ',Lu
  call Abend()
end subroutine Error

end subroutine DaFile_checkarg

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
! Copyright (C) 1991, Per-Olof Widmark                                 *
!               1993,1996,1997, Markus P. Fuelscher                    *
!               1996, Luis Serrano-Andres                              *
!               2001, Roland Lindh                                     *
!               2016, Steven Vancoillie                                *
!***********************************************************************

subroutine bDaFile(Lu,iOpt,Buf,lBuf,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Control direct access I/O operations                             *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
!               If Lu={40,50,60,70,80,90} we are possibly concerned    *
!               with a multi file unit (c.f. allocdisk)                *
!     iOpt    : integer, input                                         *
!               option code                                            *
!               iOpt= 0 dummy write                                    *
!               iOpt= 99 dummy read (return buf(1)=1 in success        *
!               iOpt= 1 synchronous write                              *
!               iOpt= 2 synchronous read                               *
!               iOpt= 5 synchronous rewind                             *
!               iOpt= 6 asynchronous write                             *
!               iOpt= 7 asynchronous read                              *
!               iOpt=10 asynchronous rewind                            *
!               Note: At present the asynchronous modes are not        *
!                     supported and work identically the synchronous   *
!                     modes                                            *
!     Buf     : array of characters, input/output                      *
!               Buffer carrying/receiving the data to write/read       *
!     lBuf    : integer, input                                         *
!               length of the buffer Buf in bytes!!!                   *
!     iDisk   : integer, input/output                                  *
!               disk address as byte offset!!!                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!     R. Lindh, University of Lund, Sweden, 2001                       *
!     Steven Vancoillie: use of byte lengths/offests (2016)            *
!                                                                      *
!***********************************************************************

use Fast_IO, only: Addr, FSCB, LuName, MaxFileSize, Multi_File, Trace
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf
character, intent(inout) :: Buf(*)
integer(kind=iwp), intent(inout) :: iDisk
integer(kind=iwp) :: iRC, kDisk, jDisk
logical(kind=iwp) :: Skip
integer(kind=iwp), external :: AixFsz

if (Trace) then
  write(u6,*) ' >>> Enter bDaFile <<<'
  write(u6,*) ' unit      :',Lu
  write(u6,*) ' name      :',LuName(Lu)
  write(u6,*) ' option    :',iOpt
  write(u6,*) ' length    :',lBuf
  write(u6,*) ' disk adr. :',iDisk
end if

Skip = .false.
if ((iOpt == 5) .or. (iOpt == 10)) then
  Addr(Lu) = 0
  iDisk = 0
  Skip = .true.
! Get filesize in 8-byte units
else if (iOpt == 0) then
  ! Dummy write. No I/O is made. Disk address is updated.
  Addr(Lu) = iDisk+lBuf
  iDisk = Addr(Lu)
  Skip = .true.
else if (iOpt == 8) then
  iRc = AixFsz(FSCB(Lu))
  iDisk = iRc
  Skip = .true.
end if

if (.not. Skip) then
  if (Multi_File(Lu) .and. MaxFileSize /= 0) then
    jDisk = iDisk
    kDisk = iDisk
    call MpDaFile(Lu,MaxFileSize,iOpt,Buf,lBuf,kDisk)
    Addr(Lu) = jDisk+lBuf
    iDisk = Addr(Lu)
  else
    call ChDaFile(Lu,iOpt,Buf,lBuf,iDisk)
  end if
end if

if (Trace) write(u6,*) ' >>> Exit bDaFile <<<'

return

end subroutine bDaFile

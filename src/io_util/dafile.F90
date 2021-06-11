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
!               2012, Victor P. Vysotskiy                              *
!               2016, Steven Vancoillie                                *
!***********************************************************************

subroutine DaFile(Lu,iOpt,Buf,lBuf,iDisk)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Control direct access I/O operations                             *
!                                                                      *
!     calling arguments:                                               *
!     Lu      : integer, input                                         *
!               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
!     iOpt    : integer, input                                         *
!               option code                                            *
!               iOpt= 0 Dummy write. No I/O is made. Disk address is   *
!               updated.                                               *
!               iOpt= 99 dummy read (return buf(1)=1 in success)       *
!               iOpt= 1 synchronous write                              *
!               iOpt= 2 synchronous read                               *
!               iOpt= 5 synchronous rewind                             *
!               iOpt= 6 asynchronous write                             *
!               iOpt= 7 asynchronous read                              *
!               iOpt= 8 position at the end of file (i.e. filesize)    *
!               iOpt=10 asynchronous rewind                            *
!               Note: At present the asynchronous modes are not        *
!                     supported and work identically the synchronous   *
!                     modes                                            *
!     Buf     : array of integers, input/output                        *
!               Buffer carrying/receiving the data to write/read       *
!     lBuf    : integer, input                                         *
!               length of the buffer Buf in bytes!!                    *
!     iDisk   : integer, input/output                                  *
!               disk address                                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark, IBM Sweden, 1991                                   *
!     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
!     L. Serrano-Andres, University of Lund, Sweden, 1996              *
!     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
!----------------------------------------------------------------------*
!                                                                      *
! History:                                                             *
!     Steven Vancoillie: use of byte lengths/offests (2016)            *
!                                                                      *
!***********************************************************************

use Definitions, only: iwp, u6
#ifdef _OLD_IO_STAT_
use Definitions, only: wp
#endif

implicit none
integer(kind=iwp), intent(in) :: Lu, iOpt, lBuf
integer(kind=iwp), intent(inout) :: Buf(*), iDisk
integer(kind=iwp) :: iRc = 0, lDisk
character(len=80) :: Text, HeadErr
integer(kind=iwp), external :: AixErr, AixRd, AixWr
#include "fio.fh"
#include "warnings.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#endif

call DaFile_checkarg(Lu,iOpt,lBuf,iDisk)
!****************** REAL I/O IS HERE **************
lDisk = iDisk
! Write to  disk
if ((iOpt == 1) .or. (iOpt == 6)) then
  HeadErr = 'Premature abort while writing buffer to disk'
! lock
# if defined (_HAVE_EXTRA_) && ! defined (_GA_)
  if (isFiM(Lu) == 0) then
# endif
    iRc = AixWr(FSCB(Lu),Buf,lBuf,lDisk)
# if defined (_HAVE_EXTRA_) && ! defined (_GA_)
  else
    iRc = FimWr(FSCB(Lu),Buf,lBuf,lDisk)
    if (iRc == -10) then
      isFiM(Lu) = 0
      iRc = AixWr(FSCB(Lu),Buf,lBuf,lDisk)
    end if
  end if
# endif
! unlock
! Read from disk
else if ((iOpt == 2) .or. (iOpt == 7) .or. (iOpt == 99)) then
  HeadErr = 'Premature abort while reading buffer from disk'
  !10 continue
# if defined (_HAVE_EXTRA_) && ! defined (_GA_)
  if (isFiM(Lu) == 0) then
# endif
    if (iOpt /= 99) then
      iRc = AixRd(FSCB(Lu),Buf,lBuf,lDisk,0)
    else
      iRc = AixRd(FSCB(Lu),Buf,lBuf,lDisk,1)
      Buf(1) = 0
      if (iRc == 0) Buf(1) = 1
      return
    end if
# if defined (_HAVE_EXTRA_) && ! defined (_GA_)
  else
    iRc = FimRd(FSCB(Lu),Buf,lBuf,lDisk)
  end if
# endif
end if

if (iRc /= 0) goto 1200

#ifdef _OLD_IO_STAT_
MxAddr(Lu) = max(MxAddr(Lu),Addr(Lu))
if (iOpt /= 0) then
  count(3,Lu) = count(3,Lu)+1
  count(4,Lu) = count(4,Lu)+real(lBuf,kind=wp)/1024.0_wp
end if
#endif
iDisk = iDisk+lBuf
Addr(Lu) = iDisk
if (Trace) write(u6,*) ' >>> Exit DaFile <<<'

return

1200 continue
iRc = AixErr(Text)
write(u6,*) HeadErr
write(u6,*) Text
write(u6,*) ' Unit      :',Lu
write(u6,*) ' Option    :',iOpt
write(u6,*) ' Buffer    :',lBuf
write(u6,*) ' Address   :',iDisk
call quit(_RC_IO_ERROR_WRITE_)

end subroutine DaFile

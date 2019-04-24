************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Per-Olof Widmark                                 *
*               1993,1996,1997, Markus P. Fuelscher                    *
*               1996, Luis Serrano-Andres                              *
*               2001, Roland Lindh                                     *
*               2016, Steven Vancoillie                                *
************************************************************************
      Subroutine bDaFile(Lu,iOpt,Buf,lBuf,iDisk)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Control direct access I/O operations                             *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*               If Lu={40,50,60,70,80,90} we are possibly concerned    *
*               with a multi file unit (c.f. allocdisk)                *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 0 dummy write                                    *
*               iOpt= 99 dummy read (return buf(1)=1 in success        *
*               iOpt= 1 synchronous write                              *
*               iOpt= 2 synchronous read                               *
*               iOpt= 5 synchronous rewind                             *
*               iOpt= 6 asynchronous write                             *
*               iOpt= 7 asynchronous read                              *
*               iOpt=10 asynchronous rewind                            *
*               Note: At present the asynchronous modes are not        *
*                     supported and work identically the synchronous   *
*                     modes                                            *
*     Buf     : array of characters, input/output                      *
*               Buffer carrying/receiving the data to write/read       *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf in bytes!!!                   *
*     iDisk   : integer, input/output                                  *
*               disk address as byte offset!!!                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*     R. Lindh, University of Lund, Sweden, 2001                       *
*     Steven Vancoillie: use of byte lengths/offests (2016)            *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

#include "SysDef.fh"
#include "blksize.fh"
#include "fio.fh"

      Character Buf(*)

      If ( Query ) Call qEnter('bDaFile')

      If ( Trace ) then
        Write (6,*) ' >>> Enter bDaFile <<<'
        Write (6,*) ' unit      :',Lu
        Write (6,*) ' name      :',LuName(Lu)
        Write (6,*) ' option    :',iOpt
        Write (6,*) ' length    :',lBuf
        Write (6,*) ' disk adr. :',iDisk
      End If

      If((iOpt.eq.5).or.(iOpt.eq.10)) then
        Addr(Lu)   = 0
        iDisk      = 0
        goto 1100
C*     Get filesize in 8-byte units
      Else If (iOpt.eq.0) Then
C*     Dummy write. No I/O is made. Disk address is updated.
        Addr(Lu)   = iDisk+lBuf
        iDisk      = Addr(Lu)
        goto 1100
      Else If ( iOpt.eq.8 ) then
        iRc=AixFsz(FSCB(Lu))
        iDisk=iRc
        goto 1100
      End If

      If ( Multi_File(Lu) .and. MaxFileSize.ne.0 ) then
         jDisk       = iDisk
         kDisk       = iDisk
         Call MpDaFile(Lu,MaxFileSize,iOpt,Buf,lBuf,kDisk)
         Addr(Lu)    = jDisk+lBuf
         iDisk       = Addr(Lu)
      Else
         Call ChDaFile(Lu,iOpt,Buf,lBuf,iDisk)
      End If

1100  Continue
      If ( Trace ) Write (6,*) ' >>> Exit bDaFile <<<'
      If ( Query ) Call qExit('bDaFile')
      Return
      End

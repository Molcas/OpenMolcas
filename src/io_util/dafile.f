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
*               2012, Victor P. Vysotskiy                              *
*               2016, Steven Vancoillie                                *
************************************************************************
#define MOLCASWRITE AixWr
#define MOLCASREAD  AixRd
      Subroutine DaFile(Lu,iOpt,Buf,lBuf,iDisk)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Control direct access I/O operations                             *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*     iOpt    : integer, input                                         *
*               option code                                            *
*               iOpt= 0 Dummy write. No I/O is made. Disk address is   *
*               updated.                                               *
*               iOpt= 99 dummy read (return buf(1)=1 in success)       *
*               iOpt= 1 synchronous write                              *
*               iOpt= 2 synchronous read                               *
*               iOpt= 5 synchronous rewind                             *
*               iOpt= 6 asynchronous write                             *
*               iOpt= 7 asynchronous read                              *
*               iOpt= 8 position at the end of file (i.e. filesize)    *
*               iOpt=10 asynchronous rewind                            *
*               Note: At present the asynchronous modes are not        *
*                     supported and work identically the synchronous   *
*                     modes                                            *
*     Buf     : array of integers, input/output                        *
*               Buffer carrying/receiving the data to write/read       *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf in bytes!!                    *
*     iDisk   : integer, input/output                                  *
*               disk address                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, IBM Sweden, 1991                                   *
*     M.P. Fuelscher, University of Lund, Sweden, 1993, 1996, 1997     *
*     L. Serrano-Andres, University of Lund, Sweden, 1996              *
*     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*     Steven Vancoillie: use of byte lengths/offests (2016)            *
*                                                                      *
************************************************************************

      Implicit Integer (A-Z)

#include "SysDef.fh"
#include "blksize.fh"
#include "fio.fh"
#ifdef _OLD_IO_STAT_
#include "ofio.fh"
#endif
#include "warnings.fh"
      Dimension Buf(*)
      Character*80 Text,HeadErr

      Character*16 TheName
      Data TheName/'DaFile'/
      Integer iRc
      Data    iRc/0/

      If ( Query ) Call qEnter(TheName)

      Call DaFile_checkarg(Lu,iOpt,lBuf,iDisk)
C****************** REAL I/O IS HERE **************
      lDisk    = iDisk
*     Write to  disk
      If ( (iOpt.eq.1) .or. (iOpt.eq.6) ) then
        HeadErr='Premature abort while writing buffer to disk'
*lock
#ifdef _GA_
        iRc = MOLCASWRITE(FSCB(Lu),Buf,lBuf,lDisk)
#else
        If(isFiM(Lu).eq.0) then
          iRc = MOLCASWRITE(FSCB(Lu),Buf,lBuf,lDisk)
        Else
          iRc = FimWr(FSCB(Lu),Buf,lBuf,lDisk)
          If(iRc.eq.-10)Then
            isFiM(Lu)=0
            iRc = MOLCASWRITE(FSCB(Lu),Buf,lBuf,lDisk)
          End If
        End If
#endif
*unlock
*     Read from disk
      Else If ( (iOpt.eq.2) .or. (iOpt.eq.7) .or. (iOpt.eq.99)) then
        HeadErr='Premature abort while reading buffer from disk'
C10    Continue
#ifdef _GA_
        if(iOpt.ne.99) then
           iRc = MOLCASREAD(FSCB(Lu),Buf,lBuf,lDisk,0)
        else
           iRc = MOLCASREAD(FSCB(Lu),Buf,lBuf,lDisk,1)
           Buf(1)=0
           if(iRc.eq.0) Buf(1)=1
           return
        endif
#else
        If (isFiM(Lu).eq.0) then
           if(iOpt.ne.99) then
             iRc = MOLCASREAD(FSCB(Lu),Buf,lBuf,lDisk,0)
           else
             iRc = MOLCASREAD(FSCB(Lu),Buf,lBuf,lDisk,1)
             Buf(1)=0
             if(iRc.eq.0) Buf(1)=1
             return
           endif
        Else
           iRc = FimRd(FSCB(Lu),Buf,lBuf,lDisk)
        End If
#endif
      End if

      If ( iRc.ne.0 ) goto 1200

#ifdef _OLD_IO_STAT_
        MxAddr(Lu)  = Max(MxAddr(Lu),Addr(Lu))
        if(iOpt.ne.0) Then
           Count(3,Lu) = Count(3,Lu)+1
           Count(4,Lu) = Count(4,Lu)+dble(lBuf)/1024.0d0
        End if
#endif
      iDisk       = iDisk+lBuf
      Addr(Lu)    = iDisk
      If ( Trace ) Write (6,*) ' >>> Exit DaFile <<<'
      If ( Query ) Call qExit(TheName)

      Return

1200  iRc = AixErr(Text)
      write (6,*) HeadErr
      write (6,*) Text
      write (6,*) ' Unit      :',Lu
      write (6,*) ' Option    :',iOpt
      write (6,*) ' Buffer    :',lBuf
      write (6,*) ' Address   :',iDisk
      call quit(_RC_IO_ERROR_WRITE_)
      End


      Subroutine DaFile_checkarg(Lu,iOpt,lBuf,iDisk)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Check arguments to iDafile                                       *
*                                                                      *
*     calling arguments:                                               *
*     Lu      : integer, input                                         *
*               logical unit number (Lu={1,2,...40,50,60,70,80,90}     *
*     iOpt    : integer, input                                         *
*               valid values: 0,1,2,5,6,7,8,10,99                      *
*     lBuf    : integer, input                                         *
*               length of the buffer Buf:   0<lBuf<2**29(2**60)        *
*     iDisk   : integer, input/output                                  *
*               disk address  >=0                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     V.P. Vysotskiy, University of Lund, Sweden, 2012                 *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
#include "fio.fh"
*2012
*VpV: a lot of checking is here.
*     Check arguments
      If ( (Lu.le.0) .or. (Lu.gt.MxFile) )
     * Call SysFileMsg(TheName,'MSG: unit', Lu,' ')

      If ( isOpen(Lu).eq.0 )
     * Call SysFileMsg(TheName,'MSG: not opened', Lu,' ')

      If ( lBuf.lt.0) then
          write (6,*) 'Invalid buffer size ',lBuf
          goto 1000
      End If

      If ( iDisk.lt.0 ) then
         write (6,*) 'Invalid disk address ',iDisk
         goto 1000
      endif

      If ( (iOpt.lt.0) .or. (iOpt.gt.10.and.iOpt.ne.99) ) then
         write (6,*) 'Invalid action code ',iOpt
         goto 1000
      endif

      If (iOpt.eq.3.or.iOpt.eq.4.or.iOpt.eq.9) then
        Write (6,*) 'DaFile: GSlist option is not in operation!'
        goto 1000
      End If

      return

1000  write (6,*) 'I/O error in ',TheName
      write (6,*) 'Unit = ',Lu
      Call QTrace()
      Call Abend()

      End

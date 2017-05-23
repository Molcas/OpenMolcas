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
* Copyright (C) 2012,2013, Victor P. Vysotskiy                         *
************************************************************************
************************************************************************
*                                                                      *
*                     Thread safe                                      *
*                     POSIX - I/O                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixPWr(Handle,Buf,nBuf,iDisk)                                     *
*                                                                      *
* A buffer is written to a file associated with the file handle. The   *
* operation is synchronous, but can be done in parallel by multiple    *
* threads                                                              *
*                                                                      *
* Input:  Handle   - This is the unique file identifier associated     *
*                    with the file. It is created by AixOpn, and must  *
*                    be used on subsequent references to the file.     *
*         Buf      - The buffer that is to be written to disk.         *
*         nBuf     - Length of the buffer in words.                    *
*         iDisk    - External disk address.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Victor P. Vysotskiy                                         *
*          University of Lund, Sweden, 2012-2013                       *
* Written: 2012-2013                                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*                                                                      *
************************************************************************
      Integer Function AixPWr(handle,Buf,nBuf,iDisk)
      Implicit Integer (a-z)

#include "SysDef.fh"

#include "blksize.fh"
#include "switch.fh"
#include "ctl.fh"
      Dimension Buf(*)
      Character*80 ErrTxt
      Character*16 TheName
#ifndef _OLD_IO_STAT_
#include "pfio.fh"
      Real*8   CPUA,CPUE,TIOA,TIOE
#endif
#include "warnings.fh"
      Data TheName/'AixPWr'/
*----------------------------------------------------------------------*
* Entry to AixWr                                                       *
*----------------------------------------------------------------------*
      AixPWr=0
      rc=0
*----------------------------------------------------------------------*
* Check if file is opened.                                             *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pHndle,n).ne.handle) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixPWr=eNtOpn
            Return
         End If
         Go To 100
      End If
      nFile=n
      desc=CtlBlk(pDesc,nFile)
#ifndef _OLD_IO_STAT_
      Call FSCB2UNIT(handle,Lu)
      Call Timing(CPUA,CPUE,TIOA,TIOE)
#endif
*----------------------------------------------------------------------*
* Position file pointer                                                *
*----------------------------------------------------------------------*
      pDisk=pHeadOffset+iDisk
#ifndef _OLD_IO_STAT_
      If(CtlBlk(pWhere,nFile).ne.pDisk) Then
         ProfData(7,Lu)=ProfData(7,Lu)+1
      End If
#endif
*----------------------------------------------------------------------*
* Write to file                                                        *
*----------------------------------------------------------------------*
      CtlBlk(pWhere,nFile)=pDisk+nBuf
      if(nBuf.gt.0) rc=c_pwrite(desc,Buf,nBuf,pDisk)
      If(rc.lt.0) Then
         idum=iPrintLevel(3)
         Call FASTIO('STATUS')
         AixPWr=AixErr(ErrTxt)
            Call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,
     *            TheName,FCtlBlk(nFile),
     *     'Premature abort while writing buffer to disk', ErrTxt)
      Else If(rc.ne.nBuf) Then
         idum=iPrintLevel(3)
         Call FASTIO('STATUS')
         AixPWr=eEof
            Call SysQuitFileMsg(_RC_IO_ERROR_WRITE_,
     *            TheName,FCtlBlk(nFile),
     *            'Premature abort while writing buffer to disk:',
     *      'Disk full? ')
      End If
#ifndef _OLD_IO_STAT_
      Call Timing(CPUA,CPUE,TIOA,TIOE)
      ProfData(1,Lu)=ProfData(1,Lu)+1
      ProfData(2,Lu)=ProfData(2,Lu)+nBuf
      ProfData(3,Lu)=ProfData(3,Lu)+TIOE
#endif
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End

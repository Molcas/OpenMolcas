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
* Copyright (C) 1990, Per-Olof Widmark                                 *
*               2012,2013, Victor P. Vysotskiy                         *
************************************************************************
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixRd(Handle,Buf,nBuf,iDisk,iErrSkip)                             *
*                                                                      *
* A buffer is read from a file associated with the file handle. The    *
* operation is asynchronous, and must be followed by a call to AixWt   *
* to ensure that data is in memory.                                    *
*                                                                      *
* Input:  Handle   - This is the unique file identifier associated     *
*                    with the file. It is created by AixOpn, and must  *
*                    be used on subsequent references to the file.     *
*         Buf      - The buffer that is to be written to disk.         *
*         nBuf     - Length of the buffer in words.                    *
*         iDisk    - External disk address.                            *
*         iErrSkip - if 0 : stop in the error case                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          S&TC, ACIS, IBM Sweden                                      *
* Written: November 1990                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*          V.P. Vysotskiy,University of Lund, Sweden, 2012-2013        *
*                                                                      *
************************************************************************
      Integer Function AixRd(handle,Buf,nBuf,iDisk,iErrSkip)
      Implicit Integer (a-z)

#include "SysDef.fh"

#include "blksize.fh"
#include "switch.fh"
#include "ctl.fh"
      Dimension Buf(*)
      Character*80 ErrTxt
      Character*16 TheName
      Character*64 Temp
#ifndef _OLD_IO_STAT_
#include "pfio.fh"
      Real*8   CPUA,CPUE,TIOA,TIOE
#endif
#include "warnings.fh"
      Data TheName/'AixRd'/
*----------------------------------------------------------------------*
* Entry to AixRd                                                       *
*----------------------------------------------------------------------*
      AixRd=0
        Temp='Premature abort while reading buffer from disk'
*----------------------------------------------------------------------*
* Check if file is opened.                                             *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pHndle,n).ne.handle) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixRd=eNtOpn
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
      If(CtlBlk(pWhere,nFile).ne.pDisk) Then
         rc=c_lseek(desc,pDisk)
#ifndef _OLD_IO_STAT_
         ProfData(8,Lu)=ProfData(8,Lu)+1
#endif
         If(rc.lt.0) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
            idum=iPrintLevel(3)
            Call FASTIO('STATUS')
            AixRd=AixErr(ErrTxt)
            Call SysWarnFileMsg(TheName, FCtlBlk(nFile),
     *                                'MSG: seek', ErrTxt)
            Call SysCondMsg('rc < 0', rc, '<', 0)
         Else If(rc.ne.pDisk) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
            idum=iPrintLevel(3)
            Call FASTIO('STATUS')
            AixRd=eInErr
            Call SysWarnFileMsg(TheName, FCtlBlk(nFile),
     *                                'MSG: seek', ' ')
            Call SysCondMsg('rc != pDisk', rc, '!=', pDisk)
         End If
      End If
      CtlBlk(pWhere,nFile)=pDisk
*----------------------------------------------------------------------*
* Read from file                                                       *
*----------------------------------------------------------------------*
      rc=c_read(desc,Buf,nBuf)
      If(rc.lt.0) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
         idum=iPrintLevel(3)
         Call FASTIO('STATUS')
         AixRd=AixErr(ErrTxt)


        Call SysQuitFileMsg(_RC_IO_ERROR_READ_,
     &    TheName,FCtlBlk(nFile),Temp, ErrTxt)

      Else If(rc.ne.nBuf) Then
            if(iErrSkip.eq.1) then
             AixRd=99
             return
            endif
         idum=iPrintLevel(3)
         Call FASTIO('STATUS')
         AixRd=eEof
            Call SysQuitFileMsg(_RC_IO_ERROR_READ_,
     *            TheName,FCtlBlk(nFile),
     *            Temp, '\nEnd of file reached ')
      End If
      CtlBlk(pWhere,nFile)=CtlBlk(pWhere,nFile)+nBuf
      iDisk=iDisk+nBuf
#ifndef _OLD_IO_STAT_
      Call Timing(CPUA,CPUE,TIOA,TIOE)
      ProfData(4,Lu)=ProfData(4,Lu)+1
      ProfData(5,Lu)=ProfData(5,Lu)+nBuf
      ProfData(6,Lu)=ProfData(6,Lu)+TIOE
#endif

*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End

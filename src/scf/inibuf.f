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
* Copyright (C) 1998, Roland Lindh                                     *
************************************************************************
      SubRoutine IniBuf(nDisc,nCore)
************************************************************************
*                                                                      *
*  Object: Initiate I/O buffer for semi-direct SCF                     *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, Sweden. October '98                  *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
#include "IOBuf.fh"
      External AllocDisk
      Integer AllocDisk
*
*     Open file for semi-direct implementation
*     nDisc in units of MByte
*     nCore in units of kByte
*
#ifdef _DEBUG_
      Call qEnter('IniBuf')
#endif
*
*     The maximum number of bytes on disk. The file size limit times
*     the number of multi files.
*
      DiskMx_MByte=DBLE(AllocDisk())*10.0D0
      DiskMx_Byte=DBLE(AllocDisk())*10.0D0*2.00**20
*
      nBuf=-99
      If (nDisc.eq.0.and.nCore.eq.0) Then
         OnDisk=.False.
         Incore=.False.
      Else If (nDisc*1024.gt.nCore) Then
         OnDisk=.True.
         Incore=.False.
         LuTMp = 32
         Call new_EAFOpen(LuTmp,'SMDINT  ')
         nBuf=2
      Else
         OnDisk=.False.
         InCore=.True.
         nBuf=1
      End If
*
*---- Adjust buffer size and allocate memory for the buffer.
*
      ipBuf = -99
      If (OnDisk.or.InCore) Then
         MemMin_Seward=1024**2 ! Real*8
         Call GetMem('IniBuf','Max','Real',iDum,MaxMem)
*        lBuf in units of Real*8 per buffer
         lBuf=(1024*nCore)/(8*nBuf)
         if(InCore) then
         MemReq=MemMin_Seward + lBuf*nBuf
*        Write (6,*) 'MemReq,MaxMem=',MemReq,MaxMem,'  lbuf=',lbuf
         If (MemReq.gt.MaxMem) Then
            lBuf=(MaxMem-MemMin_Seward)/nBuf
            If (lBuf.lt.0) Then
               nCore=(((MaxMem*3)/4)*8)/1024
            Else
               nCore=(lBuf*8)/1024
            End If
         Else
            nCore=(lBuf*8)/1024
         End If
         nCore=((nCore+7)/8)*8
         lBuf=(1024*nCore)/(8*nBuf)
         end if
C        Write (6,*) 'OnDisk=',OnDisk
C        Write (6,*) 'Incore=',Incore
C        Write (6,*) 'nBuf=',nBuf
C        Write (6,*) 'IniBuf: nDisc=',nDisc,'MByte'
C        Write (6,*) 'nCore=',nCore,'kByte'
C        Write (6,*) 'lBuf=',lBuf
*
*------- Allocate I/O Buffer
         Call GetMem('IOBuf','Allo','Real',ipBuf,lBuf*nBuf)
*debug
c        call dcopy_(lBuf*nBuf,Work(ipBuf),1,Zero,0)
      End If
*     Write (6,*) 'ipBuf=',ipBuf
*
#ifdef _DEBUG_
      Call QExit('IniBuf')
#endif
      Return
      End

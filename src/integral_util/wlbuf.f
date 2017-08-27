************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine WLBuf
      Implicit Real*8 (a-h,o-z)
#include "IOBuf.fh"
#include "SysDef.fh"
#include "WrkSpc.fh"
*
      jpBuf(i,j)=(j-1)*lBuf+i-1+ipBuf
*
      If (iStatIO.eq.Mode_Read) Then
*        Write (6,*) 'In WLbuf'
         If (OnDisk) Call new_EAFWait(LuTmp,id)
         Return
      End If
c     Disk_Save=Disk
*     Write (6,*) 'Enter WLBuf: Disk,iPos,iBuf=',Disk,iPos,iBuf
      If (InCore.and.iBuf.eq.2) Then
         Call WarningMessage(2,
     &               'Error in in-core semi-direct implementation')
         Call Abend()
      End If
*
*---- If any data in buffer write buffer to disk.
*
      If (OnDisk) Then
*        Write (6,*) 'In WLbuf'
         Call new_EAFWait(LuTmp,id)
      End If
      If (iPos.ne.1) Then
         temp=Disk+DBLE(lBuf*RtoB)
*        Write (6,*) 'temp,DiskMx_Byte=',temp,DiskMx_Byte
         If (temp.le.DiskMx_Byte) Then
            Disk_2 = Disk_1
            Disk_1 = Disk
*           If (OnDisk) Write (*,*) 'Disk=',Disk,' lBuf*RtoI=',lBuf*RtoI
c           Write (6,*) 'WLBuf write on disk @',Disk,'iBuf=',iBuf
            If (OnDisk) Call new_EAFWrite(LuTmp,Work(jpBuf(1,iBuf)),
     &                                lBuf*RtoI,Disk)
*---------- Put a dummy record at the end
            temp=Disk+DBLE(lBuf*RtoB)
*           Write (6,*) 'temp,DiskMx_Byte=',temp,DiskMx_Byte
            If (temp.le.DiskMx_Byte.and.OnDisk) Then
c              Write (6,*) 'WLBuf write on disk @',Disk,'iBuf=',iBuf
               iZero=0
               Call ICopy(lBuf*RtoI,iZero,0,Work(jpBuf(1,iBuf)),1)
               Call
     &new_EAFWrite(LuTmp,Work(jpBuf(1,iBuf)),lBuf*RtoI,Disk)
            End If
         Else
            Call WarningMessage(2,'WLBuf: Disc is full!')
            Write (6,*) 'temp           =',temp
            Write (6,*) 'DiskMx_Byte    =',DiskMx_Byte
            Call FastIO('STATUS')
            Call Abend()
         End If
      End If
      iPos = 1
*
c     If (Disk_save.ne.Disk) Then
c        Write (6,*) 'Enter WLBuf: Disk @:',Disk_Save
c        Write (6,*) 'Exit  WLBuf: Disk @:',Disk
c     End If
*     Write (*,*) 'Exit WLBuf: Disk,iPos,iBuf=',Disk,iPos,iBuf
      Return
      End

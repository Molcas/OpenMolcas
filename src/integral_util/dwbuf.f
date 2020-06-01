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
      Subroutine dWBuf(Array,nArray)
      Use dEAF
      Use IOBUF
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
      Real*8 Array(nArray)
*
*     Write (6,*) 'Enter WBuf: iPos @',iPos,' iBuf,lBuf=',iBuf,lBuf
      If (InCore.and.iBuf.eq.2) Then
         Call WarningMessage(2,
     &               'Error in in-core semi-direct implementation')
         Call Abend()
      End If
      IODone=.True.
      iArray = 1
      mArray=nArray
10       Left = lBuf-iPos+1
         If (mArray.gt.Left) Then
            Call dCopy_(Left,Array(iArray),1,Buffer(iPos,iBuf),1)
            iArray=iArray+Left
            mArray=mArray-Left
            iPos=1
*
*---------- Wait for previous buffer to complete I/O.
*           Disk=32 after writing the control!
*
            If (Disk.ne.32.0D0.and.OnDisk) Then
*              Write (6,*) 'In dwbuf'
               Call EAFWait(LuTmp,id)
            End If
*
*---------- Put current buffer on disk and change buffer.
*
            temp=Disk+DBLE(lBuf*RtoB)
*           Write (6,*) 'temp=',temp
            If (temp.le.DiskMx_Byte) Then
               Disk_2 = Disk_1
               Disk_1 = Disk
c              Write (6,*) 'WBuf write on disk @',Disk,'iBuf=',iBuf
               If (OnDisk) Then
*                 Write (6,*) 'In dwbuf'
                  Call dEAFAWrite(LuTmp,Buffer(1,iBuf),
     &                            lBuf*RtoI,Disk,id)
               End If
               If (iBuf.eq.1) Then
                  iBuf = 2
               Else
                  iBuf = 1
               End If
            Else
               Call WarningMessage(2,'WBuf: Disc is full!!')
               Call Abend()
            End If
         Else
*           Write (6,*) ' Add ',mArray,'elements to buffer',iPos,ibuf
            Call dCopy_(mArray,Array(iArray),1,Buffer(iPos,iBuf),1)
            iPos=iPos+mArray
            mArray=0
         End If
      If(mArray.gt.0) goto 10
*
*     Write (6,*) 'Exit WBuf: iPos @',iPos,'iBuf=',iBuf
*     Call GetMem('WBuf','Check','Real',iDum,iDum)
      Return
      End

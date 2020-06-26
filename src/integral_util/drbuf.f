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
      Subroutine dRBuf(Array,nArray,Copy)
      Use dEAF
      use IOBUF
      Implicit Real*8 (a-h,o-z)
      Logical Copy
#include "SysDef.fh"
      Real*8 Array(nArray)
*
*     Write (6,*) 'Enter RBuf: iPos @',iPos,'iBuf=,lBuf',iBuf,lBuf
      If (InCore.and.iBuf.eq.2) Then
         Call WarningMessage(2,
     &                    'Error in in-core semi-direct implementation')
         Call Abend()
      End If
      iArray = 1
      mArray=nArray
chjw  Do While (mArray.ge.1)
10       If (iPos.eq.1) Then
*
*---------- Wait for pending buffer.
*
            If (OnDisk) Then
*              Write (6,*) 'In drbuf.'
               Call EAFWait(LuTmp,id)
            End If
*
*---------- Get the next buffer, make sure that request is not beyond
*           the disc limitation.
*
            temp=Disk+DBLE(lBuf*RtoB)
*           Write (6,*) 'temp=',temp
            If (temp.le.DiskMx_Byte) Then
               If (iBuf.eq.1) Then
                  jBuf=2
               Else
                  jBuf=1
               End If
               Disk_2 = Disk_1
               Disk_1 = Disk
c              Write (6,*) 'RBuf aread on disk @',Disk,'jBuf=',jBuf
               If (OnDisk) Then
*                 Write (6,*) 'In drbuf.'
                  Call dEAFARead(LuTmp,Buffer(1,jBuf),lBuf*RtoI,Disk,id)
               End If
            End If
         End If
         Left = lBuf-iPos+1
         If (mArray.gt.Left) Then
            If (Copy) Call dCopy_(Left,Buffer(iPos,iBuf),1,
     &                            Array(iArray),1)
            iArray=iArray+Left
            mArray=mArray-Left
            iPos=1
*           Write (6,*) 'LuTmp,Disk=',LuTmp,Disk
*---------- Swap buffer
            If (iBuf.eq.1) Then
               iBuf=2
            Else
               iBuf=1
            End If
         Else
*           Write (6,*) ' Copy ',mArray,'elements from buffer',iPos
            If (Copy) Call dCopy_(mArray,Buffer(iPos,iBuf),1,
     &                            Array(iArray),1)
            iPos=iPos+mArray
            mArray=0
         End If
      If(mArray.gt.0) goto 10
chjw  End Do
*
C     Write (6,*) 'Exit RBuf: iPos @',iPos,'iBuf=',iBuf
      Return
      End

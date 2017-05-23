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
      SubRoutine GenBMp(irc,X,n,m,Lunit,nStp,StpSiz,Color)
      Implicit None
      Integer irc, n, m, Lunit, nStp
      Real*8  X(n,m), StpSiz
      Character*1 Color
#include "WrkSpc.fh"

      Character*6 SecNam
      Parameter (SecNam = 'GenBMp')

      Integer i, j, iBin, nBin, ipBin, ipBMp, nCh
      Real*8  Step, absX
      Character*1 myColor

      Integer nBin_Def, iUpLim, i0, i255
      Real*8  Step_Def
      Character*1 Color_Def

      Integer  iRnge
      External iRnge

C     Defaults and limits.
C     --------------------

      nBin_Def  = 5
      Step_Def  = 1.0d-2
      Color_Def = 'R'
      iUpLim    = 999999
      i0        = 0
      i255      = 255

C     Check input.
C     ------------

      irc = 0
      If (n.lt.1 .or. m.lt.1) Then
         Return
      End If
      If (n.gt.iUpLim .or. m.gt.iUpLim) Then
         irc = 1
         Return
      End If
      If (Lunit .lt. 1) Then
         irc = 2
         Return
      End If

      If (nStp.lt.2 .or. nStp.gt.256) Then
         nBin = nBin_Def
      Else
         nBin = nStp
      End If
      If (StpSiz .le. 0.0d0) Then
         Step = Step_Def
      Else
         Step = StpSiz
      End If
      myColor = Color
      Call UpCase(myColor)
      If (myColor.ne.'R' .and. myColor.ne.'G' .and. myColor.ne.'B') Then
         myColor = Color_Def
      End If

C     Set up bins.
C     ------------

      nCh = 255/(nBin-1)

      Call GetMem('Bins','Allo','Real',ipBin,nBin)
      Call GetMem('iBMp','Allo','Inte',ipBMp,nBin)

      Work(ipBin) = 1.0d0
      Do iBin = 2,nBin-1
         Work(ipBin-1+iBin)  = Work(ipBin+iBin-2)*Step
      End Do
      Work(ipBin-1+nBin) = -1.0d0

      iWork(ipBMp-1+nBin) = 255
      Do iBin = nBin-1,1,-1
         iWork(ipBMp-1+iBin) = iWork(ipBMp+iBin) - nCh
      End Do

C     Generate bitmap file.
C     Note the special loop structure.
C     --------------------------------

      Write(Lunit,'(2(1X,I6))') m,n ! #col,#row
      If (myColor .eq. 'R') Then ! red
         Do i = n,1,-1
            Do j = 1,m
               absX = abs(X(i,j))
               iBin = iRnge(absX,Work(ipBin),nBin)
               If (iWork(ipBmp+iBin-1) .eq. 255) Then
                  Write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
               Else
                  Write(Lunit,'(4(1X,I3))') iWork(ipBmp+iBin-1),i0,i0,i0
               End If
            End Do
         End Do
      Else If (myColor .eq. 'G') Then ! green
         Do i = n,1,-1
            Do j = 1,m
               absX = abs(X(i,j))
               iBin = iRnge(absX,Work(ipBin),nBin)
               If (iWork(ipBmp+iBin-1) .eq. 0) Then
                  Write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
               Else
                  Write(Lunit,'(4(1X,I3))') i0,iWork(ipBmp+iBin-1),i0,i0
               End If
            End Do
         End Do
      Else If (myColor .eq. 'B') Then ! blue
         Do i = n,1,-1
            Do j = 1,m
               absX = abs(X(i,j))
               iBin = iRnge(absX,Work(ipBin),nBin)
               If (iWork(ipBmp+iBin-1) .eq. 0) Then
                  Write(Lunit,'(4(1X,I3))') i255,i255,i255,i0
               Else
                  Write(Lunit,'(4(1X,I3))') i0,i0,iWork(ipBmp+iBin-1),i0
               End If
            End Do
         End Do
      Else
         Call SysAbendMsg(SecNam,'Logical error!',
     &                    '(Should never happen)')
      End If

C     De-allocations.
C     ---------------

      Call GetMem('iBMp','Free','Inte',ipBMp,nBin)
      Call GetMem('Bins','Free','Real',ipBin,nBin)

      End

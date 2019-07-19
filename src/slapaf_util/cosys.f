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
      Subroutine CoSys(Cent,R,xyz)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 Cent(3,3), r(3), xyz(3,2)
      Integer iComp(3)
      Logical Linear
*
      iRout=221
      iPrint=nPrint(iRout)
      Lu=6
      If (iPrint.ge.99) Then
         Call RecPrt('CoSys: Cent',' ',Cent,3,3)
      End If
*
*---- Check if linear
*
      Co=Zero
      Crap=Zero
      RR1=Zero
      RR2=Zero
      ThrAcos=1.0D-6
      Do i = 1, 3
         Co = Co + (Cent(i,1)-Cent(i,2))*(Cent(i,3)-Cent(i,2))
         Crap = Crap + ( (Cent(i,3)-Cent(i,2))
     &                  +(Cent(i,1)-Cent(i,2)) )**2
         RR1=RR1+(Cent(i,1)-Cent(i,2))**2
         RR2=RR2+(Cent(i,3)-Cent(i,2))**2
      End Do
      Co=Co/Sqrt(RR1*RR2)
      Crap=Crap/Sqrt(RR1*RR2)
      If (iPrint.ge.99) Then
         Write (6,*) 'Co=',Co
         Write (6,*) 'Crap,Sqrt(Crap)=',Crap,Sqrt(Crap)
      End If
      If (Sqrt(Crap).lt.1.0D-6) Then
         Fi=Pi-ArSin(Sqrt(Crap))
         Si=Sqrt(Crap)
      Else
       if(Co.gt.One.and.Co.lt.One+ThrAcos) Co=One
       if(Co.lt.-One.and.Co.gt.-One-ThrAcos) Co=-One
       if(Co.gt.One+ThrAcos.or.Co.lt.-One-ThrAcos) then
         Call WarningMessage(2,'Error in CoSys')
         write(6,*) 'Error in cosys: arcos(',Co,')'
         Call Abend
       endif
         Fi=ArCos(Co)
         Si=Sqrt(One-Co**2)
      End If
      If (iPrint.ge.99) Then
         Write (6,*) 'Fi,Pi=',Fi,Pi
         Write (6,*) 'Pi-Fi=',Pi-Fi
      End If
*
      Linear=Abs(Pi-Fi).lt.1.0D-13
*
*---- Form reference axis
*
      RR=Zero
      Do i = 1, 3
         R(i)=Cent(i,3)-Cent(i,1)
         RR=RR+R(i)**2
      End Do
*
      Call DScal_(3,One/Sqrt(RR),R,1)
      If (iPrint.ge.99) Then
         Write (6,*) 'Linear=',Linear
         Write (6,*) 'RR=',RR
         Call RecPrt('R',' ',R,3,1)
      End If
*
 99   Continue
      If (Linear) Then
*
      nComp=0
      Do i = 1, 3
         If (R(i).eq.Zero) Then
            nComp=nComp+1
            iComp(nComp)=i
         End If
      End Do
*
*---- Compute the WDC B-matrix
*
      call dcopy_(6,[Zero],0,xyz,1)
      If (nComp.eq.0) Then
*        Write (*,*) ' Case nComp.eq.0'
*
         xyz(1,1)=R(1)
         xyz(2,1)=R(2)
         xyz(3,1)=-R(3)
         r12=R(1)**2+R(2)**2-R(3)**2
         r11=R(1)**2+R(2)**2+R(3)**2
         xyz(1,1)=xyz(1,1)-(r12/r11)*R(1)
         xyz(2,1)=xyz(2,1)-(r12/r11)*R(2)
         xyz(3,1)=xyz(3,1)-(r12/r11)*R(3)
         r2=Sqrt(xyz(1,1)**2+xyz(2,1)**2+xyz(3,1)**2)
         xyz(1,1)=xyz(1,1)/r2
         xyz(2,1)=xyz(2,1)/r2
         xyz(3,1)=xyz(3,1)/r2
         xyz(1,2)=R(2)*xyz(3,1)-R(3)*xyz(2,1)
         xyz(2,2)=R(3)*xyz(1,1)-R(1)*xyz(3,1)
         xyz(3,2)=R(1)*xyz(2,1)-R(2)*xyz(1,1)
         r2=Sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
         xyz(1,2)=xyz(1,2)/r2
         xyz(2,2)=xyz(2,2)/r2
         xyz(3,2)=xyz(3,2)/r2
*
      Else If (nComp.eq.1) Then
*        Write (*,*) ' Case nComp.eq.1'
*
         i=iComp(1)
         xyz(i,1)=One
         xyz(1,2)=R(2)*xyz(3,1)-R(3)*xyz(2,1)
         xyz(2,2)=R(3)*xyz(1,1)-R(1)*xyz(3,1)
         xyz(3,2)=R(1)*xyz(2,1)-R(2)*xyz(1,1)
         r2=Sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
         xyz(1,2)=xyz(1,2)/r2
         xyz(2,2)=xyz(2,2)/r2
         xyz(3,2)=xyz(3,2)/r2
*
      Else If (nComp.eq.2) Then
*        Write (*,*) ' Case nComp.eq.2'
*
         i=iComp(1)
         xyz(i,1)=One
         i=iComp(2)
         xyz(i,2)=One
*
      Else
         Call WarningMessage(2,'Error in CoSys')
         Write (Lu,*) ' CoSys: nComp.eq.3'
         Call Abend()
      End If
*
      Else     ! Non-linear
*
*        Form the cross product R12xR32
         RR=Zero
         Do i = 1, 3
            If (i.eq.1) Then
               j=2
               k=3
            Else If (i.eq.2) Then
               j=3
               k=1
            Else
               j=1
               k=2
            End If
C           j=i+1
C           If (j.gt.3) j=j-3
C           k=j+1
C           If (k.gt.3) k=k-3
            R21j=Cent(j,1)-Cent(j,2)
            R21k=Cent(k,1)-Cent(k,2)
            R23j=Cent(j,3)-Cent(j,2)
            R23k=Cent(k,3)-Cent(k,2)
            xyz(i,2) = R21j*R23k-R21k*R23j
            RR = RR + xyz(i,2)**2
         End Do
         If (RR.eq.Zero) Then
            Linear=.True.
            Go To 99
         End If
         Call DScal_(3,One/Sqrt(RR),xyz(1,2),1)
         If (iPrint.ge.99) Write (6,*) 'RR=',RR
*
         RR=Zero
         Do i = 1, 3
            If (i.eq.1) Then
               j=2
               k=3
            Else If (i.eq.2) Then
               j=3
               k=1
            Else
               j=1
               k=2
            End If
C           j=i+1
C           If (j.gt.3) j=j-3
C           k=j+1
C           If (k.gt.3) k=k-3
            xyz(i,1) = xyz(j,2)*R(k)-xyz(k,2)*R(j)
            RR = RR + xyz(i,1)**2
         End Do
         Call DScal_(3,One/Sqrt(RR),xyz(1,1),1)
         If (iPrint.ge.99) Write (6,*) 'RR=',RR
         If (iPrint.ge.99) Then
            Call RecPrt('xyz',' ',xyz,3,2)
         End If
*
      End if
      If (iPrint.ge.99) Then
         Call RecPrt(' Reference Axis',' ',R,3,1)
         Call RecPrt(' Perpendicular Axes',' ',xyz,3,2)
      End If
*
      Return
      End

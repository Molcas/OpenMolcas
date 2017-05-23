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
      Function HSR(xyz,nCent)
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "info_slapaf.fh"
      Real*8   xyz(3,nCent)
*
      xyz0(i,j)=Work(ipRef+(j-1)*3+i-1)
*
*     Call QEnter('HSR')
      iRout=54
      iPrint=nPrint(iRout)
*
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the radius of the hypersphere
*
C     Call RecPrt('HSR: xyz',' ',xyz,3,nCent)
C     Call RecPrt('HSR: xyz0',' ',Work(ipRef),3,nCent)
      Lu=6
      HSR=Zero
      TMass=Zero
      Do iCent = 1, nCent
         Fact=DBLE(iDeg(xyz(1,iCent),iOper,nSym))
         xMass=Fact*Work(ipCM+iCent-1)
         TMass=TMass+xmass
C        Write (*,*) 'xMass=',xMass
         Do ixyz = 1, 3
C           Write (*,*)xyz(ixyz,iCent),xyz0(ixyz,iCent)
            temp=xyz(ixyz,iCent)-xyz0(ixyz,iCent)
            HSR=HSR + xMass*temp**2
         End Do
      End Do
      If (HSR.ne.0d0) HSR=Sqrt(HSR/TMass)
*
      If (iPrint.ge.5) Then
         Write (Lu,*)
         Write (Lu,'(16X,A)') '*************************'//
     &                        '*************************'
         Write (Lu,'(16X,A)') '* Radius of hypersphere /'//
     &                        ' au*amu**(1/2)/amu(1/2) *'
         Write (Lu,'(16X,A)') '*************************'//
     &                        '*************************'
         Write (Lu,'(25X,F18.4,A)') HSR
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Call QExit('HSR')
      Return
      End

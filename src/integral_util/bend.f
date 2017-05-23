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
      Subroutine Bend(xyz,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB)
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
      Real*8   Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent),
     &        BRij(3,2), dBRij(3,2,3,2),
     &        BRjk(3,2), dBRjk(3,2,3,2)
      Logical lWrite, ldB, lWarn
      Character*8 Label
*                                                                      *
************************************************************************
*                                                                      *
*define _TIME_
*                                                                      *
************************************************************************
*                                                                      *
*
#ifdef _TIME_
      Call QEnter('Bend')
#endif
*
      mCent=2
      Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
      Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
      Co=Zero
      Crap=Zero
*     BRij and BRjk should be normalized
      Do i = 1, 3
         Co=Co+BRij(i,1)*BRjk(i,2)
         Crap=Crap+(BRjk(i,2)+BRij(i,1))**2
      End Do
      Crap = Sqrt(Crap)
*
*.... Special care for cases close to linearity
*
      If (Crap.lt.1.0D-4) Then
         Si=Half*Crap
         Fir=Pi-ArSin(Si)
      Else
         Fir=ArCos(Co)
         If (Abs(Co).gt.One) Co=Sign(One,Co)
         Si=Sqrt(One-Co**2)
      End If
*
      If (Fir.lt.1.0d-13) Then
         Fir=Zero
#ifdef _TIME_
         Call QExit('Bend')
#endif
         Return
      Else If (Abs(Fir-Pi).lt.1.0d-13) Then
         Fir=Pi
#ifdef _TIME_
         Call QExit('Bend')
#endif
         Return
      End If
      dFir=180.0D0*Fir/Pi
      If ((Abs(dFir).gt.177.5 .or. Abs(dFir).lt.2.5).and.lWarn)
     &   Write (6,*) ' Valence angle close to end in '//
     &               'range of definition'
      If (lWrite) Write (6,'(1X,A,A,F10.4,A,F10.6,A)') Label,
     &            ' : Angle=    ', dFir,'   / Degree  ',Fir,' / rad'
*
*---- Compute the WDC B-matrix
*
      If (Si.eq.Zero) Then
*------- Dummy assignment for a linear system!
         call dcopy_(3*nCent,0.0D0,0,Bf,1)
      Else
         Do i = 1, 3
            Bf(i,1)= (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
            Bf(i,3)= (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
*.......    Utilize translational invariance.
            Bf(i,2)=-(Bf(i,1)+Bf(i,3))
         End Do
      End If
*     Call RecPrt('Bf',' ',Bf,9,1)
*
*---- Compute the cartesian derivative of the B-Matrix.
*
      If (ldB) Then
*
*        dBf=-11.11111
         If (Si.eq.Zero) Then
            Call WarningMessage(2,'Bend: Si.eq.0.0D')
            Call Abend()
         End If
         Do i = 1, 3
            Do j = 1, i
               dBf(i,1,j,1)=( -Si*Bf(i,1)*BRij(j,1)
     &                        +Co*dBRij(i,1,j,1)
     &                        -Bf(j,1)*(Co*Bf(i,1)*Rij1
     &                        +Si*BRij(i,1)) ) / (Si*Rij1)
               dBf(i,1,j,3)=(-Si*Bf(i,1)*BRjk(j,2)
     &                       +dBRij(i,1,j,2)
     &                       -Bf(j,3)*Co*Bf(i,1)*Rjk1)
     &                       / (Si*Rjk1)
*              Write (*,*) '13',dBf(i,1,j,3), i, j
               dBf(i,3,j,1)=(-Si*Bf(i,3)*BRij(j,1)
     &                       +dBRjk(i,2,j,1)
     &                       -Bf(j,1)*Co*Bf(i,3)*Rij1)
     &                       / (Si*Rij1)
               dBf(i,3,j,3)=( -Si*Bf(i,3)*BRjk(j,2)
     &                        +Co*dBRjk(i,2,j,2)
     &                        -Bf(j,3)*(Co*Bf(i,3)*Rjk1
     &                        +Si*BRjk(i,2)) ) / (Si*Rjk1)
*
               dBf(j,1,i,1)=dBf(i,1,j,1)
               dBf(j,3,i,1)=dBf(i,1,j,3)
               dBf(j,1,i,3)=dBf(i,3,j,1)
               dBf(j,3,i,3)=dBf(i,3,j,3)
*
               dBf(i,1,j,2)=-(dBf(i,1,j,1)+dBf(i,1,j,3))
               dBf(j,2,i,1)=dBf(i,1,j,2)
               dBf(j,1,i,2)=-(dBf(j,1,i,1)+dBf(j,1,i,3))
               dBf(i,2,j,1)=dBf(j,1,i,2)
               dBf(i,3,j,2)=-(dBf(i,3,j,1)+dBf(i,3,j,3))
               dBf(j,2,i,3)=dBf(i,3,j,2)
               dBf(j,3,i,2)=-(dBf(j,3,i,1)+dBf(j,3,i,3))
               dBf(i,2,j,3)=dBf(j,3,i,2)
*
               dBf(i,2,j,2)=-(dBf(i,2,j,1)+dBf(i,2,j,3))
               dBf(j,2,i,2)=dBf(i,2,j,2)
*
            End Do
         End Do
*        Call RecPrt('dBf','(9F9.1)',dBf,9,9)
*
      End If
*
#ifdef _TIME_
      Call QExit('Bend')
#endif
      Return
      End

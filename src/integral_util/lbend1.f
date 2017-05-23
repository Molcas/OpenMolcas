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
      Subroutine LBend1(Cent,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB,
     &                  Axis,Perp_Axis1)
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8  Bf(3,nCent), xxx(3,3), dBf(3,nCent,3,nCent),
     &        BRij(3,2), dBRij(3,2,3,2),
     &        BRjk(3,2), dBRjk(3,2,3,2), Fir,
     &        uMtrx(3,3), uVec(3,3), Scr1(3,3), Scr2(3,3),
     &        Axis(3),Perp_Axis1(3), Cent(3,3)
      Logical lWrite, ldB, lWarn, Linear
      Character*8 Label
*
      iRout=220
      iPrint=nPrint(iRout)
*
      Lu=6
*
*     Call QEnter('LBend1')
*
      If (iPrint.ge.99) Then
         Call RecPrt('LBend1: Axis',' ',Axis,3,1)
         Call RecPrt('LBend1: Perp_Axis1',' ',Perp_Axis1,3,1)
      End If
*
      call dcopy_(3,Axis,      1,uVec(1,1),1)
      call dcopy_(3,Perp_Axis1,1,uVec(1,2),1)
      call dcopy_(3,Zero,0,uVec(1,3),1)
*
*---- Project the coordinates to the plane
*
      Call DGEMM_('T','N',
     &            3,3,3,
     &            1.0d0,uVec,3,
     &            Cent,3,
     &            0.0d0,xxx,3)
      xxx(3,1)=Zero
      xxx(3,2)=Zero
      xxx(3,3)=Zero
      If (iPrint.ge.99) Then
         Call RecPrt('Original coordinates','(3F24.12)',Cent,3,3)
         Call RecPrt('uVec',' ',uVec,3,3)
         Call RecPrt('Projected coordinates','(3F24.12)',xxx,3,3)
      End If
*
      mCent=2
      Call Strtch(xxx(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
      Call Strtch(xxx(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
*
*---- We better be very careful here in order not to lose accuracy!
*
      Co=Zero
      Crap=Zero
      Do i = 1, 3
         Co=Co+BRij(i,1)*BRjk(i,2)
         Crap = Crap + (BRjk(i,2)+BRij(i,1))**2
      End Do
      Linear=.True.
      If (iPrint.ge.99) Then
         Call RecPrt('BRij','(3F24.12)',BRij,3,2)
         Call RecPrt('BRjk','(3F24.12)',BRjk,3,2)
         Write (6,*) ' Diff=',Abs(ArCos(Co)-Pi)
         Write (6,'(A,F24.16)') ' Co=',Co
         Write (6,'(A,F24.16)') ' Sqrt(Crap)=',Sqrt(Crap)
      End If
*
*.... Special care for cases close to linearity
*
      If (Sqrt(Crap).lt.1.0D-4) Then
         Fir=Pi-ArSin(Sqrt(Crap))
         Si=Sqrt(Crap)
      Else
         Fir=ArCos(Co)
         Si=Sqrt(One-Co**2)
      End If
*
      If (Abs(Fir-Pi).gt.1.0D-13) Then
         If (iPrint.ge.99)
     &      Write (Lu,*) ' LBend: Use none linear formulae'
         Linear=.False.
      Else
         If (iPrint.ge.99) Write (Lu,*) ' LBend: Use linear formulae'
      End If
*
      dFir=180.0D0*Fir/Pi
      If (lWrite) Then
         Write (Lu,'(1X,A,A,F10.6,A,F12.8,A)') Label,
     &            ' : Projected Angle=', dFir,'/degree, ',Fir,'/rad'
      End If
*
      call dcopy_(9,Zero,0,uMtrx,1)
      If (Linear) Then
         uMtrx(2,1)=-One/Rij1
      Else
         Do i = 1, 3
            uMtrx(i,1)=(Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
         End Do
      End If
      Call DGEMM_('N','N',
     &            3,1,3,
     &            1.0d0,uVec,3,
     &            uMtrx,3,
     &            0.0d0,Scr2,3)
      Bf(1,1)=Scr2(1,1)
      Bf(2,1)=Scr2(2,1)
      Bf(3,1)=Scr2(3,1)
*
      call dcopy_(9,Zero,0,uMtrx,1)
      If (Linear) Then
         uMtrx(2,1)=-One/Rjk1
      Else
         Do i = 1, 3
            uMtrx(i,1)=(Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
         End Do
      End If
      Call DGEMM_('N','N',
     &            3,1,3,
     &            1.0d0,uVec,3,
     &            uMtrx,3,
     &            0.0d0,Scr2,3)
      Bf(1,3)=Scr2(1,1)
      Bf(2,3)=Scr2(2,1)
      Bf(3,3)=Scr2(3,1)
*
      Do i = 1, 3
*....... Utilize translational invariance.
         Bf(i,2)=-(Bf(i,1)+Bf(i,3))
      End Do
*
*---- Compute the cartesian derivative of the B-Matrix.
*
      If (ldB) Then
*
*........ 1,1 Block
*
         call dcopy_(9,Zero,0,uMtrx,1)
         If (Linear) Then
            uMtrx(1,2)=-One/Rij1**2
            uMtrx(2,1)=uMtrx(1,2)
         Else
            Do i = 1, 2
               Bfi1=(Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
               Do j = 1, 2
                  Bfj1=(Co*BRij(j,1)-BRjk(j,2))/(Si*Rij1)
                  uMtrx(i,j)=( -Si*Bfi1*BRij(j,1)
     &                      +Co*dBRij(i,1,j,1)
     &                      -Bfj1*(Co*Bfi1*Rij1
     &                      +Si*BRij(i,1)) ) / (Si*Rij1)
               End Do
            End Do
         End If
         Call DGEMM_('N','T',
     &               3,3,3,
     &               1.0d0,uMtrx,3,
     &               uVec,3,
     &               0.0d0,Scr1,3)
         Call DGEMM_('N','N',
     &               3,3,3,
     &               1.0d0,uVec,3,
     &               Scr1,3,
     &               0.0d0,Scr2,3)
         dBf(1,1,1,1)=Scr2(1,1)
         dBf(2,1,1,1)=Scr2(2,1)
         dBf(2,1,2,1)=Scr2(2,2)
         dBf(3,1,1,1)=Scr2(3,1)
         dBf(3,1,2,1)=Scr2(3,2)
         dBf(3,1,3,1)=Scr2(3,3)
*
*....... 1,3 Block
*
         If (Linear) Then
            dBf(1,1,1,3)=Zero
            dBf(2,1,1,3)=Zero
            dBf(2,1,2,3)=Zero
            dBf(3,1,1,3)=Zero
            dBf(3,1,2,3)=Zero
            dBf(3,1,3,3)=Zero
         Else
            call dcopy_(9,Zero,0,uMtrx,1)
            Do i = 1, 2
               Bfi1=(Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
               Do j = 1, 2
                  Bfj3= (Co*BRjk(j,2)-BRij(j,1))/(Si*Rjk1)
                  uMtrx(i,j)=(-Si*Bfi1*BRjk(j,2)
     &                      +dBRij(i,1,j,2)
     &                      -Bfj3*Co*Bfi1*Rjk1)
     &                      / (Si*Rjk1)
               End Do
            End Do
            Call DGEMM_('N','T',
     &                  3,3,3,
     &                  1.0d0,uMtrx,3,
     &                  uVec,3,
     &                  0.0d0,Scr1,3)
            Call DGEMM_('N','N',
     &                  3,3,3,
     &                  1.0d0,uVec,3,
     &                  Scr1,3,
     &                  0.0d0,Scr2,3)
            dBf(1,1,1,3)=Scr2(1,1)
            dBf(2,1,1,3)=Scr2(2,1)
            dBf(2,1,2,3)=Scr2(2,2)
            dBf(3,1,1,3)=Scr2(3,1)
            dBf(3,1,2,3)=Scr2(3,2)
            dBf(3,1,3,3)=Scr2(3,3)
         End If
*
*....... 3,1 Block
*
         If (Linear) Then
            dBf(1,3,1,1)=Zero
            dBf(2,3,1,1)=Zero
            dBf(2,3,2,1)=Zero
            dBf(3,3,1,1)=Zero
            dBf(3,3,2,1)=Zero
            dBf(3,3,3,1)=Zero
         Else
            call dcopy_(9,Zero,0,uMtrx,1)
            Do i = 1, 2
               Bfi3= (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
               Do j = 1, 2
                  Bfj1=(Co*BRij(j,1)-BRjk(j,2))/(Si*Rij1)
                  uMtrx(i,j)=(-Si*Bfi3*BRij(j,1)
     &                      +dBRjk(i,2,j,1)
     &                      -Bfj1*Co*Bfi3*Rij1)
     &                      / (Si*Rij1)
               End Do
            End Do
            Call DGEMM_('N','T',
     &                  3,3,3,
     &                  1.0d0,uMtrx,3,
     &                  uVec,3,
     &                  0.0d0,Scr1,3)
            Call DGEMM_('N','N',
     &                  3,3,3,
     &                  1.0d0,uVec,3,
     &                  Scr1,3,
     &                  0.0d0,Scr2,3)
            dBf(1,3,1,1)=Scr2(1,1)
            dBf(2,3,1,1)=Scr2(2,1)
            dBf(2,3,2,1)=Scr2(2,2)
            dBf(3,3,1,1)=Scr2(3,1)
            dBf(3,3,2,1)=Scr2(3,2)
            dBf(3,3,3,1)=Scr2(3,3)
         End If
*
*....... 3,3 Block
*
         call dcopy_(9,Zero,0,uMtrx,1)
         If (Linear) Then
            uMtrx(1,2)= One/Rjk1**2
            uMtrx(2,1)= uMtrx(1,2)
         Else
            Do i = 1, 2
               Bfi3= (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
               Do j = 1, 2
                  Bfj3= (Co*BRjk(j,2)-BRij(j,1))/(Si*Rjk1)
                  uMtrx(i,j)=( -Si*Bfi3*BRjk(j,2)
     &                      +Co*dBRjk(i,2,j,2)
     &                      -Bfj3*(Co*Bfi3*Rjk1
     &                      +Si*BRjk(i,2)) ) / (Si*Rjk1)
               End Do
            End Do
         End If
         Call DGEMM_('N','T',
     &               3,3,3,
     &               1.0d0,uMtrx,3,
     &               uVec,3,
     &               0.0d0,Scr1,3)
         Call DGEMM_('N','N',
     &               3,3,3,
     &               1.0d0,uVec,3,
     &               Scr1,3,
     &               0.0d0,Scr2,3)
         dBf(1,3,1,3)=Scr2(1,1)
         dBf(2,3,1,3)=Scr2(2,1)
         dBf(2,3,2,3)=Scr2(2,2)
         dBf(3,3,1,3)=Scr2(3,1)
         dBf(3,3,2,3)=Scr2(3,2)
         dBf(3,3,3,3)=Scr2(3,3)
*
         Do i = 1, 3
            Do j = 1, i
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
*
      End If
*
      If (iPrint.ge.99) Then
         Call RecPrt('Bf',' ',Bf,3,nCent)
         If (ldB) Then
            Call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)
         End If
      End If
*     Call QExit('LBend1')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(lWarn)
      End

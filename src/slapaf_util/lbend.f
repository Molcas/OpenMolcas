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
      Subroutine LBend(Cent,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB,
     &                  Axis,Perp_Axis1,Force)
      Implicit Real*8  (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8  Bf(3,nCent), xxx(3,3), dBf(3,nCent,3,nCent),
     &        BRij(3,2), dBRij(3,2,3,2),
     &        BRjk(3,2), dBRjk(3,2,3,2), Fir,
     &        uMtrx(3,3), uVec(3,3), Scr1(3,3), Scr2(3,3),
     &        Axis(3),Perp_Axis1(3), Cent(3,3)
      Logical lWrite, ldB, lWarn, Linear, Force
      Character*8 Label
*
      iRout=220
      iPrint=nPrint(iRout)
*
      Lu=6
*
*     Call QEnter('LBend')
*
      If (iPrint.ge.99) Then
         Write(6,*) 'LBend: Force ',Force
         Call RecPrt('LBend: Axis',' ',Axis,3,1)
         Call RecPrt('LBend: Perp_Axis1',' ',Perp_Axis1,3,1)
      End If
*
      call dcopy_(3,Axis,      1,uVec(1,1),1)
      call dcopy_(3,Perp_Axis1,1,uVec(1,2),1)
      call dcopy_(3,[Zero],0,uVec(1,3),1)
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
*.... Swap atoms to ensure the complementary angle is always Pi
*
      Middle=2
      If (Force) Then
         R1=(xxx(1,1)-xxx(1,2))**2+(xxx(2,1)-xxx(2,2))**2
         R2=(xxx(1,2)-xxx(1,3))**2+(xxx(2,2)-xxx(2,3))**2
         R3=(xxx(1,3)-xxx(1,1))**2+(xxx(2,3)-xxx(2,1))**2
         If ((R1.ge.R3).and.(R1.ge.R2)) Then
           Middle=3
         Else If (R2.ge.R3) Then
           Middle=1
         End If
      End If
      If (Middle.ne.2) Then
         Call DSwap_(3,xxx(1,2),1,xxx(1,Middle),1)
         If (iPrint.ge.99) Then
            Call RecPrt('Swapped coordinates','(3F24.12)',xxx,3,3)
         End If
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
      End Do
      Do i = 1, 3
         Crap = Crap + (BRjk(i,2)-Sign(One,Co)*BRij(i,1))**2
      End Do
      Crap = Sqrt(Crap)
      Linear=.True.
      If (iPrint.ge.99) Then
         Call RecPrt('BRij','(3F24.12)',BRij,3,2)
         Call RecPrt('BRjk','(3F24.12)',BRjk,3,2)
         Write (6,*) ' Rij1=',Rij1
         Write (6,*) ' Rjk1=',Rjk1
         Write (6,*) ' Diff=',Abs(ArCos(Co)-Pi)
         Write (6,'(A,F24.16)') ' Co=',Co
         Write (6,'(A,F24.16)') ' Crap=',Crap
      End If
*
*.... Special care for cases close to linearity
*
      If (Crap.lt.1.0D-4) Then
         Si=Crap
         If (Co.lt.Zero) Then
            Fir=Pi-ArSin(Si)
         Else
            Fir=ArSin(Si)
         End If
      Else
         If (Abs(Co).gt.One) Co=Sign(One,Co)
         Fir=ArCos(Co)
         Si=Sqrt(One-Co**2)
      End If
*
*     If (Abs(Fir-Pi).gt.1.0D-13) Then
      If (Abs(Si).gt.1.0D-13) Then
         If (iPrint.ge.99)
     &      Write (Lu,*) ' LBend: Use nonlinear formulae'
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
      call dcopy_(9,[Zero],0,uMtrx,1)
      If (Linear) Then
         uMtrx(2,1)=-Sign(One,Co)/Rij1
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
      call dcopy_(9,[Zero],0,uMtrx,1)
      If (Linear) Then
         uMtrx(2,1)=One/Rjk1
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
*....... 1,1 Block
*
         call dcopy_(9,[Zero],0,uMtrx,1)
         If (Linear) Then
            If (Co.gt.Zero) Then
               uMtrx(1,2)=One/Rij1**2
               If (Rij1.lt.Rjk1) uMtrx(1,2)=Two*uMtrx(1,2)
            Else
               uMtrx(1,2)=Two/Rij1**2-One/(Rij1**2+Rij1*Rjk1)
            End If
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
         call dcopy_(9,[Zero],0,uMtrx,1)
         If (Linear) Then
            If (Co.gt.Zero) Then
               If (Rij1.lt.Rjk1) Then
                  uMtrx(1,2)=-One/(Rij1*Rjk1)
                  uMtrx(2,1)=Zero
               Else
                  uMtrx(1,2)=Zero
                  uMtrx(2,1)=One/(Rij1*Rjk1)
               End If
            Else
               uMtrx(1,2)=One/(Rij1**2+Rij1*Rjk1)
               uMtrx(2,1)=-One/(Rjk1**2+Rjk1*Rij1)
            End If
         Else
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
         dBf(1,1,1,3)=Scr2(1,1)
         dBf(2,1,1,3)=Scr2(2,1)
         dBf(2,1,2,3)=Scr2(2,2)
         dBf(3,1,1,3)=Scr2(3,1)
         dBf(3,1,2,3)=Scr2(3,2)
         dBf(3,1,3,3)=Scr2(3,3)
*
*....... 3,1 Block
*
         call dcopy_(9,[Zero],0,uMtrx,1)
         If (Linear) Then
            If (Co.gt.Zero) Then
               If (Rjk1.lt.Rij1) Then
                  uMtrx(1,2)=One/(Rjk1*Rij1)
                  uMtrx(2,1)=Zero
               Else
                  uMtrx(1,2)=Zero
                  uMtrx(2,1)=-One/(Rjk1*Rij1)
               End If
            Else
               uMtrx(1,2)=-One/(Rjk1**2+Rjk1*Rij1)
               uMtrx(2,1)=One/(Rij1**2+Rij1*Rjk1)
            End If
         Else
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
         dBf(1,3,1,1)=Scr2(1,1)
         dBf(2,3,1,1)=Scr2(2,1)
         dBf(2,3,2,1)=Scr2(2,2)
         dBf(3,3,1,1)=Scr2(3,1)
         dBf(3,3,2,1)=Scr2(3,2)
         dBf(3,3,3,1)=Scr2(3,3)
*
*....... 3,3 Block
*
         call dcopy_(9,[Zero],0,uMtrx,1)
         If (Linear) Then
            If (Co.gt.Zero) Then
               uMtrx(1,2)=-One/Rjk1**2
               If (Rjk1.lt.Rij1) uMtrx(1,2)=Two*uMtrx(1,2)
            Else
               uMtrx(1,2)=-Two/Rjk1**2+One/(Rjk1**2+Rjk1*Rij1)
            End If
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
*.... Swap atoms back
*
      If (Middle.ne.2) Then
         Call DSwap_(3,Bf(1,2),1,Bf(1,Middle),1)
         If (ldB) Then
           Call DSwap_(3*nCent*3,dBf(1,1,1,2),1,dBf(1,1,1,Middle),1)
           Call DSwap_(3*nCent,dBf(1,2,1,1),3*nCent,
     &                         dBf(1,Middle,1,1),3*nCent)
           Call DSwap_(3*nCent,dBf(2,2,1,1),3*nCent,
     &                         dBf(2,Middle,1,1),3*nCent)
           Call DSwap_(3*nCent,dBf(3,2,1,1),3*nCent,
     &                         dBf(3,Middle,1,1),3*nCent)
         End If
      End If
*
      If (iPrint.ge.99) Then
         Call RecPrt('Bf',' ',Bf,3,nCent)
         If (ldB) Then
            Call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)
         End If
      End If
*     Call QExit('LBend')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_logical(lWarn)
      End

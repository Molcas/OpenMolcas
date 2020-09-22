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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RigRot(CoorIn,rM,nAtm)

************************************************************************
*                                                                      *
* Object: to compute the rotational constants and the rotational       *
*         spectrum within the rigid-rotor model.                       *
*                                                                      *
*    Reference: P.W. Atkins, in "Molecular Quantum Mechanics",         *
*               Oxford University Press, Oxford/New York,              *
*               ch. 11.2, pp. 289.                                     *
*               Michael Tinkham, in "Group Theory and Quantum          *
*               Mechanics", McGraw-Hill Book Company, New York,        *
*               ch. 7-14, pp. 250.                                     *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              DCopy  (ESSL)                                           *
*              Jacob                                                   *
*              DSwap  (ESSL)                                           *
*              Order                                                   *
*              TriPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      use Sizes_of_Seward, only: S
      use Real_Info, only: TMass
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 CoorIn(3,nAtm), rM(nAtm), XI(3)
      Real*8, Dimension(:,:), Allocatable :: Coor, Vec
      Real*8, Dimension(:), Allocatable :: Hess, En
      Character*80 Label
      Logical Linear, RR_Show
#include "constants2.fh"
*
      iRout = 117
      iPrint = nPrint(iRout)
      RR_Show = iPrint.ge.6
      if(iprintlevel(-1).lt.3) RR_Show=.false.
      Call qEnter('RigRot')
*
      If (RR_Show) Then
         Write (6,*)
         Call CollapseOutput(1,'   Rigid rotor info:')
         Write (6,'(3X,A)')    '   -----------------'
         Write (6,*)
      End If
*
      call dcopy_(9,[Zero],0,PAx,1)
      call dcopy_(3,[One],0,PAx,4)
      If (TMass.eq.Zero) Go To 99
      Linear=.False.
      If (iPrint.ge.99) Then
         Call RecPrt(' In RigRot: CoorIn',' ',CoorIn,3,nAtm)
         Call RecPrt(' In RigRot: Mass',' ',rM,1,nAtm)
      End If
      If (RR_Show) Then
      Write (6,*)
      Write (6,'(19X,A,F10.5)') ' Total mass (a) :', TMass/uToau
      Write (6,*)
      Write (6,'(19X,A)')       ' Center of mass '
      Write (6,'(19X,3A11)')    '       X   ','       Y   ',
     &                          '       Z   '
      Write (6,'(19X,3F11.5)') (CoM(i),i=1,3)
      Write (6,*)
      End If
*
*     Translate coordinate system to center of mass.
*
      Call mma_allocate(Coor,3,nAtm)
      Do 30 iCar = 1, 3
         Do 31 iAtom = 1, nAtm
            Coor(iCar,iAtom)=CoorIn(iCar,iAtom)-CoM(iCar)
 31      Continue
 30   Continue
*
      If (RR_Show) Then
      Write (6,'(19X,A)') ' Reference system based on center of mass'
      Write (6,'(19X,A)')
     &              ' Coordinates and Masses of Atoms, in au and A'
      Write (6,'(19X,A)')
     &      '       X          Y          Z        Mass'
      Write (6,'(19X,4F11.5)') ((Coor(j,i),j=1,3),rM(i)/utoau,i=1,nAtm)
      Write (6,*)
      End If
      If (nAtm.eq.1) Then
         Call mma_deallocate(Coor)
         Go To 99
      End If
*
*     Construct the moment of inertia tensor
*
      ii = 0
      Do 40 iCar = 1, 3
         Do 41 jCar = 1, iCar-1
            ii = ii + 1
            rMI(ii) = Zero
            Do 42 iAtom = 1, nAtm
               rMI(ii) = rMI(ii) -
     &                 Coor(iCar,iAtom) * Coor(jCar,iAtom) * rM(iAtom)
 42         Continue
 41      Continue
         ii = ii + 1
         rMI(ii) = Zero
         Do 43 iAtom = 1, nAtm
            rMI(ii) = rMI(ii) + rM(iAtom) * (
     &                Coor(1,iAtom)**2 +
     &                Coor(2,iAtom)**2 +
     &                Coor(3,iAtom)**2 -
     &                Coor(iCar,iAtom)**2)
 43      Continue
 40   Continue
      Call mma_deallocate(Coor)
      If (RR_Show) Then
      Write (6,'(19X,A)') ' The Moment of Inertia Tensor / au'
      Write (6,'(19X,14X,3A)')'    X     ','     Y    ','    Z     '
      Write (6,'(19X,A,12X,3(E11.4))') ' X',rMI(1)
      Write (6,'(19X,A,12X,3(E11.4))') ' Y',rMI(2), rMI(3)
      Write (6,'(19X,A,12X,3(E11.4))') ' Z',rMI(4), rMI(5), rMI(6)
      Write (6,*)
      End If
*
*     Diagonalize and find principle axis
*
      Call mma_Allocate(Hess,6)
      call dcopy_(6,rMI,1,Hess,1)
      Call Jacob(Hess,Pax,3,3)
      Prin(1) = Hess(1)
      Prin(2) = Hess(3)
      Prin(3) = Hess(6)
*
*     Sort the prinipal axis such that z' is the one with the lowest
*     eigenvalue.
*
      Do 100 i = 1 , 2
         Do 101 j = i+1, 3
            If (Prin(i).lt.Prin(j)) Then
               Save=Prin(i)
               Prin(i)=Prin(j)
               Prin(j)=Save
               Call DSwap_(3,PAx(1+(i-1)*3),1,PAx(1+(j-1)*3),1)
            End If
 101     Continue
 100  Continue
      Call mma_deallocate(Hess)
      If (RR_Show) Then
      Write (6,'(19X,A)')
     &         ' The Principal Axes and Moments of Inertia (au)'
      Write (6,'(19X,A,3(E11.4))') ' Eigenvalues :', (Prin(i),i=1,3)
      Write (6,'(19X,14X,3A)')'    X''    ','     Y''   ','    Z''    '
      Write (6,'(19X,A)')        ' Eigenvectors:'
      Write (6,'(19X,A,3(E11.4))')
     &      ' X            ',Pax(1),Pax(4),Pax(7)
      Write (6,'(19X,A,3(E11.4))')
     &      ' Y            ',Pax(2),Pax(5),Pax(8)
      Write (6,'(19X,A,3(E11.4))')
     &      ' Z            ',Pax(3),Pax(6),Pax(9)
      Write (6,*)
C     Call Put_dArray('PAX',Pax,9)
      Write (6,'(19X,A)') ' The Rotational Constants'
      Write (6,'(19X,A)') '         (cm-1)            (GHz)'
      If (Prin(1).ge.1.D-3)
     &    Write (6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(1),
     &                          1.0D-9*auToHz*Half/Prin(1)
      If (Prin(2).ge.1.D-3)
     &    Write (6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(2),
     &                          1.0D-9*auToHz*Half/Prin(2)
      If (Prin(3).ge.1D-3)
     &    Write (6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(3),
     &                          1.0D-9*auToHz*Half/Prin(3)
      End If
      If (Prin(1).eq.Zero .and.
     &    Prin(2).eq.Zero .and.
     &    Prin(3).eq.Zero ) Go To 99
      If (RR_Show) Then
      Write (6,*)
      Write (6,*)
      Write (6,'(19X,A)') ' *******************************************'
      Write (6,'(19X,A)') ' *                                         *'
      Write (6,'(19X,A)') ' * R I G I D - R O T O R   A N A L Y S I S *'
      Write (6,'(19X,A)') ' *                                         *'
      Write (6,'(19X,A)') ' *******************************************'
      Write (6,*)
      Write (6,'(19X,A,I3)') ' j(Max):', S%jMax
      Write (6,*)
      End If
*
*     Order the three principal moments of inertia, Ia<=Ib<=Ic
*
      call dcopy_(3,Prin,1,XI,1)
      Call Order_Axis(XI,3)
*
*     Asymmetry parameter, Tinkham formula 7-96
*
      B = Half / XI(2)
      C = Half / XI(3)
      If (XI(1).le.1.D-3) Then
         Linear = .True.
         rKappa = -One
      Else
         A = Half / XI(1)
         If (Abs(A-C).le.1.D-10) Then
            rKappa = 0.0D00
         Else
            rKappa = (Two*B-A-C) / (A-C)
         End If
      End If
*
*     Change order if molecule is oblate, Ia=Ib<=Ic
*
      If (Abs(XI(1)-XI(2)).lt.0.000001D0 .and.
     &    Abs(XI(1)-XI(3)).gt.0.000001D0) Then
         Save = XI(1)
         XI(1) = XI(3)
         XI(3) = Save
      End If
*
*     Construct Hessian Matrix
*     Iz=XI(1), Ix=XI(2), and Iy=XI(3)
*
*     Set up constants, see Tinkham fomula 7-89a and 7-89b.
*
      If (Abs(XI(1)).gt.1.0D-3) Then
         A = Half / XI(1)
      Else
         A = Zero
      End If
      B = (One/XI(2) + One/XI(3))/Four
      C =(One/XI(2) - One/XI(3))/Eight
      Do 60 i = 1, 80
         Label(i:i) = ' '
 60   Continue
*     Iz=0 linear rotor
      If (Abs(A).lt.1.0d-3) Then
         Label(1:13)=' Linear Rotor'
*     Ix=/=Iy asymmetric top
      Else If (C .gt. 1.0D-3) Then
         Label(1:15)=' Asymmetric Top'
*     Iz=Ix=Iy
      Else If (Abs(A-B).lt.1.0D-10) Then
         Label(1:14)=' Spherical Top'
      Else If(A.gt.B) Then
         Label(1:25)=' Symmetric Top, Prolate'
      Else
         Label(1:24)=' Symmetric Top, Oblate'
      End If
      If (RR_Show) Then
      Write (6,'(19X,A,A)') ' Rotor Type:',Label(1:mylen(Label))
      Write (6,'(19X,A,F7.3)') ' Asymmetry parameter:',rKappa
      Write (6,'(19X,A)')       ' Prolate = -1'
      Write (6,'(19X,A)')       ' Oblate  =  1'
      Write (6,*)
      End If
*
      nEn = (S%jMax+1)*(S%jMax+2)*(S%jMax+3)/6
      Call mma_Allocate(En,nEn)
      iEn = 1
      nHess = (2*S%jMax+1)*(2*S%jMax+2)/2
      Call mma_allocate(Hess,nHess)
      Call mma_Allocate(Vec,2*S%jMax+1,2*S%jMax+1)
      Do 80 j = 0, S%jMax
         mDim = 2*j+1
         nTri = mDim*(mDim+1)/2
         call dcopy_(nTri,[Zero],0,Hess,1)
         call dcopy_(mDim**2,[Zero],0,Vec,1)
         call dcopy_(mDim,[One],0,Vec,mDim+1)
         If (iPrint.ge.99)
     &       Call RecPrt(' Vec',' ',Vec,mDim,mDim)
         k1 = 1
         Do 81 k = -j, j
            kk = k1*(k1+1)/2
*
*           Formula 7-89a, (j,k|H|j,k), diagonal term
*
            Hess(kk) = B*DBLE(j*(j+1)) + (A-B)*DBLE(k**2)
            If (Linear) Hess(kk)=B*DBLE(j*(j+1))
            If (k+2.gt.j .or. Linear) Go to 82
            k2 = k1 + 2
            kk2 = k2*(k2-1)/2 + k1
*
*           Formula 7-89b, (j,k+2|H|j,k), the only off diagonal term.
*
            Hess(kk2) = C*Sqrt(DBLE( (j*(j+1)-k*(k+1)) *
     &                  (j*(j+1)-(k+1)*(k+2)) ) )
 82         k1 = k1 + 1
 81      Continue
         If (iPrint.ge.99)
     &       Call TriPrt(' Hessian',' ',Hess,mDim)
         Call Jacob(Hess,Vec,mDim,mDim)
         If (iPrint.ge.99)
     &      Call TriPrt(' Hessian',' ',Hess,mDim)
         Do 83 i = 1, mDim
            En(iEn+I-1) = Hess(i*(i+1)/2) * auTocm
 83      Continue
         Call Order_Axis(En(iEn),mDim)
         iEn = iEn + mDim
 80   Continue
      Call mma_deallocate(Vec)
      Call mma_deallocate(Hess)
*
*     Output
*
      If (RR_Show) Then
      Write (6,*)
      Write (6,'(19X,A)') ' Rotational energies / cm-1'
      iEn = 1
      Do 90 j = 0, S%jMax
         Write (6,*)
         If (Linear) Then
             Write (6,'(19X,A,I2,A,F8.3)')
     &           ' E(J=',J,') = ',En(iEn)
             iEn = iEn + (2*j+1)
         Else
            Do 91 kappa = -j, j
               Write (6,'(19X,A,I2,A,I2,A,F12.3)')
     &           ' E(J=',J,',kappa=',kappa,') = ',En(iEn)
               iEn = iEn + 1
 91         Continue
         End If
 90   Continue
      End If
      Call mma_deallocate(En)
*
 99   Continue
      If (RR_Show) Then
         Call CollapseOutput(0,'   Rigid rotor info:')
         Write (6,*)
      End If
      Call qExit('RigRot')
      Return
      End

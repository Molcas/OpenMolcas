!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Setup1(ExpA,nPrim,ExpB,mPrim,A,B,rKappa,Pcoor,ZInv)
!***********************************************************************
!                                                                      *
!     Object : to compute some data which is needed for the one-       *
!              electron integrals.                                     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!***********************************************************************
      use Constants, only: Zero, One
      Implicit None
      Integer nPrim, mPrim
      Real*8, Intent(In):: ExpA(nPrim), ExpB(mPrim), ZInv(nPrim,mPrim),
     &                     A(3), B(3)
      Real*8, Intent(Out):: rKappa(nPrim,mPrim), Pcoor(nPrim,mPrim,3)

      Real*8 ab
      Integer iPrim, jPrim
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' *** ExpA ***',' ',ExpA,1,nPrim)
      Call RecPrt(' *** ExpB ***',' ',ExpB,1,mPrim)
      Call RecPrt(' *** ZInv ***',' ',ZInv,nPrim,mPrim)
      Write (6,*) 'A(:)=',A(:)
      Write (6,*) 'B(:)=',B(:)
#endif
      ab  = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2

      If (ab.ne.Zero) Then
      Do jPrim = 1, mPrim
         Do iPrim = 1, nPrim
            rKappa(iPrim,jPrim) = Exp(- ExpA(iPrim) * ExpB(jPrim) * ab *
     &                            ZInv(iPrim,jPrim))
            Pcoor(iPrim,jPrim,1)=(ExpA(iPrim)*A(1)+ExpB(jPrim)*B(1)) *
     &                            ZInv(iPrim,jPrim)
            Pcoor(iPrim,jPrim,2)=(ExpA(iPrim)*A(2)+ExpB(jPrim)*B(2)) *
     &                            ZInv(iPrim,jPrim)
            Pcoor(iPrim,jPrim,3)=(ExpA(iPrim)*A(3)+ExpB(jPrim)*B(3)) *
     &                            ZInv(iPrim,jPrim)
         End Do
      End Do
      Else
        rKappa(:,:)=One
        PCoor(:,:,1)=A(1)
        PCoor(:,:,2)=A(2)
        PCoor(:,:,3)=A(3)
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt(' *** Kappa ***',' ',rKappa, nPrim, mPrim)
      Call RecPrt(' ***   Px  ***',' ',Pcoor(1,1,1),nPrim,mPrim)
      Call RecPrt(' ***   Py  ***',' ',Pcoor(1,1,2),nPrim,mPrim)
      Call RecPrt(' ***   Pz  ***',' ',Pcoor(1,1,3),nPrim,mPrim)
#endif
!
      End SubRoutine Setup1

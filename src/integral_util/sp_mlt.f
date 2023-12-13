!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Sp_Mlt(W_In,ne,W_out,nVec,C,nab)
      use Constants, only: Zero
      Implicit None
      Integer ne, nVec, nab
      Real*8 W_In(ne,nVec), W_Out(nVec,nab), C(ne,nab)

      Integer, Parameter:: meMax=10
      Integer iAux(meMax+1)
      Integer iab, me, iVec, ie
!
      Do iab = 1, nab
         me=0
         Do ie = 1, ne
            If (C(ie,iab).ne.Zero) Then
               me=me+1
               iAux(me)=ie
               If (me.gt.meMax) Exit
            End If
         End Do
!
         Select Case(me)

         Case (1)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
            End Do
         Case (2)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
            End Do
         Case (3)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
            End Do
         Case (4)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
            End Do
         Case (5)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
            End Do
         Case (6)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
            End Do
         Case (7)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
     &                        +C(iAux(7),iab)*W_In(iAux(7),iVec)
            End Do
         Case (8)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
     &                        +C(iAux(7),iab)*W_In(iAux(7),iVec)
     &                        +C(iAux(8),iab)*W_In(iAux(8),iVec)
            End Do
         Case (9)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
     &                        +C(iAux(7),iab)*W_In(iAux(7),iVec)
     &                        +C(iAux(8),iab)*W_In(iAux(8),iVec)
     &                        +C(iAux(9),iab)*W_In(iAux(9),iVec)
            End Do
         Case (10)
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
     &                        +C(iAux(7),iab)*W_In(iAux(7),iVec)
     &                        +C(iAux(8),iab)*W_In(iAux(8),iVec)
     &                        +C(iAux(9),iab)*W_In(iAux(9),iVec)
     &                        +C(iAux(10),iab)*W_In(iAux(10),iVec)
            End Do
         Case Default
            Do iVec = 1, nVec
               W_Out(iVec,iab)=C(iAux(1),iab)*W_In(iAux(1),iVec)
     &                        +C(iAux(2),iab)*W_In(iAux(2),iVec)
     &                        +C(iAux(3),iab)*W_In(iAux(3),iVec)
     &                        +C(iAux(4),iab)*W_In(iAux(4),iVec)
     &                        +C(iAux(5),iab)*W_In(iAux(5),iVec)
     &                        +C(iAux(6),iab)*W_In(iAux(6),iVec)
     &                        +C(iAux(7),iab)*W_In(iAux(7),iVec)
     &                        +C(iAux(8),iab)*W_In(iAux(8),iVec)
     &                        +C(iAux(9),iab)*W_In(iAux(9),iVec)
     &                        +C(iAux(10),iab)*W_In(iAux(10),iVec)
            End Do
            Do ie = iAux(meMax+1), ne
               If (C(ie,iab).ne.Zero)
     &            Call DaXpY_(nVec,C(ie,iab),W_In(ie,1),ne,
     &                                      W_Out(1,iab),1)
            End Do
         End Select
      End Do
!
      Return
      End Subroutine Sp_Mlt

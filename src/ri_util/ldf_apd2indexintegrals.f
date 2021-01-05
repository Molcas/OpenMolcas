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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_APD2IndexIntegrals_11(A,B,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: get integrals (J|K) using atom pair-driven code.
C
C     NOTE: this code is for debugging only and requires that for each
C           atom C,D at least one atom pair can be found for which
C           (a) it contains atom A (C)
C           (b) it has no linear dependence in the linear dependence in
C               the one-center aux space.
C
      Implicit None
      Integer A
      Integer B
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

      Character*25 SecNam
      Parameter (SecNam='LDF_APD2IndexIntegrals_11')

      Integer  LDF_nBasAux_Atom, LDF_nBasAux_Pair
      External LDF_nBasAux_Atom, LDF_nBasAux_Pair

      Integer AC, AC_, iAC, BD, BD_, iBD
      Integer l, nAC, nBD
      Integer ip_Int, l_Int
      Integer nA, iRow0
      Integer nB, iCol0
      Integer K, JK, ip, ip0

      Integer i, j
      Integer A2AP
      Integer AP_Atoms
      Integer AP_1CLinDep
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)

      ! Get an atom pair that contains atom A and another containing C
      Call LDF_SetA2AP()
      AC=0
      l=A2AP(1,A)
      If (l.gt.0) Then
         ip=A2AP(2,A)-1
         iAC=0
         Do While (iAC.lt.l .and. AC.eq.0)
            iAC=iAC+1
            AC_=iWork(ip+iAC)
            If (AP_1CLinDep(1,AC_).eq.0) Then
               AC=AC_
            End If
         End Do
      End If
      If (AC.eq.0) Then
         Write(6,'(A,A,I5,A)')
     &   SecNam,': no atom pair containing atom ',A,
     &   ' has linearly independent one-center functions only!'
         Call WarningMessage(2,SecNam//': Suitable AP not found!')
         Call LDF_Quit(1)
      End If
      BD=0
      l=A2AP(1,B)
      If (l.gt.0) Then
         ip=A2AP(2,B)-1
         iBD=0
         Do While (iBD.lt.l .and. BD.eq.0)
            iBD=iBD+1
            BD_=iWork(ip+iBD)
            If (AP_1CLinDep(1,BD_).eq.0) Then
               BD=BD_
            End If
         End Do
      End If
      If (BD.eq.0) Then
         Write(6,'(A,A,I5,A)')
     &   SecNam,': no atom pair containing atom ',B,
     &   ' has linearly independent one-center functions only!'
         Call WarningMessage(2,SecNam//': Suitable AP not found!')
         Call LDF_Quit(1)
      End If
      Call LDF_UnsetA2AP()

      ! Allocate integral array
      nAC=LDF_nBasAux_Pair(AC)
      nBD=LDF_nBasAux_Pair(BD)
      l_Int=nAC*nBD
      Call GetMem('APD2I_11','Allo','Real',ip_Int,l_Int)

      ! Compute integrals
      Call LDF_ComputeIntegrals_JK_2P(AC,BD,l_Int,Work(ip_Int))

      ! Copy integrals to return
      nA=LDF_nBasAux_Atom(A)
      If (AP_Atoms(1,AC).eq.A) Then
         iRow0=0
      Else If (AP_Atoms(2,AC).eq.A) Then
         iRow0=LDF_nBasAux_Atom(AP_Atoms(1,AC))
      Else
         Call WarningMessage(2,SecNam//': Logical error! [1]')
         Call LDF_Quit(1)
         iRow0=0 ! avoid compiler warning
      End If
      nB=LDF_nBasAux_Atom(B)
      If (AP_Atoms(1,BD).eq.B) Then
         iCol0=0
      Else If (AP_Atoms(2,BD).eq.B) Then
         iCol0=LDF_nBasAux_Atom(AP_Atoms(1,BD))
      Else
         Call WarningMessage(2,SecNam//': Logical error! [2]')
         Call LDF_Quit(1)
         iCol0=0 ! avoid compiler warning
      End If
      l=nA*nB
      If (l.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      Else
         JK=0
         ip0=ip_Int-1
         Do K=iCol0+1,iCol0+nB
            ip=ip0+nAC*(K-1)
            Do J=iRow0+1,iRow0+nA
               JK=JK+1
               xInt(JK)=Work(ip+J)
            End Do
         End Do
      End If

      ! Deallocate integral array
      Call GetMem('APD2I_11','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_APD2IndexIntegrals_12(A,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: get integrals (J|K) using atom pair-driven code.
C
C     NOTE: this code is for debugging only and requires that for
C           atom A at least one atom pair can be found for which
C           (a) it contains atom A
C           (b) it has no linear dependence in the linear dependence in
C               the one-center aux space.
C
      Implicit None
      Integer A
      Integer CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

      Character*25 SecNam
      Parameter (SecNam='LDF_APD2IndexIntegrals_12')

      Integer  LDF_nBasAux_Atom, LDF_nBasAux_Pair
      External LDF_nBasAux_Atom, LDF_nBasAux_Pair

      Integer AB, AB_, iAB, C, D
      Integer ip0, ip, l
      Integer nAB, nCD
      Integer ip_Int, l_Int
      Integer nA, K, JK
      Integer iRow0, iCol0

      Integer i, j
      Integer A2AP
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Quick return
      If (AP_2CFunctions(1,CD).lt.1) Return

      ! Get an atom pair that contains atom A and another containing C
      Call LDF_SetA2AP()
      AB=0
      l=A2AP(1,A)
      If (l.gt.0) Then
         ip=A2AP(2,A)-1
         iAB=0
         Do While (iAB.lt.l .and. AB.eq.0)
            iAB=iAB+1
            AB_=iWork(ip+iAB)
            If (AP_1CLinDep(1,AB_).eq.0) Then
               AB=AB_
            End If
         End Do
      End If
      If (AB.eq.0) Then
         Write(6,'(A,A,I5,A)')
     &   SecNam,': no atom pair containing atom ',A,
     &   ' has linearly independent one-center functions only!'
         Call WarningMessage(2,SecNam//': Suitable AP not found!')
         Call LDF_Quit(1)
      End If
      Call LDF_UnsetA2AP()

      ! Allocate integral array
      nAB=LDF_nBasAux_Pair(AB)
      nCD=LDF_nBasAux_Pair(CD)
      l_Int=nAB*nCD
      Call GetMem('APD2I_12','Allo','Real',ip_Int,l_Int)

      ! Compute integrals
      Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_Int,Work(ip_Int))

      ! Copy integrals to return
      nA=LDF_nBasAux_Atom(A)
      l=nA*AP_2CFunctions(1,CD)
      If (l.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      Else
         If (AP_Atoms(1,AB).eq.A) Then
            iRow0=0
         Else If (AP_Atoms(2,AB).eq.A) Then
            iRow0=LDF_nBasAux_Atom(AP_Atoms(1,AB))
         Else
            Call WarningMessage(2,SecNam//': Logical error!')
            Call LDF_Quit(1)
            iRow0=0 ! avoid compiler warning
         End If
         C=AP_Atoms(1,CD)
         D=AP_Atoms(2,CD)
         If (C.eq.D) Then
            iCol0=LDF_nBasAux_Atom(C)-AP_1CLinDep(1,CD)
         Else
            iCol0=LDF_nBasAux_Atom(C)+LDF_nBasAux_Atom(D)
     &           -AP_1CLinDep(1,CD)
         End If
         JK=0
         ip0=ip_Int-1
         Do K=iCol0+1,iCol0+AP_2CFunctions(1,CD)
            ip=ip0+nAB*(K-1)
            Do J=iRow0+1,iRow0+nA
               JK=JK+1
               xInt(JK)=Work(ip+J)
            End Do
         End Do
      End If

      ! Deallocate integral array
      Call GetMem('APD2I_12','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_APD2IndexIntegrals_22(AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: get integrals (J|K) using atom pair-driven code.
C
C     NOTE: this code is for debugging only.
C
      Implicit None
      Integer AB
      Integer CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*25 SecNam
      Parameter (SecNam='LDF_APD2IndexIntegrals_22')

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Integer nAB, nCD
      Integer ip_Int, l_Int
      Integer l, ip, ip0
      Integer JK, K
      Integer iRow0, iCol0

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Quick return
      If (AP_2CFunctions(1,AB).lt.1 .or. AP_2CFunctions(1,CD).lt.1) Then
         Return
      End If

      ! Allocate integral array
      nAB=LDF_nBasAux_Pair(AB)
      nCD=LDF_nBasAux_Pair(CD)
      l_Int=nAB*nCD
      Call GetMem('APD2I_22','Allo','Real',ip_Int,l_Int)

      ! Compute integrals
      Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_Int,Work(ip_Int))

      ! Copy integrals to return
      l=AP_2CFunctions(1,AB)*AP_2CFunctions(1,CD)
      If (l.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      Else
         iRow0=LDF_nBasAux_Pair(AB)-AP_2CFunctions(1,AB)
         iCol0=LDF_nBasAux_Pair(CD)-AP_2CFunctions(1,CD)
         JK=0
         ip0=ip_Int-1
         Do K=iCol0+1,iCol0+AP_2CFunctions(1,CD)
            ip=ip0+nAB*(K-1)
            Do J=iRow0+1,iRow0+AP_2CFunctions(1,AB)
               JK=JK+1
               xInt(JK)=Work(ip+J)
            End Do
         End Do
      End If

      ! Deallocate integral array
      Call GetMem('APD2I_22','Free','Real',ip_Int,l_Int)

      End

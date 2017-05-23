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
      Subroutine LDF_APD3IndexIntegrals_1(AB,C,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: get integrals (uv|J) using atom pair-driven code
C
C     NOTE: this code is for debugging only and requires that at least
C           one atom pair can be found for which
C           (a) it contains atom C
C           (b) it has no linear dependence in the one-center aux space
C
      Implicit None
      Integer AB
      Integer C
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

      Character*24 SecNam
      Parameter (SecNam='LDF_APD3IndexIntegrals_1')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair

      Integer iCD, CD_, CD
      Integer ip, l
      Integer A, B
      Integer nuv, M
      Integer ip_Int, l_Int

      Integer i, j
      Integer A2AP
      Integer AP_Atoms
      Integer AP_1CLinDep
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)

      ! Get an atom pair that contains atom C
      Call LDF_SetA2AP()
      CD=0
      l=A2AP(1,C)
      If (l.gt.0) Then
         ip=A2AP(2,C)-1
         iCD=0
         Do While (iCD.lt.l .and. CD.eq.0)
            iCD=iCD+1
            CD_=iWork(ip+iCD)
            If (AP_1CLinDep(1,CD_).eq.0) Then
               CD=CD_
            End If
         End Do
      End If
      If (CD.eq.0) Then
         Write(6,'(A,A,I5,A)')
     &   SecNam,': no atom pair containing atom ',C,
     &   ' has linearly independent one-center functions only!'
         Call WarningMessage(2,SecNam//': Suitable AP not found!')
         Call LDF_Quit(1)
      End If
      Call LDF_UnsetA2AP()

      ! Allocate integral array
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nuv=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      M=LDF_nBasAux_Pair(CD)
      l_Int=nuv*M
      Call GetMem('APD3I_1','Allo','Real',ip_Int,l_Int)

      ! Compute integrals using atom pair-driven code
      Call LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_Int,Work(ip_Int))

      ! Copy integrals to return
      l=nuv*LDF_nBasAux_Atom(C)
      If (l.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      Else
         If (AP_Atoms(1,CD).eq.C) Then
            Call dCopy_(l,Work(ip_Int),1,xInt,1)
         Else If (AP_Atoms(2,CD).eq.C) Then
            ip=ip_Int+nuv*LDF_nBasAux_Atom(AP_Atoms(1,CD))
            Call dCopy_(l,Work(ip),1,xInt,1)
         Else
            Call WarningMessage(2,SecNam//': Logical error!')
            Call LDF_Quit(1)
         End If
      End If

      ! Deallocate integral array
      Call GetMem('APD3I_1','Free','Real',ip_Int,l_Int)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_APD3IndexIntegrals_2(AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: get integrals (uv|J) using atom pair-driven code
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

      Character*24 SecNam
      Parameter (SecNam='LDF_APD3IndexIntegrals_2')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair, LDF_nBasAux_Atom
      External LDF_nBas_Atom, LDF_nBasAux_Pair, LDF_nBasAux_Atom

      Integer A, B, C, D
      Integer nuv, M
      Integer ip_Int, l_Int
      Integer ip, l

      Integer i, j
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Quick return
      If (AP_2CFunctions(1,CD).lt.1) Return

      ! Allocate memory for integrals
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nuv=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      M=LDF_nBasAux_Pair(CD)
      l_Int=nuv*M
      Call GetMem('APD3I_2','Allo','Real',ip_Int,l_Int)

      ! Compute integrals
      Call LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_Int,Work(ip_Int))

      ! Copy integrals to return
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)
      If (C.eq.D) Then
         ip=ip_Int+nuv*(LDF_nBasAux_Atom(C)-AP_1CLinDep(1,CD))
      Else
         ip=ip_Int+nuv*(LDF_nBasAux_Atom(C)+LDF_nBasAux_Atom(D)
     &                 -AP_1CLinDep(1,CD))
      End If
      l=nuv*AP_2CFunctions(1,CD)
      If (l.gt.l_xInt_) Then
         Call WarningMessage(2,
     &               SecNam//': Insufficient integral array dimension!')
         Call LDF_Quit(1)
      Else
         Call dCopy_(l,Work(ip),1,xInt,1)
      End If

      ! Deallocate integrals
      Call GetMem('APD3I_2','Free','Real',ip_Int,l_Int)

      End

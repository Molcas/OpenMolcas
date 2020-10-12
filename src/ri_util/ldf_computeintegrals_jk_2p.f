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
      Subroutine LDF_ComputeIntegrals_JK_2P(AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: compute integrals (J_AB|K_CD).
C
      Implicit None
      Integer AB, CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int2.fh"
#include "ldf_atom_pair_info.fh"

      Character*26 SecNam
      Parameter (SecNam='LDF_ComputeIntegrals_JK_2P')

#if defined (_DEBUGPRINT_)
      Logical  isSymmetric
      External isSymmetric
#endif

      Integer  LDF_nBasAux_Pair
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Pair
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Integer MAB, MCD
      Integer l_xInt
      Integer C, D
      Integer dShell
      Integer nAuxShell_C, ipC
      Integer nAuxShell_D, ipD
      Integer lS, kShell, lShell
      Integer nnShl, klS
      Integer ip_SewWrk, l_SewWrk

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Get row dimension of integrals
      MAB=LDF_nBasAux_Pair(AB)

      If (AB.eq.CD) Then
         ! Check integral dimension
         l_xInt=MAB**2
         If (l_xInt.ne.l_xInt_) Then
            Call WarningMessage(2,
     &                      SecNam//': integral dimension problem! [0]')
            Call LDF_Quit(1)
         End If
         ! Use same code as fitting coefficient calculation
         ! (where (J|K) was called the G matrix)
         Call LDF_ComputeGMat(AB,MAB,xInt)
      Else
         ! Get column dimension of integrals
         MCD=LDF_nBasAux_Pair(CD)
         ! Check integral dimension
         l_xInt=MAB*MCD
         If (l_xInt.gt.l_xInt_) Then
            Call WarningMessage(2,
     &                      SecNam//': integral dimension problem! [1]')
            Call LDF_Quit(1)
         End If
         ! Init integral array
         Call Cho_dZero(xInt,l_xInt)
         ! Set indices
         Call LDF_SetIndx_JK_2P(AB,CD)
         ! Allocate Seward memory
         Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
         Call xSetMem_Ints(l_SewWrk)
         ! Compute integrals
         dShell=nShell_Valence+nShell_Auxiliary+1
         C=AP_Atoms(1,CD)
         D=AP_Atoms(2,CD)
         nAuxShell_C=LDF_nAuxShell_Atom(C)
         ipC=LDF_lAuxShell_Atom(C)-1
         kShell=dShell
         Do lS=1,nAuxShell_C
            lShell=iWork(ipC+lS)
            Call LDF_CJK2P_Col(AB,kShell,lShell,l_xInt,xInt)
         End Do
         If (D.ne.C) Then
            nAuxShell_D=LDF_nAuxShell_Atom(D)
            ipD=LDF_lAuxShell_Atom(D)-1
            Do lS=1,nAuxShell_D
               lShell=iWork(ipD+lS)
               Call LDF_CJK2P_Col(AB,kShell,lShell,l_xInt,xInt)
            End Do
         End If
         If (AP_2CFunctions(1,CD).gt.0) Then
            nnShl=l_CD_2CList_2
            Do klS=0,nnShl-1
               kShell=iWork(ip_CD_2CList+3*klS)
               lShell=iWork(ip_CD_2CList+3*klS+1)
               SPCD=iWork(ip_CD_2CList+3*klS+2)
               Call LDF_CJK2P_Col(AB,kShell,lShell,l_xInt,xInt)
            End Do
         End If
         ! Deallocate Seward memory
         Call xRlsMem_Ints()
         ! Unset indices
         Call LDF_UnsetIndx_JK_2P()
      End If

#if defined (_DEBUGPRINT_)
      ! Check symmetry
      If (AB.eq.CD) Then
         If (.not.isSymmetric(xInt,MAB,1.0d-14)) Then
            Call WarningMessage(2,SecNam//': (J|K) != (K|J)')
            Write(6,'(A,I9,1X,I9)') 'AB,CD=',AB,CD
            Call LDF_Quit(1)
         End If
      End If
#endif

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CJK2P_Col(AB,kShell,lShell,l_xInt,xInt)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Compute integrals (J_AB | kShell lShell)
C
      Implicit None
      Integer AB
      Integer kShell, lShell
      Integer l_xInt
      Real*8  xInt(l_xInt)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int2.fh"
#include "ldf_atom_pair_info.fh"

      External Int_LDF_JK_2P

      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Integer dShell
      Integer A, B
      Integer nAuxShell_A, nAuxShell_B
      Integer ipA, ipB
      Integer iShell, jShell
      Integer nnShl
      Integer jS, ijS

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      dShell=nShell_Valence+nShell_Auxiliary+1

      SHC=kShell
      SHD=lShell

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nAuxShell_A=LDF_nAuxShell_Atom(A)
      nAuxShell_B=LDF_nAuxShell_Atom(B)
      ipA=LDF_lAuxShell_Atom(A)-1
      ipB=LDF_lAuxShell_Atom(B)-1

      iShell=dShell
      SHA=iShell
      Do jS=1,nAuxShell_A
         jShell=iWork(ipA+jS)
         SHB=jShell
         Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                  xInt,l_xInt,Int_LDF_JK_2P)
      End Do
      If (B.ne.A) Then
         Do jS=1,nAuxShell_B
            jShell=iWork(ipB+jS)
            SHB=jShell
            Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                     xInt,l_xInt,Int_LDF_JK_2P)
         End Do
      End If
      If (AP_2CFunctions(1,AB).gt.0) Then
         nnShl=l_AB_2CList_2
         Do ijS=0,nnShl-1
            iShell=iWork(ip_AB_2CList+3*ijS)
            jShell=iWork(ip_AB_2CList+3*ijS+1)
            SPAB=iWork(ip_AB_2CList+3*ijS+2)
            SHA=iShell
            SHB=jShell
            Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                     xInt,l_xInt,Int_LDF_JK_2P)
         End Do
      End If

      End

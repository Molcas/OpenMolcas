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
* Copyright (C) 2010,2012, Thomas Bondo Pedersen                       *
************************************************************************
      Subroutine LDF_ComputeApproximateIntegrals(Mode,Add,AB,CD,
     &                                           l_xInt_,xInt)
C
C=======================================================================
C          WARNING: THIS ROUTINE MAY REQUIRE A LOT OF MEMORY !!!
C=======================================================================
C
C     Thomas Bondo Pedersen, October 2010, January 2012
C
C     Purpose: compute integrals using local DF coefficients.
C
C     Input:
C       Mode   --- 1: robust, 2: nonrobust, 3: half-and-half
C       Add    --- .True. if result should be added to xInt
C       AB, CD --- atom pairs for which the integrals are to be computed
C
C     LDF info must be properly set up before calling this routine.
C
C     Robust representation:
C
C     (u_A v_B | k_C l_D) = sum_K_CD (u_A v_A | K_CD) * C[k_C l_D,K_CD]
C                         + sum_J_AB C[u_A v_B,J_AB] *
C                           { (k_C l_D | J_AB)
C                           - sum_K_CD (J_AB | K_CD) * C[k_C l_D,K_CD] }
C
C     Non-robust representation:
C
C     (u_A v_B | k_C l_D) = sum_J_AB C[u_A v_B,J_AB] *
C                           sum_K_CD (J_AB | K_CD) * C[k_C l_D,K_CD]
C
C     Half-and-half representation:
C
C     (u_A v_B | k_C l_D) ={sum_K_CD (u_A v_A | K_CD) * C[k_C l_D,K_CD]
C                          +sum_J_AB C[u_A v_B,J_AB] * (k_C l_D | J_AB)}
C                         *1/2
C
C     Although the three representations are equivalent when AB=CD, the
C     robust and half-and-half representation is always used if
C     requested.
C
      Implicit None
      Integer Mode
      Logical Add
      Integer AB
      Integer CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#if defined (_DEBUG_)
#include "localdf_bas.fh"
#endif

      Character*31 SecNam
      Parameter (SecNam='LDF_ComputeApproximateIntegrals')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair
#if defined (_DEBUG_)
      Logical  isSymmetric
      External isSymmetric
      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom
#endif

      Integer A, B, C, D
      Integer nAB, MAB
      Integer nCD, MCD
      Integer l_xInt
      Integer ip_3Int, l_3Int
      Integer ip_2Int, l_2Int
      Integer ip_C, l_C, l_C_CD, l_C_AB
      Integer ip, l
#if defined (_DEBUG_)
      Integer ip_iAB, l_iAB
      Integer ip_iCD, l_iCD
      Integer ipA, ipB, ipC, ipD
      Integer k
      Integer ij, ji, kl, lk
      Integer iS, jS, kS, lS
      Integer iShell, jShell, kShell, lShell
      Integer nShell_A, nShell_B, nShell_C, nShell_D
      Integer n
#endif

      Integer i, j
      Integer AP_Atoms
#if defined (_DEBUG_)
      Integer nBasSh, iAB, iCD
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iAB(i,j)=iWork(ip_iAB-1+nShell_A*(j-1)+i)
      iCD(i,j)=iWork(ip_iCD-1+nShell_C*(j-1)+i)
#endif
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Get atoms of pairs
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)

      ! Get dimensions
      nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      MAB=LDF_nBasAux_Pair(AB)
      nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      MCD=LDF_nBasAux_Pair(CD)

      ! Return if nothing to do
      If (nAB.lt.1 .or. nCD.lt.1 .or. MAB.lt.1 .or. MCD.lt.1) Return

      ! Check dimension of integral array
      l_xInt=nAB*nCD
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! If integrals are not to be updated, initialize
      If (.not.Add) Call Cho_dZero(xInt,l_xInt)

      ! Separate codes for different integral representations
      If (Mode.eq.1 .or. Mode.eq.3) Then
         ! Allocate 3-index integrals
         l_3Int=max(nAB*MCD,nCD*MAB)
         Call GetMem('CAI3Int','Allo','Real',ip_3Int,l_3Int)
         ! Compute (u_A v_B | K_CD)
         Call LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_3Int,Work(ip_3Int))
         ! Allocate coefficients
         l_C_AB=nAB*MAB
         l_C_CD=nCD*MCD
         l_C=max(l_C_AB,l_C_CD)
         Call GetMem('CAIC','Allo','Real',ip_C,l_C)
         ! Read coefficients for CD
         Call LDF_CIO_ReadC(CD,Work(ip_C),l_C_CD)
         ! Compute integral contribution
         ! (u_A v_B | k_C l_D) +=
         ! sum_K_CD (u_A v_B | K_CD)*C[k_C l_D,K_CD]
         Call dGeMM_('N','T',nAB,nCD,MCD,
     &               1.0d0,Work(ip_3Int),nAB,Work(ip_C),nCD,
     &               1.0d0,xInt,nAB)
         ! Compute (k_C l_D | J_AB), if needed
         If (CD.ne.AB) Then
            Call LDF_ComputeIntegrals_uvJ_2P(CD,AB,l_3Int,Work(ip_3Int))
         End If
         If (Mode.eq.1) Then ! robust, 2-index integrals needed
            ! Allocate (J_AB | K_CD)
            l_2Int=MAB*MCD
            Call GetMem('CAI2Int','Allo','Real',ip_2Int,l_2Int)
            ! Compute (J_AB | K_CD)
            Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_2Int,Work(ip_2Int))
            ! Update:
            ! (k_C l_D | J_AB) -= sum_K_CD C[k_C l_D,K_CD]*(J_AB | K_CD)
            Call dGeMM_('N','T',nCD,MAB,MCD,
     &                  -1.0d0,Work(ip_C),nCD,Work(ip_2Int),MAB,
     &                  1.0d0,Work(ip_3Int),nCD)
            ! Deallocate (J_AB | K_CD)
            Call GetMem('CAI2Int','Allo','Real',ip_2Int,l_2Int)
         End If
         ! Read coefficients for AB, if different from CD
         If (CD.ne.AB) Then
            Call LDF_CIO_ReadC(AB,Work(ip_C),l_C_AB)
         End If
         ! Compute integral contribution
         ! (u_A v_B | k_C l_D) +=
         ! sum_J_AB C[u_A v_B,J_AB]*(k_C l_D | J_AB)
         Call dGeMM_('N','T',nAB,nCD,MAB,
     &               1.0d0,Work(ip_C),nAB,Work(ip_3Int),nCD,
     &               1.0d0,xInt,nAB)
         If (Mode.eq.3) Then ! half-and-half, scale by 1/2
            Call dScal_(nAB*nCD,0.5d0,xInt,1)
         End If
         ! Deallocate coefficients
         Call GetMem('CAIC','Free','Real',ip_C,l_C)
         ! Deallocate 3-index integrals
         Call GetMem('CAI3Int','Free','Real',ip_3Int,l_3Int)
      Else If (Mode.eq.2) Then ! non-robust representation
         ! Allocate (J_AB | K_CD)
         l_2Int=MAB*MCD
         Call GetMem('CAI2Int','Allo','Real',ip_2Int,l_2Int)
         ! Compute (J_AB | K_CD)
         Call LDF_ComputeIntegrals_JK_2P(AB,CD,l_2Int,Work(ip_2Int))
         ! Allocate intermediate
         l=MAB*nCD
         Call GetMem('CAIM','Allo','Real',ip,l)
         ! Allocate coefficients
         l_C_AB=nAB*MAB
         l_C_CD=nCD*MCD
         l_C=max(l_C_AB,l_C_CD)
         Call GetMem('CAIC','Allo','Real',ip_C,l_C)
         ! Read coefficients for CD
         Call LDF_CIO_ReadC(CD,Work(ip_C),l_C_CD)
         ! Compute intermediate
         ! M[[J_AB,k_C l_D] = sum_K_CD (J_AB | K_CD)*C[k_C l_D,K_CD]
         Call dGeMM_('N','T',MAB,nCD,MCD,
     &               1.0d0,Work(ip_2Int),MAB,Work(ip_C),nCD,
     &               0.0d0,Work(ip),MAB)
         ! Read coefficients for AB (if different from CD)
         If (AB.ne.CD) Then
            Call LDF_CIO_ReadC(AB,Work(ip_C),l_C_AB)
         End If
         ! Compute integrals: sum_J_AB C[u_A v_B,J_AB]*M[J_AB,k_C l_D]
         Call dGeMM_('N','N',nAB,nCD,MAB,
     &               1.0d0,Work(ip_C),nAB,Work(ip),MAB,
     &               1.0d0,xInt,nAB)
         ! Deallocate coefficients
         Call GetMem('CAIC','Free','Real',ip_C,l_C)
         ! Deallocate intermediate
         Call GetMem('CAIM','Free','Real',ip,l)
         ! Deallocate 2-index integrals
         Call GetMem('CAI2Int','Free','Real',ip_2Int,l_2Int)
      Else
         Call WarningMessage(2,SecNam//': illegal Mode')
         Call LDF_Quit(1)
      End If

#if defined (_DEBUG_)
      If (AB.eq.CD) Then
         If (.not. isSymmetric(xInt,nAB,1.0d-14)) Then
            Call WarningMessage(2,SecNam//': (AB|CD) != (CD|AB)')
            Call LDF_Quit(1)
         End If
      End If
      If (A.eq.B) Then
         nShell_A=LDF_nShell_Atom(A)
         nShell_B=nShell_A
         ipA=LDF_lShell_Atom(A)-1
         ipB=ipA
         l_iAB=nShell_A*nShell_B
         Call GetMem('DBGiAB','Allo','Inte',ip_iAB,l_iAB)
         n=0
         Do jS=1,nShell_B
            jShell=iWork(ipB+jS)
            Do iS=1,nShell_A
               iShell=iWork(ipA+iS)
               iWork(ip_iAB-1+nShell_A*(jS-1)+iS)=n
               n=n+nBasSh(iShell)*nBasSh(jShell)
            End Do
         End Do
         Do kl=1,nCD
            Do jS=1,nShell_B
               jShell=iWork(ipB+jS)
               iS=jS
               iShell=jShell
               Do j=1,nBasSh(jShell)
                  Do i=j+1,nBasSh(iShell)
                     ij=nAB*(kl-1)+iAB(iS,jS)+nBasSh(iShell)*(j-1)+i
                     ji=nAB*(kl-1)+iAB(jS,iS)+nBasSh(jShell)*(i-1)+j
                     If (abs(xInt(ij)-xInt(ji)).gt.1.0d-14) Then
                        Call WarningMessage(2,
     &                               SecNam//': [1] (AB|CD) != (BA|CD)')
                        Call LDF_Quit(1)
                     End If
                  End Do
               End Do
               Do iS=jS+1,nShell_A
                  iShell=iWork(ipA+iS)
                  Do j=1,nBasSh(jShell)
                     Do i=1,nBasSh(iShell)
                        ij=nAB*(kl-1)+iAB(iS,jS)+nBasSh(iShell)*(j-1)+i
                        ji=nAB*(kl-1)+iAB(jS,iS)+nBasSh(jShell)*(i-1)+j
                        If (abs(xint(ij)-xint(ji)).gt.1.0d-14) Then
                           Call WarningMessage(2,
     &                               SecNam//': [2] (AB|CD) != (BA|CD)')
                           Call LDF_Quit(1)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Call GetMem('DBGiAB','Free','Inte',ip_iAB,l_iAB)
      End If
      If (C.eq.D) Then
         nShell_C=LDF_nShell_Atom(C)
         nShell_D=nShell_C
         ipC=LDF_lShell_Atom(C)-1
         ipD=ipC
         l_iCD=nShell_C*nShell_D
         Call GetMem('DBGiCD','Allo','Inte',ip_iCD,l_iCD)
         n=0
         Do jS=1,nShell_D
            jShell=iWork(ipD+jS)
            Do iS=1,nShell_C
               iShell=iWork(ipC+iS)
               iWork(ip_iCD-1+nShell_C*(jS-1)+iS)=n
               n=n+nBasSh(iShell)*nBasSh(jShell)
            End Do
         End Do
         Do ij=1,nAB
            Do lS=1,nShell_D
               lShell=iWork(ipD+lS)
               kS=lS
               kShell=lShell
               Do l=1,nBasSh(lShell)
                  do k=l+1,nBasSh(kShell)
                     kl=iCD(kS,lS)+nBasSh(kShell)*(l-1)+k
                     lk=iCD(lS,kS)+nBasSh(lShell)*(k-1)+l
                     kl=nAB*(kl-1)+ij
                     lk=nAB*(lk-1)+ij
                     If (abs(xInt(kl)-xInt(lk)).gt.1.0d-14) Then
                        Call WarningMessage(2,
     &                               SecNam//': [1] (AB|CD) != (AB|DC)')
                        Call LDF_Quit(1)
                     End If
                  End Do
               End Do
               Do kS=lS+1,nShell_C
                  kShell=iWork(ipC+kS)
                  Do l=1,nBasSh(lShell)
                     Do k=1,nBasSh(kShell)
                        kl=iCD(kS,lS)+nBasSh(kShell)*(l-1)+k
                        lk=iCD(lS,kS)+nBasSh(lShell)*(k-1)+l
                        kl=nAB*(kl-1)+ij
                        lk=nAB*(lk-1)+ij
                        If (abs(xInt(kl)-xInt(lk)).gt.1.0d-14) Then
                           Call WarningMessage(2,
     &                               SecNam//': [2] (AB|CD) != (AB|DC)')
                           Call LDF_Quit(1)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
         Call GetMem('DBGiCD','Free','Inte',ip_iCD,l_iCD)
      End If
#endif

      End

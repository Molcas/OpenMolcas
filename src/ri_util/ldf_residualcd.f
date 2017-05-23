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
      Subroutine LDF_ResidualCD(iAtomPair,ip_CBar,l_CBar,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Minimize the number of two-center functions in the
C     auxiliary basis for atom pair iAtomPair. This is done by Cholesky
C     decomposition of the residual matrix
C
C     R[uv,kl] = (uv|kl) - (uv_|kl_)
C
C     where uv and kl are two-center products included in the auxiliary
C     basis, and
C
C     (uv_|kl_) = sum[J=1,M1] CBar[uv_,J]*CBar[kl_,J]
C
C     is the current DF approximation to those integrals using
C     one-center auxiliary functions only.
C
C     Two-center info for the atom pair is updated in this routine.
C
C     Return code irc is zero on successful completion.
C
C     Note: the underlying assumption is that few two-center functions
C     are needed so that the integrals can be stored in core.
C
      Implicit None
      Integer iAtomPair
      Integer ip_CBar, l_CBar
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"

      External Integral_WrOut_LDF_G

#if defined (_DEBUG_)
      Logical  isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
      External isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
#endif

      Integer  LDF_nBasAux_Pair, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBasAux_Pair, LDF_nShell_Atom, LDF_lShell_Atom

      Integer M, M1
      Integer ijS, ij
      Integer ip_Int, l_Int
      Integer uvS, klS
      Integer ip_SewWrk, l_SewWrk
      Integer iAtom, jAtom
      Integer nShell_iAtom, nShell_jAtom
      Integer iS, jS
      Integer jShell
      Integer ipi, ipj
      Integer nuv
      Integer ip_kOff, l_kOff
      Integer ipk, ip0
      Integer ip_CB, l_CB
      Integer ipCB, ipCBar
      Integer K, K2
      Integer uS, vS
      Integer u, v
      Integer ip_Vec, l_Vec
      Integer ip_ID, l_ID
      Integer nVec

      Real*8 Thr

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer IndxG2
      Integer L2C
      Integer kOff
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)
      L2C(i,j)=iWork(ip_2CList-1+l_2CList_1*(j-1)+i)
      kOff(i,j)=iWork(ip_kOff-1+nShell_iAtom*(j-1)+i)

      ! Init return code
      irc=0

      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         ! Set index array for G matrix
         ! (This includes the one-center functions)
         Call LDF_SetIndxG(iAtomPair)
         ! Get number of auxiliary functions
         M=LDF_nBasAux_Pair(iAtomPair)
         ! Get number of one-center functions
         M1=M-AP_2CFunctions(1,iAtomPair)
         ! Reset indices to exclude one-center functions
         Do ijS=1,l_IndxG2_2
            Do ij=1,l_IndxG2_1
               If (IndxG2(ij,ijS).gt.0) Then
                  iWork(ip_IndxG2-1+l_IndxG2_1*(ijS-1)+ij)=
     &            iWork(ip_IndxG2-1+l_IndxG2_1*(ijS-1)+ij)-M1
               End If
            End Do
         End Do
         nRow_G=AP_2CFunctions(1,iAtomPair)
         ! Allocate memory for integrals
         l_Int=nRow_G**2
         Call GetMem('ResidG','Allo','Real',ip_Int,l_Int)
         ! Compute integrals
         Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
         Call xSetMem_Ints(l_SewWrk)
         Call Cho_dZero(Work(ip_Int),l_Int)
         Do klS=1,l_2CList_2
            SHC=L2C(1,klS)
            SHD=L2C(2,klS)
            SPCD=L2C(3,klS)
            Do uvS=klS,l_2CList_2
               SHA=L2C(1,uvS)
               SHB=L2C(2,uvS)
               SPAB=L2C(3,uvS)
               Call Eval_IJKL(SHA,SHB,SHC,SHD,Work(ip_Int),l_Int,
     &                        Integral_WrOut_LDF_G)
            End Do
         End Do
         Call xRlsMem_Ints()
#if defined (_DEBUG_)
         ! Check integral matrix for symmetry and non-negative diagonals
         ! and Cauchy-Schwarz inequality
         If (.not.isSymmetric(Work(ip_Int),nRow_G,1.0d-15)) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: Integral matrix is not symmetric!'
            irc=irc+1
         End If
         If (.not.hasNonnegativeDiagonal(Work(ip_Int),nRow_G)) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: Integral matrix has negative diagonals!'
            irc=irc+1
         End If
         If (.not.obeysCauchySchwarz(Work(ip_Int),nRow_G,1.0d-14)) Then
            Write(6,'(A)')
     &   'LDF_ResidualCD: Integral matrix does not obey C-S inequality!'
            irc=irc+1
         End If
         If (irc.ne.0) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: => Integral matrix not PSD!!'
            Call LDF_Quit(1)
         End If
#endif
         ! Set up offset array to shell blocks of CBar
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         nShell_iAtom=LDF_nShell_Atom(iAtom)
         nShell_jAtom=LDF_nShell_Atom(jAtom)
         l_kOff=nShell_iAtom*nShell_jAtom
         Call GetMem('kOff','Allo','Inte',ip_kOff,l_kOff)
         ipi=LDF_lShell_Atom(iAtom)-1
         ipj=LDF_lShell_Atom(jAtom)-1
         nuv=0
         Do jS=1,nShell_jAtom
            jShell=iWork(ipj+jS)
            ipk=ip_kOff-1+nShell_iAtom*(jS-1)
            Do iS=1,nShell_iAtom
               iWork(ipk+iS)=nuv
               nuv=nuv+nBasSh(iWork(ipi+iS))*nBasSh(jShell)
            End Do
         End Do
         ! Allocate tmp array for CBar elements corresponding to
         ! two-center auxiliary functions
         l_CB=nRow_G*M1
         Call GetMem('CB','Allo','Real',ip_CB,l_CB)
         ! Copy subblocks of CBar
         Do K=0,M1-1
            ipCB=ip_CB+nRow_G*K
            ipCBar=ip_CBar-1+nuv*K
            Do K2=0,nRow_G-1
               ip0=AP_2CFunctions(2,iAtomPair)+4*K2
               uS=iWork(ip0)
               u=iWork(ip0+1)
               vS=iWork(ip0+2)
               v=iWork(ip0+3)
               Work(ipCB+K2)=
     &            Work(ipCBar+kOff(uS,vS)+nBasSh(iWork(ipi+uS))*(v-1)+u)
            End Do
         End Do
         ! Compute residual matrix
         Call dGeMM_('N','T',nRow_G,nRow_G,M1,
     &               -1.0d0,Work(ip_CB),nRow_G,Work(ip_CB),nRow_G,
     &               1.0d0,Work(ip_Int),nRow_G)
         ! Deallocate tmp CBar array
         Call GetMem('CB','Free','Real',ip_CB,l_CB)
         ! Deallocate offset array
         Call GetMem('kOff','Free','Inte',ip_kOff,l_kOff)
         ! Allocate ID array
         l_ID=nRow_G
         Call GetMem('ID','Allo','Inte',ip_ID,l_ID)
         ! Allocate Cholesky vector array
         l_Vec=nRow_G**2
         Call GetMem('Vec','Allo','Real',ip_Vec,l_Vec)
         ! Incomplete CD of residual matrix
         Thr=Thr_Accuracy
         nVec=0
#if defined (_DEBUG_)
         ! Check residual matrix for symmetry and non-negative diagonals
         ! and Cauchy-Schwarz inequality
         If (.not.isSymmetric(Work(ip_Int),nRow_G,1.0d-15)) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: Residual matrix is not symmetric!'
            irc=irc+1
         End If
         If (.not.hasNonnegativeDiagonal(Work(ip_Int),nRow_G)) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: Residual matrix has negative diagonals!'
            irc=irc+1
         End If
         If (.not.obeysCauchySchwarz(Work(ip_Int),nRow_G,1.0d-14)) Then
            Write(6,'(A)')
     &   'LDF_ResidualCD: Residual matrix does not obey C-S inequality!'
            irc=irc+1
         End If
         If (irc.ne.0) Then
            Write(6,'(A)')
     &      'LDF_ResidualCD: => Residual matrix not PSD!!'
            Call LDF_Quit(1)
         End If
#endif
         Call CD_InCore_P(Work(ip_Int),nRow_G,Work(ip_Vec),nRow_G,
     &                    iWork(ip_ID),nVec,Thr,irc)
         If (irc.ne.0) Then
            Write(6,'(A,I8)')
     &      'LDF_ResidualCD: CD_InCore_P returned code',irc
            Call GetMem('Vec','Free','Real',ip_Vec,l_Vec)
            Call GetMem('ID','Free','Inte',ip_ID,l_ID)
            Call GetMem('ResidG','Free','Real',ip_Int,l_Int)
            Call LDF_UnsetIndxG()
            irc=1
            Return
         End If
         ! Deallocate Cholesky vector array
         Call GetMem('Vec','Free','Real',ip_Vec,l_Vec)
         ! Deallocate integral array
         Call GetMem('ResidG','Free','Real',ip_Int,l_Int)
         ! Reset two-center auxiliary basis info
         Call LDF_Reset2CF(iAtomPair,iWork(ip_ID),nRow_G,nVec)
         ! Deallocate ID array
         Call GetMem('ID','Free','Inte',ip_ID,l_ID)
         ! Unset index arrays for G matrix
         Call LDF_UnsetIndxG()
      End If

c Avoid unused argument warnings
      If (.False.) Call Unused_integer(l_CBar)
      End

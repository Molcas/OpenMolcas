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
      Subroutine LDF_ComputeGMat(iAtomPair,M,G)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Compute G matrix for given atom pair:
C
C        G[JK] = (J|K)
C
C     where J and K are auxiliary functions of atom pair iAtomPair,
C     J,K=1,2,3,....,M.
C
C     M is the total number of auxiliary functions (1C and 2C) minus the
C     number of linearly dependent functions, as computed by the
C     function LDF_nBasAux_Pair(iAtomPair).
C
      Implicit None
      Integer iAtomPair
      Integer M
      Real*8  G(M,M)
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "WrkSpc.fh"

#if defined (_DEBUGPRINT_)
      Character*15 SecNam
      Parameter (SecNam='LDF_ComputeGMat')
      Logical  isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
      External isSymmetric, hasNonnegativeDiagonal, obeysCauchySchwarz
      Integer irc
#endif

      Integer iAtom, jAtom
      Integer iS, jS, klS, uvS, nnShl
      Integer ipi, ipj
      Integer iShell, jShell, kShell, lShell, uShell, vShell, dShell
      Integer nAuxShell_i, nAuxShell_j
      Integer l_IndxG, l_IndxG2, l_2CList
      Integer ip_SewWrk, l_SewWrk
      Integer l_G

      Logical IndxG_Here

      Integer  LDF_nAuxShell_Atom
      Integer  LDF_lAuxShell_Atom
      External LDF_nAuxShell_Atom
      External LDF_lAuxShell_Atom

      External Integral_WrOut_LDF_G

      Integer i, j
      Integer AP_Atoms
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Check if G indexation needs to be done; if so, do it.
      l_IndxG=l_IndxG_1*l_IndxG_2
      l_IndxG2=l_IndxG2_1*l_IndxG2_2
      l_2CList=l_2CList_1*l_2CList_2
      IndxG_Here=l_IndxG.lt.1.and.l_IndxG2.lt.1.and.l_2CList.lt.1
      If (IndxG_Here) Then
         Call LDF_SetIndxG(iAtomPair)
      End If

      ! Allocate memory for Seward
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Get atoms of pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Compute G matrix
      nRow_G=M
      l_G=M*M
      Call Cho_dZero(G,l_G)
      dShell=nShell_Valence+nShell_Auxiliary+1
      SHA=dShell
      SHC=dShell
      ipi=LDF_lAuxShell_Atom(iAtom)-1
      nAuxShell_i=LDF_nAuxShell_Atom(iAtom)
      Do jS=1,nAuxShell_i
         jShell=iWork(ipi+jS)
         SHD=jShell
         Do iS=jS,nAuxShell_i
            iShell=iWork(ipi+iS)
            SHB=iShell
            Call Eval_IJKL(dShell,iShell,dShell,jShell,G,l_G,
     &                     Integral_WrOut_LDF_G)
         End Do
      End Do
      If (jAtom.ne.iAtom) Then
         ipj=LDF_lAuxShell_Atom(jAtom)-1
         nAuxShell_j=LDF_nAuxShell_Atom(jAtom)
         Do jS=1,nAuxShell_j
            jShell=iWork(ipj+jS)
            SHD=jShell
            Do iS=1,nAuxShell_i
               iShell=iWork(ipi+iS)
               SHB=iShell
               Call Eval_IJKL(dShell,iShell,dShell,jShell,G,l_G,
     &                        Integral_WrOut_LDF_G)
            End Do
         End Do
         Do jS=1,nAuxShell_j
            jShell=iWork(ipj+jS)
            SHD=jShell
            Do iS=jS,nAuxShell_j
               iShell=iWork(ipj+iS)
               SHB=iShell
               Call Eval_IJKL(dShell,iShell,dShell,jShell,G,l_G,
     &                        Integral_WrOut_LDF_G)
            End Do
         End Do
      End If
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         nnShl=l_2CList_2
         Do uvS=0,nnShl-1
            uShell=iWork(ip_2CList+3*uvS)
            vShell=iWork(ip_2CList+3*uvS+1)
            SPCD=iWork(ip_2CList+3*uvS+2)
            SHC=uShell
            SHD=vShell
            Do iS=1,nAuxShell_i
               iShell=iWork(ipi+iS)
               SHB=iShell
               Call Eval_IJKL(dShell,iShell,uShell,vShell,G,l_G,
     &                        Integral_WrOut_LDF_G)
            End Do
         End Do
         If (jAtom.ne.iAtom) Then
            ipj=LDF_lAuxShell_Atom(jAtom)-1
            nAuxShell_j=LDF_nAuxShell_Atom(jAtom)
            Do uvS=0,nnShl-1
               uShell=iWork(ip_2CList+3*uvS)
               vShell=iWork(ip_2CList+3*uvS+1)
               SPCD=iWork(ip_2CList+3*uvS+2)
               SHC=uShell
               SHD=vShell
               Do jS=1,nAuxShell_j
                  jShell=iWork(ipj+jS)
                  SHB=jShell
                  Call Eval_IJKL(dShell,jShell,uShell,vShell,G,l_G,
     &                           Integral_WrOut_LDF_G)
               End Do
            End Do
         End If
         Do klS=0,nnShl-1
            kShell=iWork(ip_2CList+3*klS)
            lShell=iWork(ip_2CList+3*klS+1)
            SPCD=iWork(ip_2CList+3*klS+2)
            SHC=kShell
            SHD=lShell
            Do uvS=klS,nnShl-1
               uShell=iWork(ip_2CList+3*uvS)
               vShell=iWork(ip_2CList+3*uvS+1)
               SPAB=iWork(ip_2CList+3*uvS+2)
               SHA=uShell
               SHB=vShell
               Call Eval_IJKL(uShell,vShell,kShell,lShell,G,l_G,
     &                        Integral_WrOut_LDF_G)
            End Do
         End Do
      End If

      ! Release Seward memory
      Call xRlsMem_Ints()

      ! Unset G indexation (if it was set here)
      If (IndxG_Here) Then
         Call LDF_UnsetIndxG()
      End If

#if defined (_DEBUGPRINT_)
      ! Check G matrix for symmetry and non-negative diagonals
      ! and Cauchy-Schwarz inequality
      irc=0
      If (.not.isSymmetric(G,M,1.0d-15)) Then
         Write(6,'(A,A)')
     &   SecNam,': G matrix is not symmetric!'
         irc=irc+1
      End If
      If (.not.hasNonnegativeDiagonal(G,M)) Then
         Write(6,'(A,A)')
     &   SecNam,': G matrix has negative diagonals!'
         irc=irc+1
      End If
      If (.not.obeysCauchySchwarz(G,M,1.0d-14)) Then
         Write(6,'(A,A)')
     &   SecNam,
     &   ': G matrix does not obey the Cauchy-Schwarz inequality!'
         irc=irc+1
      End If
      If (irc.ne.0) Then
         Write(6,'(A,A)')
     &   SecNam,': => G matrix not PSD!!'
         Call LDF_Quit(1)
      End If
#endif

      End

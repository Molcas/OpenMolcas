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
      Subroutine LDF_ComputeIntegrals_uvJ_2P(AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Compute integrals of type (uv|J) where (uv| belongs to atom pair
C     AB and |J) belongs to atom pair CD.
C
      Implicit None
      Integer AB, CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_atom_pair_info.fh"

      Character*27 SecNam
      Parameter (SecNam='LDF_ComputeIntegrals_uvJ_2P')

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      Integer  LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Pair

      Integer l_xInt
      Integer M
      Integer dShell
      Integer A, B, C, D
      Integer iS, jS, iShell, jShell, nShell_A, nShell_B
      Integer ipA, ipB, ipC, ipD
      Integer n, ip0
      Integer nnShl, ijS
      Integer K, ij, ji, ij0, ji0
      Integer ip_SewWrk, l_SewWrk

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer i2CList
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      i2CList(i,j)=iWork(ip_2CList-1+l_2CList_1*(j-1)+i)

      ! Set index arrays for auxiliary functions on atom pair CD
      Call LDF_SetIndxG(CD)

      ! Set dummy shell
      dShell=nShell_Valence+nShell_Auxiliary+1

      ! Get atoms of atom pairs
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)

      ! Set row dimension
      nRow_uvJ=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      ! Set column dimension
      M=LDF_nBasAux_Pair(CD)
      ! Set integral dimension
      l_xInt=nRow_uvJ*M
      ! Check integral array dimension
      If (l_xInt_.lt.l_xInt) Then
         Call WarningMessage(2,SecNam//': integral dimension problem')
         Write(6,'(A,I9,1X,I9)')
     &   'AB,CD...............',AB,CD
         Write(6,'(A,I9,1X,I9,1X,I9,1X,I9)')
     &   'A,B,C,D.............',A,B,C,D
         Write(6,'(A,I9,1X,I9)')
     &   'nRow_uvJ,M..........',nRow_uvJ,M
         Write(6,'(A,I9,1X,I9)')
     &   'nRow_uvJ*M,l_xInt_..',l_xInt,l_xInt_
         Call LDF_Quit(1)
      End If
      ! Init integral array
      Call Cho_dZero(xInt,l_xInt)

      ! Set number of shells on atoms A and B
      nShell_A=LDF_nShell_Atom(A)
      nShell_B=LDF_nShell_Atom(B)

      ! Allocate offset array for rows of (uv|J)
      l_iOff=nShell_A*nShell_B
      Call GetMem('iOff','Allo','Inte',ip_iOff,l_iOff)

      ! Set row offset array
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1
      n=0
      Do jS=1,nShell_B
         jShell=iWork(ipB+jS)
         ip0=ip_iOff-1+nShell_A*(jS-1)
         Do iS=1,nShell_A
            iWork(ip0+iS)=n
            iShell=iWork(ipA+iS)
            n=n+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do
#if defined (_DEBUG_)
      If (n.ne.nRow_uvJ) Then
         Call WarningMessage(2,SecNam//': row dimension problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate memory for Seward
      Call GetMem('GetMax','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals (uv|J)
      ipC=LDF_lAuxShell_Atom(C)-1
      Do iS=1,LDF_nAuxShell_Atom(C)
         iShell=iWork(ipC+iS)
         Call LDF_CI_uvJ(A,B,dShell,iShell,l_xInt,xInt)
      End Do
      If (D.ne.C) Then
         ipD=LDF_lAuxShell_Atom(D)-1
         Do jS=1,LDF_nAuxShell_Atom(D)
            jShell=iWork(ipD+jS)
            Call LDF_CI_uvJ(A,B,dShell,jShell,l_xInt,xInt)
         End Do
      End If
      If (AP_2CFunctions(1,CD).gt.0) Then
         nnShl=l_2CList_2
         Do ijS=1,nnShl
            iShell=i2CList(1,ijS)
            jShell=i2CList(2,ijS)
            SPAB=i2CList(3,ijS) ! SPAB stored in localdf_int.fh
            Call LDF_CI_uvJ(A,B,iShell,jShell,l_xInt,xInt)
         End Do
      End If

      ! Release Seward memory
      Call xRlsMem_Ints()

      ! Postprocessing for A=B:
      ! Only lower (block) triangle computed => transpose to get upper
      If (A.eq.B) Then
         ipA=LDF_lShell_Atom(A)-1
         ipB=ipA
         Do K=1,M
            ip0=nRow_uvJ*(K-1)
            Do jS=2,nShell_B
               jShell=iWork(ipB+jS)
               Do iS=1,jS-1
                  iShell=iWork(ipA+iS)
                  ij0=ip0+iWork(ip_iOff-1+nShell_A*(jS-1)+iS)
                  ji0=ip0+iWork(ip_iOff-1+nShell_B*(iS-1)+jS)
                  Do j=1,nBasSh(jShell)
                     Do i=1,nBasSh(iShell)
                        ij=ij0+nBasSh(iShell)*(j-1)+i
                        ji=ji0+nBasSh(jShell)*(i-1)+j
                        xInt(ij)=xInt(ji)
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If

      ! Deallocate offset array for columns of (uv|J)
      Call GetMem('iOff','Free','Inte',ip_iOff,l_iOff)
      l_iOff=0
      ip_iOff=0
      nRow_uvJ=0
      iOffuv=0

      ! Unset index arrays for auxiliary functions on atom pair CD
      Call LDF_UnsetIndxG()

      End

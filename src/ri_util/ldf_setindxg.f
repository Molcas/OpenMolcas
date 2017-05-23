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
      Subroutine LDF_SetIndxG(iAtomPair)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: set index arrays for auxiliary functions.
C
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam='LDF_SetIndxG')

      Integer iAtom, jAtom
      Integer iS, jS, ijS, nnShl
      Integer maxAuxSh, nShell, maxAuxSP
      Integer ip, ip0, ipIG0, ipi, ipj
      Integer k, n, nn, ij, i2CF
      Integer iShell, jShell
      Integer nShell_i, nShell_j
      Integer ipLD, iLD
      Integer iLD_Atom, iLD_Shell, iLD_Indx
      Integer l_IndxG, l_IndxG2, l_2CList
      Integer M, niS, njS

      Integer  LDF_nShell_Atom, LDF_nBasSh_Atom
      Integer  LDF_lShell_Atom
      Integer  LDF_nAuxShell_Atom, LDF_nBasAuxSh_Atom
      Integer  LDF_lAuxShell_Atom, LDF_nBasAux_Pair
      External LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_lShell_Atom
      External LDF_nAuxShell_Atom, LDF_nBasAuxSh_Atom
      External LDF_lAuxShell_Atom, LDF_nBasAux_Pair

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Check if index arrays are already set
      ! Deallocate them if so....
      l_IndxG=l_IndxG_1*l_IndxG_2
      If (l_IndxG.gt.0) Then
         Call WarningMessage(1,SecNam//': IndxG is already allocated!')
         Call GetMem('IndxG','Free','Inte',ip_IndxG,l_IndxG)
         ip_IndxG=0
         l_IndxG=0
         l_IndxG_1=0
         l_IndxG_2=0
      End If
      l_IndxG2=l_IndxG2_1*l_IndxG2_2
      If (l_IndxG2.gt.0) Then
         Call WarningMessage(1,SecNam//': IndxG2 is already allocated!')
         Call GetMem('IndxG2','Free','Inte',ip_IndxG2,l_IndxG2)
         ip_IndxG2=0
         l_IndxG2=0
         l_IndxG2_1=0
         l_IndxG2_2=0
      End If
      l_2CList=l_2CList_1*l_2CList_2
      If (l_2CList.gt.0) Then
         Call WarningMessage(1,SecNam//': 2CList is already allocated!')
         Call GetMem('G_2CList','Free','Inte',ip_2CList,l_2CList)
         ip_2CList=0
         l_2CList_1=0
         l_2CList_2=0
      End If

      ! Get atoms of pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Get number of auxiliary functions
      M=LDF_nBasAux_Pair(iAtomPair)

      ! Find largest auxiliary shell of atom pair
      maxAuxSh=0
      Do iS=1,LDF_nAuxShell_Atom(iAtom)
         maxAuxSh=max(maxAuxSh,LDF_nBasAuxSh_Atom(iS,iAtom))
      End Do
      If (iAtom.ne.jAtom) Then
         Do jS=1,LDF_nAuxShell_Atom(jAtom)
            maxAuxSh=max(maxAuxSh,LDF_nBasAuxSh_Atom(jS,jAtom))
         End Do
      End If

      ! Allocate index arrays for G matrix
      nShell=nShell_Valence+nShell_Auxiliary+1
      l_IndxG_1=maxAuxSh
      l_IndxG_2=nShell
      l_IndxG=maxAuxSh*nShell
      Call GetMem('IndxG','Allo','Inte',ip_IndxG,l_IndxG)
      l_IndxG2_1=0
      l_IndxG2_2=0
      l_IndxG2=0
      ip_IndxG2=0
      nnShl=0
      maxAuxSP=0
      nShell_i=0
      nShell_j=0
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         nShell_i=LDF_nShell_Atom(iAtom)
         nShell_j=LDF_nShell_Atom(jAtom)
         Do jS=1,nShell_j
            n=LDF_nBasSh_Atom(jS,jAtom)
            Do iS=1,nShell_i
               maxAuxSP=max(maxAuxSP,n*LDF_nBasSh_Atom(iS,iAtom))
            End Do
         End Do
         nnShl=nShell_i*nShell_j
         l_IndxG2_1=maxAuxSP
         l_IndxG2_2=nnShl
         l_IndxG2=maxAuxSP*nnShl
         Call GetMem('IndxG2','Allo','Inte',ip_IndxG2,l_IndxG2)
      End If

      ! Set index arrays for G matrix
      ! IndxG(j,jShl): row/column index in G matrix for auxiliary basis
      !                function j in auxiliary shell jShl of atom iAtom
      !                or jAtom.
      ! IndxG2(ij,ijShl): row/column index in G matrix for 2-center
      !                   auxiliary function ij of shell pair ijShl of
      !                   atom pair iAtomPair.
      ! Note: row/column indices corresponding to linearly dependent
      ! auxiliary functions are zero.
      Call iZero(iWork(ip_IndxG),l_IndxG)
      ip=LDF_lAuxShell_Atom(iAtom)-1
      k=0
      Do iS=1,LDF_nAuxShell_Atom(iAtom)
         iShell=iWork(ip+iS)
         ipIG0=ip_IndxG-1+maxAuxSh*(iShell-1)
         Do i=1,nBasSh(iShell)
            k=k+1
            iWork(ipIG0+i)=k
         End Do
      End Do
      If (jAtom.ne.iAtom) Then
         ip=LDF_lAuxShell_Atom(jAtom)-1
         Do jS=1,LDF_nAuxShell_Atom(jAtom)
            jShell=iWork(ip+jS)
            ipIG0=ip_IndxG-1+maxAuxSh*(jShell-1)
            Do i=1,nBasSh(jShell)
               k=k+1
               iWork(ipIG0+i)=k
            End Do
         End Do
      End If
      If (AP_1CLinDep(1,iAtomPair).gt.0) Then
         ! zero entries corresponding to lin dep functions
         ipLD=AP_1CLinDep(2,iAtomPair)-1
         Do iLD=0,AP_1CLinDep(1,iAtomPair)-1
            iLD_Atom=iWork(ipLD+3*iLD+1)
            iLD_Shell=iWork(ipLD+3*iLD+2)
            iLD_Indx=iWork(ipLD+3*iLD+3)
            ip=LDF_lAuxShell_Atom(iLD_Atom)-1
            iShell=iWork(ip+iLD_Shell)
            ipIG0=ip_IndxG-1+maxAuxSh*(iShell-1)
            iWork(ipIG0+iLD_Indx)=0
         End Do
         ! reset row/column indices
         ip=LDF_lAuxShell_Atom(iAtom)-1
         k=0
         Do iS=1,LDF_nAuxShell_Atom(iAtom)
            iShell=iWork(ip+iS)
            ipIG0=ip_IndxG-1+maxAuxSh*(iShell-1)
            Do i=1,nBasSh(iShell)
               If (iWork(ipIG0+i).gt.0) Then
                  k=k+1
                  iWork(ipIG0+i)=k
               End If
            End Do
         End Do
         If (jAtom.ne.iAtom) Then
            ip=LDF_lAuxShell_Atom(jAtom)-1
            Do jS=1,LDF_nAuxShell_Atom(jAtom)
               jShell=iWork(ip+jS)
               ipIG0=ip_IndxG-1+maxAuxSh*(jShell-1)
               Do i=1,nBasSh(jShell)
                  If (iWork(ipIG0+i).gt.0) Then
                     k=k+1
                     iWork(ipIG0+i)=k
                  End If
               End Do
            End Do
         End If
      End If
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         ! Include 2-center functions
         Call iZero(iWork(ip_IndxG2),l_IndxG2)
         ip0=AP_2CFunctions(2,iAtomPair)-1
         Do i2CF=0,AP_2CFunctions(1,iAtomPair)-1
            ip=ip0+4*i2CF
            iS=iWork(ip+1)
            i=iWork(ip+2)
            jS=iWork(ip+3)
            j=iWork(ip+4)
            ijS=nShell_i*(jS-1)+iS
            ij=LDF_nBasSh_Atom(iS,iAtom)*(j-1)+i
            k=k+1
            iWork(ip_IndxG2-1+maxAuxSP*(ijS-1)+ij)=k
         End Do
      End If
#if defined (_DEBUG_)
      If (k.ne.M) Then
         Write(6,'(A,I8,1X,I8)') SecNam//': k,M=',k,M
         Call WarningMessage(2,SecNam//': G matrix dimension problem!')
         Call LDF_Quit(1)
      End If
#endif

      ! Set up shell pair list for 2-center functions
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         ! count shell pairs
         n=0
         Do ijS=0,nnShl-1
            nn=0
            Do ij=0,maxAuxSP-1
               nn=nn+iWork(ip_IndxG2+maxAuxSP*ijS+ij)
            End Do
            n=n+min(nn,1)
         End Do
         ! Allocate list
         l_2CList_1=3
         l_2CList_2=n
         l_2CList=l_2CList_1*l_2CList_2
         Call GetMem('G_2CList','Allo','Inte',ip_2CList,l_2CList)
         ! Set list
         niS=LDF_nShell_Atom(iAtom)
         ipi=LDF_lShell_Atom(iAtom)-1
         njS=LDF_nShell_Atom(jAtom)
         ipj=LDF_lShell_Atom(jAtom)-1
         n=0
         Do ijS=1,nnShl
            nn=0
            Do ij=0,maxAuxSP-1
               nn=nn+iWork(ip_IndxG2+maxAuxSP*(ijS-1)+ij)
            End Do
            If (nn.gt.0) Then
               jS=(ijS-1)/niS+1
               iS=ijS-niS*(jS-1)
               iWork(ip_2CList+3*n)=iWork(ipi+iS)
               iWork(ip_2CList+3*n+1)=iWork(ipj+jS)
               iWork(ip_2CList+3*n+2)=ijS
               n=n+1
            End If
         End Do
      End If

      ! Set nRow_G
      nRow_G=M

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_UnsetIndxG()
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Unset index arrays for G matrix.
C
      Implicit None
#include "localdf_int.fh"

      Integer l_IndxG, l_IndxG2, l_2CList

      l_IndxG=l_IndxG_1*l_IndxG_2
      If (l_IndxG.gt.0) Then
         Call GetMem('IndxG','Free','Inte',ip_IndxG,l_IndxG)
         ip_IndxG=0
         l_IndxG=0
         l_IndxG_1=0
         l_IndxG_2=0
      End If
      l_IndxG2=l_IndxG2_1*l_IndxG2_2
      If (l_IndxG2.gt.0) Then
         Call GetMem('IndxG2','Free','Inte',ip_IndxG2,l_IndxG2)
         ip_IndxG2=0
         l_IndxG2=0
         l_IndxG2_1=0
         l_IndxG2_2=0
      End If
      l_2CList=l_2CList_1*l_2CList_2
      If (l_2CList.gt.0) Then
         Call GetMem('G_2CList','Free','Inte',ip_2CList,l_2CList)
         ip_2CList=0
         l_2CList_1=0
         l_2CList_2=0
      End If
      If (l_iOff.gt.0) Then
         Call GetMem('iOff','Free','Inte',ip_iOff,l_iOff)
         ip_iOff=0
         l_iOff=0
      End If

      nRow_G=0
      nRow_uvJ=0
      iOffuv=0

      End

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
      Subroutine LDF_ComputeIntegrals_uvJ(iAtomPair,l_xInt,xInt)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Compute integrals of type (uv|J) for atom pair iAtomPair.
C
      Implicit None
      Integer iAtomPair
      Integer l_xInt
      Real*8  xInt(l_xInt)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_atom_pair_info.fh"

#if defined (_DEBUG_)
      Character*24 SecNam
      Parameter (SecNam='LDF_ComputeIntegrals_uvJ')
#endif

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      Integer  LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Pair

      Integer M
      Integer dShell
      Integer iAtom, jAtom
      Integer iS, jS, iShell, jShell, nShell_iAtom, nShell_jAtom
      Integer ipi, ipj
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

      ! Init integral array
      Call Cho_dZero(xInt,l_xInt)

      ! Set dummy shell
      dShell=nShell_Valence+nShell_Auxiliary+1

      ! Get atoms of atom pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Set row dimension
      nRow_uvJ=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
#if defined (_DEBUG_)
      If (l_xInt.ne.nRow_uvJ*LDF_nBasAux_Pair(iAtomPair)) Then
         Call WarningMessage(2,SecNam//': integral dimension problem')
         Write(6,'(A,I9)')
     &   'iAtomPair...........',iAtomPair
         Write(6,'(A,I9,1X,I9)')
     &   'iAtom,jAtom.........',iAtom,jAtom
         Write(6,'(A,I9,1X,I9)')
     &   'nRow_uvJ,M..........',nRow_uvJ,LDF_nBasAux_Pair(iAtomPair)
         Write(6,'(A,I9,1X,I9)')
     &   'nRow_uvJ*M,l_xInt...',nRow_uvJ*LDF_nBasAux_Pair(iAtomPair),
     &                          l_xInt
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate offset array for rows of (uv|J)
      nShell_iAtom=LDF_nShell_Atom(iAtom)
      nShell_jAtom=LDF_nShell_Atom(jAtom)
      l_iOff=nShell_iAtom*nShell_jAtom
      Call GetMem('iOff','Allo','Inte',ip_iOff,l_iOff)

      ! Set row offset array
      ipi=LDF_lShell_Atom(iAtom)-1
      ipj=LDF_lShell_Atom(jAtom)-1
      n=0
      Do jS=1,nShell_jAtom
         jShell=iWork(ipj+jS)
         ip0=ip_iOff-1+nShell_iAtom*(jS-1)
         Do iS=1,nShell_iAtom
            iWork(ip0+iS)=n
            iShell=iWork(ipi+iS)
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
      ipi=LDF_lAuxShell_Atom(iAtom)-1
      Do iS=1,LDF_nAuxShell_Atom(iAtom)
         iShell=iWork(ipi+iS)
         Call LDF_CI_uvJ(iAtom,jAtom,dShell,iShell,l_xInt,xInt)
      End Do
      If (jAtom.ne.iAtom) Then
         ipj=LDF_lAuxShell_Atom(jAtom)-1
         Do jS=1,LDF_nAuxShell_Atom(jAtom)
            jShell=iWork(ipj+jS)
            Call LDF_CI_uvJ(iAtom,jAtom,dShell,jShell,l_xInt,xInt)
         End Do
      End If
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         nnShl=l_2CList_2
         Do ijS=1,nnShl
            iShell=i2CList(1,ijS)
            jShell=i2CList(2,ijS)
            SPAB=i2CList(3,ijS) ! SPAB stored in localdf_int.fh
            Call LDF_CI_uvJ(iAtom,jAtom,iShell,jShell,l_xInt,xInt)
         End Do
      End If

      ! Release Seward memory
      Call xRlsMem_Ints()

      ! Postprocessing for iAtom=jAtom:
      ! Only lower (block) triangle computed => transpose to get upper
      If (iAtom.eq.jAtom) Then
         M=LDF_nBasAux_Pair(iAtomPair)
         ipi=LDF_lShell_Atom(iAtom)-1
         ipj=ipi
         Do K=1,M
            ip0=nRow_uvJ*(K-1)
            Do jS=2,nShell_jAtom
               jShell=iWork(ipj+jS)
               Do iS=1,jS-1
                  iShell=iWork(ipi+iS)
                  ij0=ip0+iWork(ip_iOff-1+nShell_iAtom*(jS-1)+iS)
                  ji0=ip0+iWork(ip_iOff-1+nShell_jAtom*(iS-1)+jS)
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

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_CI_uvJ(kAtom,lAtom,iShell,jShell,l_xInt,xInt)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Compute integrals (uv|J) where J is an auxiliary function in
C     shell pair iShell,jShell [iShell=dummy-shell if 1-center].
C
      Implicit None
      Integer kAtom, lAtom
      Integer iShell, jShell
      Integer l_xInt
      Real*8  xInt(l_xInt)
#include "WrkSpc.fh"
#include "localdf_int.fh"

      Character*10 SecNam
      Parameter (SecNam='LDF_CI_uvJ')

      External Int_LDF_uvJ

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer nShell_kAtom, nShell_lAtom
      Integer kS, lS
      Integer kShell, lShell
      Integer ipk, ipl

      nShell_kAtom=LDF_nShell_Atom(kAtom)
      nShell_lAtom=LDF_nShell_Atom(lAtom)
      ipk=LDF_lShell_Atom(kAtom)-1
      ipl=LDF_lShell_Atom(lAtom)-1

      SHA=iShell
      SHB=jShell

      If (kAtom.eq.lAtom) Then
         ! NOTE: only lower (block) triangular part is computed!!!
         Do lS=1,nShell_lAtom
            lShell=iWork(ipl+lS)
            SHD=lShell
            Do kS=lS,nShell_kAtom
               kShell=iWork(ipk+kS)
               SHC=kShell
               ! iOffuv = row offset
               iOffuv=iWork(ip_iOff-1+nShell_kAtom*(lS-1)+kS)
               Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,l_xInt,
     &                        Int_LDF_uvJ)
            End Do
         End Do
      Else If (kAtom.gt.lAtom) Then
         Do lS=1,nShell_lAtom
            lShell=iWork(ipl+lS)
            SHD=lShell
            Do kS=1,nShell_kAtom
               kShell=iWork(ipk+kS)
               SHC=kShell
               ! iOffuv = row offset
               iOffuv=iWork(ip_iOff-1+nShell_kAtom*(lS-1)+kS)
               Call Eval_IJKL(iShell,jShell,kShell,lShell,xInt,l_xInt,
     &                        Int_LDF_uvJ)
            End Do
         End Do
      Else
         Call WarningMessage(2,SecNam//': kAtom<lAtom')
         Call LDF_Quit(1)
      End If

      End

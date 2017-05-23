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
      Subroutine LDF_RemoveLinDep(iAtomPair,Z,ID,M,nZ)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Handle linear dependence: remove linearly dependent rows of Z and
C     update atom pair info.
C
      Implicit None
      Integer iAtomPair
      Integer M, nZ
      Real*8  Z(M,nZ)
      Integer ID(M)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"

      Character*8 Label

      Integer ip_Included, l_Included
      Integer ip, l
      Integer ip0
      Integer iAtom, jAtom
      Integer n1CLinDep, n2CFunctions
      Integer iS, jS, iShell, jShell
      Integer K, ij, ijS, i2CF
      Integer nShell_i

      Integer  LDF_nShell_Atom, LDF_nBasSh_Atom
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Integer i, j
      Integer nBasSh, Included, AP_Atoms, AP_1CLinDep, AP_2CFunctions
      Integer IndxG, IndxG2
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      Included(i)=iWork(ip_Included-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      IndxG(i,j)=iWork(ip_IndxG-1+l_IndxG_1*(j-1)+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)

      ! Return if nothing to do
      If (nZ.ge.M) Return

      ! Set up Included array
      ! Included(J)=0: aux. function J is not included
      ! Included(J)=1: aux. function J is included
      l_Included=M
      Call GetMem('Incl','Allo','Inte',ip_Included,l_Included)
      Call iZero(iWork(ip_Included),l_Included)
      ip0=ip_Included-1
      Do J=1,nZ
         iWork(ip0+ID(J))=1
      End Do

      ! Remove linearly dependent rows of Z vectors
      ! = the leading nZ x nZ part of Z will contain the
      !   Z vectors corresponding to the linearly independent
      !   auxiliary basis functions
      ! Note that no reordering is done, implying that the Z vectors
      ! reflect the ordering of the G matrix with linearly dependent
      ! columns AND rows removed.
      l=nZ
      Call GetMem('ZTmp','Allo','Real',ip,l)
      Do J=1,nZ
         K=0
         Do I=1,M
            If (Included(I).eq.1) Then
               Work(ip+K)=Z(I,J)
               K=K+1
            End If
         End Do
         Call dCopy_(nZ,Work(ip),1,Z(1,J),1)
      End Do
      Call GetMem('ZTmp','Free','Real',ip,l)

      ! Get atoms of atom pair
      iAtom=AP_Atoms(1,iAtomPair)
      jAtom=AP_Atoms(2,iAtomPair)

      ! Count number of linearly dependent 1-center functions
      n1CLinDep=AP_1CLinDep(1,iAtomPair)
      ip0=LDF_lAuxShell_Atom(iAtom)-1
      Do iS=1,LDF_nAuxShell_Atom(iAtom)
         iShell=iWork(ip0+iS)
         Do i=1,nBasSh(iShell)
            J=IndxG(i,iShell)
            If (J.gt.0) Then
               If (Included(J).eq.0) Then
                  n1CLinDep=n1CLinDep+1
               End If
            End If
         End Do
      End Do
      If (jAtom.ne.iAtom) Then
         ip0=LDF_lAuxShell_Atom(jAtom)-1
         Do jS=1,LDF_nAuxShell_Atom(jAtom)
            jShell=iWork(ip0+jS)
            Do i=1,nBasSh(jShell)
               J=IndxG(i,jShell)
               If (J.gt.0) Then
                  If (Included(J).eq.0) Then
                     n1CLinDep=n1CLinDep+1
                  End If
               End If
            End Do
         End Do
      End If

      ! Reallocate and reset 1C lindep info
      If (n1CLinDep.gt.0) Then
         Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
         l=3*AP_1CLinDep(1,iAtomPair)
         If (l.gt.0) Then
            ip=AP_1CLinDep(2,iAtomPair)
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
         l=3*n1CLinDep
         Call GetMem(Label,'Allo','Inte',ip,l)
         iWork(ip_AP_1CLinDep+2*(iAtompair-1))=n1CLinDep
         iWork(ip_AP_1CLinDep+2*(iAtompair-1)+1)=ip
         ip0=LDF_lAuxShell_Atom(iAtom)-1
         n1CLinDep=0
         Do iS=1,LDF_nAuxShell_Atom(iAtom)
            iShell=iWork(ip0+iS)
            Do i=1,nBasSh(iShell)
               J=IndxG(i,iShell)
               If (J.gt.0) Then
                  If (Included(J).eq.0) Then
                     iWork(ip+3*n1CLinDep)=iAtom
                     iWork(ip+3*n1CLinDep+1)=iS
                     iWork(ip+3*n1CLinDep+2)=i
                     n1CLinDep=n1CLinDep+1
                  End If
               Else
                  iWork(ip+3*n1CLinDep)=iAtom
                  iWork(ip+3*n1CLinDep+1)=iS
                  iWork(ip+3*n1CLinDep+2)=i
                  n1CLinDep=n1CLinDep+1
               End If
            End Do
         End Do
         If (jAtom.ne.iAtom) Then
            ip0=LDF_lAuxShell_Atom(jAtom)-1
            Do jS=1,LDF_nAuxShell_Atom(jAtom)
               jShell=iWork(ip0+jS)
               Do i=1,nBasSh(jShell)
                  J=IndxG(i,jShell)
                  If (J.gt.0) Then
                     If (Included(J).eq.0) Then
                        iWork(ip+3*n1CLinDep)=jAtom
                        iWork(ip+3*n1CLinDep+1)=jS
                        iWork(ip+3*n1CLinDep+2)=i
                        n1CLinDep=n1CLinDep+1
                     End If
                  Else
                     iWork(ip+3*n1CLinDep)=jAtom
                     iWork(ip+3*n1CLinDep+1)=jS
                     iWork(ip+3*n1CLinDep+2)=i
                     n1CLinDep=n1CLinDep+1
                  End If
               End Do
            End Do
         End If
      End If

      ! Remove 2C functions (if any)
      If (AP_2CFunctions(1,iAtomPair).gt.0) Then
         ! Count contributing 2C functions
         n2CFunctions=0
         nShell_i=LDF_nShell_Atom(iAtom)
         ip0=AP_2CFunctions(2,iAtomPair)-1
         Do i2CF=0,AP_2CFunctions(1,iAtomPair)-1
            ip=ip0+4*i2CF
            iS=iWork(ip+1)
            i=iWork(ip+2)
            jS=iWork(ip+3)
            j=iWork(ip+4)
            ijS=nShell_i*(jS-1)+iS
            ij=LDF_nBasSh_Atom(iS,iAtom)*(j-1)+i
            K=IndxG2(ij,ijS)
            If (K.gt.0) Then
               If (Included(K).eq.1) Then
                  n2CFunctions=n2CFunctions+1
               End If
            End If
         End Do
         ! Set data for those contributing
         ! Not needed, of course, if none has been removed by the CD
         If (n2CFunctions.lt.AP_2CFunctions(1,iAtomPair)) Then
            Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
            If (n2CFunctions.eq.0) Then
               k=4*AP_2CFunctions(1,iAtomPair)
               ip0=AP_2CFunctions(2,iAtomPair)
               Call GetMem(Label,'Free','Inte',ip0,k)
               iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=0
               iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=0
            Else
               l=4*n2CFunctions
               Call GetMem(Label,'Allo','Inte',ip,l)
               nShell_i=LDF_nShell_Atom(iAtom)
               ip0=AP_2CFunctions(2,iAtomPair)-1
               n2CFunctions=0
               Do i2CF=0,AP_2CFunctions(1,iAtomPair)-1
                  iS=iWork(ip0+4*i2CF+1)
                  i=iWork(ip0+4*i2CF+2)
                  jS=iWork(ip0+4*i2CF+3)
                  j=iWork(ip0+4*i2CF+4)
                  ijS=nShell_i*(jS-1)+iS
                  ij=LDF_nBasSh_Atom(iS,iAtom)*(j-1)+i
                  K=IndxG2(ij,ijS)
                  If (K.gt.0) Then
                     If (Included(K).eq.1) Then
                        iWork(ip+4*n2CFunctions)=iS
                        iWork(ip+4*n2CFunctions+1)=i
                        iWork(ip+4*n2CFunctions+2)=jS
                        iWork(ip+4*n2CFunctions+3)=j
                        n2CFunctions=n2CFunctions+1
                     End If
                  End If
               End Do
               k=4*AP_2CFunctions(1,iAtomPair)
               ip0=AP_2CFunctions(2,iAtomPair)
               Call GetMem(Label,'Free','Inte',ip0,k)
               iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=n2CFunctions
               iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip
            End If
         End If
      End If

      ! Deallocate Included array
      Call GetMem('Incl','Free','Inte',ip_Included,l_Included)

      End

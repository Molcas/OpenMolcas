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
      Subroutine LDF_CIO_ReadC_WithLinDep(iAtomPair,C,l_C)
      Implicit None
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
      Call LDF_CIO_ReadC_wLD(iAtomPair,C,l_C)
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CIO_ReadC_wLD(iAtomPair,C,l_C)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Read LDF fitting coefficients for atom pair iAtomPair.
C              The coefficients are read from memory (buffer) if
C              possible (though the caller need not worry about that).
C              Columns of C corresponding to linearly dependent
C              auxiliary (one-center) functions are included and contain
C              zeros. For reading without linearly dependent columns,
C              use LDF_CIO_ReadC.
C
C     NOTE: this routine may require memory corresponding to at least 1
C           column of the coefficient matrix.
C
      Implicit None
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
#include "WrkSpc.fh"
#include "ldf_cio.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"

      Logical  LDF_isLinDep
      External LDF_isLinDep

      Integer  LDF_nBas_Atom
      Integer  LDF_nBasAux_Pair
      Integer  LDF_nBasAux_Pair_WithLinDep
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBas_Atom
      External LDF_nBasAux_Pair
      External LDF_nBasAux_Pair_WithLinDep
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Integer iAddr
      Integer iAtom, jAtom
      Integer l, l_WithLinDep, l_WithoutLinDep
      Integer ii, iS, ipi
      Integer jj, jS, ipj
      Integer ipC
      Integer ip_Scr, l_Scr

      Integer i, j
      Integer nBasSh
      Integer AP_DiskC
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_DiskC(i)=iWork(ip_AP_DiskC-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      If (Lu_LDFC.lt.1) Then
         Call WarningMessage(2,'LDF_CIO_ReadC_wLD: Lu_LDFC<1')
         Call LDF_Quit(1)
      End If

      If (AP_1CLinDep(1,iAtomPair).eq.0) Then
         ! Use LDF_CIO_ReadC if no linear dependence
         Call LDF_CIO_ReadC(iAtomPair,C,l_C)
         Return
      Else If (AP_1CLinDep(1,iAtomPair).lt.0) Then
         ! This would be bug....
         Call WarningMessage(2,
     &         'LDF_CIO_ReadC_wLD: AP_1CLinDep<0 !?!')
         Call LDF_Quit(1)
      Else
         ! lin dep columns are returned as zero vectors
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         l=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
         l_WithLinDep=l*LDF_nBasAux_Pair_WithLinDep(iAtomPair)
         If (l_WithLinDep.gt.l_C) Then
            Call WarningMessage(2,
     &         'LDF_CIO_ReadC_wLD: insufficient array dimension')
            Call LDF_Quit(1)
         End If
         If (iAtomPair.gt.LastAtomPair) Then ! read from disk
            Call GetMem('GetMax','Max ','Real',ip_Scr,l_Scr)
            l_WithoutLinDep=l*LDF_nBasAux_Pair(iAtomPair)
            If (l_Scr.ge.l_WithoutLinDep) Then ! read full C
               l_Scr=l_WithoutLinDep
               Call GetMem('RdCScr1','Allo','Real',ip_Scr,l_Scr)
               Call LDF_CIO_ReadC(iAtomPair,Work(ip_Scr),l_Scr)
               iAddr=ip_Scr
               ipC=1
               ipi=LDF_lAuxShell_Atom(iAtom)-1
               Do iS=1,LDF_nAuxShell_Atom(iAtom)
                  Do ii=1,nBasSh(iWork(ipi+iS))
                     If (LDF_isLinDep(ii,iS,iAtom,iAtomPair)) Then
                        Call Cho_dZero(C(ipC),l)
                     Else
                        Call dCopy_(l,Work(iAddr),1,C(ipC),1)
                        iAddr=iAddr+l
                     End If
                     ipC=ipC+l
                  End Do
               End Do
               If (jAtom.ne.iAtom) Then
                  ipj=LDF_lAuxShell_Atom(jAtom)-1
                  Do jS=1,LDF_nAuxShell_Atom(jAtom)
                     Do jj=1,nBasSh(iWork(ipj+jS))
                        If (LDF_isLinDep(jj,jS,jAtom,iAtomPair)) Then
                           Call Cho_dZero(C(ipC),l)
                        Else
                           Call dCopy_(l,Work(iAddr),1,C(ipC),1)
                           iAddr=iAddr+l
                        End If
                        ipC=ipC+l
                     End Do
                  End Do
               End If
               If (AP_2CFunctions(1,iAtomPair).gt.0) Then
                  Call dCopy_(l*AP_2CFunctions(1,iAtomPair),
     &                       Work(iAddr),1,C(ipC),1)
               End If
               Call GetMem('RdCScr1','Free','Real',ip_Scr,l_Scr)
           Else ! read 1 column of C at a time
               If (l_Scr.lt.l) Then
                  Call WarningMessage(2,
     &                 'LDF_CIO_ReadC_wLD: Insufficient memory!')
                  Call LDF_Quit(1)
               End If
               l_Scr=l
               Call GetMem('RdCScr2','Allo','Real',ip_Scr,l_Scr)
               iAddr=AP_DiskC(iAtomPair)
               ipC=1
               ipi=LDF_lAuxShell_Atom(iAtom)-1
               Do iS=1,LDF_nAuxShell_Atom(iAtom)
                  Do ii=1,nBasSh(iWork(ipi+iS))
                     If (LDF_isLinDep(ii,iS,iAtom,iAtomPair)) Then
                        Call Cho_dZero(C(ipC),l)
                     Else
                        Call dDAFile(Lu_LDFC,2,Work(ip_Scr),l,iAddr)
                        Call dCopy_(l,Work(ip_Scr),1,C(ipC),1)
                     End If
                     ipC=ipC+l
                  End Do
               End Do
               If (jAtom.ne.iAtom) Then
                  ipj=LDF_lAuxShell_Atom(jAtom)-1
                  Do jS=1,LDF_nAuxShell_Atom(jAtom)
                     Do jj=1,nBasSh(iWork(ipj+jS))
                        If (LDF_isLinDep(jj,jS,jAtom,iAtomPair)) Then
                           Call Cho_dZero(C(ipC),l)
                        Else
                           Call dDAFile(Lu_LDFC,2,Work(ip_Scr),l,iAddr)
                           Call dCopy_(l,Work(ip_Scr),1,C(ipC),1)
                        End If
                        ipC=ipC+l
                     End Do
                  End Do
               End If
               If (AP_2CFunctions(1,iAtomPair).gt.0) Then
                  Call dDAFile(Lu_LDFC,2,
     &                         C(ipC),l*AP_2CFunctions(1,iAtomPair),
     &                         iAddr)
               End If
               Call GetMem('RdCScr2','Free','Real',ip_Scr,l_Scr)
           End If
         Else ! copy from buffer
            iAddr=iWork(ip_LDFC_Blocks-1+iAtomPair)
            ipC=1
            ipi=LDF_lAuxShell_Atom(iAtom)-1
            Do iS=1,LDF_nAuxShell_Atom(iAtom)
               Do ii=1,nBasSh(iWork(ipi+iS))
                  If (LDF_isLinDep(ii,iS,iAtom,iAtomPair)) Then
                     Call Cho_dZero(C(ipC),l)
                  Else
                     Call dCopy_(l,Work(iAddr),1,C(ipC),1)
                     iAddr=iAddr+l
                  End If
                  ipC=ipC+l
               End Do
            End Do
            If (jAtom.ne.iAtom) Then
               ipj=LDF_lAuxShell_Atom(jAtom)-1
               Do jS=1,LDF_nAuxShell_Atom(jAtom)
                  Do jj=1,nBasSh(iWork(ipj+jS))
                     If (LDF_isLinDep(jj,jS,jAtom,iAtomPair)) Then
                        Call Cho_dZero(C(ipC),l)
                     Else
                        Call dCopy_(l,Work(iAddr),1,C(ipC),1)
                        iAddr=iAddr+l
                     End If
                     ipC=ipC+l
                  End Do
               End Do
            End If
            If (AP_2CFunctions(1,iAtomPair).gt.0) Then
               Call dCopy_(l*AP_2CFunctions(1,iAtomPair),
     &                    Work(iAddr),1,C(ipC),1)
            End If
         End If
      End If

      End

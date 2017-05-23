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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_GetQuadraticDiagonal(ip)
C
C     NOTE: this code is redundant!!
C     As of March 2011, the diagonal is stored quadratically from the
C     beginning.
C
C     TODO/FIXME: remove this file.
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: get AP diagonals stored quadratically (for A=B).
C              On return, ip can be used in the same way as ip_AP_Diag
C              from ldf_atom_pair_info.fh.
C              To deallocate the quadratic diagonal,
C              call LDF_FreeQuadraticDiagonal(ip).
C
      Implicit None
      Integer ip
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "ldf_qdiag.fh"

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer AB
      Integer A, B
      Integer ip_QDAA, l_QDAA
      Integer ipQ
      Integer nSA, ipA
      Integer iSA, iSB
      Integer iShellA, iShellB
      Integer iA, iB
      Integer n
      Integer ip_iOffPk, l_iOffPk

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer ipDPk
      Integer iOffPk
      Integer iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      ipDPk(i)=iWork(ip_AP_Diag-1+i)
      iOffPk(i,j)=iWork(ip_iOffPk-1+nSA*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      Call WarningMessage(2,
     &'LDF_GetQuadraticDiagonal: this code is redundant, don''t use it')
      Call LDF_Quit(1)

      If (l_AP_QDiag.eq.NumberOfAtomPairs) Then ! already there
         ip=ip_AP_QDiag
      Else If (l_AP_QDiag.eq.0) Then ! set up quadratic diagonal blocks
         l_AP_QDiag=NumberOfAtomPairs
         Call GetMem('QDPtr','Allo','Inte',ip_AP_QDiag,l_AP_QDiag)
         ip=ip_AP_QDiag
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            If (A.eq.B) Then
               l_QDAA=LDF_nBas_Atom(A)**2
               Call GetMem('QDAA','Allo','Real',ip_QDAA,l_QDAA)
               ipQ=ip_QDAA-1
               nSA=LDF_nShell_Atom(A)
               ipA=LDF_lShell_Atom(A)-1
               l_iOffPk=nSA**2
               Call GetMem('iOffPk','Allo','Inte',ip_iOffPk,l_iOffPk)
               n=0
               Do iSA=1,nSA
                  iShellA=iWork(ipA+iSA)
                  Do iSB=1,iSA-1
                     iWork(ip_iOffPk-1+nSA*(iSB-1)+iSA)=n
                     iWork(ip_iOffPk-1+nSA*(iSA-1)+iSB)=n
                     iShellB=iWork(ipA+iSB)
                     n=n+nBasSh(iShellA)*nBasSh(iShellB)
                  End Do
                  iWork(ip_iOffPk-1+nSA*(iSA-1)+iSA)=n
                  n=n+nBasSh(iShellA)*(nBasSh(iShellA)+1)/2
               End Do
               Do iSB=1,nSA
                  iShellB=iWork(ipA+iSB)
                  Do iSA=1,nSA
                     iShellA=iWork(ipA+iSA)
                     If (iSA.eq.iSB) Then
                        Do iB=1,nBasSh(iShellB)
                           Do iA=1,nBasSh(iShellA)
                              Work(ipQ+nBasSh(iShellA)*(iB-1)+iA)=
     &                        Work(ipDPk(AB)-1+iOffPk(iSA,iSA)
     &                             +iTri(iA,iB))
                           End Do
                        End Do
                     Else If (iSA.gt.iSB) Then
                        Call dCopy_(nBasSh(iShellA)*nBasSh(iShellB),
     &                             Work(ipDPk(AB)+iOffPk(iSA,iSB)),1,
     &                             Work(ipQ+1),1)
                     Else
                        Do iA=1,nBasSh(iShellA)
                           Call dCopy_(nBasSh(iShellB),
     &                                Work(ipDPk(AB)+iOffPk(iSB,iSA)
     &                                     +nBasSh(iShellB)*(iA-1)),1,
     &                                Work(ipQ+iA),nBasSh(iShellA))
                        End Do
                     End If
                     ipQ=ipQ+nBasSh(iShellA)*nBasSh(iShellB)
                  End Do
               End Do
               Call GetMem('iOffPk','Free','Inte',ip_iOffPk,l_iOffPk)
               iWork(ip_AP_QDiag-1+AB)=ip_QDAA
            Else
               iWork(ip_AP_QDiag-1+AB)=iWork(ip_AP_Diag-1+AB)
            End If
         End Do
      Else
         Call WarningMessage(2,'QDiag management corrupted!')
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_FreeQuadraticDiagonal(ip)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: free quadratic diagonals (see also LDF_QuadraticDiagonal)
C
      Implicit None
      Integer ip
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_qdiag.fh"

      Character*25 SecNam
      Parameter (SecNam='LDF_FreeQuadraticDiagonal')

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer AB
      Integer A, B
      Integer ip_QDAA, l_QDAA

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (l_AP_QDiag.gt.0) Then
         If (ip.ne.ip_AP_QDiag) Then
            Call WarningMessage(2,SecNam//': ip mismatch!')
            Call LDF_Quit(1)
         End If
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            If (A.eq.B) Then
               l_QDAA=LDF_nBas_Atom(A)**2
               ip_QDAA=iWork(ip_AP_QDiag-1+AB)
               Call GetMem('QDAA','Free','Real',ip_QDAA,l_QDAA)
            End If
         End Do
         Call GetMem('QDPtr','Free','Inte',ip_AP_QDiag,l_AP_QDiag)
         ip_AP_QDiag=0
         l_AP_QDiag=0
         ip=0
      End If

      End

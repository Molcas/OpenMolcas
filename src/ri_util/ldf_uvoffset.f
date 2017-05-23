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
      Subroutine LDF_uvOffset(AB,nSA,nSB,iOff)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute offset array to shell blocks of, e.g.,
C              fitting coefficients. Diagonal atom pairs (A=B)
C              are treated as quadratic (i.e. no LT storage).
C
      Implicit None
      Integer AB
      Integer nSA
      Integer nSB
      Integer iOff(nSA,nSB)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*12 SecNam
      Parameter (SecNam='LDF_uvOffset')

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom
#if defined (_DEBUG_)
      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom
#endif

      Integer A, B
      Integer ipA, ipB
      Integer iSA, iSB
      Integer iShellB
      Integer iCount

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      If (nSA.ne.LDF_nShell_Atom(A) .or. nSB.ne.LDF_nShell_Atom(B)) Then
         Call WarningMessage(2,SecNam//': illegal nSA/nSB')
         Call LDF_Quit(1)
      End If
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1
      iCount=0
      Do iSB=1,nSB
         iShellB=iWork(ipB+iSB)
         Do iSA=1,nSA
            iOff(iSA,iSB)=iCount
            iCount=iCount+nBasSh(iWork(ipA+iSA))*nBasSh(iShellB)
         End Do
      End Do
#if defined (_DEBUG_)
      If (iCount.ne.LDF_nBas_Atom(A)*LDF_nBas_Atom(B)) Then
         Call WarningMessage(2,SecNam//': iCount error')
         Call LDF_Quit(1)
      End If
#endif

      End

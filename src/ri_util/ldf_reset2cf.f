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
      Subroutine LDF_Reset2CF(iAtomPair,ID,M2,M2_new)
C
C     Thomas Bondo Pedersen, July 2010.
C     - based on LDF_RemoveLinDep by T.B. Pedersen.
C
C     Purpose: Reset two-center function info.
C
      Implicit None
      Integer iAtomPair
      Integer M2, M2_new
      Integer ID(M2_new)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_int.fh"

      Integer  LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_nShell_Atom, LDF_nBasSh_Atom

      Character*8 Label

      Integer ip, l, ip_new, l_new, K
      Integer ip_Included, l_Included
      Integer nShell_iAtom, n2CF, i2CF, iS, jS, ijS, ij
      Integer iAtom

      Integer i, j
      Integer Included
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer IndxG2
      Included(i)=iWork(ip_Included-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)

      If (M2.ne.AP_2CFunctions(1,iAtomPair)) Then
         Call WarningMessage(2,'LDF_Reset2CF: Illegal M2')
         Call LDF_Quit(1)
      End If
      If (M2_new.eq.M2) Return

      Write(Label,'(A,I5.5)') '2CF',iAtomPair-1

      If (M2_new.lt.0) Then
         Call WarningMessage(2,'LDF_Reset2CF: M2_new<0')
         Call LDF_Quit(1)
      Else If (M2_new.eq.0) Then
         l=4*M2
         ip=AP_2CFunctions(2,iAtomPair)
         Call GetMem(Label,'Free','Inte',ip,l)
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=0
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=0
      Else If (M2_new.lt.M2) Then
         l_Included=M2
         Call GetMem('Included','Allo','Inte',ip_Included,l_Included)
         Call iZero(iWork(ip_Included),l_Included)
         ip=ip_Included-1
         Do K=1,M2_new
            iWork(ip+ID(K))=1
         End Do
         l_new=4*M2_new
         Call GetMem(Label,'Allo','Inte',ip_new,l_new)
         ip=AP_2CFunctions(2,iAtomPair)
         iAtom=AP_Atoms(1,iAtomPair)
         nShell_iAtom=LDF_nShell_Atom(iAtom)
         n2CF=0
         Do i2CF=0,M2-1
            iS=iWork(ip+4*i2CF)
            i=iWork(ip+4*i2CF+1)
            jS=iWork(ip+4*i2CF+2)
            j=iWork(ip+4*i2CF+3)
            ijS=nShell_iAtom*(jS-1)+iS
            ij=LDF_nBasSh_Atom(iS,iAtom)*(j-1)+i
            K=IndxG2(ij,ijS)
            If (K.gt.0) Then
               If (Included(K).eq.1) Then
                  iWork(ip_new+4*n2CF)=iS
                  iWork(ip_new+4*n2CF+1)=i
                  iWork(ip_new+4*n2CF+2)=jS
                  iWork(ip_new+4*n2CF+3)=j
                  n2CF=n2CF+1
               End If
            End If
         End Do
#if defined (_DEBUG_)
         If (n2CF.ne.M2_new) Then
            Call WarningMessage(2,'LDF_Reset2CF: n2CF != M2_new')
            Call LDF_Quit(1)
         End If
#endif
         l=4*M2
         ip=AP_2CFunctions(2,iAtomPair)
         Call GetMem(Label,'Free','Inte',ip,l)
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=M2_new
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip_new
         Call GetMem('Included','Free','Inte',ip_Included,l_Included)
      Else
         Call WarningMessage(2,'LDF_Reset2CF: M2_new>M2')
         Call LDF_Quit(1)
      End If

      End

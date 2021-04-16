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
      Subroutine LDF_Map2CF(AB,nRow,nCol,Map)
C
C     Thomas Bondo Pedersen, January 2011
C
C     Purpose: set up map from 2CFunctions of atom pair AB to the
C              product functions of AB
C              Note that, if nCol>1 and A=B, two product functions are
C              mapped to: one with u>v and one with u<v.
C
      Implicit None
      Integer AB
      Integer nRow, nCol
      Integer Map(nRow,nCol)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*10 SecNam
      Parameter (SecNam='LDF_Map2CF')

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer A, B
      Integer nSA, nSB
      Integer ipA, ipB
      Integer iSA, iSB
      Integer iShellA, iShellB
      Integer ip_kOff, l_kOff
      Integer i2C
      Integer iA, iB, iAB, iBA

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer List2C
      Integer kOff
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      List2C(i,j)=iWork(AP_2CFunctions(2,AB)-1+4*(j-1)+i)
      kOff(i,j)=iWork(ip_kOff-1+nSA*(j-1)+i)

      If (AP_2CFunctions(1,AB).gt.0) Then
         If (nRow.lt.AP_2CFunctions(1,AB)) Then
            Call WarningMessage(2,
     &                           SecNam//': insufficient row dimension')
            Call LDF_Quit(1)
         End If
         If (nCol.lt.1) Then
            Call WarningMessage(2,
     &                           SecNam//': insufficient col dimension')
            Call LDF_Quit(1)
         End If
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nSA=LDF_nShell_Atom(A)
         nSB=LDF_nShell_Atom(B)
         ipA=LDF_lShell_Atom(A)-1
         l_kOff=nSA*nSB
         Call GetMem('kOff','Allo','Inte',ip_kOff,l_kOff)
         Call LDF_uvOffset(AB,nSA,nSB,iWork(ip_kOff))
         Do i2C=1,AP_2CFunctions(1,AB)
            iSA=List2C(1,i2C)
            iA=List2C(2,i2C)
            iSB=List2C(3,i2C)
            iB=List2C(4,i2C)
            iShellA=iWork(ipA+iSA)
            iAB=kOff(iSA,iSB)+nBasSh(iShellA)*(iB-1)+iA
            Map(i2C,1)=iAB
         End Do
         If (A.eq.B .and. nCol.gt.1) Then
            ipB=ipA
            Do i2C=1,AP_2CFunctions(1,AB)
               iSA=List2C(1,i2C)
               iA=List2C(2,i2C)
               iSB=List2C(3,i2C)
               iB=List2C(4,i2C)
               iShellB=iWork(ipB+iSB)
               iBA=kOff(iSB,iSA)+nBasSh(iShellB)*(iA-1)+iB
               Map(i2C,2)=iBA
            End Do
         End If
         Call GetMem('kOff','Free','Inte',ip_kOff,l_kOff)
      End If

      End

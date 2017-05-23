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
      Subroutine LDF_Q2LT(A,Q,LT)
C
C     Thomas Bondo Pedersen, March 2011.
C
      Implicit None
      Integer A
      Real*8  Q(*)
      Real*8  LT(*)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

      Character*8 SecNam
      Parameter (SecNam='LDF_Q2LT')

      Integer  LDF_nShell_Atom, LDF_lShell_Atom, LDF_nBas_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom, LDF_nBas_Atom

      Integer nS, ipS
      Integer ip_iOff, l_iOff
      Integer iS, jS
      Integer iShell, jShell
      Integer n, ip0, ip1
      Integer ipLT, ipQ, l
      Integer ii, jj

      Integer i, j
      Integer nBasSh
      Integer iOff
      Integer iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iOff(i,j)=iWork(ip_iOff-1+nS*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      nS=LDF_nShell_Atom(A)
      ipS=LDF_lShell_Atom(A)-1
      l_iOff=nS**2
      Call GetMem('iOff','Allo','Inte',ip_iOff,l_iOff)
      n=0
      ip0=ip_iOff-1
      Do jS=1,nS
         jShell=iWork(ipS+jS)
         ip1=ip0+nS*(jS-1)
         Do iS=1,nS
            iShell=iWork(ipS+iS)
            iWork(ip1+iS)=n
            n=n+nBasSh(iShell)*nBasSh(jShell)
         End Do
      End Do
      If (n.ne.LDF_nBas_Atom(A)**2) Then
         Call WarningMessage(2,SecNam//': dimension error')
         Call LDF_Quit(1)
      End If
      ipLT=1
      Do iS=1,nS
         iShell=iWork(ipS+iS)
         Do jS=1,iS-1
            jShell=iWork(ipS+jS)
            l=nBasSh(iShell)*nBasSh(jShell)
            ipQ=iOff(iS,jS)+1
            Call dCopy_(l,Q(ipQ),1,LT(ipLT),1)
            ipLT=ipLT+l
         End Do
         ipLT=ipLT-1
         Do jj=1,nBasSh(iShell)
            ipQ=iOff(iS,iS)+nBasSh(iShell)*(jj-1)
            Do ii=jj,nBasSh(iShell)
               LT(ipLT+iTri(ii,jj))=Q(ipQ+ii)
            End Do
         End Do
         ipLT=ipLT+1+nBasSh(iShell)*(nBasSh(iShell)+1)/2
      End Do
      ipLT=ipLT-1
      l=LDF_nBas_Atom(A)*(LDF_nBas_Atom(A)+1)/2
      If (ipLT.ne.l) Then
         Call WarningMessage(2,SecNam//': ipLT != l')
         Call LDF_Quit(1)
      End If
      Call GetMem('iOff','Free','Inte',ip_iOff,l_iOff)

      End

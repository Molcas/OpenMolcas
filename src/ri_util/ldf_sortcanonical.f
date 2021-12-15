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
* Copyright (C) 2012, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SortCanonical(AB,l_Scr,Scr,
     &                             OffsetDefined,
     &                             l_iOff1,l_iOff2,iOff,
     &                             l_X,X)
C
C     Thomas Bondo Pedersen, September 2012.
C
C     Change ordering from shell-pair blocks to 'canonical' order
C     [nbas(A) x nbas(B)].
C
      Implicit None
      Integer AB
      Integer l_Scr
      Real*8  Scr(l_Scr)
      Logical OffSetDefined
      Integer l_iOff1, l_iOff2
      Integer iOff(l_iOff1,l_iOff2)
      Integer l_X
      Real*8  X(l_X)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"

      Character*17 SecNam
      Parameter (SecNam='LDF_SortCanonical')

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer A, B
      Integer nA, nB, nAB
      Integer nSA, nSB
      Integer ipA, ipB
      Integer iS, jS
      Integer iShell, jShell
      Integer ni, nj
      Integer i, j
      Integer ij
      Integer ip0

      Logical TestIt
#if defined (_DEBUGPRINT_)
      Parameter (TestIt=.True.)
#else
      Parameter (TestIt=.False.)
#endif
      Real*8   ChkSum(2), ChkNorm(2), diff(2)
      real*8 ddot_
      external ddot_
      Real*8   Tol
      Parameter (Tol=1.0d-14)

      Integer m, n
      Integer nBasSh
      Integer AP_Atoms
      nBasSh(m)=iWork(ip_nBasSh-1+m)
      AP_Atoms(m,n)=iWork(ip_AP_Atoms-1+2*(n-1)+m)

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nA=LDF_nBas_Atom(A)
      nB=LDF_nBas_Atom(B)
      nAB=nA*nB
      nSA=LDF_nShell_Atom(A)
      nSB=LDF_nShell_Atom(B)
      If (l_X.lt.nAB .or. l_Scr.lt.nAB .or. l_iOff1.lt.nSA .or.
     &    l_iOff2.lt.nSB) Then
         Call WarningMessage(2,
     &                SecNam//': input arrays not properly dimensioned')
         Write(6,'(A,7I10)') 'l_X,l_Scr,nAB,l_iOff1,nSA,l_iOff2,nSB=',
     &   l_X,l_Scr,nAB,l_iOff1,nSA,l_iOff2,nSB
         Call LDF_Quit(1)
      End If
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1

      If (.not.OffsetDefined) Then
         ij=0
         Do jS=1,nSB
            jShell=iWork(ipB+jS)
            nj=nBasSh(jShell)
            Do iS=1,nSA
               iShell=iWork(ipA+iS)
               ni=nBasSh(iShell)
               iOff(iS,jS)=ij
               ij=ij+ni*nj
            End Do
         End Do
         OffsetDefined=.True.
      End If

      Call dCopy_(nAB,X,1,Scr,1)
      ij=0
      Do jS=1,nSB
         jShell=iWork(ipB+jS)
         nj=nBasSh(jShell)
         Do j=1,nj
            Do iS=1,nSA
               iShell=iWork(ipA+iS)
               ni=nBasSh(iShell)
               ip0=iOff(iS,jS)+ni*(j-1)
               Do i=1,ni
                  ij=ij+1
                  X(ij)=Scr(ip0+i)
               End Do
            End Do
         End Do
      End Do

      ! Self test
      If (TestIt) Then
         ChkSum(1)=0.0d0
         Do i=1,nAB
            ChkSum(1)=ChkSum(1)+Scr(i)
         End Do
         ChkSum(2)=0.0d0
         Do i=1,nAB
            ChkSum(2)=ChkSum(2)+X(i)
         End Do
         ChkNorm(1)=sqrt(dDot_(nAB,Scr,1,Scr,1))
         ChkNorm(2)=sqrt(dDot_(nAB,X,1,X,1))
         diff(1)=ChkSum(1)-ChkSum(2)
         diff(2)=ChkNorm(1)-ChkNorm(2)
         If (abs(diff(1)).gt.Tol .or. abs(diff(2)).gt.Tol) Then
            Call WarningMessage(2,SecNam//': sorting failed!')
            Write(6,'(A,1P,3D20.12)') 'Sum(src,trg,diff)= ',
     &      ChkSum(1),ChkSum(2),diff(1)
            Write(6,'(A,1P,3D20.12)') 'Norm(src,trg,diff)=',
     &      ChkNorm(1),ChkNorm(2),diff(2)
            Call LDF_Quit(1)
         End If
      End If

      End

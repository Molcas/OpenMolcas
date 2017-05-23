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
* Copyright (C) 2010,2011, Thomas Bondo Pedersen                       *
************************************************************************
      Subroutine LDF_ComputeC(iAtomPair,ip_C,l_C,ip_Z,l_Z,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Input:
C        iAtomPair - index of atom pair in atom pair list
C        ip_C, l_C - pointer and dimension of array
C                    containing CBar on entry and C on exit
C        ip_Z, l_Z - pointer and dimension of array containing Z
C
C     Output:
C        irc - return code
C
C     Return codes:
C        irc=-1: input error
C        irc=0:  all is fine
C        irc=1:  internal error
C
C     Purpose:
C
C     Compute
C
C     C[uv,J] = { CBar[uv,J] - sum(I=J+1,M) C[uv,I]*Z[I,J] }/Z[J,J]
C
C     where CBar is contained in C on entry. (For Z vectors, see
C     LDF_ComputeCBar.)
C
C     If any of the columns of C has norm less than 1.0d-16, the column
C     and corresponding auxiliary basis function will be removed.
C
      Implicit None
      Integer iAtomPair
      Integer ip_C, l_C
      Integer ip_Z, l_Z
      Integer irc
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      Character*12 SecNam
      Parameter (SecNam='LDF_ComputeC')
#endif

      Real*8 Small
      Parameter (Small=1.0d-24)

      Integer  LDF_nBasAux_Pair, LDF_nBas_Atom
      External LDF_nBasAux_Pair, LDF_nBas_Atom

      Integer nuv
      Integer M

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Init return code
      irc=0
#if defined (_DEBUG_)
      If (iAtomPair.lt.1 .or. iAtomPair.gt.NumberOfAtomPairs) Then
         Write(6,'(A,A,I8)')
     &   SecNam,': iAtomPair out of bounds:',iAtomPair
         irc=-1
         Return
      End If
#endif

      ! Get number of AO products
      nuv=LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &   *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))

      ! Get number of auxiliary functions
      M=LDF_nBasAux_Pair(iAtomPair)

#if defined (_DEBUG_)
      If (l_C.ne.nuv*M .or. l_Z.ne.M*(M+1)/2) Then
         Write(6,'(A,A)') SecNam,': array dimension error'
         Write(6,'(A,2I10)') 'l_C,nuv*M:',l_C,nuv*M
         Write(6,'(A,2I10)') 'l_Z,M*(M+1)/2:',l_Z,M*(M+1)/2
         irc=-1
         Return
      End If
#endif

      ! Compute C
      Call LDF_ComputeC_(Work(ip_C),Work(ip_Z),nuv,M,irc)

      ! Clean rows corresponding to 2C functions, i.e.
      ! C(uv,J)=delta(uv,J) where J is a 2C function
      Call LDF_CleanC(iAtomPair,Work(ip_C),nuv,M)

      ! Remove columns with norm <=Small, and update lindep/2CFun info
      ! accordingly.
      Call LDF_RemoveColumnsOfC(iAtomPair,Work(ip_C),nuv,M,Small)

#ifndef _DEBUG_
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(l_C)
         Call Unused_integer(l_Z)
      End If
#endif
      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_ComputeC_(C,Z,nuv,M,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Compute C (see LDF_ComputeC)
C
      Implicit None
      Integer nuv, M
      Real*8  C(nuv,M)
      Real*8  Z(*)
      Integer irc

      Real*8 ZJJ_inv
#if defined (_DEBUG_)
      Real*8 TooSmall
      Parameter (TooSmall=1.0d-20)
#endif

      Integer i, j
      Integer iTri
      iTri(i,j)=i*(i-1)/2+j

      irc=0
#if defined (_DEBUG_)
      Do J=1,M
         If (Z(iTri(J,J)).lt.TooSmall) Then
            irc=irc+1
         End If
      End Do
      If (irc.gt.0) Then
         Write(6,'(I9,A,I9,A,1P,D16.6)')
     &   irc,' of',M,' Z diagonals are smaller than',TooSmall
         irc=1
         Return
      End If
#endif

      Do J=M,1,-1
         ZJJ_inv=1.0d0/Z(iTri(J,J))
         Call dScal_(nuv,ZJJ_inv,C(1,J),1)
         Do I=J-1,1,-1
            Call dAXPY_(nuv,-Z(iTri(J,I)),C(1,J),1,C(1,I),1)
         End Do
      End Do

      End
************************************************************************
************************************************************************
************************************************************************
      Subroutine LDF_RemoveColumnsOfC(iAtomPair,C,nuv,M,Thr)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Remove columns (and corresponding auxiliary functions)
C     with 2-norm less than Thr.
C
      Implicit None
      Integer iAtomPair
      Integer nuv, M
      Real*8  C(nuv,M)
      Real*8  Thr
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"

      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      Integer  LDF_nShell_Atom, LDF_nBasSh_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nShell_Atom, LDF_nBasSh_Atom

      Character*8 Label

      Integer ip_TOC, l_TOC, ip, MM
      Integer I, J, K, l
      Integer iAtom, jAtom
      Integer iS, jS, iShell, jShell
      Integer ip0, n1CLinDep, n2CFunctions
      Integer nShell_i, ij, ijS
      Integer i2CF

      real*8 ddot_
      external ddot_

      Integer nBasSh
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      Integer IndxG
      Integer IndxG2
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      IndxG(i,j)=iWork(ip_IndxG-1+l_IndxG_1*(j-1)+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)

      l_TOC=M
      Call GetMem('RCCTOC','Allo','Inte',ip_TOC,l_TOC)

      MM=0
      ip=ip_TOC-1
      Do I=1,M
         If (sqrt(dDot_(nuv,C(1,I),1,C(1,I),1)).gt.Thr) Then
            MM=MM+1
            iWork(ip+I)=MM
         Else
            iWork(ip+I)=0
         End If
      End Do

      If (MM.lt.M) Then
         Do I=1,M
            J=iWork(ip+I)
            If (J.gt.0 .and. J.lt.I) Then
               Call dCopy_(nuv,C(1,I),1,C(1,J),1)
            End If
         End Do
         !Call Cho_dZero(C(1,MM+1),nuv*(M-MM))
         Call LDF_SetIndxG(iAtomPair)
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         n1CLinDep=AP_1CLinDep(1,iAtomPair)
         ip0=LDF_lAuxShell_Atom(iAtom)-1
         Do iS=1,LDF_nAuxShell_Atom(iAtom)
            iShell=iWork(ip0+iS)
            Do i=1,nBasSh(iShell)
               J=IndxG(i,iShell)
               If (J.gt.0) Then
                  If (iWork(ip+J).eq.0) Then
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
                     If (iWork(ip+J).eq.0) Then
                        n1CLinDep=n1CLinDep+1
                     End If
                  End If
               End Do
            End Do
         End If
         If (n1CLinDep.gt.0) Then
            Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
            l=3*AP_1CLinDep(1,iAtomPair)
            If (l.gt.0) Then
               ip=AP_1CLinDep(2,iAtomPair)
               Call GetMem(Label,'Free','Inte',ip,l)
            End If
            l=3*n1CLinDep
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_AP_1CLinDep+2*(iAtomPair-1))=n1CLinDep
            iWork(ip_AP_1CLinDep+2*(iAtomPair-1)+1)=ip
            ip0=LDF_lAuxShell_Atom(iAtom)-1
            n1CLinDep=0
            Do iS=1,LDF_nAuxShell_Atom(iAtom)
               iShell=iWork(ip0+iS)
               Do i=1,nBasSh(iShell)
                  J=IndxG(i,iShell)
                  If (J.gt.0) Then
                     If (iWork(ip_TOC-1+J).eq.0) Then
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
                        If (iWork(ip_TOC-1+J).eq.0) Then
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
         If (AP_2CFunctions(1,iAtomPair).gt.0) Then
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
                  If (iWork(ip_TOC-1+K).gt.0) Then
                     n2CFunctions=n2CFunctions+1
                  End If
               End If
            End Do
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
                        If (iWork(ip_TOC-1+K).gt.0) Then
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
         Call LDF_UnsetIndxG()
      End If

      Call GetMem('RCCTOC','Free','Inte',ip_TOC,l_TOC)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CleanC(AB,C,nuv,M)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: clean C such that C(uv,J)=delta(uv,J) for 2C functions J.
C              This is only approximately true without this clean-up.
C
      Implicit None
      Integer AB
      Integer nuv, M
      Real*8  C(nuv,M)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

#if defined (_DEBUG_)
      Character*10 SecNam
      Parameter (SecNam='LDF_CleanC')
#endif

      Integer  LDF_nBasAux_Atom
      External LDF_nBasAux_Atom

      Integer ip_Map, l_Map
      Integer A, B
      Integer M1
      Integer i2C, j2C
      Integer nRow, nCol

      Integer i, j
      Integer Map
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      Map(i,j)=iWork(ip_Map-1+AP_2CFunctions(1,AB)*(j-1)+i)

      If (AP_2CFunctions(1,AB).lt.1) Return ! nothing to do
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nRow=AP_2CFunctions(1,AB)
      If (A.eq.B) Then
         nCol=2
      Else
         nCol=1
      End If
      l_Map=nRow*nCol
      Call GetMem('2CMap','Allo','Inte',ip_Map,l_Map)
      Call LDF_Map2CF(AB,nRow,nCol,iWork(ip_Map))
      M1=LDF_nBasAux_Atom(A)
      If (B.ne.A) Then
         M1=M1+LDF_nBasAux_Atom(B)
      End If
      M1=M1-AP_1CLinDep(1,AB)
#if defined (_DEBUG_)
      If ((M1+AP_2CFunctions(1,AB)).gt.M) Then
         Call WarningMessage(2,SecNam//': M1+M2 != M')
         Call LDF_Quit(1)
      End If
#endif
      Do J=1,M1
         Do i2C=1,AP_2CFunctions(1,AB)
            C(Map(i2C,1),J)=0.0d0
         End Do
      End Do
      Do j2C=1,AP_2CFunctions(1,AB)
         J=M1+j2C
         Do i2C=1,j2C-1
            C(Map(i2C,1),J)=0.0d0
         End Do
         C(Map(j2C,1),J)=1.0d0
         Do i2C=j2C+1,AP_2CFunctions(1,AB)
            C(Map(i2C,1),J)=0.0d0
         End Do
      End Do
      If (A.eq.B) Then
         Do J=1,M1
            Do i2C=1,AP_2CFunctions(1,AB)
               C(Map(i2C,2),J)=0.0d0
            End Do
         End Do
         Do j2C=1,AP_2CFunctions(1,AB)
            J=M1+j2C
            Do i2C=1,j2C-1
               C(Map(i2C,2),J)=0.0d0
            End Do
            C(Map(j2C,1),J)=1.0d0
            Do i2C=j2C+1,AP_2CFunctions(1,AB)
               C(Map(i2C,2),J)=0.0d0
            End Do
         End Do
      End If
      Call GetMem('2CMap','Free','Inte',ip_Map,l_Map)

      End

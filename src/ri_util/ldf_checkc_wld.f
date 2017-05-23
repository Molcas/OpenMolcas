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
      Subroutine LDF_CheckAllC_wLD(DoPrint)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: check all coefficients.
C
      Implicit None
      Logical DoPrint
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nBas_Atom, LDF_nBasAux_Pair_wLD

      Real*8 NormA, NormB, NormAB, NormTot

      Integer AB
      Integer ip_C, l_C, l_thisC

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      l_C=0
      Do AB=1,NumberOfAtomPairs
         l_C=max(l_C,LDF_nBas_Atom(AP_Atoms(1,AB))
     &              *LDF_nBas_Atom(AP_Atoms(2,AB))
     &              *LDF_nBasAux_Pair_wLD(AB))
      End Do
      Call GetMem('CCHE_C','Allo','Real',ip_C,l_C)
      Do AB=1,NumberOfAtomPairs
         l_thisC=LDF_nBas_Atom(AP_Atoms(1,AB))
     &          *LDF_nBas_Atom(AP_Atoms(2,AB))
     &          *LDF_nBasAux_Pair_wLD(AB)
         Call LDF_CIO_ReadC_wLD(AB,Work(ip_C),l_thisC)
         Call LDF_CheckC_wLD(DoPrint,AB,l_thisC,Work(ip_C),
     &                       NormA,NormB,NormAB,NormTot)
      End Do
      Call GetMem('CCHE_C','Free','Real',ip_C,l_C)

      End
      Subroutine LDF_CheckC_wLD(DoPrint,AB,l_C,C,
     &                          NormA,NormB,NormAB,NormTot)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute and optionally print norm of fitting coefficients
C              with linearly dependent functions retained.
C              It is checked that coefficients are symmetric when A=B,
C              i.e. C(uv,J)=C(vu,J).
C              It is checked that coefficients corresponding to LinDep
C              functions are zero.
C              It is checked that 2CFunction coefficients are one for
C              the corresponding product functions.
C
      Implicit None
      Logical DoPrint
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Real*8 NormA
      Real*8 NormB
      Real*8 NormAB
      Real*8 NormTot
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Character*14 SecNam
      Parameter (SecNam='LDF_CheckC_wLD')

      Real*8 Small
      Parameter (Small=1.0d-20)
      Real*8 Tol
      Parameter (Tol=1.0d-12)

      Logical  isConstant
      External isConstant

      Integer  LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nBasAux_Atom, LDF_nBasAux_Pair_wLD
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      real*8 ddot_
      external ddot_

      Integer A, B
      Integer n, M
      Integer ip, l
      Integer ip_Col, l_Col, ipCol
      Integer iCount
      Integer nASA, nASB
      Integer ipA, ipB
      Integer kAtom, kS, k
      Integer i1CLinDep
      Integer kShell
      Integer nErr, mErr
      Integer ip_Map, l_Map
      Integer i2CF, j2CF
      Integer iRow, iCol
      Integer ip0
      Integer ip_iOff, l_iOff
      Integer nSA, iSA, jSA, iA, jA, iShell, jShell, ijA, jiA

      Real*8 NormCol, Tst

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_1CLinDep
      Integer AP_2CFunctions
      Integer iOff
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      iOff(i,j)=iWork(ip_iOff-1+nSA*(j-1)+i)

      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      n=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      M=LDF_nBasAux_Pair_wLD(AB)
      If (n.lt.1 .or. M.lt.1) Then
         If (DoPrint) Then
            Write(6,'(A,A,2I10)')
     &      SecNam,': C(n,M) has zero dimension: n,M=',n,M
         End If
         Return
      End If

      l=n*M
      If (l_C.lt.l) Then
         Call WarningMessage(2,SecNam//': l_C<l')
         Write(6,'(A,2I10)') 'l_C,l=',l_C,l
         Call LDF_Quit(1)
      End If

      ip=1
      l=n*LDF_nBasAux_Atom(A)
      NormA=dDot_(l,C(ip),1,C(ip),1)
      NormTot=NormA
      NormA=sqrt(NormA)
      ip=ip+l
      If (B.ne.A) Then
         l=n*LDF_nBasAux_Atom(B)
         NormB=dDot_(l,C(ip),1,C(ip),1)
         NormTot=NormTot+NormB
         NormB=sqrt(NormB)
         ip=ip+l
      Else
         NormB=0.0d0
      End If
      If (AP_2CFunctions(1,AB).gt.0) Then
         l=n*AP_2CFunctions(1,AB)
         NormAB=dDot_(l,C(ip),1,C(ip),1)
         NormTot=NormTot+NormAB
         NormAB=sqrt(NormAB)
         ip=ip+l
      Else
         NormAB=0.0d0
      End If
      NormTot=sqrt(NormTot)

      If (DoPrint) Then
         Write(6,'(A,I10,1P,4(1X,A,D20.10))')
     &   'AB=',AB,'||CA||=',NormA,'||CB||=',NormB,'||CAB||=',NormAB,
     &   '||CTot||=',NormTot
         Call xFlush(6)
      End If

      ! Check symmetry when A=B
      If (A.eq.B) Then
         nSA=LDF_nShell_Atom(A)
         ipA=LDF_lShell_Atom(A)-1
         l_iOff=nSA**2
         Call GetMem('CiOff','Allo','Inte',ip_iOff,l_iOff)
         Call LDF_uvOffset(AB,nSA,nSA,iWork(ip_iOff))
         nErr=0
         ip=0
         Do iCol=1,M
            mErr=0
            Do jSA=1,nSA
               jShell=iWork(ipA+jSA)
               Do jA=1,nBasSh(jShell)-1
                  Do iA=jA+1,nBasSh(jShell)
                     ijA=iOff(jSA,jSA)+nBasSh(jShell)*(jA-1)+iA
                     jiA=iOff(jSA,jSA)+nBasSh(jShell)*(iA-1)+jA
                     If (abs(C(ip+ijA)-C(ip+jiA)).gt.Tol) Then
                        mErr=mErr+1
                     End If
                  End Do
               End Do
               Do iSA=jSA+1,nSA
                  iShell=iWork(ipA+iSA)
                  Do jA=1,nBasSh(jShell)
                     Do iA=1,nBasSh(iShell)
                        ijA=iOff(iSA,jSA)+nBasSh(iShell)*(jA-1)+iA
                        jiA=iOff(jSA,iSA)+nBasSh(jShell)*(iA-1)+jA
                        If (abs(C(ip+ijA)-C(ip+jiA)).gt.Tol) Then
                           mErr=mErr+1
                        End If
                     End Do
                  End Do
               End Do
            End Do
            If (mErr.ne.0) Then
               Write(6,'(A,I10,A,I10,A,I10,A,1P,D20.10,A)')
     &         'AB=',AB,' J=',iCol,' #Err=',mErr,
     &         ': C(uv,J) != C(vu,J)  (Tol=',Tol,')'
               nErr=nErr+mErr
            End If
            ip=ip+n
         End Do
         If (nErr.gt.0) Then
            Call WarningMessage(2,SecNam//': C not symmetric')
            Call LDF_Quit(1)
         End If
         Call GetMem('CiOff','Free','Inte',ip_iOff,l_iOff)
      End If

      ! Check that LinDep cols are zero
      If (AP_1CLinDep(1,AB).gt.0) Then
         l_Col=LDF_nAuxShell_Atom(A)
         If (B.ne.A) Then
            l_Col=l_Col+LDF_nAuxShell_Atom(B)
         End If
         Call GetMem('CCol','Allo','Inte',ip_Col,l_Col)
         nASA=LDF_nAuxShell_Atom(A)
         ipA=LDF_lAuxShell_Atom(A)-1
         ip=ip_Col-1
         iCount=0
         Do kS=1,nASA
            iWork(ip+kS)=iCount
            iCount=iCount+n*nBasSh(iWork(ipA+kS))
         End Do
         If (B.ne.A) Then
            nASB=LDF_nAuxShell_Atom(B)
            ipB=LDF_lAuxShell_Atom(B)-1
            ip=ip+nASA
            Do kS=1,nASB
               iWork(ip+kS)=iCount
               iCount=iCount+n*nBasSh(iWork(ipB+kS))
            End Do
         Else
            nASB=nASA
            ipB=ipA
         End If
         If (iCount.gt.(l_C-n*AP_2CFunctions(1,AB))) Then
            Call WarningMessage(2,SecNam//': Logical error [1]')
            Call LDF_Quit(1)
         End If
         nErr=0
         ip=AP_1CLinDep(2,AB)-1
         Do i1CLinDep=1,AP_1CLinDep(1,AB)
            kAtom=iWork(ip+3*(i1CLinDep-1)+1)
            kS=iWork(ip+3*(i1CLinDep-1)+2)
            k=iWork(ip+3*(i1CLinDep-1)+3)
            If (kAtom.eq.A) Then
               If (kS.lt.1 .or. kS.gt.nASA) Then
                  Call WarningMessage(2,SecNam//': Logical error [3]')
                  Call LDF_Quit(1)
               End If
               kShell=iWork(ipA+kS)
               If (kShell.le.nShell_Valence .or.
     &             kShell.gt.(nShell_Valence+nShell_Auxiliary)) Then
                  Call WarningMessage(2,SecNam//': Logical error [5]')
                  Call LDF_Quit(1)
               End If
               If (k.lt.1 .or. k.gt.nBasSh(kShell)) Then
                  Call WarningMessage(2,SecNam//': Logical error [7]')
                  Call LDF_Quit(1)
               End If
               ipCol=iWork(ip_Col-1+kS)+n*(k-1)+1
               NormCol=sqrt(dDot_(n,C(ipCol),1,C(ipCol),1))
               If (NormCol.gt.Small) Then
                  nErr=nErr+1
                  Write(6,'(A,I10,1X,A,I10,1X,A,1P,D20.10)')
     &            'AB=',AB,'Norm of LinDep col',i1CLinDep,
     &            'is non-zero:',NormCol
                  Write(6,'(3X,A,3I10)')
     &            '- info about LinDep col: k,kS,kAtom=',k,kS,kAtom
               End If
            Else If (kAtom.eq.B) Then
               If (kS.lt.1 .or. kS.gt.nASB) Then
                  Call WarningMessage(2,SecNam//': Logical error [4]')
                  Call LDF_Quit(1)
               End If
               kShell=iWork(ipB+kS)
               If (kShell.le.nShell_Valence .or.
     &             kShell.gt.(nShell_Valence+nShell_Auxiliary)) Then
                  Call WarningMessage(2,SecNam//': Logical error [6]')
                  Call LDF_Quit(1)
               End If
               If (k.lt.1 .or. k.gt.nBasSh(kShell)) Then
                  Call WarningMessage(2,SecNam//': Logical error [8]')
                  Call LDF_Quit(1)
               End If
               ipCol=iWork(ip_Col-1+nASA+kS)+n*(k-1)+1
               NormCol=sqrt(dDot_(n,C(ipCol),1,C(ipCol),1))
               If (NormCol.gt.Small) Then
                  nErr=nErr+1
                  Write(6,'(A,I10,1X,A,I10,1X,A,1P,D20.10)')
     &            'AB=',AB,'Norm of LinDep col',i1CLinDep,
     &            'is non-zero:',NormCol
                  Write(6,'(3X,A,3I10)')
     &            '- info about LinDep col: k,kS,kAtom=',k,kS,kAtom
               End If
            Else
               Call WarningMessage(2,SecNam//': Logical error [2]')
               Call LDF_Quit(1)
            End If
         End Do
         Call GetMem('CCol','Free','Inte',ip_Col,l_Col)
         If (nErr.ne.0) Then
            Call WarningMessage(2,SecNam//': non-zero LinDep col found')
            Write(6,'(A,I10,A,2I10)') 'Atom pair=',AB,'  Atoms=',A,B
            Write(6,'(A,I10)') 'Number of non-zero LinDep cols=',nErr
            Write(6,'(A,I10)') 'Total number of LinDep cols=',
     &                         AP_1CLinDep(1,AB)
            Call LDF_Quit(1)
         End If
      End If

      ! Check 2C parts
      If (AP_2CFunctions(1,AB).gt.0) Then
         l_Map=AP_2CFunctions(1,AB)
         Call GetMem('2CMap','Allo','Inte',ip_Map,l_Map)
         Call LDF_Map2CF(AB,l_Map,1,iWork(ip_Map))
         l=AP_2CFunctions(1,AB)**2
         Call GetMem('2CBlock','Allo','Real',ip,l)
         iCol=LDF_nBasAux_Atom(A)
         If (B.ne.A) Then
            iCol=iCol+LDF_nBasAux_Atom(B)
         End If
         nErr=0
         Do j2CF=1,AP_2CFunctions(1,AB)
            ip0=ip-1+AP_2CFunctions(1,AB)*(j2CF-1)
            Do i2CF=1,AP_2CFunctions(1,AB)
               iRow=iWork(ip_Map-1+i2CF)
               Work(ip0+i2CF)=C(n*iCol+iRow)
            End Do
            iRow=iWork(ip_Map-1+j2CF)
            Tst=sqrt(dDot_(M,C(iRow),n,C(iRow),n))
            If (abs(Tst-1.0d0).gt.Tol) Then
               nErr=nErr+1
               Write(6,'(2X,A,I10,1X,A,I10,1X,A,1P,D20.10)')
     &         'AB=',AB,'2CFunction=',j2CF,'norm of C row=',Tst
               Write(6,'(2X,A)') 'C row:'
               Call Cho_Output(C,iRow,iRow,1,M,n,M,1,6)
            End If
            iCol=iCol+1
         End Do
         Do i2CF=1,AP_2CFunctions(1,AB)
            Work(ip-1+AP_2CFunctions(1,AB)*(i2CF-1)+i2CF)=
     &               Work(ip-1+AP_2CFunctions(1,AB)*(i2CF-1)+i2CF)-1.0d0
         End Do
         If (.not.isConstant(Work(ip),l,0.0d0,Tol)) Then
            Call WarningMessage(2,
     &                     SecNam//': 2C subblock minus unit not zero!')
            Call Cho_Head(SecNam//': 2C subblock minus unit','-',80,6)
            Call Cho_Output(Work(ip),1,AP_2CFunctions(1,AB),
     &                      1,AP_2CFunctions(1,AB),AP_2CFunctions(1,AB),
     &                      AP_2CFunctions(1,AB),1,6)
            Write(6,'(/,2X,A,1P,D20.10)')
     &      'RMS=',sqrt(dDot_(l,Work(ip),1,Work(ip),1)/dble(l))
            If (nErr.ne.0) Then
               Call WarningMessage(2,
     &                              SecNam//': non-unit norm of 2C row')
               Write(6,'(A,I10,A,2I10)') 'Atom pair=',AB,'  Atoms=',A,B
               Write(6,'(A,I10)') 'Number of non-unit 2C rows=',nErr
               Write(6,'(A,I10)') 'Total number of 2C rows=',
     &                            AP_2CFunctions(1,AB)
            End If
            Call LDF_Quit(1)
         End If
         If (nErr.ne.0) Then
            Call WarningMessage(2,SecNam//': non-unit norm of 2C row')
            Write(6,'(A,I10,A,2I10)') 'Atom pair=',AB,'  Atoms=',A,B
            Write(6,'(A,I10)') 'Number of non-unit 2C rows=',nErr
            Write(6,'(A,I10)') 'Total number of 2C rows=',
     &                         AP_2CFunctions(1,AB)
            Call LDF_Quit(1)
         End If
         Call GetMem('2CBlock','Free','Real',ip,l)
         Call GetMem('2CMap','Free','Inte',ip_Map,l_Map)
      End If

      End

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
      Subroutine LDF_ComputeValenceIntegrals(AB,CD,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: Compute valence integrals (u_A v_B | k_C l_D) for atom
C              pairs AB and CD.
C
      Implicit None
      Integer AB
      Integer CD
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_atom_pair_info.fh"

      Character*27 SecNam
      Parameter (SecNam='LDF_ComputeValenceIntegrals')

#if defined (_DEBUG_)
      Logical  isSymmetric
      External isSymmetric
      Integer  ji, lk
#endif

      External Int_LDF_SQ

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Integer A, B, C, D
      Integer nAB, nCD
      Integer l_xInt
      Integer nShell_A, nShell_B, nShell_C, nShell_D
      Integer ipA, ipB, ipC, ipD
      Integer maxAB, maxCD
      Integer n, n_
      Integer iS, jS, kS, lS
      Integer iShell, jShell, kShell, lShell
      Integer iS1, kS1, nij, nkl, nijkl
      Integer ip_iAB, l_iAB, ip_iCD, l_iCD
      Integer ip_SQ, l_SQ
      Integer ip_SewWrk, l_SewWrk
      Integer ij, ij0, kl, kl0, ijkl, ijkl0
      Integer iCol, iRow, ip, ip0
      Integer k, l

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer iAB, iCD
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      iAB(i,j)=iWork(ip_iAB-1+nShell_A*(j-1)+i)
      iCD(i,j)=iWork(ip_iCD-1+nShell_C*(j-1)+i)

      ! Get atoms of pairs
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)

      ! Get dimensions
      nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      nShell_A=LDF_nShell_Atom(A)
      nShell_B=LDF_nShell_Atom(B)
      nShell_C=LDF_nShell_Atom(C)
      nShell_D=LDF_nShell_Atom(D)
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1
      ipC=LDF_lShell_Atom(C)-1
      ipD=LDF_lShell_Atom(D)-1

      ! Check integral dimension
      l_xInt=nAB*nCD
      If (l_xInt.gt.l_xInt_) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Initialize integral array
      Call Cho_dZero(xInt,l_xInt)

      ! Allocate and set index array to shell rows and cols of integrals
      l_iAB=nShell_A*nShell_B
      Call GetMem('CVIiAB','Allo','Inte',ip_iAB,l_iAB)
      maxAB=0
      n=0
      Do jS=1,nShell_B
         jShell=iWork(ipB+jS)
         Do iS=1,nShell_A
            iShell=iWork(ipA+iS)
            iWork(ip_iAB-1+nShell_A*(jS-1)+iS)=n
            n_=nBasSh(iShell)*nBasSh(jShell)
            n=n+n_
            maxAB=max(maxAB,n_)
         End Do
      End Do
#if defined (_DEBUG_)
      If (n.ne.nAB) Then
         Call WarningMessage(2,SecNam//': n != nAB')
         Call LDF_Quit(1)
      End If
#endif
      If (CD.eq.AB) Then
         l_iCD=0
         ip_iCD=ip_iAB
         maxCD=maxAB
      Else
         l_iCD=nShell_C*nShell_D
         Call GetMem('CVIiCD','Allo','Inte',ip_iCD,l_iCD)
         maxCD=0
         n=0
         Do jS=1,nShell_D
            jShell=iWork(ipD+jS)
            Do iS=1,nShell_C
               iShell=iWork(ipC+iS)
               iWork(ip_iCD-1+nShell_C*(jS-1)+iS)=n
               n_=nBasSh(iShell)*nBasSh(jShell)
               n=n+n_
               maxCD=max(maxCD,n_)
            End Do
         End Do
      End If
#if defined (_DEBUG_)
      If (n.ne.nCD) Then
         Call WarningMessage(2,SecNam//': n != nCD')
         Call LDF_Quit(1)
      End If
#endif

      ! Allocate memory for largest shell quadruple
      l_SQ=maxAB*maxCD
      Call GetMem('CVISQ','Allo','Real',ip_SQ,l_SQ)

      ! Allocate memory for Seward
      Call GetMem('Max','Max ','Real',ip_SewWrk,l_SewWrk)
      Call xSetMem_Ints(l_SewWrk)

      ! Compute integrals
      Do lS=1,nShell_D
         lShell=iWork(ipD+lS)
         If (C.eq.D) Then
            kS1=lS
         Else
            kS1=1
         End If
         Do kS=kS1,nShell_C
            kShell=iWork(ipC+kS)
            nkl=nBasSh(kShell)*nBasSh(lShell)
            Do jS=1,nShell_B
               jShell=iWork(ipB+jS)
               If (A.eq.B) Then
                  iS1=jS
               Else
                  iS1=1
               End If
               Do iS=iS1,nShell_A
                  iShell=iWork(ipA+iS)
                  nij=nBasSh(iShell)*nBasSh(jShell)
                  nijkl=nij*nkl
                  Call Cho_dZero(Work(ip_SQ),nijkl)
                  SHA=iShell
                  SHB=jShell
                  SHC=kShell
                  SHD=lShell
                  Call Eval_IJKL(iShell,jShell,kShell,lShell,
     &                           Work(ip_SQ),nijkl,
     &                           Int_LDF_SQ)
                  Do l=1,nBasSh(lShell)
                     kl0=nBasSh(kShell)*(l-1)
                     Do k=1,nBasSh(kShell)
                        kl=kl0+k
                        iCol=iCD(kS,lS)+kl
                        ip0=nAB*(iCol-1)
                        ijkl0=ip_SQ-1+nij*(kl-1)
                        Do j=1,nBasSh(jShell)
                           ij0=nBasSh(iShell)*(j-1)
                           Do i=1,nBasSh(iShell)
                              ij=ij0+i
                              iRow=iAB(iS,jS)+ij
                              ip=ip0+iRow
                              ijkl=ijkl0+ij
                              xInt(ip)=Work(ijkl)
                              If (A.eq.B) Then
                                 iRow=iAB(jS,iS)
     &                               +nBasSh(jShell)*(i-1)+j
                                 ip=ip0+iRow
                                 xInt(ip)=Work(ijkl)
                              End If
                           End Do
                        End Do
                        If (D.eq.C) Then
                           iCol=iCD(lS,kS)
     &                         +nBasSh(lShell)*(k-1)+l
                           ip0=nAB*(iCol-1)
                           Do j=1,nBasSh(jShell)
                              ij0=nBasSh(iShell)*(j-1)
                              Do i=1,nBasSh(iShell)
                                 ij=ij0+i
                                 iRow=iAB(iS,jS)+ij
                                 ip=ip0+iRow
                                 ijkl=ijkl0+ij
                                 xInt(ip)=Work(ijkl)
                                 If (A.eq.B) Then
                                    iRow=iAB(jS,iS)
     &                                  +nBasSh(jShell)*(i-1)+j
                                    ip=ip0+iRow
                                    xInt(ip)=Work(ijkl)
                                 End If
                              End Do
                           End Do
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End Do

#if defined (_DEBUG_)
      If (AB.eq.CD) Then
         If (.not. isSymmetric(xInt,nAB,1.0d-14)) Then
            Call WarningMessage(2,SecNam//': (AB|CD) != (CD|AB)')
            Call LDF_Quit(1)
         End If
      End If
      If (A.eq.B) Then
         Do kl=1,nCD
            Do jS=1,nShell_B
               jShell=iWork(ipB+jS)
               iS=jS
               iShell=jShell
               Do j=1,nBasSh(jShell)
                  Do i=j+1,nBasSh(iShell)
                     ij=nAB*(kl-1)+iAB(iS,jS)+nBasSh(iShell)*(j-1)+i
                     ji=nAB*(kl-1)+iAB(jS,iS)+nBasSh(jShell)*(i-1)+j
                     If (abs(xInt(ij)-xInt(ji)).gt.1.0d-14) Then
                        Call WarningMessage(2,
     &                               SecNam//': [1] (AB|CD) != (BA|CD)')
                        Call LDF_Quit(1)
                     End If
                  End Do
               End Do
               Do iS=jS+1,nShell_A
                  iShell=iWork(ipA+iS)
                  Do j=1,nBasSh(jShell)
                     Do i=1,nBasSh(iShell)
                        ij=nAB*(kl-1)+iAB(iS,jS)+nBasSh(iShell)*(j-1)+i
                        ji=nAB*(kl-1)+iAB(jS,iS)+nBasSh(jShell)*(i-1)+j
                        If (abs(xint(ij)-xint(ji)).gt.1.0d-14) Then
                           Call WarningMessage(2,
     &                               SecNam//': [2] (AB|CD) != (BA|CD)')
                           Call LDF_Quit(1)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
      If (C.eq.D) Then
         Do ij=1,nAB
            Do lS=1,nShell_D
               lShell=iWork(ipD+lS)
               kS=lS
               kShell=lShell
               Do l=1,nBasSh(lShell)
                  do k=l+1,nBasSh(kShell)
                     kl=iCD(kS,lS)+nBasSh(kShell)*(l-1)+k
                     lk=iCD(lS,kS)+nBasSh(lShell)*(k-1)+l
                     kl=nAB*(kl-1)+ij
                     lk=nAB*(lk-1)+ij
                     If (abs(xInt(kl)-xInt(lk)).gt.1.0d-14) Then
                        Call WarningMessage(2,
     &                               SecNam//': [1] (AB|CD) != (AB|DC)')
                        Call LDF_Quit(1)
                     End If
                  End Do
               End Do
               Do kS=lS+1,nShell_C
                  kShell=iWork(ipC+kS)
                  Do l=1,nBasSh(lShell)
                     Do k=1,nBasSh(kShell)
                        kl=iCD(kS,lS)+nBasSh(kShell)*(l-1)+k
                        lk=iCD(lS,kS)+nBasSh(lShell)*(k-1)+l
                        kl=nAB*(kl-1)+ij
                        lk=nAB*(lk-1)+ij
                        If (abs(xInt(kl)-xInt(lk)).gt.1.0d-14) Then
                           Call WarningMessage(2,
     &                               SecNam//': [2] (AB|CD) != (AB|DC)')
                           Call LDF_Quit(1)
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      End If
#endif

      ! Deallocation
      Call xRlsMem_Ints()
      Call GetMem('CVISQ','Free','Real',ip_SQ,l_SQ)
      If (l_iCD.gt.0) Then
         Call GetMem('CVIiCD','Free','Inte',ip_iCD,l_iCD)
      End If
      Call GetMem('CVIiAB','Free','Inte',ip_iAB,l_iAB)

      End

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
      Subroutine LDF_SetIntegralPrescreeningInfo()
C
C     Thomas Bondo Pedersen, November 2010.
C
C     Purpose: Allocate and Set LDF integral prescreening info.
C              To deallocate, call LDF_UnsetIntegralPrescreeningInfo().
C
C     LDF info must be properly set up before calling this routine.
C
      Implicit None
#include "WrkSpc.fh"
#include "localdf.fh"
#include "localdf_int.fh"
#include "localdf_bas.fh"
#include "ldf_integral_prescreening_info.fh"
#include "ldf_atom_pair_info.fh"

#if defined (_DEBUG_)
      Character*31 SecNam
      Parameter (SecNam='LDF_SetIntegralPrescreeningInfo')
#endif

      Integer  LDF_nAtom, LDF_nAuxShell_Atom, LDF_nShell_Atom
      Integer  LDF_lShell_Atom, LDF_nBas_Atom
      External LDF_nAtom, LDF_nAuxShell_Atom, LDF_nShell_Atom
      External LDF_lShell_Atom, LDF_nBas_Atom

      Integer LDF_GlobalToAtomicShell

      Character*8 Label

      Integer nAtom
      Integer A, B, AB
      Integer ip, l
      Integer ip_myOffset, l_myOffset, l_myOffset_1
      Integer nnShl
      Integer iS, jS, ijS
      Integer iShell, jShell
      Integer ii, jj, ij0, ij
      Integer n1, n2
      Integer nSA, nSB, ipA, ipB
      Integer ip_D, l_D
      Integer ipD

      Real*8  Tmax, Tsum
      Real*8  CutInt_, CutInt

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
      Integer AP_2CFunctions
      Integer IndxG2
      Integer myOffset
      Integer iTri
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)
      myOffset(i,j)=iWork(ip_myOffset-1+l_myOffset_1*(j-1)+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      ! Get number of atoms
      nAtom=LDF_nAtom()

      ! Allocate GDiag 1C arrays
      l_GDiag_1C=2*nAtom
      Call GetMem('GD1C','Allo','Inte',ip_GDiag_1C,l_GDiag_1C)
      Do A=1,nAtom
         l=LDF_nAuxShell_Atom(A)
         If (l.gt.0) Then
            Write(Label,'(A,I5.5)') 'GD1',A-1
            Call GetMem(Label,'Allo','Real',ip,l)
            iWork(ip_GDiag_1C+2*(A-1))=l
            iWork(ip_GDiag_1C+2*(A-1)+1)=ip
         Else
            iWork(ip_GDiag_1C+2*(A-1))=0
            iWork(ip_GDiag_1C+2*(A-1)+1)=0
         End If
      End Do
      l_GDiag_1C_Mx=nAtom
      Call GetMem('GD1CMx','Allo','Real',ip_GDiag_1C_Mx,l_GDiag_1C_Mx)
      l_GDiag_1C_Sm=nAtom
      Call GetMem('GD1CSm','Allo','Real',ip_GDiag_1C_Sm,l_GDiag_1C_Sm)

      ! Set GDiag 1C info
      Call LDF_GetCutInt(CutInt_)
      CutInt=1.0d-99
      Call LDF_SetCutInt(CutInt)
      Call GetMem('GetMax','Max ','Real',ip,l)
      Call xSetMem_Ints(l)
      Do A=1,nAtom
         l=iWork(ip_GDiag_1C+2*(A-1))
         If (l.gt.0) Then
            ip=iWork(ip_GDiag_1C+2*(A-1)+1)
            Call LDF_SIPI_G1C(A,l,Work(ip),Work(ip_GDiag_1C_Mx-1+A),
     &                                     Work(ip_GDiag_1C_Sm-1+A))
#if defined (_DEBUG_)
            If (Work(ip_GDiag_1C_Sm-1+A).lt.0.0d0) Then
               Call WarningMessage(2,SecNam//': [1] sum < 0')
               Call LDF_Quit(1)
            End If
#endif
            Work(ip_GDiag_1C_Sm-1+A)=sqrt(Work(ip_GDiag_1C_Sm-1+A))
         Else
            Work(ip_GDiag_1C_Mx-1+A)=0.0d0
            Work(ip_GDiag_1C_Sm-1+A)=0.0d0
         End If
      End Do
      Call xRlsMem_Ints()
      Call LDF_SetCutInt(CutInt_) ! restore CutInt

      ! Allocate and set GDiag 2C arrays
      If (LDF2) Then
         l_GDiag_2C_Mx=NumberOfAtompairs
         Call GetMem('GD2CMx','Allo','Real',ip_GDiag_2C_Mx,
     &                                       l_GDiag_2C_Mx)
         l_GDiag_2C_Sm=NumberOfAtompairs
         Call GetMem('GD2CSm','Allo','Real',ip_GDiag_2C_Sm,
     &                                       l_GDiag_2C_Sm)
         l_GDiag_2C=2*NumberOfAtomPairs
         Call GetMem('GD2C','Allo','Inte',ip_GDiag_2C,l_GDiag_2C)
         Do AB=1,NumberOfAtomPairs
            If (AP_Atoms(1,AB).eq.AP_Atoms(2,AB)) Then
               l_D=LDF_nBas_Atom(AP_Atoms(1,AB))
               l_D=l_D*(l_D+1)/2
               Call GetMem('DiaLT','Allo','Real',ip_D,l_D)
               Call LDF_Q2LT(AP_Atoms(1,AB),
     &                       Work(iWork(ip_AP_DiagBak-1+AB)),Work(ip_D))
            Else
               l_D=0
               ip_D=iWork(ip_AP_DiagBak-1+AB)
            End If
            If (AP_2CFunctions(1,AB).gt.0) Then
               Call LDF_SetIndxG(AB)
               nnShl=l_2CList_2
               l=nnShl
               Write(Label,'(A,I5.5)') 'GD2',AB-1
               Call GetMem(Label,'Allo','Real',ip,l)
               iWork(ip_GDiag_2C+2*(AB-1))=l
               iWork(ip_GDiag_2C+2*(AB-1)+1)=ip
               A=AP_Atoms(1,AB)
               B=AP_Atoms(2,AB)
               nSA=LDF_nShell_Atom(A)
               ipA=LDF_lShell_Atom(A)-1
               nSB=LDF_nShell_Atom(B)
               ipB=LDF_lShell_Atom(B)-1
               l_myOffset_1=nSA
               l_myOffset=l_myOffset_1*nSB
               Call GetMem('myOffset','Allo','Inte',ip_myOffset,
     &                                               l_myOffset)
               Tsum=0.0d0
               If (A.eq.B) Then
                  ii=0
                  Do iS=1,l_myOffset_1
                     n1=nBasSh(iWork(ipA+iS))
                     Do jS=1,iS-1
                        iWork(ip_myOffset-1+l_myOffset_1*(jS-1)+iS)=ii
                        iWork(ip_myOffset-1+l_myOffset_1*(iS-1)+jS)=ii
                        ii=ii+n1*nBasSh(iWork(ipA+jS))
                     End Do
                     iWork(ip_myOffset-1+l_myOffset_1*(iS-1)+iS)=ii
                     ii=ii+n1*(n1+1)/2
                  End Do
                  Do ijS=0,nnShl-1
                     iShell=iWork(ip_2CList+3*ijS)
                     iS=LDF_GlobalToAtomicShell(A,iShell)
                     jShell=iWork(ip_2CList+3*ijS+1)
                     jS=LDF_GlobalToAtomicShell(B,jShell)
                     SPAB=iWork(ip_2CList+3*ijS+2)
                     Tmax=0.0d0
                     If (iS.eq.jS) Then
                        Do jj=1,nBasSh(jShell)
                           ij0=nBasSh(iShell)*(jj-1)
                           Do ii=1,nBasSh(iShell)
                              ij=ij0+ii
                              If (IndxG2(ij,SPAB).gt.0) Then
                                 ipD=ip_D-1+myOffset(iS,jS)
     &                              +iTri(ii,jj)
                                 Tmax=max(Tmax,Work(ipD))
                                 Tsum=Tsum+Work(ipD)
                              End If
                           End Do
                        End Do
                     Else If (iS.gt.jS) Then
                        Do jj=1,nBasSh(jShell)
                           ij0=nBasSh(iShell)*(jj-1)
                           Do ii=1,nBasSh(iShell)
                              ij=ij0+ii
                              If (IndxG2(ij,SPAB).gt.0) Then
                                 ipD=ip_D-1+myOffset(iS,jS)+ij
                                 Tmax=max(Tmax,Work(ipD))
                                 Tsum=Tsum+Work(ipD)
                              End If
                           End Do
                        End Do
                     Else ! iS < jS
                        Do jj=1,nBasSh(iShell)
                           ij0=nBasSh(jShell)*(jj-1)
                           Do ii=1,nBasSh(jShell)
                              ij=ij0+ii
                              If (IndxG2(ij,SPAB).gt.0) Then
                                 ipD=ip_D-1+myOffset(jS,iS)+ij
                                 Tmax=max(Tmax,Work(ipD))
                                 Tsum=Tsum+Work(ipD)
                              End If
                           End Do
                        End Do
                     End If
                     Work(ip+ijS)=Tmax
                  End Do
               Else
                  ii=0
                  Do jS=1,nSB
                     n2=nBasSh(iWork(ipB+jS))
                     Do iS=1,nSA
                        iWork(ip_myOffset-1+l_myOffset_1*(jS-1)+iS)=ii
                        ii=ii+nBasSh(iWork(ipA+iS))*n2
                     End Do
                  End Do
                  Do ijS=0,nnShl-1
                     iShell=iWork(ip_2CList+3*ijS)
                     iS=LDF_GlobalToAtomicShell(A,iShell)
                     jShell=iWork(ip_2CList+3*ijS+1)
                     jS=LDF_GlobalToAtomicShell(B,jShell)
                     SPAB=iWork(ip_2CList+3*ijS+2)
                     Tmax=0.0d0
                     Do jj=1,nBasSh(jShell)
                        ij0=nBasSh(iShell)*(jj-1)
                        Do ii=1,nBasSh(iShell)
                           ij=ij0+ii
                           If (IndxG2(ij,SPAB).gt.0) Then
                              ipD=ip_D-1+myOffset(iS,jS)+ij
                              Tmax=max(Tmax,Work(ipD))
                              Tsum=Tsum+Work(ipD)
                           End If
                        End Do
                     End Do
                     Work(ip+ijS)=Tmax
                  End Do
               End If
               Call GetMem('myOffset','Free','Inte',ip_myOffset,
     &                                               l_myOffset)
               l_myOffset=0
               l_myOffset_1=0
               Call LDF_UnsetIndxG()
               Tmax=Work(ip)
               Do ijS=1,nnShl-1
                  Tmax=max(Tmax,Work(ip+ijS))
               End Do
               Work(ip_GDiag_2C_Mx-1+AB)=Tmax
#if defined (_DEBUG_)
               If (Tsum.lt.0.0d0) Then
                  Call WarningMessage(2,SecNam//': [2] sum < 0')
                  Call LDF_Quit(1)
               End If
#endif
               Work(ip_GDiag_2C_Sm-1+AB)=sqrt(Tsum)
            Else
               iWork(ip_GDiag_2C+2*(AB-1))=0
               iWork(ip_GDiag_2C+2*(AB-1)+1)=0
               Work(ip_GDiag_2C_Mx-1+AB)=0.0d0
               Work(ip_GDiag_2C_Sm-1+AB)=0.0d0
            End If
            If (l_D.gt.0) Then
               Call GetMem('DiaLT','Free','Real',ip_D,l_D)
            End If
         End Do
      Else
         ip_GDiag_2C=0
         l_GDiag_2C=0
         ip_GDiag_2C_Mx=0
         l_GDiag_2C_Mx=0
         ip_GDiag_2C_Sm=0
         l_GDiag_2C_Sm=0
      End If

      ! Allocate and set IDiag arrays
      l_IDiag_Mx=NumberOfAtomPairs
      Call GetMem('IDiag_Mx','Allo','Real',ip_IDiag_Mx,l_IDiag_Mx)
      l_IDiag_Sm=NumberOfAtomPairs
      Call GetMem('IDiag_Sm','Allo','Real',ip_IDiag_Sm,l_IDiag_Sm)
      l_IDiag=2*NumberOfAtomPairs
      Call GetMem('IDiag','Allo','Inte',ip_IDiag,l_IDiag)
      Do AB=1,NumberOfAtomPairs
         l=LDF_nShell_Atom(AP_Atoms(1,AB))
     &    *LDF_nShell_Atom(AP_Atoms(2,AB))
         If (l.gt.0) Then
            Write(Label,'(A,I5.5)') 'IDI',AB-1
            Call GetMem(Label,'Allo','Real',ip,l)
            iWork(ip_IDiag+2*(AB-1))=l
            iWork(ip_IDiag+2*(AB-1)+1)=ip
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            If (A.eq.B) Then
               l_D=LDF_nBas_Atom(A)
               l_D=l_D*(l_D+1)/2
               Call GetMem('DiaLT','Allo','Real',ip_D,l_D)
               Call LDF_Q2LT(AP_Atoms(1,AB),
     &                       Work(iWork(ip_AP_DiagBak-1+AB)),Work(ip_D))
            Else
               l_D=0
               ip_D=iWork(ip_AP_DiagBak-1+AB)
            End If
            nSA=LDF_nShell_Atom(A)
            ipA=LDF_lShell_Atom(A)-1
            nSB=LDF_nShell_Atom(B)
            ipB=LDF_lShell_Atom(B)-1
            Tsum=0.0d0
            If (A.eq.B) Then
               ipD=ip_D-1
               Do iS=1,nSA
                  Do jS=1,iS-1
                     Tmax=0.0d0
                     n1=nBasSh(iWork(ipA+iS))*nBasSh(iWork(ipB+jS))
                     Do ij=1,n1
                        Tmax=max(Tmax,Work(ipD+ij))
                        Tsum=Tsum+Work(ipD+ij)
                     End Do
                     ipD=ipD+n1
                     Work(ip-1+nSA*(jS-1)+iS)=Tmax
                     Work(ip-1+nSA*(iS-1)+jS)=Tmax
                  End Do
                  Tmax=0.0d0
                  n1=nBasSh(iWork(ipA+iS))*(nBasSh(iWork(ipA+iS))+1)/2
                  Do ij=1,n1
                     Tmax=max(Tmax,Work(ipD+ij))
                     Tsum=Tsum+Work(ipD+ij)
                  End Do
                  ipD=ipD+n1
                  Work(ip-1+nSA*(iS-1)+iS)=Tmax
               End Do
            Else
               ipD=ip_D-1
               Do jS=1,nSB
                  Do iS=1,nSA
                     Tmax=0.0d0
                     n1=nBasSh(iWork(ipA+iS))*nBasSh(iWork(ipB+jS))
                     Do ij=1,n1
                        Tmax=max(Tmax,Work(ipD+ij))
                        Tsum=Tsum+Work(ipD+ij)
                     End Do
                     ipD=ipD+n1
                     Work(ip-1+nSA*(jS-1)+iS)=Tmax
                  End Do
               End Do
            End If
            Tmax=Work(ip)
            Do ijS=1,nSA*nSB-1
               Tmax=max(Tmax,Work(ip+ijS))
            End Do
            Work(ip_IDiag_Mx-1+AB)=Tmax
#if defined (_DEBUG_)
            If (Tsum.lt.0.0d0) Then
               Call WarningMessage(2,SecNam//': [3] sum < 0')
               Call LDF_Quit(1)
            End If
#endif
            Work(ip_IDiag_Sm-1+AB)=sqrt(Tsum)
            If (l_D.gt.0) Then
               Call GetMem('DiaLT','Free','Real',ip_D,l_D)
            End If
         Else
            iWork(ip_IDiag+2*(AB-1))=0
            iWork(ip_IDiag+2*(AB-1)+1)=0
            Work(ip_IDiag_Mx-1+AB)=0.0d0
            Work(ip_IDiag_Sm-1+AB)=0.0d0
         End If
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Integer Function LDF_GlobalToAtomicShell(A,iShell)
      Implicit None
      Integer A, iShell
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nShell_Atom, LDF_lShell_Atom

      Integer iS, iS_, nS, ip

      nS=LDF_nShell_Atom(A)
      ip=LDF_lShell_Atom(A)-1
      iS=0
      iS_=0
      Do While (iS_.lt.nS)
         iS_=iS_+1
         If (iWork(ip+iS_).eq.iShell) Then
            iS=iS_
            iS_=nS+1 ! break loop
         End If
      End Do
      If (iS.eq.0) Then
         Call WarningMessage(2,
     &   'LDF_GlobalToAtomicShell: shell not found!')
         Call LDF_Quit(1)
      End If

      LDF_GlobalToAtomicShell=iS

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SIPI_G1C(A,l,Gmax_S,Gmax,Gsum)
      Implicit None
      Integer A
      Integer l
      Real*8  Gmax_S(l)
      Real*8  Gmax
      Real*8  Gsum
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int3.fh"

      External Int_LDF_Gmax_S

      Integer  LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Integer l_Res
      Parameter (l_Res=2)
      Real*8 Res(l_Res)

      Integer nAuxShellA, ipA
      Integer iS, iShell, dShell

      nAuxShellA=LDF_nAuxShell_Atom(A)
      If (l.ne.nAuxShellA) Then
         Call WarningMessage(2,'LDF_SIPI_G1C: dimension error!')
         Call LDF_Quit(1)
      End If
      ipA=LDF_lAuxShell_Atom(A)-1

      Gmax=0.0d0
      Gsum=0.0d0

      dShell=nShell_Valence+nShell_Auxiliary+1
      SHA=dShell
      SHC=dShell
      Do iS=1,nAuxShellA
         iShell=iWork(ipA+iS)
         SHB=iShell
         SHD=iShell
         Call Eval_IJKL(dShell,iShell,dShell,iShell,Res,l_Res,
     &                  Int_LDF_Gmax_S)
         Gmax_S(iS)=Res(1)
         Gmax=max(Gmax,Res(1))
         Gsum=Gsum+Res(2)
      End Do

      SHA=0
      SHB=0
      SHC=0
      SHD=0

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetIntegralPrescreeningInfo()
C
C     Thomas Bondo Pedersen, November 2010.
C
C     Purpose: Deallocate LDF integral prescreening info.
C
      Implicit None
#include "WrkSpc.fh"
#include "ldf_integral_prescreening_info.fh"

      Integer  LDF_nAtomPair, LDF_nAtom
      External LDF_nAtomPair, LDF_nAtom

      Character*8 Label

      Integer A, AB
      Integer ip, l

      If (l_GDiag_1C.gt.0) Then
         Do A=1,LDF_nAtom()
            l=iWork(ip_GDiag_1C+2*(A-1))
            If (l.gt.0) Then
               ip=iWork(ip_GDiag_1C+2*(A-1)+1)
               Write(Label,'(A,I5.5)') 'GD1',A-1
               Call GetMem(Label,'Free','Real',ip,l)
            End If
         End Do
         Call GetMem('GD1C','Free','Inte',ip_GDiag_1C,l_GDiag_1C)
         ip_GDiag_1C=0
         l_GDiag_1C=0
      End If
      If (l_GDiag_1C_Mx.gt.0) Then
         Call GetMem('GD1CMx','Free','Real',ip_GDiag_1C_Mx,
     &                                       l_GDiag_1C_Mx)
         ip_GDiag_1C_Mx=0
         l_GDiag_1C_Mx=0
      End If
      If (l_GDiag_1C_Sm.gt.0) Then
         Call GetMem('GD1CSm','Free','Real',ip_GDiag_1C_Sm,
     &                                       l_GDiag_1C_Sm)
         ip_GDiag_1C_Sm=0
         l_GDiag_1C_Sm=0
      End If

      If (l_GDiag_2C.gt.0) Then
         Do AB=1,LDF_nAtomPair()
            l=iWork(ip_GDiag_2C+2*(AB-1))
            If (l.gt.0) Then
               ip=iWork(ip_GDiag_2C+2*(AB-1)+1)
               Write(Label,'(A,I5.5)') 'GD2',AB-1
               Call GetMem(Label,'Free','Real',ip,l)
            End If
         End Do
         Call GetMem('GD2C','Free','Inte',ip_GDiag_2C,l_GDiag_2C)
         ip_GDiag_2C=0
         l_GDiag_2C=0
      End If
      If (l_GDiag_2C_Mx.gt.0) Then
         Call GetMem('GD2CMx','Free','Real',ip_GDiag_2C_Mx,
     &                                       l_GDiag_2C_Mx)
         ip_GDiag_2C_Mx=0
         l_GDiag_2C_Mx=0
      End If
      If (l_GDiag_2C_Sm.gt.0) Then
         Call GetMem('GD2CSm','Free','Real',ip_GDiag_2C_Sm,
     &                                       l_GDiag_2C_Sm)
         ip_GDiag_2C_Sm=0
         l_GDiag_2C_Sm=0
      End If

      If (l_IDiag.gt.0) Then
         Do AB=1,LDF_nAtomPair()
            l=iWork(ip_IDiag+2*(AB-1))
            If (l.gt.0) Then
               ip=iWork(ip_IDiag+2*(AB-1)+1)
               Write(Label,'(A,I5.5)') 'IDI',AB-1
               Call GetMem(Label,'Free','Real',ip,l)
            End If
         End Do
         Call GetMem('IDiag','Free','Inte',ip_IDiag,l_IDiag)
         ip_IDiag=0
         l_IDiag=0
      End If
      If (l_IDiag_Mx.gt.0) Then
         Call GetMem('IDiag_Mx','Free','Real',ip_IDiag_Mx,l_IDiag_Mx)
         ip_IDiag_Mx=0
         l_IDiag_Mx=0
      End If
      If (l_IDiag_Sm.gt.0) Then
         Call GetMem('IDiag_Sm','Free','Real',ip_IDiag_Sm,l_IDiag_Sm)
         ip_IDiag_Sm=0
         l_IDiag_Sm=0
      End If

      End

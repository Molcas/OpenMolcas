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
      Subroutine LDF_ComputeOverlapBlock(AB,l_S,S)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute overlap integral S(uv) for uv on AB.
C              It is assumed that the proper initialization of one-el
C              integral data is set before calling this routine
C              (by calling LDF_SetOneEl with argument 'Mltpl  0').
C
      use iSD_data
      Implicit None
      Integer AB
      Integer l_S
      Real*8  S(l_S)
      Integer iPrim, jPrim, iAng, jAng, nOrder, lFinal, lScrt1,
     &        lScrt2, MemKrn, ixyz, MemKer, mFinal, mScrt1, mScrt2,
     &        nElem
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "localdf_bas.fh"
#include "rmat_option.fh"
#include "ldf_oneel.fh"
#include "nsd.fh"
#include "setup.fh"
#include "ldf_atom_pair_info.fh"
#include "property_label.fh"

      Real*8, Dimension(:), Allocatable :: Final, Scrtch,
     &                                     ScrSph, Kern
      Character*23 SecNam
      Parameter (SecNam='LDF_ComputeOverlapBlock')

      External MltInt, MltMem

      Integer  LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom
      External LDF_nBas_Atom, LDF_nShell_Atom, LDF_lShell_Atom

      Character*8 Label

      Logical Do_PGamma

      Integer A, B
      Integer nSA, nSB
      Integer iSA, iSB
      Integer iShellA, iShellB
      Integer ipA, ipB
      Integer lS
      Integer iPrint
      Integer nGrid
      Integer iAddPot
      Integer nOrdOp
      Integer ip_SOInt, l_SOInt
      Integer ip, l
      Integer iBas, iCmp, iAO
      Integer jBas, jCmp, jAO

      Real*8  PtChrg(1)

      Integer i, j
      Integer nBasSh
      Integer AP_Atoms
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)
*                                                                      *
************************************************************************
*                                                                      *
      ! Check operator label
      If (OperatorLabel.ne.'Mltpl  0') Then
         Call WarningMessage(2,SecNam//': illegal operator label')
         Write(6,'(A,A)') 'OperatorLabel=',OperatorLabel
         Call LDF_Quit(1)
      End If

      ! Get atoms
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      nSA=LDF_nShell_Atom(A)
      nSB=LDF_nShell_Atom(B)
      ipA=LDF_lShell_Atom(A)-1
      ipB=LDF_lShell_Atom(B)-1

      ! Check dimension
      lS=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      If (lS.lt.1) Return
      If (lS.gt.l_S) Then
         Call WarningMessage(2,
     &                SecNam//': insufficient integral array dimension')
         Call LDF_Quit(1)
      End If

      ! Set arguments for call to OneEl_IJ
      RMat_type_integrals=.False. ! stored in common block...
      Do_PGamma=.True.
      Label=OperatorLabel
      PLabel=' '
      iPrint=0
      nGrid=1
      PtChrg(1)=0.0d0
      iAddPot=0
      nOrdOp=0

      ! Allocate temporary integral array
      l_SOInt=0
      Do iSB=1,nSB
         iShellB=iWork(ipB+iSB)
         Do iSA=1,nSA
            iShellA=iWork(ipA+iSA)
            l_SOInt=max(l_SOInt,nBasSh(iShellA)*nBasSh(iShellB))
         End Do
      End Do
      Call GetMem('SBlock','Allo','Real',ip_SOInt,l_SOInt)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate scratch for the evaluation of the integrals             *
*                                                                      *
      lFinal=1
      lScrt1=1
      lScrt2=1
      MemKrn=1
      Do iSB=1,nSB
         iShellB=iWork(ipB+iSB)
         Do iSA=1,nSA
            iShellA=iWork(ipA+iSA)
*
            iPrim=iSD(5,iShellA)
            jPrim=iSD(5,iShellB)
            iBas=iSD(3,iShellA)
            jBas=iSD(3,iShellB)
            iAng=iSD(1,iShellA)
            jAng=iSD(1,iShellB)
*
            mFinal=nIC*iPrim*jPrim*nElem(iAng)*nElem(jAng)
            lFinal=Max(lFinal,mFinal)
*
            mScrt1=nIC*Max(iPrim,jBas)*Max(iBas,jPrim)
     &         *nElem(iAng)*nElem(jAng)
            lScrt1=Max(mScrt1,lScrt1)
*
            mScrt2=nIC*iBas*jBas*nElem(iAng)*nElem(jAng)
            lScrt2=Max(mScrt2,lScrt2)
*
            Call MltMem(nOrder,MemKer,iAng,jAng,nOrdOp)
            MemKrn=Max(MemKer*iPrim*jPrim,MemKrn)
         End Do
      End Do
*
      Call mma_allocate(Final,lFinal,label='Final')
      Call mma_allocate(Scrtch,lScrt1,label='Scrtch')
      Call mma_allocate(ScrSph,lScrt2,label='ScrSph')
      Call mma_allocate(Kern,MemKrn,label='Kern')
*                                                                      *
************************************************************************
*                                                                      *
      ! Compute integrals
      Call Cho_dZero(S,lS)
      ip=1
      Do iSB=1,nSB
         iShellB=iWork(ipB+iSB)
         Do iSA=1,nSA
            iShellA=iWork(ipA+iSA)
            l=nBasSh(iShellA)*nBasSh(iShellB)
            Call Cho_dZero(Work(ip_SOInt),l)
            Call OneEl_IJ(iShellA,iShellB,iPrint,Do_PGamma,
     &                    Work(ip_xZeta),Work(ip_xZI),Work(ip_xKappa),
     &                    Work(ip_xPCoor),
     &                    MltInt,MltMem,Label,iWork(ip_lOper),nComp,
     &                    Work(ip_CCoor),
     &                    nOrdOp,iWork(ip_kOper),
     &                    iStabO,nStabO,nIC,
     &                    PtChrg,nGrid,iAddPot,Work(ip_SOInt),l,
     &                    Final,lFinal,Scrtch,lScrt1,
     &                    ScrSph,lScrt2,Kern,MemKrn)

            iCmp=iSD(2,iShellA)
            iBas=iSD(3,iShellA)
            iAO=iSD(7,iShellA)
            jCmp=iSD(2,iShellB)
            jBas=iSD(3,iShellB)
            jAO=iSD(7,iShellB)
            Call LDF_SortOverlapBlock(Work(ip_SOInt),iBas,jBas,
     &                                iCmp,jCmp,iAO,jAO,S(ip),l)
            ip=ip+l
         End Do
      End Do
*
      ! Deallocation
      Call mma_deallocate(Final)
      Call mma_deallocate(Scrtch)
      Call mma_deallocate(ScrSph)
      Call mma_deallocate(Kern)
*
      Call GetMem('SBlock','Free','Real',ip_SOInt,l_SOInt)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SortOverlapBlock(SOInt,iBas,jBas,iCmp,jCmp,iAO,jAO,
     &                                S,lS)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: extract block of overlap matrix.
C
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (a-h,o-z)
      Integer iBas, jBas
      Real*8  SOInt(iBas*jBas,*)
      Integer iCmp, jCmp
      Integer iAO, jAO
      Integer lS
      Real*8  S(lS)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

#if defined (_DEBUG_)
      Character*20 SecNam
      Parameter (SecNam='LDF_SortOverlapBlock')
#endif

      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      iSOShl(i)=iWork(ip_iSOShl-1+i)

      iShell=iSOShl(iAOtSO(iAO+1,0))
      jShell=iSOShl(iAOtSO(jAO+1,0))
#if defined (_DEBUG_)
      If (nBasSh(iShell)*nBasSh(jShell).gt.lS) Then
         Call WarningMessage(2,SecNam//': array dimension problem')
         Call LDF_Quit(1)
      End If
#endif
      lSO=0
      If (iShell.eq.jShell) Then
         Do iComp=1,iCmp
            iSO0=iAOtSO(iAO+iComp,0)-1
            Do jComp=1,iComp-1
               jSO0=iAOtSO(jAO+jComp,0)-1
               lSO=lSO+1
               Do jj=1,jBas
                  jSOj=iShlSO(jSO0+jj)
#if defined (_DEBUG_)
                  If (jSOj.lt.1 .or. jSOj.gt.nBasSh(jShell)) Then
                     Call WarningMessage(2,
     &                                    SecNam//' jSOj out of bounds')
                     Call LDF_Quit(1)
                  End If
#endif
                  i_j_0=nBasSh(iShell)*(jSOj-1)
                  ij0=iBas*(jj-1)
                  Do ii=1,iBas
                     iSOi=iShlSO(iSO0+ii)
#if defined (_DEBUG_)
                     If (iSOi.lt.1 .or. iSOi.gt.nBasSh(iShell)) Then
                        Call WarningMessage(2,
     &                                    SecNam//' iSOi out of bounds')
                        Call LDF_Quit(1)
                     End If
                     If ((i_j_0+iSOi).lt.1 .or. (i_j_0+iSOi).gt.lS) Then
                        Call WarningMessage(2,
     &                                     SecNam//' i_j out of bounds')
                        Call LDF_Quit(1)
                     End If
#endif
                     S(i_j_0+iSOi)=SOInt(ij0+ii,lSO)
                  End Do
               End Do
            End Do
            jComp=iComp
            jSO0=iAOtSO(jAO+jComp,0)-1
            lSO=lSO+1
            Do jj=1,jBas
               jSOj=iShlSO(jSO0+jj)
#if defined (_DEBUG_)
               If (jSOj.lt.1 .or. jSOj.gt.nBasSh(jShell)) Then
                  Call WarningMessage(2,SecNam//' jSOj out of bounds')
                  Call LDF_Quit(1)
               End If
#endif
               i_j_0=nBasSh(iShell)*(jSOj-1)
               ij0=iBas*(jj-1)
               Do ii=jj,iBas
                  iSOi=iShlSO(iSO0+ii)
#if defined (_DEBUG_)
                  If (iSOi.lt.1 .or. iSOi.gt.nBasSh(iShell)) Then
                     Call WarningMessage(2,
     &                                    SecNam//' iSOi out of bounds')
                     Call LDF_Quit(1)
                  End If
                  If ((i_j_0+iSOi).lt.1 .or. (i_j_0+iSOi).gt.lS) Then
                     Call WarningMessage(2,SecNam//' i_j out of bounds')
                     Call LDF_Quit(1)
                  End If
#endif
                  S(i_j_0+iSOi)=SOInt(ij0+ii,lSO)
                  S(nBasSh(jShell)*(iSOi-1)+jSOj)=SOInt(ij0+ii,lSO)
               End Do
            End Do
         End Do
      Else
         Do iComp=1,iCmp
            iSO0=iAOtSO(iAO+iComp,0)-1
            Do jComp=1,jCmp
               jSO0=iAOtSO(jAO+jComp,0)-1
               lSO=lSO+1
               Do jj=1,jBas
                  jSOj=iShlSO(jSO0+jj)
#if defined (_DEBUG_)
                  If (jSOj.lt.1 .or. jSOj.gt.nBasSh(jShell)) Then
                     Call WarningMessage(2,
     &                                    SecNam//' jSOj out of bounds')
                     Call LDF_Quit(1)
                  End If
#endif
                  i_j_0=nBasSh(iShell)*(jSOj-1)
                  ij0=iBas*(jj-1)
                  Do ii=1,iBas
                     iSOi=iShlSO(iSO0+ii)
#if defined (_DEBUG_)
                     If (iSOi.lt.1 .or. iSOi.gt.nBasSh(iShell)) Then
                        Call WarningMessage(2,
     &                                    SecNam//' iSOi out of bounds')
                        Call LDF_Quit(1)
                     End If
                     If ((i_j_0+iSOi).lt.1 .or. (i_j_0+iSOi).gt.lS) Then
                        Call WarningMessage(2,
     &                                     SecNam//' i_j out of bounds')
                        Call LDF_Quit(1)
                     End If
#endif
                     S(i_j_0+iSOi)=SOInt(ij0+ii,lSO)
                  End Do
               End Do
            End Do
         End Do
      End If

      End

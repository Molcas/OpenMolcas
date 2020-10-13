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
      Subroutine LDF_ComputeAuxInt(ip_ABV)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute integrals over the auxiliary functions (incl. 2C
C              functions, if present):
C
C                 \int J(r)dr
C
C              Argument ip_ABV (input) must point to an auxiliary basis
C              vector allocated with subroutine LDF_AllocateAuxBasVector
C              such that
C                 iWork(ip_ABV-1+A)
C              points to the block of integrals for functions centered
C              on atom A and
C                 iWork(ip_ABV-1+LDF_nAtom()+AB)
C              points to the block of integrals for 2C functions on atom
C              pair AB.
C
      Implicit None
      Integer ip_ABV
#include "WrkSpc.fh"
#include "localdf.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Character*8 Label

      Integer A, AB
      Integer ip, l
      Integer ip_ABV0

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Set up one-electron integral info
      Label='Mltpl  0'
      Call LDF_SetOneEl(Label)

      ! Compute one-center integrals
      ip_ABV0=ip_ABV-1
      Do A=1,LDF_nAtom()
         l=LDF_nBasAux_Atom(A)
         ip=iWork(ip_ABV0+A)
         Call LDF_ComputeAuxInt_1(A,l,Work(ip))
      End Do

      ! Compute two-center integrals
      If (LDF2) Then
         ip_ABV0=ip_ABV0+LDF_nAtom()
         Do AB=1,NumberOfAtomPairs
            l=AP_2CFunctions(1,AB)
            If (l.gt.0) Then
               ip=iWork(ip_ABV0+AB)
               Call LDF_ComputeAuxInt_2(AB,l,Work(ip))
            End If
         End Do
      End If

      ! Unset one-electron integral info
      Call LDF_UnsetOneEl(Label)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ComputeAuxInt_1(A,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute \int J(r)dr for auxiliary functions J(r) on atom
C              A.
C
      use iSD_data
      Implicit None
      Integer A
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
*
      Integer iPrim, jPrim, iAng, jAng, nOrder, lFinal, lScrt1,
     &        lScrt2, MemKrn, ixyz, iBas, MemKer, nElem, mFinal,
     &        mScrt1, mScrt2
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "localdf_bas.fh"
#include "rmat_option.fh"
#include "ldf_oneel.fh"
#include "nsd.fh"
#include "setup.fh"
#include "property_label.fh"
      Real*8, Dimension(:), Allocatable :: Final, Scrtch,
     &                                     ScrSph, Kern
      Character*19 SecNam
      Parameter (SecNam='LDF_ComputeAuxInt_1')

      External MltInt, MltMem

      Integer  LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom
      External LDF_nBasAux_Atom, LDF_nAuxShell_Atom, LDF_lAuxShell_Atom

      Character*8 Label

      Logical Do_PGamma

      Integer l_xInt
      Integer ipInt
      Integer jS, ipS
      Integer iShell, jShell
      Integer iPrint
      Integer nGrid
      Integer iAddPot
      Integer nOrdOp
      Integer ip_SOInt, l_SOInt
      Integer jBas, jCmp, jAO

      Real*8  PtChrg(1)

      Integer i
      Integer nBasSh
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nBasSh(i)=iWork(ip_nBasSh-1+i)
*                                                                      *
************************************************************************
*                                                                      *
      ! Check operator label
      If (OperatorLabel.ne.'Mltpl  0') Then
         Call WarningMessage(2,SecNam//': illegal operator label')
         Call LDF_Quit(1)
      End If

      ! Check dimension
      l_xInt=LDF_nBasAux_Atom(A)
      If (l_xInt.lt.1) Return
      If (l_xInt.gt.l_xInt_) Then
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
      ipS=LDF_lAuxShell_Atom(A)-1
      Do jS=1,LDF_nAuxShell_Atom(A)
         jShell=iWork(ipS+jS)
         l_SOInt=max(l_SOInt,nBasSh(jShell))
      End Do
      Call GetMem(' SO ','Allo','Real',ip_SOInt,l_SOInt)
*
      iShell=nShell_Valence+nShell_Auxiliary+1
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate scratch for the evaluation of the integrals             *
*                                                                      *
      lFinal=1
      lScrt1=1
      lScrt2=1
      MemKrn=1
      Do jS=1,LDF_nAuxShell_Atom(A)
         jShell=iWork(ipS+jS)
*
         iPrim=iSD(5,iShell)
         jPrim=iSD(5,jShell)
         iBas=iSD(3,iShell)
         jBas=iSD(3,jShell)
         iAng=iSD(1,iShell)
         jAng=iSD(1,jShell)
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
*
      Call mma_allocate(Final,lFinal,label='Final')
      Call mma_allocate(Scrtch,lScrt1,label='Scrtch')
      Call mma_allocate(ScrSph,lScrt2,label='ScrSph')
      Call mma_allocate(Kern,MemKrn,label='Kern')
*                                                                      *
************************************************************************
*                                                                      *

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      ipInt=1
      Do jS=1,LDF_nAuxShell_Atom(A)
         jShell=iWork(ipS+jS)
         Call Cho_dZero(Work(ip_SOInt),nBasSh(jShell))
         Call OneEl_IJ(iShell,jShell,iPrint,Do_PGamma,
     &                 Work(ip_xZeta),Work(ip_xZI),Work(ip_xKappa),
     &                 Work(ip_xPCoor),
     &                 MltInt,MltMem,Label,iWork(ip_lOper),nComp,
     &                 Work(ip_CCoor),
     &                 nOrdOp,iWork(ip_kOper),
     &                 iStabO,nStabO,nIC,
     &                 PtChrg,nGrid,iAddPot,Work(ip_SOInt),l_SOInt,
     &                 Final,lFinal,Scrtch,lScrt1,
     &                 ScrSph,lScrt2,Kern,MemKrn)
         jCmp=iSD(2,jShell)
         jBas=iSD(3,jShell)
         jAO=iSD(7,jShell)
#if defined (_DEBUGPRINT_)
         If (jBas*jCmp.ne.nBasSh(jShell)) Then
            Call WarningMessage(2,SecNam//': jBas*jCmp != nBasSh')
            Write(6,'(A,5(1X,I10))')
     &      'A,jBas,jCmp,jShell,nBasSh(jShell)=',
     &      A,jBas,jCmp,jShell,nBasSh(jShell)
            Call LDF_Quit(1)
         End If
#endif
         Call LDF_SortAuxInt_1(Work(ip_SOInt),jBas,jCmp,jAO,
     &                         xInt(ipInt),nBasSh(jShell))
         ipInt=ipInt+nBasSh(jShell)
      End Do
*
      Call mma_deallocate(Final)
      Call mma_deallocate(Scrtch)
      Call mma_deallocate(ScrSph)
      Call mma_deallocate(Kern)
*
#if defined (_DEBUGPRINT_)
      If ((ipInt-1).ne.l_xInt) Then
         Call WarningMessage(2,
     &                    SecNam//': integral array dimension problem!')
         Call LDF_Quit(1)
      End If
#endif
      Call GetMem(' SO ','Free','Real',ip_SOInt,l_SOInt)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SortAuxInt_1(SOInt,jBas,jCmp,jAO,PrpInt,nPrp)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: sort the integrals of auxiliary functions.
C
C     NOTE: symmetry not implemented!
C
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (a-h,o-z)
      Integer jBas
      Integer jCmp
      Real*8  SOInt(jBas,jCmp)
      Integer jAO
      Integer nPrp
      Real*8  PrpInt(nPrp)
#include "WrkSpc.fh"
#include "localdf_bas.fh"

#if defined (_DEBUGPRINT_)
      Character*16 SecNam
      Parameter (SecNam='LDF_SortAuxInt_1')
#endif

      Integer jComp, jSO0, j

      Integer i
      Integer iShlSO
      iShlSO(i)=iWork(ip_iShlSO-1+i)

      Do jComp=1,jCmp
         jSO0=iAOtSO(jAO+jComp,0)-1
         Do j=1,jBas
#if defined (_DEBUGPRINT_)
            If (iShlSO(jSO0+j).lt.1 .or. iShlSO(jSO0+j).gt.nPrp) Then
               Call WarningMessage(2,
     &                          SecNam//': PrpInt index out of bounds!')
               Call LDF_Quit(1)
            End If
#endif
            PrpInt(iShlSO(jSO0+j))=SOInt(j,jComp)
         End Do
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ComputeAuxInt_2(AB,l_xInt_,xInt)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: compute \int J(r)dr for 2C auxiliary functions J(r) on
C              atom pair AB.
C
      use iSD_data
      Implicit None
      Integer AB
      Integer l_xInt_
      Real*8  xInt(l_xInt_)
      Integer iPrim, jPrim, iAng, jAng, nOrder, lFinal, lScrt1,
     &        lScrt2, MemKrn, ixyz, MemKer, nElem, mFinal, mScrt1,
     &        mScrt2
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "ldf_atom_pair_info.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"
#include "ldf_oneel.fh"
#include "rmat_option.fh"
#include "nsd.fh"
#include "setup.fh"
#include "property_label.fh"
      Real*8, Dimension(:), Allocatable :: Final, Scrtch,
     &                                     ScrSph, Kern
      Character*19 SecNam
      Parameter (SecNam='LDF_ComputeAuxInt_2')

      External MltInt, MltMem

      Integer  LDF_nBasAux_Pair
      External LDF_nBasAux_Pair

      Character*8 Label

      Logical Do_PGamma

      Integer l_xInt
      Integer ijS
      Integer iShell, jShell
      Integer iPrint
      Integer nGrid
      Integer iAddPot
      Integer nOrdOp
      Integer ip_SOInt, l_SOInt
      Integer iBas, iCmp, iAO
      Integer jBas, jCmp, jAO
      Integer iCount, jCount

      Real*8 PtChrg(1)

      Integer i, j
      Integer nBasSh
      Integer AP_2CFunctions
      Integer AP_2CList
*
*     Statement functions
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      AP_2CList(i,j)=iWork(ip_2CList-1+l_2CList_1*(j-1)+i)
*                                                                      *
************************************************************************
*                                                                      *

      ! Check operator label
      If (OperatorLabel.ne.'Mltpl  0') Then
         Call WarningMessage(2,SecNam//': illegal operator label')
         Call LDF_Quit(1)
      End If

      ! Check integral array dimension
      l_xInt=AP_2CFunctions(1,AB)
      If (l_xInt.lt.1) Return
      If (l_xInt.gt.l_xInt_) Then
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

      ! Set info for 2C functions
      Call LDF_SetIndxG(AB)
      Call LDF_ColMod(LDF_nBasAux_Pair(AB)-AP_2CFunctions(1,AB))

      ! Allocate temporary integral array
      l_SOInt=0
      Do ijS=1,l_2CList_2
         iShell=AP_2CList(1,ijS)
         jShell=AP_2CList(2,ijS)
         l_SOInt=max(l_SOInt,nBasSh(iShell)*nBasSh(jShell))
      End Do
      Call GetMem(' SO ','Allo','Real',ip_SOInt,l_SOInt)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate scratch for the evaluation of the integrals             *
*                                                                      *
      lFinal=1
      lScrt1=1
      lScrt2=1
      MemKrn=1
      Do ijS=1,l_2CList_2
         iShell=AP_2CList(1,ijS)
         jShell=AP_2CList(2,ijS)
*
         iPrim=iSD(5,iShell)
         jPrim=iSD(5,jShell)
         iBas=iSD(3,iShell)
         jBas=iSD(3,jShell)
         iAng=iSD(1,iShell)
         jAng=iSD(1,jShell)
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
*
      Call mma_Allocate(Final,lFinal,label='Final')
      Call mma_allocate(Scrtch,lScrt1,label='Scrtch')
      Call mma_allocate(ScrSph,lScrt2,label='ScrSph')
      Call mma_allocate(Kern,MemKrn,label='Kern')
*                                                                      *
************************************************************************
*                                                                      *

      ! Compute integrals
      Call Cho_dZero(xInt,l_xInt)
      iCount=0
      Do ijS=1,l_2CList_2
         iShell=AP_2CList(1,ijS)
         jShell=AP_2CList(2,ijS)
         SPAB=AP_2CList(3,ijS)
         SHA=iShell
         SHB=jShell
         Call Cho_dZero(Work(ip_SOInt),nBasSh(iShell)*nBasSh(jShell))
         Call OneEl_IJ(iShell,jShell,iPrint,Do_PGamma,
     &                 Work(ip_xZeta),Work(ip_xZI),Work(ip_xKappa),
     &                 Work(ip_xPCoor),
     &                 MltInt,MltMem,Label,iWork(ip_lOper),nComp,
     &                 Work(ip_CCoor),
     &                 nOrdOp,iWork(ip_kOper),
     &                 iStabO,nStabO,nIC,
     &                 PtChrg,nGrid,iAddPot,Work(ip_SOInt),l_SOInt,
     &                 Final,lFinal,Scrtch,lScrt1,
     &                 ScrSph,lScrt2,Kern,MemKrn)
         iCmp=iSD(2,iShell)
         iBas=iSD(3,iShell)
         iAO=iSD(7,iShell)
         jCmp=iSD(2,jShell)
         jBas=iSD(3,jShell)
         jAO=iSD(7,jShell)
         jCount=0
         Call LDF_SortAuxInt_2(Work(ip_SOInt),iBas,jBas,iCmp,jCmp,
     &                         iAO,jAO,jCount,xInt,l_xInt)
         iCount=iCount+jCount
      End Do
*
      Call mma_deallocate(Final)
      Call mma_deallocate(Scrtch)
      Call mma_deallocate(ScrSph)
      Call mma_deallocate(Kern)
*
      If (iCount.ne.l_xInt) Then
         Call WarningMessage(2,SecNam//': missing integrals!!')
         Call LDF_Quit(1)
      End If

      ! Deallocate temp int array
      Call GetMem(' SO ','Free','Real',ip_SOInt,l_SOInt)

      ! Unset info for 2C functions
      Call LDF_UnsetIndxG()

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SortAuxInt_2(SOInt,iBas,jBas,iCmp,jCmp,
     &                            iAO,jAO,iCount,PrpInt,nPrp)
C
C     Thomas Bondo Pedersen, February 2011.
C     - based on SOSctt by R. Lindh.
C
C     Purpose: sort the integrals of 2C auxiliary functions.
C
C     NOTE: symmetry not implemented!
C
      use SOAO_Info, only: iAOtSO
      Implicit Real*8 (a-h,o-z)
      Integer iBas, jBas
      Real*8  SOInt(iBas*jBas,*)
      Integer iCmp, jCmp
      Integer iAO, jAO
      Integer iCount
      Integer nPrp
      Real*8  PrpInt(nPrp)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "localdf_int.fh"

#if defined (_DEBUGPRINT_)
      Character*16 SecNam
      Parameter (SecNam='LDF_SortAuxInt_2')
      Character*53 Str
      Parameter (Str=
     &          ': extracting to the same array position several times')
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      IndxG2(i,j)=iWork(ip_IndxG2-1+l_IndxG2_1*(j-1)+i)

#if defined (_DEBUGPRINT_)
      nErr=0
      Do iComp=1,iCmp
         Do ii=0,iBas-1
            If (iSOShl(iAOtSO(iAO+iComp,0)+ii).ne.SHA) Then
               nErr=nErr+1
            End If
         End Do
      End Do
      If (nErr.ne.0) Then
         Call WarningMessage(2,SecNam//': shell problem [SHA]')
         Call LDF_Quit(1)
      End If
      nErr=0
      Do jComp=1,jCmp
         Do jj=0,jBas-1
            If (iSOShl(iAOtSO(jAO+jComp,0)+jj).ne.SHB) Then
               nErr=nErr+1
            End If
         End Do
      End Do
      If (nErr.ne.0) Then
         Call WarningMessage(2,SecNam//': shell problem [SHB]')
         Call LDF_Quit(1)
      End If
      If (SHA.eq.SHB) Then
         If (iBas.ne.jBas .or. iCmp.ne.jCmp .or. iAO.ne.jAO) Then
            Call WarningMessage(2,SecNam//': shell problem [SHA=SHB]')
            Call LDF_Quit(1)
         End If
      End If
      lPrpInt=max(nPrp,1)
      Call GetMem('PrpSav','Allo','Real',ipPrpInt,lPrpInt)
      Call dCopy_(nPrp,PrpInt,1,Work(ipPrpInt),1)
#endif

      iCount=0
      lSO=0
      If (SHA.eq.SHB) Then
         Do iComp=1,iCmp
            iSO0=iAOtSO(iAO+iComp,0)-1
            Do jComp=1,iComp-1
               jSO0=iAOtSO(jAO+jComp,0)-1
               lSO=lSO+1
               Do jj=1,jBas
                  jSOj=iShlSO(jSO0+jj)
                  i_j_0=nBasSh(SHA)*(jSOj-1)
                  Do ii=1,iBas
                     iSOi=iShlSO(iSO0+ii)
                     i_j=i_j_0+iSOi
                     If (IndxG2(i_j,SPAB).gt.0) Then
#if defined (_DEBUGPRINT_)
                        If (IndxG2(i_j,SPAB).gt.nPrp) Then
                           Call WarningMessage(2,
     &                              SecNam//': index out of bounds [1]')
                           Call LDF_Quit(1)
                        End If
                        Tst=PrpInt(IndxG2(i_j,SPAB))
     &                     -Work(ipPrpInt-1+IndxG2(i_j,SPAB))
                        If (abs(Tst).gt.1.0d-14) Then
                           Call WarningMessage(2,SecNam//Str//' [1]')
                           Call LDF_Quit(1)
                        End If
#endif
                        ij=iBas*(jj-1)+ii
                        PrpInt(IndxG2(i_j,SPAB))=SOInt(ij,lSO)
                        iCount=iCount+1
                     End If
                  End Do
               End Do
            End Do
            jSO0=iSO0
            lSO=lSO+1
            Do jj=1,jBas
               jSOj=iShlSO(jSO0+jj)
               i_j_0=nBasSh(SHA)*(jSOj-1)
               Do ii=jj,iBas
                  iSOi=iShlSO(iSO0+ii)
                  i_j=i_j_0+iSOi
                  If (IndxG2(i_j,SPAB).gt.0) Then
#if defined (_DEBUGPRINT_)
                     If (IndxG2(i_j,SPAB).gt.nPrp) Then
                        Call WarningMessage(2,
     &                            SecNam//': index out of bounds [1.2]')
                        Call LDF_Quit(1)
                     End If
                     Tst=PrpInt(IndxG2(i_j,SPAB))
     &                  -Work(ipPrpInt-1+IndxG2(i_j,SPAB))
                     If (abs(Tst).gt.1.0d-14) Then
                        Call WarningMessage(2,SecNam//Str//' [1.2]')
                        Call LDF_Quit(1)
                     End If
#endif
                     ij=iBas*(jj-1)+ii
                     PrpInt(IndxG2(i_j,SPAB))=SOInt(ij,lSO)
                     iCount=iCount+1
                  End If
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
                  i_j_0=nBasSh(SHA)*(jSOj-1)
                  Do ii=1,iBas
                     iSOi=iShlSO(iSO0+ii)
                     i_j=i_j_0+iSOi
                     If (IndxG2(i_j,SPAB).gt.0) Then
#if defined (_DEBUGPRINT_)
                        If (IndxG2(i_j,SPAB).gt.nPrp) Then
                           Call WarningMessage(2,
     &                              SecNam//': index out of bounds [2]')
                           Call LDF_Quit(1)
                        End If
                        Tst=PrpInt(IndxG2(i_j,SPAB))
     &                     -Work(ipPrpInt-1+IndxG2(i_j,SPAB))
                        If (abs(Tst).gt.1.0d-14) Then
                           Call WarningMessage(2,SecNam//Str//' [2]')
                           Call LDF_Quit(1)
                        End If
#endif
                        ij=iBas*(jj-1)+ii
                        PrpInt(IndxG2(i_j,SPAB))=SOInt(ij,lSO)
                        iCount=iCount+1
                     End If
                  End Do
               End Do
            End Do
         End Do
      End If

#if defined (_DEBUGPRINT_)
      Call GetMem('PrpSav','Free','Real',ipPrpInt,lPrpInt)
#endif

      End

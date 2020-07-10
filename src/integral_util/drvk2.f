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
* Copyright (C) 1990,1991,1993, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drvk2(Cmpct,DoFock,DoGrad)
************************************************************************
*                                                                      *
*  Object: to precompute all pair entites as zeta, kappa, P and the    *
*          integral prescreening vector to be used with the Schwartz   *
*          inequlity.                                                  *
*                                                                      *
* Called from: Drv2El or Server (DP case)                              *
*                                                                      *
* Calling    : QEnter                                                  *
*              mHrr                                                    *
*              DCopy   (ESSL)                                          *
*              MemRys                                                  *
*              PSOAO0                                                  *
*              DCR                                                     *
*              k2Loop                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN.                              *
*             June '91, modified for k2 loop.                          *
*             Modified for direct SCF, January '93                     *
************************************************************************
      use Real_Spherical
      use k2_setup
      use iSD_data
      use k2_arrays
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
      External Cmpct
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "lundio.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
*     Local arrays
      Real*8  Coor(3,4)
      Integer   iAngV(4), iCmpV(4), iDCRR(0:7), iShllV(2)
      Logical DoFock, force_part_save, DoGrad, ReOrder, Rls
      Character*100 Get_ProgName, ProgName
      Character*8 Method
      Real*8, Allocatable:: HRRMtrx(:,:), Scr(:,:)
      Real*8, Allocatable:: Knew(:), Lnew(:), Pnew(:), Qnew(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      nElem(i)=(i+1)*(i+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      If (k2_Status.eq.Produced) Return
*
      ReOrder=.False.
      ProgName=Get_ProgName()
      If (Index(ProgName,'scf').ne.0) Then
         Call Get_cArray('Relax Method',Method,8)
         ReOrder=Method.eq.'KS-DFT'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 240
      iPrint = nPrint(iRout)
*     iPrint = 99
      Call QEnter('Drvk2')
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      DoGrad_=DoGrad
      DoHess_=.False.
      la_=iAngMx
      mabMin_=nabSz(Max(la_,la_)-1)+1
      mabMax_=nabSz(la_+la_)
      ne_=(mabMax_-mabMin_+1)
      nHrrMtrx=ne_*nElem(la_)*nElem(la_)
      call mma_allocate(HRRMtrx,nHRRMtrx,2,Label='HRRMatrix')
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for k2 data. Observe that the call to Allok2
*     can be done elsewhere and then this call will simply result in
*     a return.
*
      Call Allok2()
*                                                                      *
************************************************************************
*                                                                      *
      nScree = 0
      mScree = 0
      jpk2 = 1
      nk2 = 0
      mk2 = 0
*
*     Allocate memory for zeta, kappa, and P.
      ipZInv = ipZeta + m2Max
      ipKab  = ipZInv + m2Max
      ipP    = ipKab  + m2Max
      ipCon  = ipP    + m2Max*3
      ipAlpha= ipCon  + m2Max
      ipBeta = ipAlpha+ m2Max
      ipInd  = ipiZet
*                                                                      *
************************************************************************
*                                                                      *
      MemTmp=0
      Do iAng = 0, iAngMx
         MemTmp=Max(MemTmp,(MaxPrm(iAng)*nElem(iAng))**2)
      End Do
      Call mma_allocate(Scr,MemTmp,3,Label='Scr')
      Call mma_allocate(Knew,m2Max,Label='Knew')
      Call mma_allocate(Lnew,m2Max,Label='Lnew')
      Call mma_allocate(Pnew,m2Max*3,Label='Pnew')
      Call mma_allocate(Qnew,m2Max*3,Label='Qnew')
*                                                                      *
************************************************************************
*                                                                      *
      If (Allocated(Sew_Scr)) Then
         Rls=.False.
         MemMax=SIZE(Sew_Scr)
C        Write (*,*) 'Drvk2: Memory already allocated:',MemMax
      Else
         Rls=.True.
         Call mma_maxDBLE(MemMax)
         Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
C        Write (*,*) 'Drvk2: Memory allocated:',MemMax
      End If
      ipMem1=1
*                                                                      *
************************************************************************
*                                                                      *
*-----Canonical double loop over shells.
*
      Do iS = 1, mSkal
         iShll  = iSD( 0,iS)
         If (AuxShell(iShll).and.iS.ne.mSkal) Go To 100
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         Coor(1:3,1)=dbsc(iCnttp)%Coor(1:3,iCnt)
*
         If (ReOrder) Call OrdExpD2C(iPrim,Shells(iShll)%Exp,iBas,
     &                                     Shells(iShll)%pCff)
*
         iAngV(1) = iAng
         iShllV(1) = iShll
         iCmpV(1) = iCmp
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            If (AuxShell(iShll).and..Not.AuxShell(jShll)) Go To 200
            If (AuxShell(jShll).and.jS.eq.mSkal) Go To 200
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            Coor(1:3,2)=dbsc(jCnttp)%Coor(1:3,jCnt)
*
            iAngV(2) = jAng
            iShllV(2) = jShll
            iCmpV(2) = jCmp
*
*           Fix for the dummy basis set
            If (AuxShell(iShll)) Coor(1:3,1)=Coor(1:3,2)
*
            Call iCopy(2,iAngV(1),1,iAngV(3),1)
            Call ICopy(2,iCmpV(1),1,iCmpV(3),1)
*
            iPrimi   = iPrim
            jPrimj   = jPrim
            nBasi    = iBas
            nBasj    = jBas
*
            iBasi = iPrimi
            jBasj = jPrimj
            kPrimk = 1
            lPriml = 1
            kBask = 1
            lBasl = 1
*
            nZeta = iPrimi * jPrimj
*
            Call ConMax(Mem_DBLE(ipCon),iPrimi,jPrimj,
     &                  Shells(iShll)%pCff,nBasi,
     &                  Shells(jShll)%pCff,nBasj)
*
            call dcopy_(6,Coor(1,1),1,Coor(1,3),1)
            If (iPrint.ge.99) Call RecPrt(' Sym. Dist. Centers',' ',
     &                                    Coor,3,4)
*
            ijS=iTri(iShell,jShell)
            If (DoFock) Then
               ipDij =ipOffD(1,ijS)
               nDCR  =ipOffD(2,ijS)
               nDij  =ipOffD(3,ijS)
            Else
               ipDij=ip_Dummy
               nDCR  =1
               nDij=1
            End If
*
            nSO = 1
*
*           Compute memory request for the primitives, i.e. how much
*           memory is needed up to the transfer equation.
*
            Call MemRys(iAngV,MemPrm)
*
*           Decide on the partioning of the shells based on
*           on the available memory and the requested memory
*
*-----------Now do a dirty trick to avoid splitting of the first
*           contracted index. Move all over on the second index.
*
            iPrims=iPrimi
            jPrims=jPrimj
            iBasi = 1
            jBasj = nZeta
            iPrimi = 1
            jPrimj = nZeta
            force_part_save=force_part_c
            force_part_c=.False.
            Call PSOAO0(nSO,MemPrm, MemMax,
     &                  iAngV, iCmpV,
     &                  iBasi,iBsInc, jBasj,jBsInc,
     &                  kBask,kBsInc, lBasl,lBsInc,
     &                  iPrimi,iPrInc,jPrimj,jPrInc,
     &                  kPrimk,kPrInc,lPriml,lPrInc,
     &                  ipMem1,ipMem2,
     &                  Mem1,  Mem2,.FALSE.)
            force_part_c=force_part_save
            ijInc = Min(jBsInc,jPrInc)
            iPrimi = iPrims
            jPrimj = jPrims
            If (iPrint.ge.59) Then
               Write (6,*) ' ************** Memory',
     &                     ' partioning **************'
               Write (6,*) ' ipMem1=',ipMem1
               Write (6,*) ' ipMem2=',ipMem2
               Write (6,*) ' Mem1=',Mem1
               Write (6,*) ' Mem2=',Mem2
               Write (6,*) ' *********************',
     &                     '**************************'
            End If
*
*           Find the Double Coset Representatives for center A and B.
*
            Call ICopy(nIrrep,iOper,1,iDCRR,1)
            nDCRR=nIrrep
*
*           Compute all pair entities (zeta, kappa, P, and [nm|nm],
*           total of six types) for all possible unique pairs of
*           centers generated for the symmetry unique centers A and B.
*
*           Write (*,*) ' Generating batches:', mk2+1,' -',
*    &                  mk2+nDCRR
*           Write (*,*) ' Shell index =', ijS, iShell, jShell
*           Write (*,*) ' jpk2=',jpk2
            nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
            nHm=nHm*nIrrep
            ijCmp=nElem(iAng)*nElem(jAng)
            If (.Not.DoGrad_) ijCmp=0
            Call k2Loop(Coor,
     &                  iAngV,iCmpV,iShllV,
     &                  iDCRR,nDCRR,Data_k2(jpk2),
     &                  Shells(iShll)%Exp,iPrimi,
     &                  Shells(jShll)%Exp,jPrimj,
     &                  Mem_DBLE(ipAlpha),Mem_DBLE(ipBeta),
     &                  Shells(iShll)%pCff,nBasi,
     &                  Shells(jShll)%pCff,nBasj,
     &                  Mem_DBLE(ipZeta),Mem_DBLE(ipZInv),
     &                  Mem_DBLE(ipKab),Mem_DBLE(ipP),Mem_INT(ipInd),
     &                  nZeta,ijInc,Mem_DBLE(ipCon),
     &                  Sew_Scr(ipMem2),Mem2,Cmpct,
     &                  nScree,mScree,mdci,mdcj,
     &                  DeDe(ipDij),nDij,nDCR  ,nHm,ijCmp,DoFock,
     &                  Scr, MemTmp,
     &                  Knew,Lnew,Pnew,Qnew,m2Max,DoGrad,
     &                  HrrMtrx,nHrrMtrx)
*
            Indk2(1,ijS) = jpk2
            Indk2(2,ijS) = nDCRR
            nData=nZeta*(nDArray+2*ijCmp)+nDScalar+nHm
            nk2 = nk2 + nData*nDCRR
            mk2 = mk2 + nDCRR
            jpk2 = 1 + nk2
*
 200        Continue
         End Do
 100     Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (Rls) Then
C        Write (6,*) 'Drvk2: Release Sew_Scr'
         Call mma_deallocate(Sew_Scr)
      End If
      Call mma_deallocate(Qnew)
      Call mma_deallocate(Pnew)
      Call mma_deallocate(Lnew)
      Call mma_deallocate(Knew)
      Call mma_deallocate(Scr)
      Call mma_deallocate(HRRMtrx)
*                                                                      *
************************************************************************
*                                                                      *
*     rScree = One -(One*mScree)/(One*nScree)
      If (iPrint.ge.19) Then
      Write (6,*)
      Write (6,*) ' *** The k2 entities has been precomputed ***'
      Write (6,'(I7,A)') mk2,' blocks of k2 data were computed and'
      Write (6,'(I7,A)') nk2,' Word(*8) of memory is used for storage.'
      If (lSchw) Then
         Write (6,*) ' Prescreening based on primitive integrals.'
      Else
         Write (6,*) ' Prescreening based on radial overlap.'
      End If
*     Write (*,'(1X,A,F7.5)') 'Pair screening ratio:',rScree
      Write (6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Call QExit('Drvk2')
      Call CWTime(TCpu2,TWall2)
      Call SavTim(2,TCpu2-TCpu1,TWall2-TWall1)
      k2_Status=Produced
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

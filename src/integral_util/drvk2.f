!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Drvk2(Cmpct,DoFock,DoGrad)
!***********************************************************************
!                                                                      *
!  Object: to precompute all pair entites as zeta, kappa, P and the    *
!          integral prescreening vector to be used with the Schwartz   *
!          inequlity.                                                  *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN.                              *
!             June '91, modified for k2 loop.                          *
!             Modified for direct SCF, January '93                     *
!***********************************************************************
      use setup, only: mSkal
      use iSD_data, only: iSD
      use k2_arrays, only: DoGrad_, DoHess_, DeDe,
     &                     ipOffD, Sew_Scr,
     &                      Create_BraKet, Destroy_BraKet, BraKet
      use Basis_Info, only: Shells, DBSC
      use Symmetry_Info, only: nIrrep, iOper
      use Gateway_global, only: force_part_c
      use Sizes_of_Seward, only: S
      use k2_structure, only: k2_Processed, k2Data, Indk2
#ifdef _DEBUGPRINT_
      use Gateway_Info, only: lSchw
#endif
      use UnixInfo, only: ProgName
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      External Cmpct
      Logical DoFock, DoGrad

!     Local arrays
      Real*8  Coor(3,4)
      Integer   iAngV(4), iCmpV(4), iDCRR(0:7), iShllV(2)
      Logical force_part_save, ReOrder, Rls
      Character(LEN=8) Method
      Real*8, Allocatable:: HRRMtrx(:,:), Scr(:,:)
      Real*8, Allocatable:: Knew(:), Lnew(:), Pnew(:), Qnew(:)
      Integer ik2
      Integer i, j, ixyz, nElem, nabSz, iTri,
     &        iS, jS, iShll, jShll, iBas, jBas,
     &        iPrim, jPrim, la_, mabMin_, mabMax_, ne_, nHrrMtrx,
     &        nScree, mScree, mk2,
     &        MemTmp, iAng, jAng, MemMax, ipMem1, iCmp, jCmp,
     &        iShell, jShell, mdci, mdcj, iCnttp, jCnttp, iCnt, jCnt,
     &        iPrimi, jPrimj, kPrimk, lPriml, nBasi, nBasj,
     &        iBasi, jBasj, kBask, lBasl, nZeta, ijS, nDCR, nDij, nSO,
     &        ipDij, iPrimS, jPrimS, nDCRR, nHm, ijCmp, ipMem2,
     &        iBsInc, jBsInc, kBsInc, lBsInc, ijInc,
     &        iPrInc, jPrInc, kPrInc, lPrInc, Mem1, Mem2, MemPrm
      Real*8 TCPU1, TCPU2, TWALL1, TWALL2

      Interface
      SubRoutine k2Loop(Coor,
     &                  iAnga,iCmpa,iShll,
     &                  iDCRR,nDCRR,
     &                  k2data,
     &                  Alpha,nAlpha,Beta, nBeta,
     &                  Alpha_,Beta_,
     &                  Coeff1,iBasn,Coeff2,jBasn,
     &                  Zeta,ZInv,Kappab,P,IndP,nZeta,IncZZ,Con,
     &                  Wrk,nWork2,
     &                  Cmpct,nScree,mScree,iStb,jStb,
     &                  Dij,nDij,nDCR,ijCmp,DoFock,
     &                  Scr,nScr,
     &                  Knew,Lnew,Pnew,Qnew,nNew,DoGrad,HMtrx,nHrrMtrx)
      use k2_structure, only: k2_type
      Implicit None
      External Cmpct
      Integer nZeta, ijCmp,  nDCRR,
     &        nAlpha, iBasn, nBeta, jBasn, nWork2, nScree, mScree,
     &        iStb, jStb, nDij, nDCR, nScr, nNew, nHRRMtrx, IncZZ
      type(k2_type), intent(inout) :: k2data(nDCRR)
      Real*8 Coor(3,4),
     &       Alpha(nAlpha), Beta(nBeta), Alpha_(nZeta), Beta_(nZeta),
     &       Coeff1(nAlpha,iBasn), Coeff2(nBeta,jBasn),
     &       Zeta(nZeta), ZInv(nZeta), Kappab(nZeta), P(nZeta,3),
     &       Con(nZeta), Wrk(nWork2), Dij(nDij,nDCR), Scr(nScr,3),
     &       Knew(nNew), Lnew(nNew), Pnew(nNew*3), Qnew(nNew*3),
     &       HMtrx(nHrrMtrx,2)
      Logical DoFock, DoGrad
      Integer iAnga(4), iCmpa(4), iShll(2), iDCRR(0:7), IndP(nZeta)
      End SubRoutine k2Loop
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement functions
!
      nElem(i)=(i+1)*(i+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!                                                                      *
!***********************************************************************
!                                                                      *
      If (k2_processed) Return
!
      ReOrder=.False.
      If (Index(ProgName,'scf').ne.0) Then
         Call Get_cArray('Relax Method',Method,8)
         ReOrder=Method.eq.'KS-DFT'
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
      Call CWTime(TCpu1,TWall1)
!                                                                      *
!***********************************************************************
!                                                                      *
      DoGrad_=DoGrad
      DoHess_=.False.
      la_=S%iAngMx
      mabMin_=nabSz(Max(la_,la_)-1)+1
      mabMax_=nabSz(la_+la_)
      ne_=(mabMax_-mabMin_+1)
      nHrrMtrx=ne_*nElem(la_)*nElem(la_)
      call mma_allocate(HRRMtrx,nHRRMtrx,2,Label='HRRMatrix')
!                                                                      *
!***********************************************************************
!                                                                      *
!     Allocate memory for k2 data. Observe that the call to Allok2
!     can be done elsewhere and then this call will simply result in
!     a return.
!
      Call Allok2()
!                                                                      *
!***********************************************************************
!                                                                      *
      nScree = 0
      mScree = 0
      mk2 = 0
!
!     Allocate memory for zeta, kappa, and P.
      Call Create_BraKet(S%m2Max,S%m2Max)
!                                                                      *
!***********************************************************************
!                                                                      *
      MemTmp=0
      Do iAng = 0, S%iAngMx
         MemTmp=Max(MemTmp,(S%MaxPrm(iAng)*nElem(iAng))**2)
      End Do
      Call mma_allocate(Scr,MemTmp,3,Label='Scr')
      Call mma_allocate(Knew,S%m2Max,Label='Knew')
      Call mma_allocate(Lnew,S%m2Max,Label='Lnew')
      Call mma_allocate(Pnew,S%m2Max*3,Label='Pnew')
      Call mma_allocate(Qnew,S%m2Max*3,Label='Qnew')
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Allocated(Sew_Scr)) Then
         Rls=.False.
         MemMax=SIZE(Sew_Scr)
!        Write (*,*) 'Drvk2: Memory already allocated:',MemMax
      Else
         Rls=.True.
         Call mma_maxDBLE(MemMax)
         If (MemMax.gt.1000) MemMax=MemMax-1000
         Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
!        Write (*,*) 'Drvk2: Memory allocated:',MemMax
      End If
      ipMem1=1
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Canonical double loop over shells.
!
      Do iS = 1, mSkal
         iShll  = iSD( 0,iS)
         If (Shells(iShll)%Aux.and.iS.ne.mSkal) Go To 100
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
         iCnttp = iSD(13,iS)
         iCnt   = iSD(14,iS)
         Coor(1:3,1)=dbsc(iCnttp)%Coor(1:3,iCnt)
!
         If (ReOrder) Call OrdExpD2C(iPrim,Shells(iShll)%Exp,iBas,
     &                                     Shells(iShll)%pCff)
!
         iAngV(1) = iAng
         iShllV(1) = iShll
         iCmpV(1) = iCmp
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            If (Shells(iShll)%Aux.and..Not.Shells(jShll)%Aux) Go To 200
            If (Shells(jShll)%Aux.and.jS.eq.mSkal) Go To 200
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
            jCnttp = iSD(13,jS)
            jCnt   = iSD(14,jS)
            Coor(1:3,2)=dbsc(jCnttp)%Coor(1:3,jCnt)
!
            iAngV(2) = jAng
            iShllV(2) = jShll
            iCmpV(2) = jCmp
!
!           Fix for the dummy basis set
            If (Shells(iShll)%Aux) Coor(1:3,1)=Coor(1:3,2)
!
            Call iCopy(2,iAngV(1),1,iAngV(3),1)
            Call ICopy(2,iCmpV(1),1,iCmpV(3),1)
!
            iPrimi   = iPrim
            jPrimj   = jPrim
            nBasi    = iBas
            nBasj    = jBas
!
            iBasi = iPrimi
            jBasj = jPrimj
            kPrimk = 1
            lPriml = 1
            kBask = 1
            lBasl = 1
!
            nZeta = iPrimi * jPrimj
!
            Call ConMax(BraKet%Eta(:),iPrimi,jPrimj,
     &                  Shells(iShll)%pCff,nBasi,
     &                  Shells(jShll)%pCff,nBasj)
!
            call dcopy_(6,Coor(1,1),1,Coor(1,3),1)
#ifdef _DEBUGPRINT_
            Call RecPrt(' Sym. Dist. Centers',' ',
     &                                    Coor,3,4)
#endif
!
            ijS=iTri(iShell,jShell)
            If (DoFock) Then
               ipDij =ipOffD(1,ijS)
               nDCR  =ipOffD(2,ijS)
               nDij  =ipOffD(3,ijS)
            Else
               ipDij= -1
               nDCR  =1
               nDij=1
            End If
!
            nSO = 1
!
!           Compute memory request for the primitives, i.e. how much
!           memory is needed up to the transfer equation.
!
            Call MemRys(iAngV,MemPrm)
!
!           Decide on the partioning of the shells based on
!           on the available memory and the requested memory
!
!-----------Now do a dirty trick to avoid splitting of the first
!           contracted index. Move all over on the second index.
!
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
#ifdef _DEBUGPRINT_
            Write (6,*) ' ************** Memory',
     &                  ' partioning **************'
            Write (6,*) ' ipMem1=',ipMem1
            Write (6,*) ' ipMem2=',ipMem2
            Write (6,*) ' Mem1=',Mem1
            Write (6,*) ' Mem2=',Mem2
            Write (6,*) ' *********************',
     &                  '**************************'
#endif
!
!           Find the Double Coset Representatives for center A and B.
!
            iDCRR(0:nIrrep-1)=iOper(0:nIrrep-1)
            nDCRR=nIrrep
!
!           Compute all pair entities (zeta, kappa, P, and [nm|nm],
!           total of six types) for all possible unique pairs of
!           centers generated for the symmetry unique centers A and B.
!
            nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
            nHm=nHm*nIrrep
            ijCmp=nElem(iAng)*nElem(jAng)
            If (.Not.DoGrad_) ijCmp=0
            ik2=Indk2(3,ijS)
            Call k2Loop(Coor,
     &                  iAngV,iCmpV,iShllV,
     &                  iDCRR,nDCRR,
     &                  k2data(:,ik2),
     &                  Shells(iShll)%Exp,iPrimi,
     &                  Shells(jShll)%Exp,jPrimj,
     &                  BraKet%xA(:),BraKet%xB(:),
     &                  Shells(iShll)%pCff,nBasi,
     &                  Shells(jShll)%pCff,nBasj,
     &                  BraKet%Zeta(:),BraKet%ZInv(:),
     &                  BraKet%KappaAB(:),BraKet%P(:,:),
     &                  BraKet%IndZet(:),
     &                  nZeta,ijInc,BraKet%Eta(:),
     &                  Sew_Scr(ipMem2),Mem2,Cmpct,
     &                  nScree,mScree,mdci,mdcj,
     &                  DeDe(ipDij),nDij,nDCR  ,ijCmp,DoFock,
     &                  Scr, MemTmp,
     &                  Knew,Lnew,Pnew,Qnew,S%m2Max,DoGrad,
     &                  HrrMtrx,nHrrMtrx)
!
            Indk2(2,ijS) = nDCRR
            mk2 = mk2 + nDCRR
!
 200        Continue
         End Do
 100     Continue
      End Do

      Call Destroy_Braket()
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Rls) Then
!        Write (6,*) 'Drvk2: Release Sew_Scr'
         Call mma_deallocate(Sew_Scr)
      End If
      Call mma_deallocate(Qnew)
      Call mma_deallocate(Pnew)
      Call mma_deallocate(Lnew)
      Call mma_deallocate(Knew)
      Call mma_deallocate(Scr)
      Call mma_deallocate(HRRMtrx)
!                                                                      *
!***********************************************************************
!                                                                      *
!     rScree = One -(One*mScree)/(One*nScree)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) ' *** The k2 entities has been precomputed ***'
      Write (6,'(I7,A)') mk2,' blocks of k2 data were computed.'
      If (lSchw) Then
         Write (6,*) ' Prescreening based on primitive integrals.'
      Else
         Write (6,*) ' Prescreening based on radial overlap.'
      End If
!     Write (*,'(1X,A,F7.5)') 'Pair screening ratio:',rScree
      Write (6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!
      k2_processed=.True.
      Call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
      Return
      End

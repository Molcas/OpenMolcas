      subroutine BP_Driver(iiS,jjS,kkS,llS,TInt,nTInt)
      use setup, only: mSkal, nAux, nSOs
      use k2_structure, only: IndK2
      use k2_arrays, only: ipDijS, Sew_Scr,Aux, DeDe, FT, iSOSym,
     &                     nDeDe, nFT, Create_BraKet, Destroy_Braket
      use iSD_data, only: iSD
      use Basis_Info, only: Shells
      use Gateway_Info, only: CutInt
      use Symmetry_Info, only: nIrrep
      use Int_Options, only: DoIntegrals, DoFock, Map4
      use Integral_interfaces, only: Int_PostProcess
#ifdef _DEBUGBREIT_
      use Breit, only: nOrdOp
      use UnixInfo, only: SuperName
#endif
      use Constants, only: Zero
      use stdalloc, only: mma_allocate
      use k2_structure, only: k2data


      implicit none
!eval_ijkl
      real(kind=wp), intent(in) :: ThrAO
      integer(kind=iwp) :: iCnttp, ijS, iOpt, iS, jCnttp, jS, kCnttp, klS, kS, lCnttp
      real(kind=wp) :: A_int, P_Eff, PP_Count, PP_Eff, PP_Eff_delta, S_Eff, ST_Eff, T
                       Twall2
      logical(kind=iwp) :: DoGrad, Indexation, Triangular
      character(len=72) :: SLine
      real(kind=wp), allocatable :: TInt(:), TMax(:,:)
      integer(kind=iwp), parameter :: nTInt = 1
      integer(kind=iwp), allocatable :: Pair_Index(:,:)
      logical(kind=iwp), external :: Rsv_GTList
!drvg1    
      integer(kind=iwp), intent(in) :: nGrad
      real(kind=wp), intent(inout) :: Grad(nGrad)
      real(kind=wp), intent(out) :: Temp(nGrad)
      #include "print.fh"
      integer(kind=iwp) :: i, iAng, iAnga(4), iAOst(4), iAOV(4), iBasAO, iBasi, iBasn
                           ijMax, ijS, ik2, iOpt, iost, ipMem1, ipMem2, iPrem, iPren,
                           iSh, iShela(4), iShlla(4), iSSDM, istabs(4), j, jAng, jBAs
                           jPrInc, jS, k2ij, k2kl, kBasAO, kBask, kBasn, kBsInc, kBtc
                           lBsInc, lPriml, lPrInc, lRealName, lS, luCMOPT2, luGamma,
                           Mem2, MemMax, MemPSO, nab, nBasI, nBasT, nBtch, ncd, nDCRR
                           nijkl, nOcc(8), nPairs, nQuad, nRys, nSkal, nSO, nZeta
      real(kind=wp) :: A_int, Cnt, Coor(3,4), P_Eff, PMax, Prem, Pren, TCpu1, TCpu2,
      logical(kind=iwp) :: ABCDeq, AeqB, CeqD, DoFock, DoGrad, EQ, Indexation, is_err
                           Skip, Triangular
      character(len=4096) :: RealName
      character(len=72) :: formt
      character(len=8) :: Method_chk
      integer(kind=iwp), allocatable :: Ind_ij(:,:), iOffAO(:)
      real(kind=wp), allocatable :: CMOPT2(:), TMax(:,:), WRK1(:), WRK2(:)
      integer(kind=iwp), save :: MemPrm
      integer(kind=iwp), external :: IsFreeUnit
      logical(kind=iwp), external :: Rsv_GTLis

!                                                                      *
#ifdef _DEBUGBREIT_
!     use the Breit option computing 1/r^3 integralas but convert to
!     conventional 1/r integrals
      If (.NOT.DoFock .and.
     &    SuperName/= 'gateway' .and.
     &    nIrrep==1) Call Set_Breit(1)
#endif
      mDCRij=1
      mDCRkl=1
      If (nIrrep==1) Then
         Do_TwoEl => TwoEl_NoSym
      Else
         Do_TwoEl => TwoEl_Sym
      End If
!     If (.NOT.Associated(Int_PostProcess)) Stop 124
!                                                                      *
!***********************************************************************
!                                                                      *
      If (.Not.Allocated(iSOSym)) Then
         Call WarningMessage(2,
     &               'Eval_Ints_: Integral environment is not set up!')
         Call Abend()
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
      NoInts=.True.
      Tmax=Zero
!                                                                      *
!***********************************************************************
!                                                                      *
!     If memory not allocated already at this point allocate!          *
!                                                                      *
      If (.Not.Allocated(Sew_Scr)) Then
!        Write (*,*) 'Eval_ints: Allocate memory'
         Call mma_MaxDBLE(MemMax)
         If (MemMax.gt.8000) MemMax=MemMax-8000
         Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      Else
!        Write (*,*) 'Eval_ints: Memory already allocated'
         MemMax=SIZE(Sew_Scr)
      End If
!     Write (*,*) 'Eval_ints: MemMax=',MemMax
      ipMem1=1
!
      Map4(1)=1
      Map4(2)=2
      Map4(3)=3
      Map4(4)=4
      iS_=Max(iiS,jjS)
      jS_=Min(iiS,jjS)
      kS_=Max(kkS,llS)
      lS_=Min(kkS,llS)
      If (iiS.ne.iS_) Then
          iTmp=Map4(1)
          Map4(1)=Map4(2)
          Map4(2)=iTmp
      End If
      If (kkS.ne.kS_) Then
          iTmp=Map4(3)
          Map4(3)=Map4(4)
          Map4(4)=iTmp
      End If
!     Write (*,*) ' -->',iS_,jS_,kS_,lS_,'<--'
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Int_Setup(iSD,mSkal,iS_,jS_,kS_,lS_,Coor,Shijij,
     &               iAngV,iCmpV,iShelV,iShllV,iAOV,iStabs)
!                                                                      *
!***********************************************************************
!                                                                      *
      iPrimi   = Shells(iShllV(1))%nExp
      jPrimj   = Shells(iShllV(2))%nExp
      kPrimk   = Shells(iShllV(3))%nExp
      lPriml   = Shells(iShllV(4))%nExp
      iBasi    = Shells(iShllV(1))%nBasis
      jBasj    = Shells(iShllV(2))%nBasis
      kBask    = Shells(iShllV(3))%nBasis
      lBasl    = Shells(iShllV(4))%nBasis
      nZeta    = iPrimi * jPrimj
      nEta     = kPrimk * lPriml
      mDij=nZeta+1 ! Dummy initialize
      mDkl=nEta+1  ! Dummy initialize
!                                                                      *
!***********************************************************************
!                                                                      *
!     partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q

      Call Create_BraKet(nZeta,nEta)
!                                                                      *
!***********************************************************************
!                                                                      *
!
!     No SO block in direct construction of the Fock matrix.
      nSO = MemSO2(iCmpV(1),iCmpV(2),iCmpV(3),iCmpV(4),
     &             iShelV(1),iShelV(2),iShelV(3),iShelV(4),
     &             iAOV(1),iAOV(2),iAOV(3),iAOV(4))
      If (nSO.eq.0) Then
        Return
      End If
!
      iS = iShelV(1)
      jS = iShelV(2)
      kS = iShelV(3)
      lS = iShelV(4)
      ijS = iTri(iS,jS)
      klS = iTri(kS,lS)
      ikS = iTri(iS,kS)
      ilS = iTri(iS,lS)
      jkS = iTri(jS,kS)
      jlS = iTri(jS,lS)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Pick up pointers to k2 entities.
!

      nDCRR = IndK2(2,ijS)
      ik2   = IndK2(3,ijS)
      nDCRS = IndK2(2,klS)
      jk2   = IndK2(3,klS)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Pick up pointers to desymmetrized 1st order density
!     matrices. Observe that the desymmetrized 1st order
!     density matrices follows the contraction index.
!
      If (DoFock) Then
         ipTmp = ipDijs
         Nr_of_D=1
         Call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,Nr_of_D)
         Call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,Nr_of_D)
         Call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,Nr_of_D)
         Call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,Nr_of_D)
         Call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,Nr_of_D)
         Call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,Nr_of_D)
!
!        Write (*,*) ' Pointers to D=',
!    &                ipDij,ipDkl,ipDik,ipDil,ipDjk,ipDjl
!
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
!     Write (6,*) ' *** Centers ***'
!     Write (6,'(3F7.3,6X,3F7.3)') ((Coor(i,j),i=1,3),j=1,2)
!     Write (6,'(3F7.3,6X,3F7.3)') ((Coor(i,j),i=1,3),j=3,4)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute memory request for the primitives, i.e.
!     how much memory is needed up to the transfer
!     equation.
      Call MemRys(iAngV,MemPrm)
!                                                                      *
!***********************************************************************
!                                                                      *
!     Decide on the partioning of the shells based on the
!     available memory and the requested memory.
      Call PSOAO0(nSO,MemPrm,MemMax,
      Call PSOAO0(nSO,MemPrm,MemMax,
     &            iAngV,iCmpV,
     &            iBasi,iBsInc,jBasj,jBsInc,
     &            kBask,kBsInc,lBasl,lBsInc,
     &            iPrimi,iPrInc,jPrimj,jPrInc,
     &            kPrimk,kPrInc,lPriml,lPrInc,
     &            ipMem1,ipMem2,
     &            Mem1,Mem2,DoFock)
#ifdef _DEBUGPRINT_
!     Write (6,*) ' ************** Memory partioning **************'
!     Write (6,*) ' ipMem1=',ipMem1
!     Write (6,*) ' ipMem2=',ipMem2
!     Write (6,*) ' Mem1=',Mem1
!     Write (6,*) ' Mem2=',Mem2
!     Write (6,*) ' iBasi,iBsInc=',iBasi,iBsInc
!     Write (6,*) ' jBasj,jBsInc=',jBasj,jBsInc
!     Write (6,*) ' kBasi,kBsInc=',kBask,kBsInc
!     Write (6,*) ' lBasl,lBsInc=',lBasl,lBsInc
!     Write (6,*) ' iPrimi,iPrInc=',iPrimi,iPrInc
!     Write (6,*) ' jPrimj,jPrInc=',jPrimj,jPrInc
!     Write (6,*) ' kPrimk,kPrInc=',kPrimk,kPrInc
!     Write (6,*) ' lPriml,lPrInc=',lPriml,lPrInc
!     Write (6,*) ' ***********************************************'
#endif
      SOInt(1:Mem1)=>Sew_Scr(ipMem1:ipMem1+Mem1-1)
      AOInt(1:Mem2)=>Sew_Scr(ipMem2:ipMem2+Mem1-1)
!                                                                      *
!***********************************************************************
!                                                                      *
      jbas_=jBasj
      lbas_=lBasl
!                                                                      *
!***********************************************************************
!                                                                      *
!     These loops will partition the contraction loops if there is not
!     enough memory to store the whole SO/AO-block simultaneously. The
!     memory partitioning is determined by PSOAO0.
!
      Do iBasAO = 1, iBasi, iBsInc
         iBasn=Min(iBsInc,iBasi-iBasAO+1)
         iAOst(1) = iBasAO-1
!
         Do jBasAO = 1, jBasj, jBsInc
            jBasn=Min(jBsInc,jBasj-jBasAO+1)
            iAOst(2) = jBasAO-1
!
!---------- Move appropiate portions of the desymmetrized 1st
!           order density matrix.
!
            If (DoFock) Then
               Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                     jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                     iCmpV(1),iCmpV(2),iShelV(1),iShelV(2),
     &                     mDCRij,ipDij,ipDDij,mDij,DeDe,nDeDe)
            End If
!
            Do kBasAO = 1, kBask, kBsInc
               kBasn=Min(kBsInc,kBask-kBasAO+1)
               iAOst(3) = kBasAO-1
!
               If (DoFock) Then
                  Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                        kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                        iCmpV(1),iCmpV(3),iShelV(1),iShelV(3),
     &                        mDCRik,ipDik,ipDDik,mDik,DeDe,nDeDe)
               End If
!
               If (DoFock) Then
                  Call Picky_(jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                        kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                        iCmpV(2),iCmpV(3),iShelV(2),iShelV(3),
     &                        mDCRjk,ipDjk,ipDDjk,mDjk,DeDe,nDeDe)
               End If
!
                Do lBasAO = 1, lBasl, lBsInc
                   lBasn=Min(lBsInc,lBasl-lBasAO+1)
                   iAOst(4) = lBasAO-1
!
                   If (DoFock) Then
                      Call Picky_(kBask,kBsInc,kPrimk,kBasAO,kBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(3),iCmpV(4),iShelV(3),iShelV(4),
     &                            mDCRkl,ipDkl,ipDDkl,mDkl,DeDe,nDeDe)
                   End If
!
                   If (DoFock) Then
                      Call Picky_(iBasi,iBsInc,iPrimi,iBasAO,iBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(1),iCmpV(4),iShelV(1),iShelV(4),
     &                            mDCRil,ipDil,ipDDil,mDil,DeDe,nDeDe)
                   End If
!
                   If (DoFock) Then
                      Call Picky_(jBasj,jBsInc,jPrimj,jBasAO,jBasn,
     &                            lBasl,lBsInc,lPriml,lBasAO,lBasn,
     &                            iCmpV(2),iCmpV(4),iShelV(2),iShelV(4),
     &                            mDCRjl,ipDjl,ipDDjl,mDjl,DeDe,nDeDe)
                   End If
!                                                                      *
!***********************************************************************
!                                                                      *
!                 Compute SO/AO-integrals
!
                  Call Do_TwoEl(iS_,jS_,kS_,lS_,Coor,
     &                          iAngV,iCmpV,iShelV,iShllV,
     &                          iAOV,iAOst,NoInts,iStabs,
     &                          iPrimi,iPrInc,jPrimj,jPrInc,
     &                          kPrimk,kPrInc,lPriml,lPrInc,
     &                          nDCRR,
     &                          nDCRS,
     &                          k2Data(:,ik2),k2Data(:,jk2),
     &                          IJeqKL,kOp,
     &                          DeDe(ipDDij),mDij,mDCRij,
     &                          DeDe(ipDDkl),mDkl,mDCRkl,
     &                          DeDe(ipDDik),mDik,mDCRik,
     &                          DeDe(ipDDil),mDil,mDCRil,
     &                          DeDe(ipDDjk),mDjk,mDCRjk,
     &                          DeDe(ipDDjl),mDjl,mDCRjl,
     &                          Shells(iShllV(1))%pCff(1,iBasAO),iBasn,
     &                          Shells(iShllV(2))%pCff(1,jBasAO),jBasn,
     &                          Shells(iShllV(3))%pCff(1,kBasAO),kBasn,
     &                          Shells(iShllV(4))%pCff(1,lBasAO),lBasn,
     &                          FT,nFT,nZeta,nEta,
     &                          SOInt,nSO,AOInt,Mem2,
     &                          Shijij,Aux,nAux)
!                                                                      *
!***********************************************************************
!                                                                      *
!
                  nijkl=iBasn*jBasn*kBasn*lBasn
#ifdef _DEBUGBREIT_
                  If (nOrdOp/=0) Then
                     If (nIrrep==1) Then
                        n=iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
                        Call ReSort_Int(AOInt,nijkl,6,n)
                     Else
                        Call ReSort_Int(SOInt,nijkl,6,nSO)
                     End If
                  End If
#endif
#ifdef _DEBUGPRINT_
                  If (nIrrep==1) Then
                     n=iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
                     Call RecPrt('AOInt',' ',AOInt,nijkl,n)
                  Else
                     Call RecPrt('SOInt',' ',SOInt,nijkl,nSO)
                  End If
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
!              Process SO/AO-integrals
!
                  If (DoIntegrals.and..Not.NoInts) Then
!                    Get max AO/SO integrals
                     If (nIrrep.eq.1) Then
                        n=nijkl*iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
                        Tmax=max(Tmax,
     &                       abs(AOInt(iDAMax_(n,AOInt,1))))
                     Else
                        n=nijkl*nSO
                        Tmax=max(Tmax,
     &                        abs(SOInt(iDAMax_(n,SOInt,1))))
                     End If
                     If (Tmax.gt.CutInt) Then
                        Call Int_PostProcess(iCmpV,iShelV,
     &                                       iBasn,jBasn,kBasn,lBasn,
     &                                       kOp,
     &                                       Shijij,iAOV,iAOst,nijkl,
     &                                       AOInt,SOInt,nSO,
     &                                       iSOSym,nSOs,
     &                                       TInt,nTInt,nIrrep)
                     Else
                        Tmax=Zero
                     End If
                  End If
!
               End Do
            End Do
         End Do
      End Do
      SOInt=>Null()
      AOInt=>Null()
      Call Destroy_BraKet()
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGBREIT_
      Call Set_Breit(0)

      Contains
      Subroutine ReSort_Int(IntRaw,nijkl,nComp,nA)
      Implicit None
      Integer, intent(In) :: nijkl, nComp, nA
      Real*8, Target :: IntRaw(nijkl*nComp*nA)

      Real*8, Pointer :: IntIn(:,:,:), IntOut(:,:,:)
      Integer iA, i_ijkl

      IntIn(1:nijkl,1:nComp,1:nA)=>IntRaw(:)
      IntOut(1:nijkl,1:1,1:nA)=>IntRaw(1:nijkl*nA)
#ifdef _DEBUGPRINT_
      Write (6,*) 'nijkl,nComp,nA=',nijkl,nComp,nA
      Call RecPrt('IntRaw',' ',IntRaw,nijkl,nComp*nA)
#endif

      Do iA = 1, nA
         Do i_ijkl = 1, nijkl
            IntOut(i_ijkl,1,iA) = IntIn(i_ijkl,1,iA)
     &                          + IntIn(i_ijkl,4,iA)
     &                          + IntIn(i_ijkl,6,iA)
         End Do
      End Do

      IntIn=>Null()
      IntOut=>Null()
      End Subroutine ReSort_Int
#endif


!-- Prepare handling of two-particle density.

call PrepP

if (Method_chk == 'CASPT2  ') then
  nBasT = 0
  do i=0,nIrrep-1
    nBasT = nBasT+nBas(i)
  end do
  nSSDM = 0

  !! The two MO indices in the half-transformed amplitude are
  !! not CASSCF but quasi-canonical orbitals.
  call mma_allocate(CMOPT2,nBasT*nBasT,Label='CMOPT2')
  LuCMOPT2 = isFreeUnit(66)
  call PrgmTranslate('CMOPT2',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',i

  do i=1,nBasT*nBasT
    read(LuCMOPT2) CMOPT2(i)
  end do
  read(LuCMOPT2) nOcc(1)
  read(LuCMOPT2) nOcc(2)
  read(LuCMOPT2) nOcc(3)
  read(LuCMOPT2) nOcc(4)
  read(LuCMOPT2) nOcc(5)
  read(LuCMOPT2) nOcc(6)
  read(LuCMOPT2) nOcc(7)
  read(LuCMOPT2) nOcc(8)
  read(LuCMOPT2) nFro(1)
  read(LuCMOPT2) nFro(2)
  read(LuCMOPT2) nFro(3)
  read(LuCMOPT2) nFro(4)
  read(LuCMOPT2) nFro(5)
  read(LuCMOPT2) nFro(6)
  read(LuCMOPT2) nFro(7)
  read(LuCMOPT2) nFro(8)
  read(LuCMOPT2) nSSDM
  if (nSSDM /= 0) then
    call mma_allocate(SSDM,nBas(0)*(nBas(0)+1)/2,2,nSSDM,Label='SSDM')
    do iSSDM=1,nSSDM
      do i=1,nBas(0)*(nBas(0)+1)/2
        read(LuCMOPT2) SSDM(i,1,iSSDM),SSDM(i,2,iSSDM)
      end do
    end do
  end if

  close(LuCMOPT2)

  write(u6,*) 'Number of Non-Frozen Occupied Orbitals = ',nOcc(1)
  write(u6,*) 'Number of     Frozen          Orbitals = ',nFro(1)

  call mma_allocate(iOffAO,nSkal+1,Label='iOffAO')
  MaxShlAO = 0
  iOffAO(1) = 0
  do iSh=1,nSkal
    nBasI = iSD(2,iSh)*iSD(3,iSh)
    if (nBasI > MaxShlAO) MaxShlAO = nBasI
    iOffAO(iSh+1) = iOffAO(iSh)+nBasI
  end do
  call mma_allocate(G_toc,MaxShlAO**4,Label='GtocCASPT2')

  LuGAMMA = isFreeUnit(65)
  call PrgmTranslate('GAMMA',RealName,lRealName)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',io

  call mma_allocate(WRK1,nOcc(1)*nOcc(1),Label='WRK1')
  call mma_allocate(WRK2,MaxShlAO*nOcc(1),Label='WRK2')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
MxPrm = 0
do iAng=0,S%iAngMx
  MxPrm = max(MxPrm,S%MaxPrm(iAng))
end do
nZeta = MxPrm*MxPrm
nEta = MxPrm*MxPrm
!
!***********************************************************************
!                                                                      *
!-- Compute entities for prescreening at shell level

call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
call Shell_MxSchwz(nSkal,TMax)
TMax_all = Zero
do iS=1,nSkal
  do jS=1,iS
    TMax_all = max(TMax_all,TMax(iS,jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,Label='Ind_ij')
nij = 0
do iS=1,nSkal
  do jS=1,iS
    if ((TMax_All*TMax(iS,jS) >= CutInt) .or. (Method_chk == 'CASPT2  ')) then
      nij = nij+1
      Ind_ij(1,nij) = iS
      Ind_ij(2,nij) = jS
    end if
  end do
end do
P_Eff = real(nij,kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
!-- Compute FLOPs for the transfer equation.

do iAng=0,S%iAngMx
  do jAng=0,iAng
    nHrrab = 0
    do i=0,iAng+1
      do j=0,jAng+1
        if (i+j <= iAng+jAng+1) then
          ijMax = min(iAng,jAng)+1
          nHrrab = nHrrab+ijMax*2+1
        end if
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
Triangular = .true.
call Init_TList(Triangular,P_Eff)
call Init_PPList
call Init_GTList
iOpt = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! In MPP case dispatch one processor to do 1-el gradients first

if ((nProcs > 1) .and. King()) then
  call Drvh1(Grad,Temp,nGrad)
  !if (nPrint(1) >= 15) call PrGrad(' Gradient excluding two-electron contribut
  Temp(:) = Zero
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_MaxDBLE(MemMax)
if (MemMax > 8000) MemMax = MemMax-8000
call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
ipMem1 = 1
ijS = 0
klS = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) exit

  ! Now do a quadruple loop over shells

  call Get_cArray('Relax Method',Method_chk,8)
  if (Method_chk /= 'CASPT2  ') then
    ijS = int((One+sqrt(Eight*TskLw-Three))/Two)
    iS = Ind_ij(1,ijS)
    jS = Ind_ij(2,ijS)
    klS = int(TskLw-real(ijS,kind=wp)*(real(ijS,kind=wp)-One)/Two)
    kS = Ind_ij(1,klS)
    lS = Ind_ij(2,klS)
  else
    iS = 1
    jS = 1
    kS = 1
    lS = 1
    !! proceed the index
    do iCnt=1,int(TskLw)-1
      call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
    end do
    Cnt = real(iCnt,kind=wp)
    !! If LoadVec is true, a new vector of the half-transformed
    !! T-amplitude is read. In the first loop, it is always true.
    !! In other loops, a new vector is read only when I- and K-th
    !! are different from the previous loop.
    !! The half back-transformation, T_{ij}^{ab} ->
    !! T_{ij}^{rho sigma}, is done somewhere in CASPT2.
    !! rho and sigma correspond to either I- or K-th shells.
    !! Occupied orbital indices (correspond to J- or L-th shells)
    !! are back-transformed on-the-fly.
    LoadVec = .true.
  end if
  Cnt = TskLw
  call CWTime(TCpu1,TWall1)
  do
    A_int = TMax(iS,jS)*TMax(kS,lS)
    Skip = .false.
    if (A_Int < CutInt) Skip = .true.
    if (.not. Skip) then
      if (iPrint >= 15) write(u6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
      call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
      if (No_batch) Skip = .true.
    end if

    if (.not. Skip) then
      call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! -------> Memory Managment <--------
      !
      ! Compute memory request for the primitives, i.e.
      ! how much memory is needed up to the transfer
      ! equation.

      call MemRys_g(iSD4,nSD,nRys,MemPrm)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ABCDeq = EQ(Coor(1,1),Coor(1,2)) .and. EQ(Coor(1,1),Coor(1,3)) .and. EQ(C
      ijklA = iSD4(1,1)+iSD4(1,2)+iSD4(1,3)+iSD4(1,4)
      if ((nIrrep == 1) .and. ABCDeq .and. (mod(ijklA,2) == 1)) Skip = .true.
    end if
    if (.not. Skip) then
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Decide on the partioning of the shells based on the
      ! available memory and the requested memory.
      !
      ! Now check if all blocks can be computed and stored at once.

      call SOAO_g(iSD4,nSD,nSO,MemPrm,MemMax,iBsInc,jBsInc,kBsInc,lBsInc,iPrInc
                  MemPSO)
      iBasi = iSD4(3,1)
      jBasj = iSD4(3,2)
      kBask = iSD4(3,3)
      lBasl = iSD4(3,4)
      !                                                                *
      !*****************************************************************
      !                                                                *
      call Int_Parm_g(iSD4,nSD,iAnga,iCmpa,iShlla,iShela,iPrimi,jPrimj,kPrimk,l
                      k2ij,ik2,nDCRR,k2kl,jk2,nDCRS,mdci,mdcj,mdck,mdcl, &
                      AeqB,CeqD,nZeta,nEta,l2DI,nab,nHmab,ncd,nHmcd,nIrrep)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Scramble arrays (follow angular index)

     do iCar=1,3
        do iSh=1,4
          JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
          if (btest(iSD4(15,iSh),iCar-1)) then
            JfGrad(iCar,iSh) = .true.
          else
            JfGrad(iCar,iSh) = .false.
          end if
        end do
      end do

      do iBasAO=1,iBasi,iBsInc
        iBasn = min(iBsInc,iBasi-iBasAO+1)
        iAOst(1) = iBasAO-1
        do jBasAO=1,jBasj,jBsInc
          jBasn = min(jBsInc,jBasj-jBasAO+1)
          iAOst(2) = jBasAO-1
          do kBasAO=1,kBask,kBsInc
            kBasn = min(kBsInc,kBask-kBasAO+1)
            iAOst(3) = kBasAO-1
            do lBasAO=1,lBasl,lBsInc
              lBasn = min(lBsInc,lBasl-lBasAO+1)
              iAOst(4) = lBasAO-1

              ! Get the 2nd order density matrix in SO basis.

              nijkl = iBasn*jBasn*kBasn*lBasn

              ! Fetch the T_i,j,kappa, lambda corresponding to
              ! kappa = k, lambda = l

#             ifdef _CD_TIMING_
              call CWTIME(Pget0CPU1,Pget0WALL1)
#             endif
              if (Method_chk == 'CASPT2  ') call CASPT2_BTAMP(LuGAMMA,iS,jS,kS,
                                                              iFnc(4)*lBasn,iOf
                                                              WRK2,G_Toc)
              call PGet0(iCmpa,iBasn,jBasn,kBasn,lBasn,Shijij,iAOV,iAOst,nijkl,
                         iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,Sew_Scr(ipMem2),Mem
              if (A_Int*PMax < CutInt) cycle
#             ifdef _CD_TIMING_
              call CWTIME(Pget0CPU2,Pget0WALL2)
              Pget_CPU = Pget_CPU+Pget0CPU2-Pget0CPU1
              Pget_Wall = Pget_Wall+Pget0WALL2-Pget0WALL1
#             endif
              if (A_Int*PMax < CutInt) cycle

              ! Compute gradients of shell quadruplet

#             ifdef _CD_TIMING_
              call CWTIME(TwoelCPU1,TwoelWall1) ! timing_cdscf
#             endif
              call TwoEl_g(Coor,iAnga,iCmpa,iShela,iShlla,iAOV,mdci,mdcj,mdck,m
                           k2data(:,ik2),k2data(:,jk2), &
                           nDCRR,nDCRS,Pren,Prem,iPrimi,iPrInc,jPrimj,jPrInc,kP
                           Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,Shells(iSD4(0
                           Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,Shells(iSD4(0
                           nZeta,nEta,Temp,nGrad,JfGrad,JndGrd,Sew_Scr(ipMem1),
                           Sew_Scr(ipMem2),Mem2,Aux,nAux,Shijij)
#             ifdef _CD_TIMING_
              call CWTIME(TwoelCPU2,TwoelWall2)
              Twoel_CPU = Twoel_CPU+TwoelCPU2-TwoelCPU1
              Twoel_Wall = Twoel_Wall+TwoelWall2-TwoelWall1
#             endif
              if (iPrint >= 15) call PrGrad(' In Drvg1: Grad',Temp,nGrad,ChDisp

            end do
          end do

        end do
      end do

      call Destroy_BraKet()

    end if
    Cnt = Cnt+One
    if (Cnt-TskHi > 1.0e-10_wp) then
      exit
    else if (Method_chk == 'CASPT2  ') then
      call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
    else
      klS = klS+1
      if (klS > ijS) then
        ijS = ijS+1
        klS = 1
      end if
      iS = Ind_ij(1,ijS)
      jS = Ind_ij(2,ijS)
      kS = Ind_ij(1,klS)
      lS = Ind_ij(2,klS)
    end if
  end do

  ! Task endpoint
  call CWTime(TCpu2,TWall2)
end do
! End of big task loop
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _MOLCAS_MPP_
call GADGOP(Temp,nGrad,'+')
#endif
!                                                                      *


      end subroutine BP_Driver

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
! Copyright (C) 1991,1993,1999,2023, Roland Lindh                      *
!               1995, Martin Schuetz                                   *
!***********************************************************************
!#define _DEBUGPRINT_
!#define _DEBUGBREIT_
      SubRoutine Eval_ijkl(iiS,jjS,kkS,llS,TInt,nTInt)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals, parallel region          *
!          contains memory partitioning and loops over uncontracted    *
!          functions...                                                *
!                                                                      *
!  Input:                                                              *
!          iiS,jjS,kkS,llS     : shell indices                         *
!          TInt                : Computed Integrals                    *
!          nTInt               : dimension of TInt                     *
!          Integ_Proc          : subroutine for post processing        *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh / Martin Schuetz,                           *
!             Dept. of Theoretical Chemistry, University of Lund,      *
!             SWEDEN.                                                  *
!             Modified for k2 loop. August '91                         *
!             Modified for direct SCF. January '93                     *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             parallel region split off in drvtwo.f, April '95         *
!             Total rehack May '99                                     *
!             Total rehack Aug '23                                     *
!***********************************************************************
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
      Implicit None
!
!     subroutine parameters
      Integer iiS,jjS,kkS,llS
      Integer nTInt
      Real*8  TInt(nTInt)
!
#include "ibas_ricd.fh"
      Integer::
     &        ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl,
     &        iBsInc,jBsInc,kBsInc,lBsInc,iPrInc,jPrInc,kPrInc,lPrInc,
     &        ipMem1, ipMem2, Mem1,Mem2

!     other local variables
      Integer iS_,jS_,kS_,lS_
      Real*8  Coor(3,4), Tmax
      Integer n
      Integer iTmp, Nr_of_D, nIJKL,
     &        ipDum , MemPrm
      Integer, External:: MemSO2
      Integer iAOst(4), iPrimi,jPrimj,kPrimk,lPriml,
     &        iBasi,jBasj,kBask,lBasl,
     &        iBasn,jBasn,kBasn,lBasn,
     &        nDCRR,nDCRS, ipTmp,
     &        mDij,mDik,mDjk,mDkl,mDil,mDjl,
     &        mDCRij,mDCRik,mDCRjk,mDCRkl,mDCRil,mDCRjl,
     &        ipDij,ipDik,ipDjk,ipDkl,ipDil,ipDjl, nZeta, nEta
      Integer nSO,iBasAO,jBasAO,kBasAO,lBasAO,
     &        iS,jS,kS,lS,ijS,klS,ikS,ilS,jkS,jlS
      Logical IJeqKL
      Integer iAngV(4),iCmpV(4), iShelV(4),iShllV(4),iAOV(4),iStabs(4),
     &        MemMax, kOp(4)
      Logical Shijij, NoInts
      Real*8, pointer:: SOInt(:), AOInt(:)
      Integer, External :: iDAMax_
      Integer ik2, jk2
!                                                                      *
!***********************************************************************
!                                                                      *
      Abstract Interface

      SubRoutine TwoEl(iS_,jS_,kS_,lS_,
     &           Coor,
     &           iAnga,iCmp,iShell,iShll,iAO,iAOst,
     &           NoInts, iStabs,
     &           nAlpha,iPrInc, nBeta,jPrInc,
     &           nGamma,kPrInc,nDelta,lPrInc,
     &           nData1,nData2,
     &           k2data1, k2data2,
     &           IJeqKL,kOp,
     &           Dij,mDij,mDCRij,Dkl,mDkl,mDCRkl,Dik,mDik,mDCRik,
     &           Dil,mDil,mDCRil,Djk,mDjk,mDCRjk,Djl,mDjl,mDCRjl,
     &           Coeff1,iBasi,Coeff2,jBasj,Coeff3,kBask,Coeff4,lBasl,
     &           FckTmp,nFT,nZeta,nEta,
     &           SOInt,nSOInt,Wrk,nWork2,
     &           Shijij,Aux,nAux)
      use k2_Structure, only: k2_type
      Implicit None
      Integer iS_,jS_,kS_,lS_
      Real*8  Coor(3,4)
      Integer iAnga(4), iCmp(4), iShell(4), iShll(4), iAO(4), iAOst(4)
      Logical NoInts
      Integer iStabs(4)
      Integer iPrInc, jPrInc, kPrInc, lPrInc
      Integer nData1, nData2
      Type(k2_type) k2data1(nData1), k2data2(nData2)
      Logical IJeqKL
      Integer kOp(4)
      Integer mDij,mDCRij,mDkl,mDCRkl,mDik,mDCRik,
     &        mDil,mDCRil,mDjk,mDCRjk,mDjl,mDCRjl
      Real*8  Dij(mDij,mDCRij),Dkl(mDkl,mDCRkl),Dik(mDik,mDCRik),
     &        Dil(mDil,mDCRil),Djk(mDjk,mDCRjk),Djl(mDjl,mDCRjl)
      Integer nAlpha, nBeta, nGamma, nDelta
      Integer iBasi, jBasj, kBask, lBasl
      Real*8 Coeff1(nAlpha,iBasi), Coeff2(nBeta,jBasj),
     &       Coeff3(nGamma,kBask), Coeff4(nDelta,lBasl)
      Integer nFT
      Real*8  FckTmp(nFT)
      Integer nZeta, nEta
      Integer nSOInt, nWork2
      Real*8  SOInt(iBasi*jBasj*kBask*lBasl,nSOInt)
      Real*8  Wrk(nWork2)
      Logical Shijij
      Integer nAux
      Real*8  Aux(nAux)


      End SubRoutine TwoEl
      End Interface

      Procedure(Twoel) :: TwoEl_NoSym, TwoEl_Sym
      Procedure(Twoel), pointer :: Do_TwoEl=>Null()
!                                                                      *
!***********************************************************************
!                                                                      *
!     Statement functions
!
      Integer i, j, iTri

      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
!                                                                      *
!***********************************************************************
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
         If (MemMax.gt.1000) MemMax=MemMax-1000
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
      End SubRoutine Eval_ijkl

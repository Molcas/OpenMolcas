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
* Copyright (C) 1990-1992,2000, Roland Lindh                           *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drvg1(Grad,Temp,nGrad)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will control the type of the two-electron integral, eg.     *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis        *
*          functions of the requested type.                            *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              SetUp_Ints                                              *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              Swap                                                    *
*              MemRg1                                                  *
*              PSOAO1                                                  *
*              PGet0                                                   *
*              TwoEl                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified for gradient calculation. January '92           *
*             Modified for SetUp_Ints. January '00                     *
************************************************************************
      use k2_setup
      use iSD_data
      use PSO_Stuff
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
      use Aces_Stuff, only: G_toc,nSSDM,SSDM
      Implicit Real*8 (A-H,O-Z)
      External Rsv_GTList
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "para_info.fh"
C#include "caspt2_grad.fh"
*     Local arrays
      Real*8  Coor(3,4), Grad(nGrad), Temp(nGrad)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4)
      Integer nHrrTb(0:iTabMx,0:iTabMx,2)
      Logical EQ, Shijij, AeqB, CeqD, lDummy,
     &        DoGrad, DoFock, Indexation,
     &        JfGrad(3,4), ABCDeq, No_Batch, Rsv_GTList,
     &        FreeK2, Verbose, Triangular
      Character Format*72
      Character*8 Method_chk
      Real*8, Allocatable:: TMax(:,:)
      Integer, Allocatable:: Ind_ij(:,:)
      Logical LoadVec ! For CASPT2 gradient
      Real*8, Allocatable :: CMOPT2(:)
      Integer, Allocatable :: iOffAO(:)
      Integer :: nOcc(8),nFro(8)
      Character*4096 RealName
      Real*8, Allocatable :: WRK1(:),WRK2(:)
************ columbus interface ****************************************
      Integer  Columbus
*
      Integer iSD4(0:nSD,4)
      save MemPrm
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      iPrint = 000000000

      iFnc(1)=0
      iFnc(2)=0
      iFnc(3)=0
      iFnc(4)=0
      PMax=Zero
      idum=0
      idum1=0
#ifdef _CD_TIMING_
      Twoel_CPU = 0.0d0
      Twoel_Wall = 0.0d0
      Pget_CPU = 0.0d0
      Pget_Wall = 0.0d0
#endif
      Call QEnter('Drvg1')
      call dcopy_(nGrad,[Zero],0,Temp,1)
*
      Call StatusLine(' Alaska:',' Computing 2-electron gradients')
*                                                                      *
************************************************************************
*                                                                      *
      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*-----Precompute k2 entities.
*
      Call Get_iScalar('Columbus',Columbus)
      Indexation=.False.
*     MP2 gradients:
      Call Get_cArray('Relax Method',Method_chk,8)
      If(Method_chk.eq.'MBPT2   ') Indexation=.True.
************ columbus interface ****************************************
* in order to access the half-sorted density matrix file
* some index arrays must be generated
      If (Columbus.eq.1) Indexation=.True.
      DoFock=.False.
      DoGrad=.True.
      ThrAO=Zero
      Call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      mSkal=nSkal
      nPairs=nSkal*(nSkal+1)/2
      nQuad =nPairs*(nPairs+1)/2
      Pren = Zero
      Prem = Zero
*                                                                      *
************************************************************************
*                                                                      *
*-----Prepare handling of two-particle density.
*
      !! integral_util/prepp.f
      Call PrepP
      If (Method_chk.eq.'CASPT2  ') Then
        nBasT = 0
        Do iIrrep = 0, nIrrep - 1
          nBasT = nBasT + nBas(iIrrep)
        End Do
        nSSDM = 0
C
        !! The two MO indices in the half-transformed amplitude are
        !! not CASSCF but quasi-canonical orbitals.
        Call mma_allocate(CMOPT2,nBasT*nBasT,Label='CMOPT2')
        Call PrgmTranslate('CMOPT2',RealName,lRealName)
        LuCMOPT2 = 61
C       Open (Unit=LuCMOPT2,
C    *        File=RealName(1:lRealName),
C    *        Status='OLD',
C    *        Form='UNFORMATTED')
C       call molcas_Open(LuCMOPT2,RealName(1:lRealName))
        Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
     &                        'DIRECT','UNFORMATTED',
     &                        iost,.FALSE.,
     &                        1,'OLD',is_error)
        Do i = 1, nBasT*nBasT
          Read (LuCMOPT2) CMOPT2(i)
        End Do
        Read (LuCMOPT2) nOcc(1)
        Read (LuCMOPT2) nOcc(2)
        Read (LuCMOPT2) nOcc(3)
        Read (LuCMOPT2) nOcc(4)
        Read (LuCMOPT2) nOcc(5)
        Read (LuCMOPT2) nOcc(6)
        Read (LuCMOPT2) nOcc(7)
        Read (LuCMOPT2) nOcc(8)
        Read (LuCMOPT2) nFro(1)
        Read (LuCMOPT2) nFro(2)
        Read (LuCMOPT2) nFro(3)
        Read (LuCMOPT2) nFro(4)
        Read (LuCMOPT2) nFro(5)
        Read (LuCMOPT2) nFro(6)
        Read (LuCMOPT2) nFro(7)
        Read (LuCMOPT2) nFro(8)
        Read (LuCMOPT2) nSSDM
C
C       read (lucmopt2) val
C       write(6,*) "val=",val
        If (nSSDM.ne.0) Then
          Call mma_allocate(SSDM,nBas(0)*(nBas(0)+1)/2,2,nSSDM,
     *                      Label='SSDM')
          Do iSSDM = 1, nSSDM
            Do i = 1, nBas(0)*(nBas(0)+1)/2
              Read (LuCMOPT2) SSDM(i,1,iSSDM),SSDM(i,2,iSSDM)
            End Do
          End Do
        End If
C       Read (LuCMOPT2) (nOcc(iSym),iSym=1,8)
C       Read (LuCMOPT2) (nFro(iSym),iSym=1,8)
        Close (LuCMOPT2)
        write(6,*) "Number of Non-Frozen Occupied Orbitals = ", nOcc(1)
        write(6,*) "Number of     Frozen          Orbitals = ", nFro(1)
C
        Call mma_allocate(iOffAO,nSkal+1,Label='iOffAO')
        MaxShlAO = 0
        iOffAO(1) = 0
        Do iSh = 1, nSkal
          nBasI = iSD(2,iSh)*iSD(3,iSh)
          If (nBasI.gt.MaxShlAO) MaxShlAO = nBasI
          iOffAO(iSh+1) = iOffAO(iSh)+nBasI
C         write(6,*) ish,nbasi
        End Do
C       write(6,*) "ioffao"
C       do i = 1, nskal
C         write(6,*) i,ioffao(i)
C       end do
        Call mma_allocate(G_toc,MaxShlAO**4,Label='GtocCASPT2')
C
C       nOcc = 10 ! nIsh(1) + nAsh(1)
        Call PrgmTranslate('GAMMA',RealName,lRealName)
        LuGamma = 60
C       Open (Unit=LuGamma,
C    *        File=RealName(1:lRealName),
C    *        Status='OLD',
C    *        Form='UNFORMATTED',
C    *        Access='DIRECT',
C    *        Recl=nOcc(1)*nOcc(1)*8)
C       call molcas_Open(LuGamma,RealName(1:lRealName))
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                        'DIRECT','UNFORMATTED',
     &                        iost,.TRUE.,
     &                        nOcc(1)*nOcc(1)*8,'OLD',is_error)
*
        Call mma_allocate(WRK1,nOcc(1)*nOcc(1),Label='WRK1')
        Call mma_allocate(WRK2,MaxShlAO*nOcc(1),Label='WRK2')
      End If

*                                                                      *
************************************************************************
*                                                                      *
      MxPrm = 0
      Do iAng = 0, iAngMx
         MxPrm = Max(MxPrm,MaxPrm(iAng))
      End Do
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
*
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax,nSkal,nSkal,Label='TMax')
      Call Shell_MxSchwz(nSkal,TMax)
      TMax_all=Zero
      Do iS = 1, nSkal
         Do jS = 1, iS
            TMax_all=Max(TMax_all,TMax(iS,jS))
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Create list of non-vanishing pairs
*
      Call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,Label='Ind_ij')
      nij=0
      Do iS = 1, nSkal
         Do jS = 1, iS
            If (TMax_All*TMax(iS,jS).ge.CutInt
     *        .or.Method_chk.eq.'CASPT2  ') Then
               nij = nij + 1
               Ind_ij(1,nij)=iS
               Ind_ij(2,nij)=jS
            End If
         End Do
      End Do
      !! For efficient back-transformation of the amplitude,
      !! CASPT2 gradient uses a different shell indexing.
C     If (Method_chk.ne.'CASPT2  ') Then
        P_Eff=Dble(nij)
C     Else
C       P_Eff=nSkal*(nSkal+1)/2
C     End If
*                                                                      *
************************************************************************
*                                                                      *
*-------Compute FLOPs for the transfer equation.
*
        Do iAng = 0, iAngMx
           Do jAng = 0, iAng
              nHrrab = 0
              Do i = 0, iAng+1
                 Do j = 0, jAng+1
                    If (i+j.le.iAng+jAng+1) Then
                       ijMax = Min(iAng,jAng)+1
                       nHrrab = nHrrab + ijMax*2+1
                    End If
                 End Do
              End Do
              nHrrTb(iAng,jAng,1)=nHrrab
              nHrrTb(jAng,iAng,1)=nHrrab
           End Do
        End Do

*                                                                      *
************************************************************************
*                                                                      *
      Triangular=.True.
      Call Alloc_TList(Triangular,P_Eff)
      Call Init_TList(Triangular,P_Eff)
      Call Init_PPList
      Call Init_GTList
      iOpt=0
*                                                                      *
************************************************************************
*                                                                      *
*     In MPP case dispatch one processor to do 1-el gradients first
*
      If (nProcs.gt.1.and.King()) Then
            Call Drvh1(Grad,Temp,nGrad)
*        If (nPrint(1).ge.15)
*    &   Call PrGrad(' Gradient excluding two-electron contribution',
*    &               Grad,lDisp(0),lIrrep,ChDisp,5)
         call dcopy_(nGrad,[Zero],0,Temp,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      ipMem1 = 1
*                                                                      *
************************************************************************
*                                                                      *
*     big loop over individual tasks, distributed over individual nodes
   10 Continue
*     make reservation of a task on global task list and get task range
*     in return. Function will be false if no more tasks to execute.
      If (.Not.Rsv_GTList(TskLw,TskHi,iOpt,lDummy)) Go To 11
*
*     Now do a quadruple loop over shells
*
      Call Get_cArray('Relax Method',Method_chk,8)
      If (Method_chk.ne.'CASPT2  ') Then
        ijS = Int((One+sqrt(Eight*TskLw-Three))/Two)
        iS = Ind_ij(1,ijS)
        jS = Ind_ij(2,ijS)
        klS = Int(TskLw-DBLE(ijS)*(DBLE(ijS)-One)/Two)
        kS = Ind_ij(1,klS)
        lS = Ind_ij(2,klS)
      Else
        iS = 1
        jS = 1
        kS = 1
        lS = 1
        !! proceed the index
        Do iCount = 1, int(TskLw)-1
          Call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
        End Do
        Count=dble(iCount)
        !! If LoadVec is true, a new vector of the half-transformed
        !! T-amplitude is read. In the first loop, it is always true.
        !! In other loops, a new vector is read only when I- and K-th
        !! are different from the previous loop.
        !! The half back-transformation, T_{ij}^{ab} ->
        !! T_{ij}^{rho sigma}, is done somewhere in CASPT2.
        !! rho and sigma correspond to either I- or K-th shells.
        !! Occupied orbital indices (correspond to J- or L-th shells)
        !! are back-transformed on-the-fly.
        LoadVec = .True.
      End If
C     End If
      Count=TskLw
      Call CWTime(TCpu1,TWall1)
  13  Continue
*
         Aint=TMax(iS,jS)*TMax(kS,lS)
         If (AInt.lt.CutInt) Go To 14
         If (iPrint.ge.15) Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
C        write(6,*) "-------------------------"
C        Write (6,'("iS,jS,kS,lS=",4i8)') iS,jS,kS,lS
*
************************************************************************
*                                                                      *
         !! integral_util/gen_isd4.f
         !! construct jQuad and iSD4
         !! the latter is copied from iSD
         !! The definition of iSD is in integral_util/def_shells.f
         Call Gen_iSD4(iS, jS, kS, lS,iSD,nSD,iSD4)
         !! alaska_util/size_sob_g.f
         !! nSO: Number of AOs of the shell
         Call Size_SO_block_g(iSD4,nSD,Petite,nSO,No_batch)
C        write(6,*) "nso = ", nso
         If (No_batch) Go To 140
*
         !! alaska_util/int_prep_g.f
         Call Int_Prep_g(iSD4,nSD,Coor,Shijij,iAOV,iStabs)

*
*                                                                      *
************************************************************************
*                                                                      *
*       --------> Memory Managment <--------
*
*        Compute memory request for the primitives, i.e.
*        how much memory is needed up to the transfer
*        equation.
*
         Call MemRys_g(iSD4,nSD,nRys,MemPrm)
*                                                                      *
************************************************************************
*                                                                      *
         !! coordinate (atom to be differentiated) is the same?
         ABCDeq=EQ(Coor(1,1),Coor(1,2)) .and.
     &          EQ(Coor(1,1),Coor(1,3)) .and.
     &          EQ(Coor(1,1),Coor(1,4))
         !! sum of the angular momentum
         ijklA=iSD4(1,1)+iSD4(1,2)
     &        +iSD4(1,3)+iSD4(1,4)
         If (nIrrep.eq.1.and.ABCDeq.and.Mod(ijklA,2).eq.1)
     &      Go To 140
*                                                                      *
************************************************************************
*                                                                      *
*        Decide on the partioning of the shells based on the
*        available memory and the requested memory.
*
*        Now check if all blocks can be computed and stored at
*        once.
*
         !! alaska_util/soao_g.f and alaska_util/psoao1.f
         Call SOAO_g(iSD4,nSD,nSO,
     &               MemPrm, MemMax,
     &               nExp,nBasis,MxShll,
     &               iBsInc,jBsInc,kBsInc,lBsInc,
     &               iPrInc,jPrInc,kPrInc,lPrInc,
     &               ipMem1,ipMem2, Mem1,  Mem2,
     &               iPrint,iFnc, MemPSO)
         !! this is the number of contracted functions
         iBasi    = iSD4(3,1)
         jBasj    = iSD4(3,2)
         kBask    = iSD4(3,3)
         lBasl    = iSD4(3,4)
C     write(6,'("iBasi,jBasj,kBask,lBasl = ",4i5)')
C    *   ibasi,jbasj,kbask,lbasl
C     write(6,'("iBasInc,jBasInc,kBasInc,lBasInc = ",4i5)')
C    *   ibasInc,jbasInc,kbasInc,lbasInc
*                                                                      *
************************************************************************
*                                                                      *
         !! iAnga : angular momentum of the shell
         !! iCmpa : number of AOs of the shell
         !! iShlla: unique shell index (?)
         !! iShl  : equal to iS,jS,kS,lS?
         Call Int_Parm_g(iSD4,nSD,iAnga,
     &                 iCmpa,iShlla,iShela,
     &                 iPrimi,jPrimj,kPrimk,lPriml,
     &                 ipCffi,jpCffj,kpCffk,lpCffl,
     &                 nExp,ipExp,ipCff,MxShll,
     &                 indij,k2ij,nDCRR,k2kl,nDCRS,
     &                 mdci,mdcj,mdck,mdcl,AeqB,CeqD,
     &                 nZeta,nEta,ipZeta,ipZI,
     &                 ipP,ipEta,ipEI,ipQ,ipiZet,ipiEta,
     &                 ipxA,ipxB,ipxG,ipxD,l2DI,nab,nHmab,ncd,nHmcd,
     &                 nIrrep)
*                                                                      *
************************************************************************
*                                                                      *
*        Scramble arrays (follow angular index)
*

         Do iCar = 1, 3
            Do iSh = 1, 4
               JndGrd(iCar,iSh) = iSD4(15+iCar,iSh)
               If (iAnd(iSD4(15,iSh),2**(iCar-1)) .eq.
     &             2**(iCar-1)) Then
                  JfGrad(iCar,iSh) = .True.
               Else
                  JfGrad(iCar,iSh) = .False.
               End If
            End Do
         End Do
*
         Do 400 iBasAO = 1, iBasi, iBsInc
           iBasn=Min(iBsInc,iBasi-iBasAO+1)
           iAOst(1) = iBasAO-1
         Do 410 jBasAO = 1, jBasj, jBsInc
           jBasn=Min(jBsInc,jBasj-jBasAO+1)
           iAOst(2) = jBasAO-1
         Do 420 kBasAO = 1, kBask, kBsInc
           kBasn=Min(kBsInc,kBask-kBasAO+1)
           iAOst(3) = kBasAO-1
         Do 430 lBasAO = 1, lBasl, lBsInc
           lBasn=Min(lBsInc,lBasl-lBasAO+1)
           iAOst(4) = lBasAO-1
*
*----------Get the 2nd order density matrix in SO basis.
*
C       write(6,'("iBasn,jBasn,kBasn,lBasn=",4i4)')
C    *    ibasn,jbasn,kbasn,lbasn
C       write(6,'("iFnc(1),iFnc(2),iFnc(3),iFnc(4)=",4i4)')
C    *    iFnc(1),iFnc(2),iFnc(3),iFnc(4)
           nijkl = iBasn*jBasn*kBasn*lBasn


*     Fetch the T_i,j,kappa, lambda corresponding to
*     kappa = k, lambda = l
*
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU1,Pget0WALL1)
#endif
           !! get the density
           !! integral_util/pget0.f
           !! PSO(nIJKL,nSO) = Sew_Scr(ipMem1) (2-body AO density?)
C          write(6,*) "nIJKL, nSO = ", nIJKL, nSO
C         write(6,'(4i3)') ifnc(1)*ibasn,ifnc(2)*jbasn,
C    *                      ifnc(3)*kbasn,ifnc(4)*lbasn
           If (Method_chk.eq.'CASPT2  ') Then
C            write(6,*) "ikS,jlS = ", ikS,jls
C            if (ikS.ge.jlS) then
C            if (jS.ne.lS) then
              Call CASPT2_BTAMP(iS,jS,kS,lS,
     *                         iFnc(1)*iBasn,iFnc(2)*jBasn,
     *                         iFnc(3)*kBasn,iFnc(4)*lBasn,
     *                         iOffAO,nBasT,
     *                         nOcc(1),CMOPT2(1+nbast*nfro(1)),
     *                         WRK1,WRK2,G_Toc)
C            else
C            Call CASPT2_BTAMP(jS,iS,lS,kS,
C    *                         iFnc(2)*jBasn,iFnc(1)*iBasn,
C    *                         iFnc(4)*lBasn,iFnc(3)*kBasn,
C    *                         iOffAO,nBasT,
C    *                         CMOPT2,WRK1,WRK2,G_Toc)
C            end if
           End If
           Call PGet0(iCmpa,iShela,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
           If (AInt*PMax.lt.CutInt) Go To 430
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU2,Pget0WALL2)
           Pget_CPU = Pget_CPU + Pget0CPU2-Pget0CPU1
           Pget_Wall = Pget_Wall + Pget0WALL2-Pget0WALL1
#endif
           If (AInt*PMax.lt.CutInt) Go To 430
*
*----------Compute gradients of shell quadruplet
*
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU1,TwoelWall1) ! timing_cdscf
#endif
C            write(6,'("shell(i,j,k,l) = ", 4i3)') is,js,ks,ls
C      call dcopy_(ngrad,0.0d+00,0,temp,1)
           Call TwoEl_g(Coor,
     &          iAnga,iCmpa,iShela,iShlla,iAOV,
     &          mdci,mdcj,mdck,mdcl,nRys,
     &          Data_k2(k2ij),nab,nHmab,nDCRR,
     &          Data_k2(k2kl),ncd,nHmcd,nDCRS,Pren,Prem,
     &          iPrimi,iPrInc,jPrimj,jPrInc,
     &          kPrimk,kPrInc,lPriml,lPrInc,
     &          Work(ipCffi+(iBasAO-1)*iPrimi),iBasn,
     &          Work(jpCffj+(jBasAO-1)*jPrimj),jBasn,
     &          Work(kpCffk+(kBasAO-1)*kPrimk),kBasn,
     &          Work(lpCffl+(lBasAO-1)*lPriml),lBasn,
     &          Mem_DBLE(ipZeta),Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,
     &          Mem_DBLE(ipEta), Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,
     &          Mem_DBLE(ipxA),Mem_DBLE(ipxB),
     &          Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,
     &          JfGrad,JndGrd,Sew_Scr(ipMem1), nSO,Sew_Scr(ipMem2),Mem2,
     &          Aux,nAux,Shijij)
C             do i = 1, 4
C             write(6,'(i3,3f20.10)') i,(temp(j+3*(i-1)),j=1,3)
C             end do
C             write(6,*) "---------------------------"
C     call cpu_time(t8)
C     ttt(6) = ttt(6) + t8-t7
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU2,TwoelWall2)
           Twoel_CPU = Twoel_CPU + TwoelCPU2-TwoelCPU1
           Twoel_Wall = Twoel_Wall + TwoelWall2-TwoelWall1
#endif
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg1: Grad',
     &                  Temp,nGrad,lIrrep,ChDisp,5)
C              Call PrGrad(' In Drvg1: Grad',
C    *                  Temp,nGrad,lIrrep,ChDisp,5)
*
 430     Continue
 420     Continue
*
 410     Continue
 400     Continue
*
 140     Continue
*
 14      Continue
         Count=Count+One
         If (Count-TskHi.gt.1.0D-10) Go To 12
         if (Method_chk.ne.'CASPT2  ') Then
           klS = klS + 1
           If (klS.gt.ijS) Then
              ijS = ijS + 1
              klS = 1
           End If
           iS = Ind_ij(1,ijS)
           jS = Ind_ij(2,ijS)
           kS = Ind_ij(1,klS)
           lS = Ind_ij(2,klS)
         Else If (Method_chk.eq.'CASPT2  ') Then
           Call CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
           ikS = Max(iS,kS)*(Max(iS,kS)-1)/2+Min(iS,kS)
           jlS = Max(jS,lS)*(Max(jS,lS)-1)/2+Min(jS,lS)
         End If
         Go To 13
*
*     Task endpoint
 12   Continue
      Call CWTime(TCpu2,TWall2)
      Call SavTim(4,TCpu2-TCpu1,TWall2-Twall1)
      Call SavStat(1,One,'+')
      Call SavStat(2,TskHi-TskLw+One,'+')
      Go To 10
 11   Continue
*     End of big task loop
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _MOLCAS_MPP_
      Call GADGOP(Temp,nGrad,'+')
#endif
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Sew_Scr)
      Call Free_GTList
      Call Free_PPList
      Call Free_TList
      Call mma_deallocate(Ind_ij)
      Call mma_deallocate(TMax)
      If (Method_chk.eq.'CASPT2  ') Then
        Close (LuGamma)
        Call mma_deallocate(iOffAO)
        Call mma_deallocate(CMOPT2)
        If (nSSDM.ne.0) Call mma_deallocate(SSDM)
        Call mma_deallocate(WRK1)
        Call mma_deallocate(WRK2)
      End If
*                                                                      *
************************************************************************
*                                                                      *
C       write(6,*) "finishing?"
      Call CloseP
C       write(6,*) "closep finished"
#ifdef _CD_TIMING_
      Drvg1_CPU = TCpu2-TCpu1
      Drvg1_Wall= TWall2-TWall1
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Prepp:'
      Write(6,*) 'Wall/CPU',Prepp_Wall, Prepp_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Pget:'
      Write(6,*) 'Wall/CPU',Pget_Wall, Pget_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Time spent in Drvg1:'
      Write(6,*) 'Wall/CPU',Drvg1_Wall, Drvg1_CPU
      Write(6,*) '-------------------------'
      Total_Dens_Wall = Prepp_Wall+Pget_Wall
      Total_Dens_CPU = Prepp_CPU+Pget_CPU
      Total_Der_Wall = Drvg1_Wall - Total_Dens_Wall
      Total_Der_CPU = Drvg1_CPU - Total_Dens_CPU
      Total_Der_Wall2 = TwoEl_Wall
      Total_Der_CPU2 = TwoEl_CPU

      Write(6,*) 'Total Time for Density:'
      Write(6,*) 'Wall/CPU',Total_Dens_Wall, Total_Dens_CPU
      Write(6,*) '-------------------------'
      Write(6,*) 'Total Time for Derivatives:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall2, Total_Der_CPU2
      Write(6,*) '-------------------------'
      Write(6,*) 'Derivative check:'
      Write(6,*) 'Wall/CPU',Total_Der_Wall, Total_Der_CPU
      Write(6,*) '-------------------------'
#endif
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*                                                                      *
************************************************************************
*                                                                      *
*
      Call Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)
*
      iPren=3+Max(1,Int(Log10(Pren+0.001D+00)))
      iPrem=3+Max(1,Int(Log10(Prem+0.001D+00)))
      Write (Format,'(A,I2,A,I2,A)') '(A,F',iPren,
     &           '.0,A,F',iPrem,'.0,A)'
      If (iPrint.ge.6) Then
      Write (6,Format)
     &   ' A total of', Pren,' entities were prescreened and',
     &                  Prem,' were kept.'
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_iSD()
      Call QExit('Drvg1')
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine CASPT2_Grad_FwdCnt(iS,jS,kS,lS,LoadVec)
*
      Implicit Real*8 (A-H,O-Z)
*
      Logical LoadVec
*
      LoadVec = .False.
      iSold = iS
      kSold = kS
*
      lS = lS + 1
      If (iS.eq.kS) Then
        If (lS.gt.jS) Then
          jS = jS + 1
          lS = 1
        End If
        If (jS.gt.iS) Then
          iS = iS + 1
          jS = 1
          kS = 1
          lS = 1
        End If
      Else
        If (lS.gt.kS) Then
          jS = jS + 1
          lS = 1
        End If
        If (jS.gt.iS) Then
          kS = kS + 1
          jS = 1
          lS = 1
        End If
        If (kS.gt.iS) Then
          iS = iS + 1
          jS = 1
          kS = 1
          lS = 1
        End If
      End IF
*
      If (iSold.ne.iS .or. kSold.ne.kS) LoadVec = .True.
*
      Return
*
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine CASPT2_BTAMP(iS,jS,kS,lS,nBasI,nBasJ,nBasK,nBasL,
     *                        iOffAO,nBasT,nOcc,
     *                        CMOPT2,WRK1,WRK2,G_Toc)
*
      Implicit Real*8 (A-H,O-Z)
*
#include "nsd.fh"
#include "etwas.fh"
#include "stdalloc.fh"
*
      Dimension iOffAO(*),CMOPT2(*),
C    *          G_toc(MaxShlAO,MaxShlAO,MaxShlAO,MaxShlAO)
     *          G_Toc(*)
      Dimension WRK1(*),WRK2(*)
*
*     Transform T_{ij}^{rho sigma} to T_{mu nu}^{rho sigma}
*     i- and k-shells correspond to rho and sigma (order?)
*     The transformed amplitude is stored as D_{j,l,i,k},
*     and it will be sorted in integral_util/pget3.f
*
*     It is possible to reduce operations. The third transformation
*     can be postponed until j-shell varies.
*
      LuGamma = 60
C     nOcc = 10 ! nIsh(1) + nAsh(1)
*
C     nBasI = iSD4(2,1)*iSD4(3,1)
C     nBasJ = iSD4(2,2)*iSD4(3,2)
C     nBasK = iSD4(2,3)*iSD4(3,3)
C     nBasL = iSD4(2,4)*iSD4(3,4)
C         write(6,'(4i3)') nbasi,nbasj,nbask,nbasl
C     write(6,'(4i3)') (isd4(2,i)*isd4(3,i),i=1,4)
*
      scal = 0.1250d+00!*0.5d+00*0.5d+00
C     scal = 0.250d+00
C     scal = 0.500d+00
C     scal = 0.1357D+00 !!! aaaaa
*
      Do kBas0 = 1, nBasK
        kBas = iOffAO(kS) + kBas0
        Do iBas0 = 1, nBasI
          iBas = iOffAO(iS) + iBas0
          iRec = iBas + nBasT*(kBas-1)
          !! Read the half-transformed-amplitude
C         write(6,*) "going to read"
          Read (Unit=LuGamma,Rec=iRec) (WRK1(i),i=1,nOcc*nOcc)
C         if (irec.eq.1) then
C         write(6,*) "the first hbt amplitude"
C          write(6,*) "ibas0,kbas0=",ibas0,kbas0
C          write(6,*) "ibas ,kbas =",ibas,kbas
C           call sqprt(wrk1,nocc)
C           call abend
C         end if
          !! do the remaining (third and fourth) transformation
          Call DGemm_('N','N',nBasJ,nOcc,nOcc,
     *                1.0D+00,CMOPT2(1+iOffAO(jS)),nBasT,WRK1,nOcc,
     *                0.0D+00,WRK2,nBasJ)
C     if (is.ne.ks)
C     call dgemm_('n','t',nbasj,nocc,nocc,
C    *            1.0d+00,cmopt2(1+ioffao(js)),nbast,wrk1,nocc,
C    *            1.0d+00,wrk2,nbasj)
          Loc = nBasJ*nBasL*(iBas0-1+nBasI*(kBas0-1))
          Call DGemm_('N','T',nBasJ,nBasL,nOcc,
     *                SCAL   ,WRK2,nBasJ,CMOPT2(1+iOffAO(lS)),nBasT,
     *                0.0D+00,G_toc(1+Loc),nBasJ)
C      if (ibas.eq.6.and.kbas.eq.1) then
C        write(6,*) "show the AO transformed for ibas,kbas=6,1"
C        call sqprt(g_toc(1+loc),nbasj)
C      end if
C         if (is.ne.ks) then
C         Call DGemm_('N','T',nBasJ,nOcc,nOcc,
C    *                1.0D+00,CMOPT2(1+iOffAO(jS)),nBasT,WRK1,nOcc,
C    *                0.0D+00,WRK2,nBasJ)
C         Call DGemm_('N','T',nBasJ,nBasL,nOcc,
C    *           SCAL*0.5d+00   ,WRK2,nBasJ,CMOPT2(1+iOffAO(lS)),nBasT,
C    *                0.5D+00,G_toc(1+Loc),nBasJ)
C         end if
C      write(6,*) "ibas,kbas = ", ibas0,kbas0
C      do i = 1, nbasj
C      do j = 1, nbasl
C      write(6,'(2i3,f20.10)') i,j,g_toc(i+nbasj*(j-1)+loc)
C      end do
C      end do
        End Do
      End Do
C     write(6,*) "after transformation"
C     do kbas0 = 1, nbask
C       kbas = ioffao(ks)+kbas0
C       do ibas0 = 1, nbasi
C         ibas = ioffao(is)+ibas0
C         do lbas0 = 1, nbasl
C           lbas = ioffao(ls)+lbas0
C           do jbas0 = 1, nbasj
C             jbas = ioffao(js)+jbas0
C     loc = jbas0-1+nbasj*(lbas0-1+nbasl*(ibas0-1+nbasi*(kbas0-1)))
C      write(6,'("asdf",4i3,f20.10)')
C    *  ibas,jbas,kbas,lbas,g_toc(1+loc)*8.0d+00
C           end do
C         end do
C       end do
C     end do
C
C     if (is.eq.ks .and. js.eq.ls) return
C     If (jS.lt.iS.and.lS.lt.iS) Then
C     If (jS.eq.lS) Then
C     call dcopy_(maxshlao**4,0.0d+00,0,g_toc,1)
      Do lBas0 = 1, nBasL
        lBas = iOffAO(lS) + lBas0
        Do jBas0 = 1, nBasJ
          jBas = iOffAO(jS) + jBas0
          iRec = jBas + nBasT*(lBas-1)
          !! Read the half-transformed-amplitude
          Read (Unit=LuGamma,Rec=iRec) (WRK1(i),i=1,nOcc*nOcc)
          !! do the remaining (third and fourth) transformation
          Call DGemm_('N','N',nBasI,nOcc,nOcc,
     *                1.0D+00,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,
     *                0.0D+00,WRK2,nBasI)
C     if (js.ne.ls)
C     call dgemm_('n','t',nbasi,nocc,nocc,
C    *            1.0d+00,cmopt2(1+ioffao(is)),nbast,wrk1,nocc,
C    *            1.0d+00,wrk2,nbasi)
          Call DGemm_('N','T',nBasI,nBasK,nOcc,
     *                SCAL   ,WRK2,nBasI,CMOPT2(1+iOffAO(kS)),nBasT,
     *                0.0D+00,WRK1,nBasI)
          Do kBas0 = 1, nBasK
          kbas = ioffao(ks)+kbas0
            Do iBas0 = 1, nBasI
          ibas = ioffao(is)+ibas0
              Loc = jBas0-1 + nBasJ*(lBas0-1
     *            + nBasL*(iBas0-1 + nBasI*(kBas0-1)))
C      write(6,'("asdf",4i3,f20.10)')
C    *  ibas,jbas,kbas,lbas,wrk1(ibas0+nbasi*(kbas0-1))*8.0d+00
              G_toc(1+Loc) = G_toc(1+Loc) + WRK1(iBas0+nBasI*(kBas0-1))
            End Do
          End Do
        End Do
      End Do
C     End If
C     call dscal_(maxshlao**4,0.5d+00,g_toc,1)
      Do lBas0 = 1, nBasL
        lBas = iOffAO(lS) + lBas0
        Do iBas0 = 1, nBasI
          iBas = iOffAO(iS) + iBas0
          iRec = iBas + nBasT*(lBas-1)
          !! Read the half-transformed-amplitude
C         write(6,*) "going to read"
          Read (Unit=LuGamma,Rec=iRec) (WRK1(i),i=1,nOcc*nOcc)
          !! do the remaining (third and fourth) transformation
          Call DGemm_('N','N',nBasJ,nOcc,nOcc,
     *                1.0D+00,CMOPT2(1+iOffAO(jS)),nBasT,WRK1,nOcc,
     *                0.0D+00,WRK2,nBasJ)
          Call DGemm_('N','T',nBasJ,nBasK,nOcc,
     *                SCAL   ,WRK2,nBasJ,CMOPT2(1+iOffAO(kS)),nBasT,
     *                0.0D+00,WRK1,nBasJ)
          Do kBas0 = 1, nBasK
            Do jBas0 = 1, nBasJ
              Loc = jBas0-1 + nBasJ*(lBas0-1
     *            + nBasL*(iBas0-1 + nBasI*(kBas0-1)))
              G_toc(1+Loc) = G_toc(1+Loc) + WRK1(jBas0+nBasJ*(kBas0-1))
            End Do
          End Do
        End Do
      End Do
*
      Do kBas0 = 1, nBasK
        kBas = iOffAO(kS) + kBas0
        Do jBas0 = 1, nBasJ
          jBas = iOffAO(jS) + jBas0
          If (jBas.ge.kBas) Then
            iRec = jBas + nBasT*(kBas-1)
          Else
            iRec = kBas + nBasT*(jBas-1)
          End If
          !! Read the half-transformed-amplitude
C         write(6,*) "going to read"
          Read (Unit=LuGamma,Rec=iRec) (WRK1(i),i=1,nOcc*nOcc)
          !! do the remaining (third and fourth) transformation
          if (jbas.ge.kbas) then
          Call DGemm_('N','N',nBasI,nOcc,nOcc,
     *                1.0D+00,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,
     *                0.0D+00,WRK2,nBasI)
          else
          Call DGemm_('N','T',nBasI,nOcc,nOcc,
     *                1.0D+00,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,
     *                0.0D+00,WRK2,nBasI)
          end if
          Call DGemm_('N','T',nBasI,nBasL,nOcc,
     *                SCAL   ,WRK2,nBasI,CMOPT2(1+iOffAO(lS)),nBasT,
     *                0.0D+00,WRK1,nBasI)
          Do lBas0 = 1, nBasL
            Do iBas0 = 1, nBasI
              Loc = jBas0-1 + nBasJ*(lBas0-1
     *            + nBasL*(iBas0-1 + nBasI*(kBas0-1)))
              G_toc(1+Loc) = G_toc(1+Loc) + WRK1(iBas0+nBasI*(lBas0-1))
            End Do
          End Do
        End Do
      End Do
*
      Return
*
      End

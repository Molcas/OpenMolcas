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
* Copyright (C) 1990,1991,1993,1998,2006,2007, Roland Lindh            *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drv2El_3Center_RI(Integral_WrOut,ThrAO)
************************************************************************
*                                                                      *
*  Object: driver for the 3 center integrals in the RI approach.       *
*                                                                      *
*          This code has three sections                                *
*          1) a 2-center section to generate the Q-vectors             *
*          2) a 3-center section to generate the R-vectors             *
*          3) a partial transpose section to generate the RI vectors   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified to minimize overhead for calculations with      *
*             small basis sets and large molecules. Sept. '93          *
*             Modified driver. Jan. '98                                *
*             Modified to 3-center ERIs for RI Jan '06                 *
*             Modified to out-of-core version Feb '07                  *
************************************************************************
      use iSD_data
      use Wrj12
      use Basis_Info, only: dbsc, nBas, nBas_Aux
      use Temporary_Parameters, only: force_out_of_core
      use Real_Info, only: CutInt
      use RICD_Info, only: LDF
      use Symmetry_Info, only: nIrrep
      use j12
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut, Rsv_Tsk
#include "Molcas.fh"
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
*
#include "lRI.fh"
#include "setup.fh"
#include "nsd.fh"
#define _no_nShs_
#include "iTOffs.fh"
*
      Character Name_R*6
      Integer iOff_3C(3,0:7), Lu_R(0:7), iAddr_R(0:7), iMax_R(2,0:7),
     &        iTtmp(0:7), NoChoVec(0:7), iOff_Rv(0:7)
      Logical Verbose, Indexation, FreeK2,
     &        DoGrad, DoFock, Out_of_Core, Rsv_Tsk, Reduce_Prt
      External Reduce_Prt

      Real*8, Allocatable :: A_Diag(:), Local_A(:,:)
      Integer, Allocatable :: SO2C(:), AB(:,:)

      Real*8, Allocatable :: Tmp(:,:), TMax_Auxiliary(:),
     &                       TMax_Valence(:,:)
      Integer, Allocatable:: TmpList(:), iRv(:), LBList(:)
      Integer, Allocatable:: Addr(:), NuMu(:,:)
      Real*8, Allocatable :: Arr_3C(:), Rv(:), Qv(:), Diag(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface

      SubRoutine Drv2El_2Center_RI(ThrAO,A_Diag,nSO_Aux,MaxCntr,SO2C)
      Real*8 ThrAO
      Real*8, Allocatable :: A_Diag(:)
      Integer, Allocatable :: SO2C(:)
      Integer nSO_Aux, MaxCntr
      End SubRoutine Drv2El_2Center_RI

      SubRoutine Post_2Center_LDF(A_Diag,AB,MaxCntr,Lu_AB,Local_A,
     &                            SO2C,nSO_Aux)
      Real*8, Allocatable :: A_Diag(:), Local_A(:,:)
      Integer, Allocatable:: SO2C(:), AB(:,:)
      Integer MaxCntr,Lu_AB,nSO_Aux
      End SubRoutine Post_2Center_LDF

      SubRoutine Post_2Center_RI(A_Diag)
      Real*8, Allocatable :: A_Diag(:)
      End SubRoutine Post_2Center_RI

      End Interface
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
*
*     Get global print level
*
      iPL=iPrintLevel(-1)
      If (iPL.eq.2) Then
         iPL=5
      Else If (iPL.eq.3) Then
         iPL=6
      Else If (iPL.eq.4) Then
         iPL=99
      Else If (iPL.eq.5) Then
         iPL=99
      End If
      nPrint(iRout)=iPL
*
*     Reduce print level if iterating
*
      If (Reduce_Prt().and.iPL.le.5) Then
         nPrint(iRout)=4
      End If
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.6) Call CWTime(TC0,TW0)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     2 - C E N T E R   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Compute the two-center integrals over the auxiliary basis
*
      Call Drv2El_2Center_RI(ThrAO,A_Diag,nSO_Aux,MaxCntr,SO2C)
*
*     Post processing to generate the Q-vectors.
*
      If (LDF) Then
*
*        Local RI
*
         Call Post_2Center_LDF(A_Diag,AB,MaxCntr,Lu_AB,Local_A,
     &                         SO2C,nSO_Aux)
*
      Else
*
*        Standard RI
*
         Call Post_2Center_RI(A_Diag)
*
      End If
*
      Call Set_Basis_Mode('Auxiliary')
      Call Nr_Shells(nSkal_Auxiliary)
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' 2-center integrals:'
         Call CWTime(TC1,TW1)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC1-TC0,' sec.',
     &                  '      Wall time:',TW1-TW0,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     3 - C E N T E R   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      Call StatusLine(' Seward:',' Computing 3-center RI integrals')
*
*     Handle both the valence and the auxiliary basis set
*
      Call Set_Basis_Mode('WithAuxiliary')
      Call SetUp_iSD
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do not generate
*     tables for indexation.
*
      Indexation = .False.
      DoGrad=.False.
      DoFock=.False.
      Call Setup_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      nSkal_Valence=nSkal-nSkal_Auxiliary
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax_Valence,nSkal_Valence,nSkal_Valence,
     &                  Label='TMax_Valence')
      Call mma_allocate(TMax_Auxiliary,nSkal_Auxiliary,
     &                  Label='TMax_Auxiliary')
*
      Call mma_allocate(Tmp,nSkal,nSkal,Label='Tmp')
      Call Shell_MxSchwz(nSkal,Tmp)
      TMax_all=Zero
      Do iS = 1, nSkal_Valence
         Do jS = 1, iS
            TMax_Valence(iS,jS)=Tmp(iS,jS)
            TMax_Valence(jS,iS)=Tmp(iS,jS)
            TMax_all=Max(TMax_all,Tmp(iS,jS))
         End Do
      End Do
      Do iS = 1, nSkal_Auxiliary-1
         iS_ = iS + nSkal_Valence
         jS_ = nSkal_Valence + nSkal_Auxiliary
         TMax_Auxiliary(iS)=Tmp(jS_,iS_)
         TMax_all=Max(TMax_all,Tmp(jS_,iS_))
      End Do
*
      Call mma_deallocate(Tmp)
*                                                                      *
************************************************************************
*                                                                      *
*     Set up indexation for Gaussian pairs.
*
*     Generate some offsets and dimensions for the J12 matrix and
*     the RI vectors.
*
      Call Setup_Aux(nIrrep,nBas,nSkal_Valence,nSkal_Auxiliary,nSO,
     &               TMax_Valence,CutInt,nSkal2,nBas_Aux,nChV,iTOffs)
*
      Call mma_Allocate(iRv,nSkal2,Label='iRv')
      iRv(:)=0
*                                                                      *
************************************************************************
*                                                                      *
*     Let us now decide on the memory partitioning
*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     Preallocate some core for Seward!
*
      Call mma_maxDBLE(MemSew)
      MemLow=Min(MemSew/2,1024*128)
      MemSew=Max(MemSew/10,MemLow)
      Call xSetMem_Ints(MemSew)
*
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     During this phase we will have three memory sections
*
*     1) the three center intergrals for a fixed {kl}
*     2) a similar block for the R-vectors
*     3) a buffer to contain subsets of the Q-vectors
*
*     Compute the max size of 1 and 2
*
      n3CMax=0
      nRvMax=0
      Call IZero(iMax_R,2*nIrrep)
      Do klS_ = 1, nSkal2
         kS = iShij(1,klS_)
         lS = iShij(2,klS_)
         nRv = nSize_Rv(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_Rv,nChV)
         nRvMax = Max (nRvMax,nRv)
         n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,nBas_Aux)
         n3CMax = Max (n3CMax,n3C)
         Do iIrrep = 0, nIrrep-1
            iMax_R(1,iIrrep)=Max(iMax_R(1,iIrrep),iOff_3C(1,iIrrep))
            iMax_R(2,iIrrep)=iMax_R(2,iIrrep)+iOff_3C(1,iIrrep)
         End Do
      End Do
*
      Call mma_allocate(Arr_3C,n3CMax,Label='Arr_3C')
      Call mma_allocate(Rv,nRvMax,Label='Rv')
*
      Call mma_maxDBLE(MaxMem)
      nQv=0
      Do iIrrep = 0, nIrrep-1
         lJ=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) lJ=lJ-1 ! remove dummy basis function
         nQv = nQv + lJ*nChV(iIrrep)
      End Do
*
*     The Q-vectors can be read in a single whole block or in chunks.
*
      If (Force_Out_of_Core) MaxMem=(8*nQv)/10
      Out_of_Core = nQv.gt.MaxMem
      nQv = Min(nQv,MaxMem)  ! note that nQv is effectively reset here
      Call mma_allocate(Qv,nQv,Label='Qv')
*                                                                      *
************************************************************************
*                                                                      *
*     In case of in-core mode read Q-vectors only once!
*
      If (.Not.Out_of_Core) Then
         mQv=1
         Do iIrrep = 0, nIrrep-1
            lJ=nBas_Aux(iIrrep)
            If (iIrrep.eq.0) lJ=lJ-1 ! remove dummy basis function
*
            If (lJ.gt.0) Then
               iAddr=0
               kQv = lJ*nChV(iIrrep)
               Call dDaFile(Lu_Q(iIrrep),2,Qv(mQv),kQv,iAddr)
               mQv = mQv + kQv
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Open files for the R-vectors.
*
      Do iIrrep = 0, nIrrep-1
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.ne.0) Then
            iSeed=55+iIrrep
            Lu_R(iIrrep)=IsFreeUnit(iSeed)
            Write(Name_R,'(A4,I2.2)') 'RVec',iIrrep
            Call DaName_MF_WA(Lu_R(iIrrep),Name_R)
         End If
         iAddr_R(iIrrep)=0
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*
      kCenter=0  ! dummy initialize
      lCenter=0  ! dummy initialize
      iS = nSkal ! point to dummy shell
*     Save this field for the time being!
      Call ICopy(nIrrep,iTOffs(3),3,iTtmp,1)
*
      Call Init_Tsk(id,nSkal2)
*
*
      klS=0
       iTask=0
C      Do klS = 1, nSkal2
 100     Continue
         If (.Not.Rsv_Tsk(id,klS)) Go To 200
C        Write (*,*) 'Processing shell-pair:',klS
         iTask=iTask+1
*
         iRv(iTask) = klS
         kS = iShij(1,klS)
         lS = iShij(2,klS)
*
*        Logic to avoid integrals with mixed muonic and electronic
*        basis.
*
         kCnttp=iSD(13,kS)
         lCnttp=iSD(13,lS)
*
         If (LDF) Then
*
*           Pick up the corresponding (K|L)^{-1} block
*
            kCenter=iSD(10,kS)
            lCenter=iSD(10,lS)
C           Write (6,*) 'kCenter, lCenter=',kCenter, lCenter
            klCenter = kCenter*(kCenter-1)/2 + lCenter
            iAdr_AB=AB(1,klCenter)
            nAB    =AB(2,klCenter)
            Call dDaFile(Lu_AB,2,Local_A(:,2),nAB**2,iAdr_AB)
C           Call RecPrt('A^-1',' ',Local_A,nAB,nAB)
*
*           Now I need some lookup tables to be used below. I need to
*           go from SO index to lO index and from a given lO index
*           back to the SO index.
*
            Call IZero(ISO2LO,2*(MaxBfn+MaxBfn_Aux))
            iLO=0
            nCase=1
            If (kCenter.ne.lCenter) nCase=2
            Do iCase = 1, nCase
               If (iCase.eq.1) Then
                  jCenter=kCenter
               Else
                  jCenter=lCenter
               End If
               Do iSO_Aux = 1, nSO_Aux
                  iCenter=SO2C(iSO_Aux)
C                 Write (6,*) 'iCenter=',iCenter
                  If (iCenter.eq.jCenter) Then
                     iLO = iLO + 1
C                    Write (6,*) 'iLO,iSO_Aux=',iLO,iSO_Aux
                     iSO2LO(1,iSO_Aux)=iLO
                     iSO2LO(2,iLO)=iSO_Aux
                  End If
               End Do
            End Do
         End If
*
         Aint_kl = TMax_Valence(kS,lS)
         If (dbsc(kCnttp)%fMass.ne.dbsc(lCnttp)%fMass) Aint_kl=0.0D0
*
         nRv = nSize_Rv(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_Rv,
     &                  nChV)
         n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,iOff_3C,
     &                  nBas_Aux)
         Arr_3C(1:n3C)=Zero
         Rv(1:nRv)=Zero
*
         Call ICopy(nIrrep,iOff_3C,3,iTOffs(3),3)
*
*        Loop over the auxiliary basis set
*
         Do jS = nSkal_Valence+1, nSkal-1
C           Write (6,*) 'jS,kS,lS=',jS,kS,lS
            If (LDF) Then
               jCenter=iSD(10,jS)
               If (jCenter.ne.kCenter .and.
     &             jCenter.ne.lCenter ) Go To 14
C              Write (6,*) 'jCenter=',jCenter
            End If
*
            Aint=Aint_kl * TMax_Auxiliary(jS-nSkal_Valence)
*
#ifdef _DEBUGPRINT_
            Write (6,*)
            Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
            Write (6,*) 'AInt,CutInt=',AInt,CutInt
            Write (6,*)
#endif
            If (AInt.lt.CutInt) Go To 14
            Call Eval_IJKL(iS,jS,kS,lS,Arr_3C,n3C,Integral_WrOut)
 14         Continue
*
*           Use a time slot to save the number of tasks and shell
*           quadrupltes process by an individual node
            Call SavStat(1,One,'+')
            Call SavStat(2,One,'+')
*
         End Do    ! jS
*                                                                      *
************************************************************************
*                                                                      *
*        Multiply the 3-center integrals with the Q-vectors
*
*        Compute HQ
*
         Call Mult_3C_Qv_S(Arr_3C,n3C,Qv,nQv,Rv,nRv,nChV,
     &                     iOff_3C,nIrrep,Out_of_Core,Lu_Q,'N')
*                                                                      *
************************************************************************
*                                                                      *
*        Write the R-vectors to disk. These will be retrieved and sort
*        afterwards in step 3.
*
         Do iIrrep = 0, nIrrep-1
            ip_R = 1 + iOff_Rv(iIrrep)
            nRv=iOff_3C(1,iIrrep)*nChV(iIrrep)
C           Write (*,*) 'iAddr_R(iIrrep)=',iAddr_R(iIrrep)
            If (nRv.gt.0) Then
               Call dDaFile(Lu_R(iIrrep),1,Rv(ip_R),nRv,iAddr_R(iIrrep))
            End If
         End Do
*                                                                      *
************************************************************************
*                                                                      *
C      End Do    ! klS

      Go To 100
 200  Continue
      nTask=iTask
*
*     Restore iTOffs(3,*)
      Call ICopy(nIrrep,iTtmp,1,iTOffs(3),3)
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(1,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
      Call Free_Tsk(id)
*
*     Set up array to store the load balance if this might be needed in
*     a gradient calculation.
*
      Call mma_allocate(TmpList,nSkal2,Label='TmpList')
      TmpList(:)=0
      Call mma_allocate(LBList,nSkal2,Label='LBList')
      LBList(:)=-1
      Do iTask = 1, nTask
         klS_ = iRv(iTask)
         TmpList(klS_) = 1
      End Do
*
      iLB=1
      Do klS_ = 1, nSkal2
         If (TmpList(klS_).eq.1) Then
            LBList(iLB)=klS_
            iLB=iLB+1
         End If
      End Do
*
*
      Call Put_iArray('LBList',LBList,nSkal2)
*
      Call mma_deallocate(LBList)
      Call mma_deallocate(TmpList)
*
      Call mma_deallocate(Rv)
      Call mma_deallocate(Arr_3C)
      Call mma_deallocate(Qv)
      Call xRlsMem_Ints()
      Call mma_deallocate(TMax_Auxiliary)
      Call mma_deallocate(TMax_Valence)
      If (LDF) Then
         Call mma_deallocate(SO2C)
         Call mma_deallocate(AB)
         Call mma_deallocate(Local_A)
         Call DaClos(Lu_AB)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Each node does now have an incomplete set of R-vectors!
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*
      Call mma_deallocate(iSSOff)
      Call mma_deallocate(ShlSO)
      Call mma_deallocate(SOShl)
      Call Free_iSD()
*
*     Let go off the Q-vectors for now!
*
      Do iIrrep = 0, nIrrep-1
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.ne.0) Call DaClos(Lu_Q(iIrrep))
      End Do
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' 3-center integrals:'
         Call CWTime(TC0,TW0)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC0-TC1,' sec.',
     &                  '      Wall time:',TW0-TW1,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     P A R T I A L   T R A N S P O S E   S E C T I O N                *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     For the interface to work fix the tables of Seward
*
      Call Set_Basis_Mode('Valence')
      Call SetUp_iSD
*                                                                      *
************************************************************************
*                                                                      *
*     Initialize for 2-electron integral evaluation. Do generate
*     tables for indexation.
*
      Indexation = .True.
      Call Setup_Ints(nSkal_Valence,Indexation,ThrAO,DoFock,DoGrad)
*
*     Initiate stuff for Cholesky style storage.
*
      Call IniCho_RI(nSkal_Valence,nChV,nIrrep,iTOffs,iShij,nSkal2)

      Call mma_allocate(Addr,nSkal2,Label='Addr') ! addr for read
      Call mma_allocate(NuMu,2,nSkal2,Label='NuMu')
*                                                                      *
************************************************************************
*                                                                      *
*     Write out the RI vectors in Cholesky format
*
*     Here we will read one chuck from the R-vector file, while we will
*     store an as large part of the RI vectors in Cholesky format.
*
      LenVec=0
      Do iIrrep = 0, nIrrep-1
         iChoVec=0
*
         nB_Aux=nBas_Aux(iIrrep)
         If (iIrrep.eq.0) nB_Aux=nB_Aux-1
         If (nB_Aux.eq.0) Go To 998
*
         iSym = iIrrep+1
*
*        NumVec: is no longer equal to the # of auxiliary functions
*
         NumVec=iTOffs(3*iIrrep+1)
         If (NumVec.eq.0) Go To 999
*
         Addr(1) = 0
         Do i=2,nTask  ! init the addr for reading vectors
               klS_ = iRv(i-1)
               kS = iShij(1,klS_)
               lS = iShij(2,klS_)
               n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               Addr(i) = Addr(i-1) + nMuNu*NumVec
         End Do
*
         LenVec_Red = iMax_R(1,iIrrep)
         n_Rv = NumVec*LenVec_Red
         Call mma_allocate(Rv,n_Rv,Label='Rv')
*
*        LenVec: # of valence Gaussian products in this irrep
*
         LenVec = iMax_R(2,iIrrep)
         Call Create_Chunk(LenVec,NumVec,IncVec)
*
         Do iVec = 1, NumVec, IncVec
            NumVec_ = Min(NumVec-iVec+1,IncVec)
*                                                                      *
************************************************************************
*                                                                      *
*           Read now the R-vectors for a fixed shell-pair and
*           irrep, but for all auxiliary functions.
*
            mMuNu=0
            Do klS_ = 1, nSkal2
               kS = iShij(1,klS_)
               lS = iShij(2,klS_)
               n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               m3C = nMuNu * NumVec_
*
               If (m3C.le.0) Go To 555
*
               MuNu_s=mMuNu+1
               MuNu_e=mMuNu+nMuNu
*
               NuMu(1,klS_) = MuNu_s
               NuMu(2,klS_) = MuNu_e
*
 555           Continue
               mMuNu = mMuNu + nMuNu
            End Do
*
            Do i = 1, nTask
               klS_ = iRv(i)
               kS = iShij(1,klS_)
               lS = iShij(2,klS_)
*
               n3C = nSize_3C(kS,lS,nBasSh,nSkal-1,nIrrep,
     &                        iOff_3C,nBas_Aux)
               nMuNu = iOff_3C(1,iIrrep)
               m3C = nMuNu * NumVec_

               If (m3C.le.0) Go To 666
*
               Call dDaFile(Lu_R(iIrrep),2,Rv,m3C,Addr(i))

*              Copy the appropriate section into the RI vectors in
*              Cholesky format.
*
               MuNu_s = NuMu(1,klS_)
               MuNu_e = NuMu(2,klS_)
               j_s=1
               j_e=NumVec_
               Call Put_Chunk(MuNu_s,MuNu_e,j_s,j_e,Rv,nMuNu,LenVec)
*
 666           Continue
            End Do
*                                                                      *
************************************************************************
*                                                                      *
*
*           Now transfer the RI vectors to disk
*
            Call Get_Chunk(LenVec,NumVec_,iChoVec,iSym,iVec)
*
         End Do   ! iVec = 1, NumVec, IncVec
*

         Call Destroy_Chunk()
         Call mma_deallocate(Rv)
*
 999     Continue
*
*        Let go off the R-vectors for good!
*
         Call DaClos(Lu_R(iIrrep))
 998     Continue
         NoChoVec(iIrrep)=iChoVec
*
      End Do    ! iIrrep
      Call mma_deallocate(NuMu)
      Call mma_deallocate(Addr)
      Call mma_deallocate(iRv)
      Call mma_deallocate(nBasSh)
      Call mma_deallocate(iShij)
*                                                                      *
************************************************************************
*                                                                      *
      iPass = 1
      iRed = 1
      Call Cho_RI_PutInfo(iPass,iRed)
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate integral environment.
*
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*
      If (iPrint.ge.6) Then
         Write (6,'(A)') ' Block-transpose:'
         Call CWTime(TC1,TW1)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC1-TC0,' sec.',
     &                  '      Wall time:',TW1-TW0,' sec.'
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     D I A G O N A L   S E C T I O N                                  *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      nDiag = 0
      Do iIrrep = 0, nIrrep-1
         nDiag = nDiag+nBas(iIrrep)
      End Do
      nDiag = nDiag*(nDiag+1)/2
      Call mma_allocate(Diag,nDiag,Label='Diag')
      Diag(:)=Zero
*
      Call Drv2El_RI_Diag(ThrAO,Diag,nDiag)
*
*     Write the diagonal to disk
*
      Call Cho_IODiag(Diag,1)
*
      Call mma_deallocate(Diag)
*
      If (iPrint.ge.6) Then
         Write (6,*) 'Diagonal vector:'
         Call CWTime(TC0,TW0)
         Write (6,'(A,F8.2,A,/,A,F8.2,A)')
     &                  '      CPU time :',TC0-TC1,' sec.',
     &                  '      Wall time:',TW0-TW1,' sec.'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Terminate Cholesky stuff here.
*
      irc = 0
      Call TermCho_RI(irc,NoChoVec,8)
      If (irc .ne. 0) Then
         Write(6,*) 'TermCho_RI returned ',irc
         Call SysAbendMsg('Drv2El_3Center_RI',
     &                    'Cholesky termination failed!',' ')
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

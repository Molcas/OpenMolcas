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
* Copyright (C) 1990,1991,1992,2000,2007, Roland Lindh                 *
*               1990, IBM                                              *
************************************************************************
      SubRoutine Drvg1_3Center_RI(Grad,Temp,nGrad,ij3,nij_Eff)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will controll the type of the two-electron integral, eg.    *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis func-  *
*          tions of the requested type.                                *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
* For RI-HF gradients read:                                            *
* "Analytical Gradients of Hartee-Fock Exchange with Density Fitting   *
* Approximations", J. Bostrom, F. Aquilante, T. B. Pedersen and R.     *
* Lindh, JCTC,  9:204-212 (2013).                                      *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified for k2 loop. August '91                         *
*             Modified for gradient calculation. January '92           *
*             Modified for SetUp_Ints. January '00                     *
*             Modified for 3-center RI gradients, March '07            *
************************************************************************
      use k2_setup
      use iSD_data
      use pso_stuff
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
      use Basis_Info
      use Sizes_of_Seward, only:S
      use Gateway_Info, only: CutInt
      use RICD_Info, only: Do_RI
      use Symmetry_Info, only: nIrrep
      use ExTerm, only: CijK, CilK, BklK, VJ
      use ExTerm, only: Ymnij, ipYmnij, nYmnij, iOff_Ymnij
      use ExTerm, only: Yij, BMP2, iMP2prpt, CMOi, DMLT
      use Data_Structures, only: Deallocate_DSBA
      Implicit Real*8 (A-H,O-Z)
      Logical, External :: Rsv_Tsk2
#include "Molcas.fh"
#include "itmax.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "cholesky.fh"
#include "setup.fh"
#include "exterm.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "bdshell.fh"
      Integer nGrad, nij_Eff
      Real*8  Grad(nGrad), Temp(nGrad)
      Integer, Allocatable:: ij3(:,:)
*     Local arrays
      Real*8  Coor(3,4)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4),
     &        nAct(0:7)
      Logical EQ, Shijij, AeqB, CeqD, DoGrad, DoFock, Indexation,
     &        JfGrad(3,4), ABCDeq, No_Batch, Found, FreeK2, Verbose
      Character Format*72, Method*8, KSDFT*16
      Character*50 CFmt
      Character(LEN=16), Parameter :: SECNAM = 'drvg1_3center_ri'
      Integer, External:: Cho_irange
*
      Integer iSD4(0:nSD,4)
      save MemPrm
      Logical FlipFlop
#include "chotime.fh"

      Real*8, Allocatable:: MaxDens(:), SDG(:), Thhalf(:)
      Integer, Allocatable:: Shij(:,:), Shij2(:,:), LBList(:)
      Real*8, Allocatable:: CVec2(:,:,:), CVec(:,:)
      Real*8, Allocatable:: Xmi(:,:,:,:)
      Real*8, Allocatable:: Tmp(:,:), TMax_Valence(:,:),
     &                      TMax_Auxiliary(:)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
#ifdef _CD_TIMING_
      Twoel3_CPU = 0.0d0
      Twoel3_Wall = 0.0d0
      Pget3_CPU = 0.0d0
      Pget3_Wall = 0.0d0
#endif
      iFnc(1)=0
      iFnc(2)=0
      iFnc(3)=0
      iFnc(4)=0
      PMax=Zero
      call dcopy_(nGrad,[Zero],0,Temp,1)
*                                                                      *
************************************************************************
*                                                                      *
      xfk=1.0D-3 ! changing this parameter tunes LK-type screening thr
*     xfk=1.0D-12! Debugging
*     xfk=0.0D-12! Debugging
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_dScalar('Cholesky Threshold',ThrCom)
      ThrCom=Max(ThrCom,1.0d-6) ! not to sacrify efficiency too much
*                                                                      *
************************************************************************
*                                                                      *
*     Handle mixed basis set
*
      If (Do_RI) Then
         Call Set_Basis_Mode('Auxiliary')
         Call Nr_Shells(nSkal_Auxiliary)
         Call Set_Basis_Mode('WithAuxiliary')
      Else
         Call Set_Basis_Mode('Valence')
         nSkal_Auxiliary=0
      End If
      Call SetUp_iSD()
*                                                                      *
************************************************************************
*                                                                      *
*-----Precompute k2 entities.
*
      Indexation=.False.
      DoFock=.False.
      DoGrad=.True.
      ThrAO=Zero
      Call SetUp_Ints(nSkal,Indexation,ThrAO,DoFock,DoGrad)
      nSkal_Valence=nSkal-nSkal_Auxiliary
      mSkal=nSkal
      nPairs=nSkal*(nSkal+1)/2
      nQuad =nPairs*(nPairs+1)/2
      Pren = Zero
      Prem = Zero
*                                                                      *
************************************************************************
*                                                                      *
      MxPrm = 0
      Do iAng = 0, S%iAngMx
         MxPrm = Max(MxPrm,S%MaxPrm(iAng))
      End Do
      nZeta = MxPrm * MxPrm
      nEta  = MxPrm * MxPrm
*                                                                      *
************************************************************************
*                                                                      *
      maxnAct=0
      If (lPSO) Then
         Call Get_iArray('nAsh',nAct,nIrrep)
         maxnnP=nnP(0)
         maxnAct=nAct(0)
         Do i=1,nIrrep-1
           maxnnP=max(maxnnP,nnP(i))
           maxnAct=max(maxnAct,nAct(i))
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---  Compute entities for prescreening at shell level
*
      Call mma_allocate(TMax_Valence,nSkal_Valence,nSkal_Valence,
     &                  Label='TMax_Valence')
      nTMax=nSkal_Auxiliary
      If (Do_RI) nTMax = Max(1,nTMax-1)
      Call mma_allocate(TMax_Auxiliary,nTMax,Label='TMax_Auxiliary')
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
      If (Do_RI) Then
         Do iS = 1, nSkal_Auxiliary-1
            iS_ = iS + nSkal_Valence
            jS_ = nSkal_Valence + nSkal_Auxiliary
            TMax_Auxiliary(iS)=Tmp(iS_,jS_)
         End Do
      End If
*
      Call mma_deallocate(Tmp)

*                                                                      *
************************************************************************
*                                                                      *
      MxInShl = 1
      Do i = 1, nSkal_Valence
         MxInShl = max(MxInShl,iSD(3,i)*iSD(2,i))
      End Do
*
*     Calculate maximum density value for each shellpair
*
      lMaxDens = nSkal_Valence*(nSkal_Valence+1)/2
      Call mma_allocate(MaxDens,lMaxDens,Label='MaxDens')
      MaxDens(:)=Zero
*
      Do iSym = 0, nSym-1
         kS=1+nSkal_Valence*iSym ! note diff wrt declaration of iBDsh
         Do j=1,nBas(iSym)
            jsh=Cho_Irange(j,iBDsh(kS),nSkal_Valence,.true.)
            Do i=1,j
               ish=Cho_Irange(i,iBDsh(kS),nSkal_Valence,.true.)
               ijS=jsh*(jsh-1)/2+ish
               Do iSO=1,nJDens
                 If (.NOT.DMLT(iSO)%Active) Cycle
                 ij=j*(j-1)/2+i
                 Dm_ij=abs(DMLT(iSO)%SB(iSym+1)%A1(ij))
                 MaxDens(ijS)=Max(MaxDens(ijS),Dm_ij)
               End Do
            End Do
         End Do
      End Do
*
      Do i = 1, 5
         If (DMLT(i)%Active) Call deallocate_DSBA(DMLT(i))
      End Do
*
*     Create list of non-vanishing pairs
*
*1)   For the valence basis set
*
*
      Call mma_allocate(Shij,2,nSkal_Valence*(nSkal_Valence+1)/2,
     &                  Label='Shij')
      nSkal2=0
      Do iS = 1, nSkal_Valence
         iiQ = iS*(iS+1)/2
         XDm_ii = MaxDens(iiQ)
         Do jS = 1, iS
            jjQ = jS*(jS+1)/2
            XDm_jj = MaxDens(jjQ)
            ijQ=iS*(iS-1)/2+jS
            XDm_ij = MaxDens(ijQ)
            XDm_max = Max(XDm_ij,XDm_ii,XDm_jj)
            Aint_ij=TMax_Valence(iS,jS)
            If (TMax_All*Aint_ij .ge. CutInt) Then
*
* --- FAQ: more aggressive screening to discard shprs formed
*          by AOs contributing mainly to the virtual MO space.
*
               If (Aint_ij*XDm_max .ge. CutInt) Then
                  nSkal2 = nSkal2 + 1
                  Shij(1,nSkal2)=iS
                  Shij(2,nSkal2)=jS
               End If
*
            End If
         End Do
      End Do
*
*2)   For the auxiliary basis set
*
      If (Do_RI) Then
         mij = nSkal_Auxiliary
         Call mma_allocate(Shij2,2,mij,Label='Shij2')
         nij=0
         Do jS = nSkal_Valence+1, nSkal-1
            If (TMax_All*TMax_Auxiliary(jS-nSkal_Valence).ge.CutInt)
     &         Then
               nij = nij + 1
               Shij2(1,nij)=nSkal
               Shij2(2,nij)=jS
            End If
         End Do
      Else
         mij = nij_Eff
         Call mma_allocate(Shij2,2,mij,Label='Shij2')
         nij=0
         Do ij = 1, nij_Eff
            iS = ij3(1,ij)
            jS = ij3(2,ij)
            If (TMax_All*TMax_Valence(iS,jS).ge.CutInt) Then
               nij = nij + 1
               Shij2(1,nij)=iS
               Shij2(2,nij)=jS
            End If
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (DoCholExch) Then
*
*        Find the largest number of contractions in any given shell of
*        auxiliary functions.
*
         MxChVInShl = 1
         If(Do_RI) Then
            Do i = nSkal_Valence+1, nSkal_Valence+nSkal_Auxiliary
               MxChVInShl = max(MxChVInShl,iSD(3,i))
            End Do
         Else
            Write (6,*) 'Not implemented for Cholesky yet!'
            Call Abend()
         End If
*
*        Find the largest set of ij. The i and j basis are due to the CD
*        of the one-particle density matrix.
*
*        nIJ1: diagonal blocks are triangularized.
*        nIJR: diagonal blocks are square.
*        nIMax: largest number of i basis in any irrep.
*
         nIJ1Max = 0
         nIJRMax = 0
         nIMax = 0
         Do iSym = 1, nIrrep
            Do iSO=1,nKDens
               nIMax = max(nIMax,nChOrb(iSym-1,iSO))
            End Do
            Do jSym = 1, nIrrep
               Do iSO=1,nKVec
                 nIJ1Max = max(nIJ1Max,nIJ1(iSym,jSym,iSO))
                 nIJRMax = max(nIJRMax,nIJR(iSym,jSym,iSO))
               End Do
            End Do
         End Do
*
*        Allocate scratch memory for step 4 (Eq. 16). This is done
*        in two steps.
*
*        First step: Sum(j) X_lj C_ij^K = C_il^K; (l=valence basis)
*        Second step: Sum(i) C_il^K X_ki = B_kl^K
*
         lCijK = nIJRMax*MxChVInShl
         lCilK = MxInShl*nIMax*MxChVInShl
         lCilK = Max(lCilK,lCijK) ! it is used as scratch in pget
         Call mma_allocate(CijK,lCilK,Label='CijK')
         If (lPSO) lCilK=Max(lCilK,maxnAct) ! used as scratch
         Call mma_allocate(CilK,lCilK,Label='CilK')
         lBklK = MxInShl**2*MxChVInShl
         Call mma_allocate(BklK,lBklK,Label='BklK')
*
         If(iMp2prpt .eq. 2) Then
            lB_mp2 = mxChVInShl*nBas(0)*nBas(0)
            Call mma_allocate(BMP2,lB_mp2,2,Label='BMP2')
         End If
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        The C_ij^K vectors are stored in triangular form. We now
*        change this to stricked rectangular/square form. Diagonal
*        elements are rescaled. In case of symmetry this is only
*        done for the blocks with isym=jSym, i.e. kSym=1
*
         kSym = 1
         nK = NumAuxVec(kSym)
*
         Do iSO=1,nKVec
           If (lSA) Go to 15
           Do iSym = 1, nSym
              jSym = iSym
*             jSym = iEor(kSym-1,jSym-1)+1
*
*             Read a whole block of C_ij^K
*
              iAdrC = iAdrCVec(kSym,iSym,iSO)
              Call mma_allocate(CVec,nIJ1(iSym,jSym,iSO),nK,
     &                          Label='CVec')
              Call dDaFile(LuCVector(kSym,iSO),2,CVec,
     &                     nIJ1(iSym,jSym,iSO)*nK,iAdrC)
*
              ni = nChOrb(iSym-1,iSO)
              Call mma_allocate(CVec2,ni,ni,nK,Label='CVec2')
*             nj=ni

              Do KAux = 1, nK
*
                 Do i = 1, ni
                    Do j = 1, i-1
                       ij = j + i*(i-1)/2
                       CVec2(i,j,KAux) = CVec(ij,KAux)
                       CVec2(j,i,KAux) = CVec(ij,KAux)
                    End Do
                    ii = i*(i+1)/2
                    CVec2(i,i,KAux) = CVec(ii,KAux)*Sqrt(Two)
                 End Do
*
              End Do
*
*             Write back to disk. Note that the file is prepared for
*             rectangular/square storage so that one can safely write
*             back the expanded set to disk without any overwrite
*             problems.
*
              iAdrC = iAdrCVec(kSym,iSym,iSO)
              Call dDaFile(LuCVector(kSym,iSO),1,CVec2,
     &                     nIJR(iSym,jSym,iSO)*nK,iAdrC)
*
              Call mma_deallocate(CVec2)
              Call mma_deallocate(CVec)
           End Do
 15        Continue
         End Do
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        Stuff used in the prescreening!
*
         MumOrb=0
         NumOrb=0
         Do iSO=1,nKDens
            Do jSym = 1, nSym
               NumOrb = NumOrb + nChOrb(jSym-1,iSO)
               MumOrb = Max( MumOrb, nChOrb(jSym-1,iSO) )
            End Do
            If (iSO.lt.nKDens) ipYmnij(iSO+1)=NumOrb
         End Do
*
*        Scratch store the index of the MOs which finds the estimate
*        according to Eq. 18 to be larger than the threshold.
*
         Call mma_allocate(Ymnij,NumOrb,Label='Ymnij')
         Ymnij(:)=0
         ipYmnij(1)=1
         Do i=2,nKDens
           ipYmnij(i)=ipYmnij(1)+ipYmnij(i)
         End Do
*
*        Make a list for each shell-pair over the largest element
*        SQRT(ABS( (mu,nu|mu,nu) ))
*
         nnSkal = nSkal_valence*(nSkal_valence+1)/2
         Call mma_allocate(SDG,nnSkal,Label='SDG')
         Call get_maxDG(SDG,nnSkal,MxBasSh)
*
*        Scratch for reduced lists of X_mi. Used in pget.
*
         nXki=MumOrb*MxBasSh*nSym
         Call mma_allocate(Yij,nXki,2,nKDens,Label='Yij')
*
*        Make a list the largest element X_mu,i for each valence shell
*        and a fixed i. X_mu,i defined in Eq. 13.
*
         Call mma_allocate(Xmi,MumOrb,nSkal_Valence,nIrrep,nKDens,
     &                     Label='Xmi')
*
         Do iSO=1,nKDens
            Call get_mXOs(iSO,Xmi(:,:,:,iSO),MumOrb,nSkal_valence,
     &                 nIrrep,nChOrb(0,iSO))
         End Do
*
      Else
*
         nXki=0
         NumOrb=0
         MumOrb=0
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     For CASSCF process the active space contribution.
*
      If (lPSO) Then
         nBas_Aux(0)=nBas_Aux(0)-1
         Call mma_allocate(Thhalf,maxnnP,Label='Thhalf')
         nThpkl=MxChVInShl*MxInShl**2
         Call mma_allocate(Thpkl,nThpkl,Label='Thpkl')
*
         Call contract_Zpk_Tpxy(Z_p_k ,nZ_p_k,
     &                          Txy   ,n_Txy,
     &                          Thhalf,maxnnP,
     &                          DMdiag ,nG1,
     &                          nnP,nBas_Aux,
     &                          nADens,nAvec,nAct,nIrrep)
*
         Call mma_deallocate(Thhalf)
         nBas_Aux(0)=nBas_Aux(0)+1
      Else
         nThpkl=1
         Call mma_allocate(Thpkl,nThpkl,Label='Thpkl')
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-------Compute FLOP's for the transfer equation.
*
      Do iAng = 0, S%iAngMx
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
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     For a parallel implementation the iterations over shell-pairs
*     are parallelized.
*
*     If only Coulombic terms are to be processed use dynamic setup.
*     Otherwise do batches exactly in the same order as Seward did the
*     2-center terms.
*
      Call Get_cArray('Relax Method',Method,8)
      If (Method.ne.'KS-DFT  ') Then
         iOpt=1
      Else
         Call Get_cArray('DFT functional',KSDFT,16)
         ExFac=Get_ExFac(KSDFT)
         iOpt=0
         If (ExFac.ne.Zero) iOpt=1
      End If
      If(.not. Do_RI) iOpt=0
*
      If (iOpt.eq.1) Then
         Call qpg_iArray('LBList',Found,nSkal2_)
         If (Found) Then
            Call mma_allocate(LBList,nSkal2_,Label='LBList')
            Call Get_iArray('LBList',LBList,nSkal2_)
         Else
            Call WarningMessage(2,'LBList not found!')
            Call Abend()
         End If
      Else
         Call mma_allocate(LBList,1,Label='LBList')
      End If
      Call Init_Tsk2(id,nSkal2,iOpt,LBList)
      Call mma_deallocate(LBList)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_MaxDBLE(MemMax)
      Call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
      ipMem1=1
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Do klS = 1, nSkal2
   10    Continue
         If (.Not.Rsv_Tsk2(id,klS)) Go To 11
*
         kS = Shij(1,klS)
         lS = Shij(2,klS)
*
         AInt_kl = TMax_Valence(kS,lS)
*
         klS_=iTri(kS,lS)
*                                                                      *
************************************************************************
*                                                                      *
*        Prescreening stuff for exchange
*
         If (DoCholExch) Then
*
*
*           For the shell-pair, (kS,lS), pick up the largest element
*           Sqrt(Abs(  (kappa,lambda|kappa,lambda) ))
*
            SDGmn=SDG(klS_)
*
*
*           Loop over the MO basis, jb and ib and approximate Y_ij
*           (Eq. 18)
*
            Do iSO=1,nKVec
              FlipFlop=.True.
              iMOleft=iSO
              iMOright=iSO
              If (lSA) iMOright=iSO+2
 20           Continue
              nj=0
              Do jSym = 1, nSym
*
                mj=0
                Do jb= 1, nChOrb(jSym-1,iMOleft)
                  Xjk=Xmi(jb,kS,jSym,iMOleft)
                  Xjl=Xmi(jb,lS,jSym,iMOleft)
*
                  jSym_s = jSym
                  if (ks.ne.ls.or.iMOright.ne.iMOleft) jSym_s=1
                  Do iSym = jsym_s, nSym
                     NumOrb_i  = nChOrb(iSym-1,iMOright)
                     If (iSym.eq.jSym.and.ks.eq.ls
     &                    .and.iMOright.eq.iMOleft) NumOrb_i = jb
*
                     Do ib = NumOrb_i, 1, -1
                       Xik=Xmi(ib,kS,iSym,iMOright)
                       Xil=Xmi(ib,lS,iSym,iMOright)
*
* ---                  Yij[mn] = (1+Pij) Xim * (mn|mn)^1/2 * Xjn
*
                       PZmnij=(Xik*Xjl+Xil*Xjk)*SDGmn
*
*                      If larger than the threshold put j in the
*                      list and exit the loop.
*
                       If ( PZmnij.ge.xfk*ThrCom ) Then
!                        orbital in the list
                         Ymnij(ipYmnij(iMOleft)+mj+nj)=jb
                         mj=mj+1
                         Go To 666
                       End If
                     End Do ! ib
                  End Do    ! iSym
 666              Continue
                End Do
*
*               The first element is to keep track on how many elements
*               that were saved.
*
!               nOrbs in the list ==> dim(ij)=nOrbs**2
                nYmnij(jSym,iMOleft)=mj
                iOff_Ymnij(jSym,iMOleft) = nj
                nj = nj + mj
*
              End Do ! jSym
              If (lSA.and.FlipFlop) Then
                FlipFlop=.False.
                itmp=iMOleft
                iMOleft=iMOright
                iMOright=itmp
                Go To 20
              EndIf
            End Do
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Do ijS = 1, nij
         iS = Shij2(1,ijS)
         jS = Shij2(2,ijS)
*
         If (Do_RI) Then
            Aint=AInt_kl*TMax_Auxiliary(jS-nSkal_Valence)
         Else
            Aint=AInt_kl*TMax_Valence(iS,jS)
         End If
         If (AInt.lt.CutInt) Go To 14
         If (iPrint.ge.15) Write (6,*) 'iS,jS,kS,lS=',iS,jS,kS,lS
*                                                                      *
************************************************************************
*                                                                      *
         Call Gen_iSD4(iS, jS, kS, lS,iSD,nSD,iSD4)
         Call Size_SO_block_g(iSD4,nSD,nSO,No_batch)
         If (No_batch) Go To 140
*
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
         ABCDeq=EQ(Coor(1,1),Coor(1,2)) .and.
     &          EQ(Coor(1,1),Coor(1,3)) .and.
     &          EQ(Coor(1,1),Coor(1,4))
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
         Call SOAO_g(iSD4,nSD,nSO,
     &               MemPrm, MemMax,
     &               iBsInc,jBsInc,kBsInc,lBsInc,
     &               iPrInc,jPrInc,kPrInc,lPrInc,
     &               ipMem1,ipMem2, Mem1,  Mem2,
     &               iFnc, MemPSO)
         iBasi    = iSD4(3,1)
         jBasj    = iSD4(3,2)
         kBask    = iSD4(3,3)
         lBasl    = iSD4(3,4)
*                                                                      *
************************************************************************
*                                                                      *
         Call Int_Parm_g(iSD4,nSD,iAnga,
     &                 iCmpa,iShlla,iShela,
     &                 iPrimi,jPrimj,kPrimk,lPriml,
     &                 k2ij,nDCRR,k2kl,nDCRS,
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
               If (iSh.eq.1.and.Do_RI) Then
                  JfGrad(iCar,iSh) = .False.
                  JndGrd(iCar,iSh) = 0
               Else If (iAnd(iSD4(15,iSh),2**(iCar-1)) .eq.
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
*
         Do 420 kBasAO = 1, kBask, kBsInc
           kBasn=Min(kBsInc,kBask-kBasAO+1)
           iAOst(3) = kBasAO-1
         Do 430 lBasAO = 1, lBasl, lBsInc
           lBasn=Min(lBsInc,lBasl-lBasAO+1)
           iAOst(4) = lBasAO-1
*
*----------Get the 2nd order density matrix in SO basis.
*
           nijkl = iBasn*jBasn*kBasn*lBasn
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU1,Pget0WALL1)
#endif
           Call PGet0(iCmpa,
     &                iBasn,jBasn,kBasn,lBasn,Shijij,
     &                iAOV,iAOst,nijkl,Sew_Scr(ipMem1),nSO,
     &                iFnc(1)*iBasn,iFnc(2)*jBasn,
     &                iFnc(3)*kBasn,iFnc(4)*lBasn,MemPSO,
     &                Sew_Scr(ipMem2),Mem2,iS,jS,kS,lS,nQuad,PMax)
#ifdef _CD_TIMING_
           CALL CWTIME(Pget0CPU2,Pget0WALL2)
           Pget3_CPU = Pget3_CPU + Pget0CPU2-Pget0CPU1
           Pget3_Wall = Pget3_Wall + Pget0WALL2-Pget0WALL1
#endif
            If (AInt*PMax.lt.CutInt) Then
               Go To 430
            End If
*
*----------Compute gradients of shell quadruplet
*
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU1,TwoelWall1)
#endif
           Call TwoEl_g(Coor,
     &          iAnga,iCmpa,iShela,iShlla,iAOV,
     &          mdci,mdcj,mdck,mdcl,nRys,
     &          Data_k2(k2ij),nab,nHmab,nDCRR,
     &          Data_k2(k2kl),ncd,nHmcd,nDCRS,Pren,Prem,
     &          iPrimi,iPrInc,jPrimj,jPrInc,
     &          kPrimk,kPrInc,lPriml,lPrInc,
     &          Shells(iSD4(0,1))%pCff(1,iBasAO),iBasn,
     &          Shells(iSD4(0,2))%pCff(1,jBasAO),jBasn,
     &          Shells(iSD4(0,3))%pCff(1,kBasAO),kBasn,
     &          Shells(iSD4(0,4))%pCff(1,lBasAO),lBasn,
     &          Mem_DBLE(ipZeta),Mem_DBLE(ipZI),Mem_DBLE(ipP),nZeta,
     &          Mem_DBLE(ipEta), Mem_DBLE(ipEI),Mem_DBLE(ipQ),nEta,
     &          Mem_DBLE(ipxA),Mem_DBLE(ipxB),
     &          Mem_DBLE(ipxG),Mem_DBLE(ipxD),Temp,nGrad,
     &          JfGrad,JndGrd,Sew_Scr(ipMem1), nSO,Sew_Scr(ipMem2),Mem2,
     &          Aux,nAux,Shijij)
#ifdef _CD_TIMING_
           Call CWTIME(TwoelCPU2,TwoelWall2)
           Twoel3_CPU = Twoel3_CPU + TwoelCPU2-TwoelCPU1
           Twoel3_Wall = Twoel3_Wall + TwoelWall2-TwoelWall1
#endif
*
            If (iPrint.ge.15)
     &         Call PrGrad(' In Drvg1_3Center_RI: Grad',
     &                  Temp,nGrad,ChDisp)
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
         End Do
*
         Go To 10
 11   Continue
*     End of big task loop
*                                                                      *
************************************************************************
*                                                                      *
*     Write Timings:
*
      If(Timings) Then
       TotCPU = tbvec(1) + tavec(1)
       TotWall = tbvec(2) + tavec(2)
       CFmt='(2x,A)'
       Write(6,*)
       Write(6,CFmt)'Cholesky Gradients timing from A and B vectors:'
       Write(6,CFmt)'-----------------------------------------------'
       Write(6,*)
       Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
       Write(6,CFmt)'                                CPU       WALL   '
       Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'

       Write(6,'(2x,A26,2f10.2)')'Density (2-center):               '
     &                           //'         ',tavec(1),tavec(2)
       Write(6,'(2x,A26,2f10.2)')'Density (3-center):               '
     &                           //'         ',tbvec(1),tbvec(2)
       Write(6,*)
       Write(6,'(2x,A26,2f10.2)')'TOTAL                             '
     &                           //'         ',TotCPU,TotWall
       Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
       Write(6,*)
*
      End If
      timings = timings_default
*                                                                      *
************************************************************************
*                                                                      *
*                         E P I L O G U E                              *
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate scratch for exchange term
*
      If (DoCholExch) Then
         Call mma_deallocate(Xmi)
         Call mma_deallocate(Yij)
         Call mma_deallocate(SDG)
         Call mma_deallocate(Ymnij)
      End If
      If (Allocated(CijK)) Call mma_deallocate(CijK)
      If (Allocated(CilK)) Call mma_deallocate(CilK)
      If (Allocated(BklK)) Call mma_deallocate(BklK)
      If (Allocated(VJ)) Call mma_deallocate(VJ)
      Do i=1,nKDens
         Call Deallocate_DSBA(CMOi(i))
      End Do
      Call mma_deallocate(MaxDens)
*
      If (Allocated(BMP2)) Call mma_deallocate(BMP2)
      If (Allocated(Thpkl)) Call mma_deallocate(Thpkl)
*
      Call mma_deallocate(Sew_Scr)
      Call Free_Tsk2(id)
      Call mma_deallocate(Shij2)
      Call mma_deallocate(Shij)
      Call mma_deallocate(TMax_Auxiliary)
      Call mma_deallocate(TMax_Valence)
*                                                                      *
************************************************************************
*                                                                      *
      Verbose = .False.
      FreeK2=.True.
      Call Term_Ints(Verbose,FreeK2)
*                                                                      *
************************************************************************
*                                                                      *
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
      Call Free_iSD()
*
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Grad)
      End

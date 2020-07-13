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
      SubRoutine Drvg1_3Center_RI(Grad,Temp,nGrad,ip_ij3,nij_Eff)
************************************************************************
*                                                                      *
*  Object: driver for two-electron integrals. The four outermost loops *
*          will controll the type of the two-electron integral, eg.    *
*          (ss|ss), (sd|pp), etc. The next four loops will generate    *
*          list of symmetry distinct centers that do have basis func-  *
*          tions of the requested type.                                *
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
*                                                                      *
************************************************************************
      use k2_setup
      use iSD_data
      use pso_stuff
      use k2_arrays, only: ipZeta, ipiZet, Mem_DBLE, Aux, Sew_Scr
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External Rsv_Tsk2
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "cholesky.fh"
#include "setup.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
      Integer iBDsh(MxShll*8)
      Common /BDshell/ iBDsh
      Integer  Cho_irange
      External Cho_irange
*     Local arrays
      Real*8  Coor(3,4), Grad(nGrad), Temp(nGrad)
      Integer iAnga(4), iCmpa(4), iShela(4),iShlla(4),
     &        iAOV(4), istabs(4), iAOst(4), JndGrd(3,4), iFnc(4),
     &        nAct(0:7)
      Integer ipXmi(5)
      Integer nHrrTb(0:iTabMx,0:iTabMx,2)
      Logical EQ, Shijij, AeqB, CeqD, DoGrad, DoFock, Indexation,
     &        JfGrad(3,4), ABCDeq, No_Batch, Rsv_Tsk2, Found,
     &        FreeK2, Verbose
      Character Format*72, Method*8, KSDFT*16
      Character*50 CFmt
      Character*16 SECNAM
      Parameter (SECNAM = 'drvg1_3center_ri')
*
      Integer iSD4(0:nSD,4)
      save MemPrm
      Logical Timings,FlipFlop
      COMMON /CHOTIME / timings
#include "ymnij.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      TMax_Valence(i,j)=Work(ipTMax-1+(j-1)*nSkal_Valence+i)
      TMax_Auxiliary(i)=Work(ipTMax-1+nSkal_Valence**2+i)
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 9
      iPrint = nPrint(iRout)
      Call QEnter('Drvg1_3Center_RI')
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
      idum=0
      idum1=0
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
      Do iAng = 0, iAngMx
         MxPrm = Max(MxPrm,MaxPrm(iAng))
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
      nTMax=nSkal_Valence**2+nSkal_Auxiliary
      If (Do_RI) nTMax = nTMax-1
      Call GetMem('TMax','Allo','Real',ipTMax,nTMax)
*
      Call Allocate_Work(ip_Tmp,nSkal**2)
      Call Shell_MxSchwz(nSkal,Work(ip_Tmp))
      TMax_all=Zero
      Do iS = 1, nSkal_Valence
         Do jS = 1, iS
            ip_Out=ip_Tmp + (jS-1)*nSkal + iS -1
            ip_In =ipTMax + (jS-1)*nSkal_Valence + iS -1
            Work(ip_In)=Work(ip_Out)
cVV: ifort 11 can't handle the code without this dummy print.
            if(iPrint.gt.100) write(6,*) ip_In, ip_Out
            ip_In =ipTMax + (iS-1)*nSkal_Valence + jS -1
            Work(ip_In)=Work(ip_Out)
            TMax_all=Max(TMax_all,Work(ip_Out))
         End Do
      End Do
      If (Do_RI) Then
         Do iS = 1, nSkal_Auxiliary-1
            iS_ = iS + nSkal_Valence
            jS_ = nSkal_Valence + nSkal_Auxiliary
            ip_Out = ip_Tmp + (iS_-1)*nSkal + jS_ -1
            ip_In  = ipTMax + nSkal_Valence**2 + iS -1
            Work(ip_In)=Work(ip_Out)
         End Do
      End If
*
      Call Free_Work(ip_Tmp)

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
      Call GetMem('MaxDensity','Allo','Real',ip_MaxDens,lMaxDens)
      Call FZero(Work(ip_MaxDens),lMaxDens)
*
      iOff=0
      Do iSym = 0, nSym-1
         iDlt=ipDMLT(1)+iOff-1
         kS=1+nSkal_Valence*iSym ! note diff wrt declaration of iBDsh
         Do j=1,nBas(iSym)
            jsh=Cho_Irange(j,iBDsh(kS),nSkal_Valence,.true.)
            Do i=1,j
               ish=Cho_Irange(i,iBDsh(kS),nSkal_Valence,.true.)
               ijS=ip_MaxDens-1+jsh*(jsh-1)/2+ish
               Do iSO=1,nJDens
                 ij=ipDMLT(iSO)+iOff-1+j*(j-1)/2+i
                 Dm_ij=abs(Work(ij))
                 Work(ijS)=Max(Work(ijS),Dm_ij)
               End Do
            End Do
         End Do
         iOff=iOff+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
*
      Call Free_Work(ipDMLT(1))
      If (nKdens.eq.2) Call Free_Work(ipDMLT(2))
*
*     Create list of non-vanishing pairs
*
*1)   For the valence basis set
*
*
      Call GetMem('ip_ij','Allo','Inte',ip_ij,
     &            nSkal_Valence*(nSkal_Valence+1))
      nSkal2=0
      Do iS = 1, nSkal_Valence
         iiQ = iS*(iS+1)/2
         XDm_ii = Work(ip_MaxDens+iiQ-1)
         Do jS = 1, iS
            jjQ = jS*(jS+1)/2
            XDm_jj = Work(ip_MaxDens+jjQ-1)
            ijQ=iS*(iS-1)/2+jS
            XDm_ij = Work(ip_MaxDens+ijQ-1)
            XDm_max = Max(XDm_ij,XDm_ii,XDm_jj)
            Aint_ij=TMax_Valence(iS,jS)
            If (TMax_All*Aint_ij .ge. CutInt) Then
*
* --- FAQ: more aggressive screening to discard shprs formed
*          by AOs contributing mainly to the virtual MO space.
*
               If (Aint_ij*XDm_max .ge. CutInt) Then
                  nSkal2 = nSkal2 + 1
                  iWork((nSkal2-1)*2+ip_ij  )=iS
                  iWork((nSkal2-1)*2+ip_ij+1)=jS
               End If
*
            End If
         End Do
      End Do
*
*2)   For the auxiliary basis set
*
      If (Do_RI) Then
         mij = nSkal_Auxiliary*2
         Call GetMem('ip_ij','Allo','Inte',ip_ij2,mij)
         nij=0
         Do jS = nSkal_Valence+1, nSkal-1
            If (TMax_All*TMax_Auxiliary(jS-nSkal_Valence).ge.CutInt)
     &         Then
               nij = nij + 1
               iWork((nij-1)*2+ip_ij2  )=nSkal
               iWork((nij-1)*2+ip_ij2+1)=jS
            End If
         End Do
      Else
         mij = 2*nij_Eff
         Call GetMem('ip_ij','Allo','Inte',ip_ij2,mij)
         nij=0
         Do ij = 1, nij_Eff
            iS = iWork(ip_ij3-1 + (ij-1)*2 +1)
            jS = iWork(ip_ij3-1 + (ij-1)*2 +2)
            If (TMax_All*TMax_Valence(iS,jS).ge.CutInt) Then
               nij = nij + 1
               iWork((nij-1)*2+ip_ij2  )=iS
               iWork((nij-1)*2+ip_ij2+1)=jS
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
         Call GetMem('CijK','Allo','Real',ip_CijK,lCilK)
         If (lPSO) lCilK=Max(lCilK,maxnAct) ! used as scratch
         Call GetMem('CilK','Allo','Real',ip_CilK,lCilK)
         lBklK = MxInShl**2*MxChVInShl
         Call GetMem('BklK','Allo','Real',ip_BklK,lBklK)
*
         If(iMp2prpt .eq. 2) Then
            lB_mp2 = mxChVInShl*nBas(0)*nBas(0)
            Do i = 1, 2
               Call GetMem('B_mp2','Allo','Real',ip_B_mp2(i),lB_mp2)
            End Do
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
         lCVec  = nIJ1Max*nK
         lCVec2 = nIJRMax*nK
         Call GetMem('C_Vector','Allo','Real',ip_CVec,  lCVec )
         Call GetMem('C_Vector2','Allo','Real',ip_CVec2,lCVec2)
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
              Call dDaFile(LuCVector(kSym,iSO),2,Work(ip_CVec),
     &                     nIJ1(iSym,jSym,iSO)*nK,iAdrC)
*
              ni = nChOrb(iSym-1,iSO)
*             nj=ni
              index1 = ip_CVec
              Do KAux = 1, nK
                 index_aux = ip_CVec2+(KAux-1)*(ni*ni)
*
                 Do i = 1, ni
                    Do j = 1, i-1
                       index2 = index_aux + i-1 + (j-1)*ni
                       index3 = index_aux + j-1 + (i-1)*ni
                       Work(index2) = Work(index1)
                       Work(index3) = Work(index1)
                       index1 = index1 + 1
                    End Do
                    index2 = index_aux + i-1 + (j-1)*ni
                    Work(index2) = Work(index1)*Sqrt(Two)
                    index1 = index1 + 1
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
              Call dDaFile(LuCVector(kSym,iSO),1,Work(ip_CVec2),
     &                     nIJR(iSym,jSym,iSO)*nK,iAdrC)
*
           End Do
 15        Continue
         End Do
*
         Call GetMem('C_Vector','Free','Real',ip_CVec,  lCVec )
         Call GetMem('C_Vector2','Free','Real',ip_CVec2,lCVec2)
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
         Call GetMem('Ymnij','Allo','Inte',ipYmnij(1),NumOrb)
         Call IZero(iWork(ipYmnij(1)),NumOrb)
         Do i=2,nKDens
           ipYmnij(i)=ipYmnij(1)+ipYmnij(i)
         End Do
*
*        Make a list for each shell-pair over the largest element
*        SQRT(ABS( (mu,nu|mu,nu) ))
*
         nnSkal = nSkal_valence*(nSkal_valence+1)/2
         Call GetMem('MaxDG','Allo','Real',ipSDG,nnSkal)
         Call get_maxDG(Work(ipSDG),nnSkal,MxBasSh)
*
*        Scratch for reduced lists of X_mi. Used in pget.
*
         nXki=MumOrb*MxBasSh*nSym
         Call GetMem('MOs_Yij','Allo','Real',jr_Xki(1),2*nKDens*nXki)
         jr_Xli(1)=jr_Xki(1)+nXki
         Do i=2,nKDens
           jr_Xki(i)=jr_Xki(i-1)+2*nXki
           jr_Xli(i)=jr_Xki(i)+nXki
         End Do
*
*        Make a list the largest element X_mu,i for each valence shell
*        and a fixed i. X_mu,i defined in Eq. 13.
*
         l_MaxXi=MumOrb*nSkal_valence*nIrrep ! X_i,ishell,isym
         Call GetMem('MaxXi','Allo','Real',ipXmi(1),l_MaxXi*nKDens)
         Do i=2,nKDens
           ipXmi(i)=ipXmi(i-1)+l_MaxXi
         End Do
         l_MaxXi=l_MaxXi*nKDens
*
         Do iSO=1,nKDens
            Call get_mXOs(iSO,Work(ipXmi(iSO)),MumOrb,nSkal_valence,
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
         Call GetMem('Thhalf','Allo','Real',ip_Thhalf,maxnnP)
         nThpkl=MxChVInShl*MxInShl**2
         Call mma_allocate(Thpkl,nThpkl,Label='Thpkl')
*
         Call contract_Zpk_Tpxy(Z_p_k ,nZ_p_k,
     &                          Txy   ,n_Txy,
     &                          Work(ip_Thhalf),maxnnP,
     &                          DMdiag ,nG1,
     &                          nnP,nBas_Aux,
     &                          nADens,nAvec,nAct,nIrrep)
*
         Call GetMem('Thhalf','Free','Real',ip_Thhalf,maxnnP)
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
      ip_LB=ip_iDummy
      If (iOpt.eq.1) Then
         Call qpg_iArray('LBList',Found,nSkal2_)
         If (Found) Then
            Call Allocate_iWork(ip_LB,nSkal2_)
            Call Get_iArray('LBList',iWork(ip_LB),nSkal2_)
         Else
            Call WarningMessage(2,'LBList not found!')
            Call Abend()
         End If
      End If
      Call Init_Tsk2(id,nSkal2,iOpt,iWork(ip_LB))
      If (iOpt.eq.1) Call Free_iWork(ip_LB)
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
         kS = iWork((klS-1)*2+ip_ij  )
         lS = iWork((klS-1)*2+ip_ij+1)
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
            SDGmn=Work(ipSDG+klS_-1)
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
                NumOrb_j  = nChOrb(jSym-1,iMOleft)-1
*
                mj=0
                Do jb=0,NumOrb_j
                  ip_Xmi=ipXmi(iMOleft) + MumOrb*nSkal_valence*(jSym-1)
                  Xjk=Work(ip_Xmi+MumOrb*(kS-1)+jb)
                  Xjl=Work(ip_Xmi+MumOrb*(lS-1)+jb)
*
                  jSym_s = jSym
                  if (ks.ne.ls.or.iMOright.ne.iMOleft) jSym_s=1
                  Do iSym = jsym_s, nSym
                     NumOrb_i  = nChOrb(iSym-1,iMOright)-1
                     If (iSym.eq.jSym.and.ks.eq.ls
     &                    .and.iMOright.eq.iMOleft) NumOrb_i = jb
*
                     Do ib = NumOrb_i, 0, -1
                       jp_Xmi=ipXmi(iMOright)+
     &                        MumOrb*nSkal_valence*(iSym-1)
                       Xik=Work(jp_Xmi+MumOrb*(kS-1)+ib)
                       Xil=Work(jp_Xmi+MumOrb*(lS-1)+ib)
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
                         iWork(ipYmnij(iMOleft)+mj+nj)=jb+1
                         mj=mj+1
                         Go To 666
                       End If
                     End Do ! ib
                  End Do    ! iSym
 666              Continue
                End Do
*
*               Trick used in pget to skip
*               If (mj.eq.NumOrb_j .and. xfk.gt.zero) mj=-mj
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
         iS = iWork((ijS-1)*2+ip_ij2  )
         jS = iWork((ijS-1)*2+ip_ij2+1)
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
         Call Size_SO_block_g(iSD4,nSD,Petite,nSO,No_batch)
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
     &               nExp,nBasis,MxShll,
     &               iBsInc,jBsInc,kBsInc,lBsInc,
     &               iPrInc,jPrInc,kPrInc,lPrInc,
     &               ipMem1,ipMem2, Mem1,  Mem2,
     &               iPrint,iFnc, MemPSO)
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
     &                 nExp,MxShll,
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
           Call PGet0(iCmpa,iShela,
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
     &                  Temp,nGrad,lIrrep,ChDisp,iPrint)
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
         Call GetMem('MaxXi','Free','Real',ipXmi(1),l_MaxXi)
         Call GetMem('MOs_Yij','Free','Real',jr_Xki(1),2*nKVec*nXki)
         Call GetMem('MaxDG','Free','Real',ipSDG,nnSkal)
         Call GetMem('Ymnij','Free','Inte',ipYmnij(1),NumOrb)
         Call GetMem('CijK','Free','Real',ip_CijK,lCijK)
         Call GetMem('CilK','Free','Real',ip_CilK,lCilK)
         Call GetMem('BklK','Free','Real',ip_BklK,lBklK)
         Call GetMem('ijList','Free','Inte',ipijList,lijList)
         Call GetMem('ijListTri','Free','Inte',ipijListTri,lijList)
         Call GetMem('JKVEC','Free','Real',ip_VJ,ljkVec)
         Do i=1,nKDens
           If (lCMOi(i).gt.0) Then
              Call GetMem('CMO_inv','FREE','Real',ip_CMOi(i), lCMOi(i))
           End If
         End DO
      End If
      Call GetMem('MaxDensity','Free','Real',ip_MaxDens,lMaxDens)
*
      If(iMp2prpt .eq. 2) Then
         Do i = 1, 2
            Call GetMem('B_mp2','Free','Real',ip_B_mp2(i),lB_mp2)
         End Do
      End If
*
      If (Allocated(Thpkl)) Call mma_deallocate(Thpkl)
*
      Call mma_deallocate(Sew_Scr)
      Call Free_Tsk2(id)
      Call GetMem('ip_ij','Free','Inte',ip_ij2,mij)
      Call GetMem('ip_ij','Free','Inte',ip_ij,nSkal*(nSkal+1))
      Call GetMem('TMax','Free','Real',ipTMax,nSkal**2)
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
      Call QExit('Drvg1_3Center_RI')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Grad)
      End

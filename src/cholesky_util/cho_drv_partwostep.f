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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_Drv_ParTwoStep(irc)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: Parallel two-step decomposition of two-electron
C              integrals.
C
      use ChoArr, only: iAtomShl
      Implicit None
      Integer irc
#include "choprint.fh"
#include "cholesky.fh"
#include "chosew.fh"
#include "chosubscr.fh"
#include "chosimri.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Integer ip_Err, l_Err
      Integer iSec, kDiag
      Integer BlockSize_Bak, iPrint_Bak, Cho_IOVec_Bak, N1_VecRD_Bak
      Integer N2_VecRd_Bak
      Integer Cho_DecAlg_Bak
      Integer nSys_Call_Bak, nDGM_Call_Bak
      Integer iSym, n, ni, nj, nnBlock, i1
      Integer ip_NVT, l_NVT
      Integer ip_nBlock, l_nBlock
      Integer nB, nB_Max
      Integer ip_Z, l_Z
      Integer ip_nVBlock, l_nVBlock
      Integer ip_iV1Block, l_iV1Block
      Integer ip_ZBlock, l_ZBLock
      Integer iBlock, jBlock, ijBlock
      Integer MinQual_Bak, MaxQual_Bak, N1_Qual_Bak, N2_Qual_Bak
      Integer MxShPr_Bak, iAlQua_Bak
      Integer N_Subtr_Bak
      Integer Mode_Screen_Bak, Cho_DecAlg_Def_Bak, ModRst_Bak
      Integer iDum

      Real*8 tCPU0,  tCPU1,  tC0, tC1
      Real*8 tWall0, tWall1, tW0, tW1
      Real*8 C0, C1, W0, W1
      Real*8 Thr_PreScreen_Bak, ThrDiag_Bak, Frac_ChVBuf_Bak, SSTau_Bak
      Real*8 Thr_SimRI_Bak, Tol_DiaChk_Bak
      Real*8 TimSec_Bak(4,nSection)
      Real*8 tInteg_Bak(2,nInteg)
      Real*8 tDecom_Bak(2,nDecom)
      Real*8 tMisc_Bak(2,nMisc)
      Real*8 Byte

      Logical lConv, Free_Z
      Logical Cho_1Center_Bak, Cho_No2Center_Bak, Cho_PreScreen_Bak
      Logical ScDiag_Bak, Cho_SScreen_Bak, Cho_SimRI_Bak, HaltIt_Bak
      Logical Did_DecDrv_Bak
      Logical Cho_UseAbs_Bak, Cho_DiaChk_Bak, Cho_Fake_Par_Bak
      Logical Cho_SimP_Bak, Cho_ReOrd_Bak, ChkOnly_Bak
      Logical Cho_IntChk_Bak, Cho_MinChk_Bak, Cho_TrcNeg_Bak
      Logical Cho_TstScreen_Bak, RstDia_Bak, RstCho_Bak
      Logical Trace_Idle_Bak

      Character(LEN=18), Parameter:: SecNam='Cho_Drv_ParTwoStep'
      Character(LEN=4), Parameter:: myName='DPTS'
      Character(LEN=2) Unt

      Real*8, Parameter:: DumTst=0.123456789d0, DumTol=1.0d-15
      Real*8, Allocatable:: Check(:)

      Integer NVT, nBlock, iTri
      Integer i, j

      NVT(i)=iWork(ip_NVT-1+i)
      nBlock(i)=iWork(ip_nBlock-1+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

C     Preliminaries.
C     ==============

#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem('Start of '//SecNam)
#endif

      ! Start overall timing
      If (iPrint .ge. Inf_Timing) Call Cho_Timer(tCPU0,tWall0)

      ! Init return code
      irc=0

      ! make a dummy allocation
      Call mma_allocate(Check,1,Label='Check')
      Check(1)=DumTst

      ! init local variables
      kDiag = 0
      lConv = .False.

C     Initialization.
C     ===============

      iSec=1
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
      End If
      Call Cho_Init(.False.,.True.)
      Call Cho_GASync()
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
         Call Cho_PrtTim('Cholesky initialization',
     &                   TimSec(2,iSec),TimSec(1,iSec),
     &                   TimSec(4,iSec),TimSec(3,iSec),
     &                   1)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After init')
#endif

C     Get diagonal.
C     =============

      iSec=2
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
         Write(LuPri,'(/,A)')
     &   '***** Starting Cholesky diagonal setup *****'
         Call Cho_Flush(LuPri)
      End If
      Call Cho_GetDiag(kDiag,lConv)
      Call Cho_GASync()
      If (lConv) Then
         ! restart is not possible, so it can not be converged!!
         Write(LuPri,'(A,A)')
     &   SecNam,': logical error: converged but not restart?!?!'
         Call Cho_Quit('Error in '//SecNam,103)
      End If
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
         Call Cho_PrtTim('Cholesky diagonal setup',
     &                   TimSec(2,iSec),TimSec(1,iSec),
     &                   TimSec(4,iSec),TimSec(3,iSec),
     &                   1)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After diagonal')
#endif

C     Cholesky decomposition in two steps.
C     ====================================

      ! Time the rest as "decomposition driver"
      Call Cho_Timer(C0,W0)

      iSec=3
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
         Write(LuPri,'(/,A)')
     &   '***** Starting Cholesky decomposition *****'
         Call Cho_Flush(LuPri)
      End If

      ! Perform first CD step: get parent diagonals and Z vectors.
      ! Z vectors contain more elements than needed at this stage!
      ! (more elements than vectors). The final Z vectors are obtained
      ! below (Cho_GetZ).
      Call Cho_P_SetAddr()
      Call Cho_DecDrv(Work(kDiag))
      Call Cho_GASync()
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC1,tW1)
         Call Cho_PrtTim('Cholesky map generation',
     &                   tC1,TimSec(1,iSec),
     &                   tW1,TimSec(3,iSec),
     &                   2)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After 1st step')
#endif

      ! Shut down and re-initialize.
      ! This is done to get rid of the vast amount of different data
      ! generated by the parallel decomposition in particular.
      ! In addition, it will make it easier to use the vector calculator
      ! externally (i.e. after seward execution)
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC0,tW0)
      End If
      Call dCopy_(4*nSection,TimSec,1,Timsec_Bak,1)
      Call dCopy_(2*nInteg,tInteg,1,tInteg_Bak,1)
      Call dCopy_(2*nDecom,tDecom,1,tDecom_Bak,1)
      Call dCopy_(2*nMisc,tMisc,1,tMisc_Bak,1)
      Trace_Idle_Bak=Trace_Idle
      Cho_UseAbs_Bak=Cho_UseAbs
      Cho_DiaChk_Bak=Cho_DiaChk
      Cho_Fake_Par_Bak=Cho_Fake_Par
      Cho_SimP_Bak=Cho_SimP
      Cho_ReOrd_Bak=Cho_ReOrd
      ChkOnly_Bak=ChkOnly
      Cho_IntChk_Bak=Cho_IntChk
      Cho_MinChk_Bak=Cho_MinChk
      Cho_TrcNeg_Bak=Cho_TrcNeg
      Cho_TstScreen_Bak=Cho_TstScreen
      RstDia_Bak=RstDia
      RstCho_Bak=RstCho
      Mode_Screen_Bak=Mode_Screen
      Cho_DecAlg_Def_Bak=Cho_DecAlg_Def
      ModRst_Bak=ModRst
      N_Subtr_Bak=N_Subtr
      Did_DecDrv_Bak=Did_DecDrv
      HaltIt_Bak=HaltIt
      Tol_DiaChk_Bak=Tol_DiaChk
      Cho_1Center_Bak=Cho_1Center
      Cho_No2Center_Bak=Cho_No2Center
      Cho_PreScreen_Bak=Cho_PreScreen
      Thr_PreScreen_Bak=Thr_PreScreen
      ThrDiag_Bak=ThrDiag
      ScDiag_Bak=ScDiag
      MinQual_Bak=MinQual
      MaxQual_Bak=MaxQual
      N1_Qual_Bak=N1_Qual
      N2_Qual_Bak=N2_Qual
      MxShPr_Bak=MxShPr
      iAlQua_Bak=iAlQua
      Cho_IOVec_Bak=Cho_IOVec
      Frac_ChVBuf_Bak=Frac_ChVBuf
      Cho_SScreen_Bak=Cho_SScreen
      SStau_Bak=SStau
      Cho_SimRI_Bak=Cho_SimRI
      Thr_SimRI_Bak=Thr_SimRI
      Cho_DecAlg_Bak=Cho_DecAlg
      BlockSize_Bak=BlockSize
      iPrint_Bak=iPrint
      Cho_IOVec_Bak=Cho_IOVec
      N1_VecRd_Bak=N1_VecRd
      N2_VecRd_Bak=N2_VecRd
      nSys_Call_Bak=nSys_Call
      nDGM_Call_Bak=nDGM_Call
      Call Cho_P_WrDiag()
      Call Cho_Final(.True.)
      Call Cho_P_OpenVR(2)
      Call Cho_X_Init(irc,0.0d0)
      If (irc .ne. 0) Then
         Write(LuPri,*) SecNam,': Cho_X_Init returned code ',irc
         irc = 1
         Go To 1 ! clear memory and return
      End If
      If (Cho_1Center) Then
         If (.NOT.Allocated(iAtomShl)) Then
            Call mma_allocate(iAtomShl,nShell,Label='iAtomShl')
            Call Cho_SetAtomShl(irc,iAtomShl,SIZE(iAtomShl))
            If (irc.ne.0) Then
               Write(LuPri,'(A,A,I8)')
     &         SecNam,': Cho_SetAtomShl returned code',irc
               irc=1
               Go To 1 ! clear memory and return
            End If
         End If
      End If
      Call dCopy_(4*nSection,TimSec_Bak,1,Timsec,1)
      Call dCopy_(2*nInteg,tInteg_Bak,1,tInteg,1)
      Call dCopy_(2*nDecom,tDecom_Bak,1,tDecom,1)
      Call dCopy_(2*nMisc,tMisc_Bak,1,tMisc,1)
      Trace_Idle=Trace_Idle_Bak
      Cho_UseAbs=Cho_UseAbs_Bak
      Cho_DiaChk=Cho_DiaChk_Bak
      Cho_Fake_Par=Cho_Fake_Par_Bak
      Cho_SimP=Cho_SimP_Bak
      Cho_ReOrd=Cho_ReOrd_Bak
      ChkOnly=ChkOnly_Bak
      Cho_IntChk=Cho_IntChk_Bak
      Cho_MinChk=Cho_MinChk_Bak
      Cho_TrcNeg=Cho_TrcNeg_Bak
      Cho_TstScreen=Cho_TstScreen_Bak
      RstDia=RstDia_Bak
      RstCho=RstCho_Bak
      Mode_Screen=Mode_Screen_Bak
      Cho_DecAlg_Def=Cho_DecAlg_Def_Bak
      ModRst=ModRst_Bak
      N_Subtr=N_Subtr_Bak
      Did_DecDrv=Did_DecDrv_Bak
      HaltIt=HaltIt_Bak
      Tol_DiaChk=Tol_DiaChk_Bak
      Cho_1Center=Cho_1Center_Bak
      Cho_No2Center=Cho_No2Center_Bak
      Cho_PreScreen=Cho_PreScreen_Bak
      Thr_PreScreen=Thr_PreScreen_Bak
      ThrDiag=ThrDiag_Bak
      ScDiag=ScDiag_Bak
      MinQual=MinQual_Bak
      MaxQual=MaxQual_Bak
      N1_Qual=N1_Qual_Bak
      N2_Qual=N2_Qual_Bak
      N1_VecRd=N1_VecRd_Bak
      N2_VecRd=N2_VecRd_Bak
      MxShPr=MxShPr_Bak
      iAlQua=iAlQua_Bak
      Cho_IOVec=Cho_IOVec_Bak
      Frac_ChVBuf=Frac_ChVBuf_Bak
      Cho_SScreen=Cho_SScreen_Bak
      SStau=SStau_Bak
      Cho_SimRI=Cho_SimRI_Bak
      Thr_SimRI=Thr_SimRI_Bak
      Cho_DecAlg=Cho_DecAlg_Bak
      BlockSize=BlockSize_Bak
      iPrint=iPrint_Bak
      Cho_IOVec=Cho_IOVec_Bak
      N1_VecRd=N1_VecRd_Bak
      N2_VecRd=N2_VecRd_Bak
      nSys_Call=nSys_Call_Bak
      nDGM_Call=nDGM_Call_Bak
      ! Allocate memory for extracting integral columns directly in
      ! reduced set from Seward. Set Seward interface to 3 (to treat
      ! columns shell pair-wise).
      IFCSEW=3
      l_iShP2RS=2*Mx2Sh
      l_iShP2Q=l_iShP2RS
      Call GetMem('ShP2RS','Allo','Inte',ip_iShP2RS,l_iShP2RS)
      Call GetMem('ShP2Q','Allo','Inte',ip_iShP2Q,l_iShP2Q)
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC1,tW1)
         Call Cho_PrtTim('Cholesky reinitialization',tC1,tC0,tW1,tW0,2)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After re-init')
#endif

      ! Get Z vectors in core
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC0,tW0)
      End If
      l_NVT=nSym
      Call GetMem('TotV','Allo','Inte',ip_NVT,l_NVT)
      Call Cho_X_GetTotV(iWork(ip_NVT),l_NVT)
      Call Cho_ZMem(irc,l_Z,iWork(ip_NVT),l_NVT,iPrint.ge.Inf_Timing,
     &              .True.)
      If (irc.ne.0) Then
         Write(LuPri,'(A,A,I6)')
     &   SecNam,': Cho_ZMem returned code',irc
         If (irc.eq.999) Then
            If (iPrint.lt.Inf_Timing) Then
               Call Cho_ZMem(irc,l_Z,iWork(ip_NVT),l_NVT,.True.,.False.)
            End If
            Call GetMem('ZMax','Max ','Real',ip_Z,l_Z)
            Call Cho_Word2Byte(l_Z,8,Byte,Unt)
            Write(LuPri,'(A,I12,A,F7.3,1X,A,A)')
     &      'Largest available memory block:',l_Z,' words (',
     &      Byte,Unt,')'
            Write(LuPri,'(A)')
     &      '=> INSUFFICIENT MEMORY FOR STORING Z VECTORS!'
            Write(LuPri,'(/,A,/,A)')
     &      'You have the following options:',
     &      '(a) Increase available memory (MOLCAS_MEM),'
            Write(LuPri,'(A)')
     &      'and/or'
            Write(LuPri,'(A,/,A,A,/,A,/,A,/,A,/,A)')
     &      '(b) -SERIAL EXECUTION:',
     &      '       Use the serial two-step algorithm by specifying ',
     &      'the keywords',
     &      '          ChoInput',
     &      '          TwoStep',
     &      '          EndChoInput',
     &      '       in Seward input.'
            Write(LuPri,'(A,/,A,A,/,A,/,A,/,A,/,A,/,A)')
     &      '    -PARALLEL EXECUTION:',
     &      '       Use the parallel one-step algorithm by specifying ',
     &      'the keywords',
     &      '          ChoInput',
     &      '          OneStep',
     &      '          Parallel',
     &      '          EndChoInput',
     &      '       in Seward input.'
            Call Cho_Quit(SecNam//': Insufficient memory for Z vectors',
     &                    101)
         End If
         irc=1
         Go To 1 ! clear memory and return
      End If
      Call GetMem('Z','Allo','Real',ip_Z,l_Z)
      l_nBlock=nSym
      Call GetMem('nBlock','Allo','Inte',ip_nBlock,l_nBlock)
      nB_Max=0
      Do iSym=1,nSym
         nB=(NVT(iSym)-1)/BlockSize+1
         nB_Max=max(nB_Max,nB)
         iWork(ip_nBlock-1+iSym)=nB
      End Do
      l_nVBlock=nB_Max*nSym
      l_iV1Block=nB_Max*nSym
      Call GetMem('nVBlock','Allo','Inte',ip_nVBlock,l_nVBlock)
      Call GetMem('iV1Block','Allo','Inte',ip_iV1Block,l_iV1Block)
      Call iZero(iWork(ip_nVBlock),l_nVBlock)
      Call iZero(iWork(ip_iV1Block),l_iV1Block)
      Do iSym=1,nSym
         i1=1
         Do iBlock=1,nBlock(iSym)-1
            iWork(ip_nVBlock-1+nB_Max*(iSym-1)+iBlock)=BlockSize
            iWork(ip_iV1Block-1+nB_Max*(iSym-1)+iBlock)=i1
            i1=i1+BlockSize
         End Do
         iWork(ip_nVBlock-1+nB_Max*(iSym-1)+nBlock(iSym))=NVT(iSym)
     &                                       -BlockSize*(nBlock(iSym)-1)
         iWork(ip_iV1Block-1+nB_Max*(iSym-1)+nBlock(iSym))=i1
      End Do
      nnBlock=nB_Max*(nB_Max+1)/2
      l_ZBlock=nnBlock*nSym
      Call GetMem('ZBlock','Allo','Inte',ip_ZBlock,l_ZBlock)
      Call iZero(iWork(ip_ZBlock),l_ZBlock)
      n=ip_Z
      Do iSym=1,nSym
         Do jBlock=1,nBlock(iSym)
            nj=iWork(ip_nVBlock-1+nB_Max*(iSym-1)+jBlock)
            ijBlock=iTri(jBlock,jBlock)
            iWork(ip_ZBlock-1+nnBlock*(iSym-1)+ijBlock)=n
            n=n+nj*(nj+1)/2
            Do iBlock=jBlock+1,nBlock(iSym)
               ni=iWork(ip_nVBlock-1+nB_Max*(iSym-1)+iBlock)
               ijBlock=iTri(iBlock,jBlock)
               iWork(ip_ZBlock-1+nnBlock*(iSym-1)+ijBlock)=n
               n=n+ni*nj
            End Do
         End Do
      End Do
      Call Cho_GetZ(irc,
     &              iWork(ip_NVT),l_NVT,
     &              iWork(ip_nBlock),l_nBlock,
     &              iWork(ip_nVBlock),nB_Max,nSym,
     &              iWork(ip_iV1Block),nB_Max,nSym,
     &              iWork(ip_ZBlock),nnBlock,nSym)
      If (irc .ne. 0) Then
         Write(LuPri,*) SecNam,': Cho_GetZ returned code ',irc
         irc = 1
         Go To 1 ! clear memory and return
      End If
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC1,tW1)
         Call Cho_PrtTim('Cholesky Z vector fetching',tC1,tC0,tW1,tW0,
     &                   2)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After GetZ')
#endif

      ! Perform second step: calculation of Cholesky vectors from Z
      ! vectors and integrals.
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC0,tW0)
      End If
      Free_Z=.True. ! let Cho_X_CompVec deallocate Z array
      Call Cho_X_CompVec(irc,
     &                   iWork(ip_NVT),l_NVT,
     &                   iWork(ip_nBlock),l_nBlock,
     &                   iWork(ip_nVBlock),nB_Max,nSym,
     &                   iWork(ip_iV1Block),nB_Max,nSym,
     &                   iWork(ip_ZBlock),nnBlock,nSym,
     &                   Free_Z)
      If (irc .ne. 0) Then
         Write(LuPri,*) SecNam,': Cho_X_CompVec returned code ',irc
         irc = 1
         Go To 1 ! clear memory and return
      End If
      ! Write restart files
      Call Cho_PTS_WrRst(irc,iWork(ip_NVT),l_NVT)
      If (irc .ne. 0) Then
         Write(LuPri,*) SecNam,': Cho_PTS_WrRst returned code ',irc
         irc = 1
         Go To 1 ! clear memory and return
      End If
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(tC1,tW1)
         Call Cho_PrtTim('Cholesky vector generation',tC1,tC0,tW1,tW0,
     &                   2)
      End If

      ! Final timing of decomposition section
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
         Call Cho_PrtTim('Cholesky decomposition',
     &                   TimSec(2,iSec),TimSec(1,iSec),
     &                   TimSec(4,iSec),TimSec(3,iSec),
     &                   1)
      End If

      ! time as "decomposition driver"
      Call Cho_Timer(C1,W1)
      tDecDrv(1)=C1-C0
      tDecDrv(2)=W1-W0
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After 2nd step.')
#endif

C     Check diagonal.
C     ===============

      iSec=4
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
         Write(LuPri,'(/,A)')
     &   '***** Starting Cholesky diagonal check *****'
         Call Cho_Flush(LuPri)
      End If
      l_Err=4
      Call GetMem('DiaErr','Allo','Real',ip_Err,l_Err)
      iPrint_Bak=iPrint
      If (iPrint.lt.Inf_Pass) Then
         iPrint=-99999999 ! suppress printing in Cho_X_CheckDiag
      End If
      Call Cho_X_CheckDiag(irc,Work(ip_Err))
      iPrint=iPrint_Bak
      If (irc .ne. 0) Then
         Write(LuPri,*) SecNam,': Cho_X_CheckDiag returned code ',irc
         irc = 1
         Go To 1 ! release memory and return
      End If
      If (Work(ip_Err+1) .gt. ThrCom) Then
         Write(LuPri,'(/,A)')
     &   'Cholesky decomposition failed!'
         Write(LuPri,'(3X,A,1P,D15.6)')
     &   'Largest integral diagonal..',Work(ip_Err+1)
         Write(LuPri,'(3X,A,1P,D15.6)')
     &   'Decomposition threshold....',ThrCom
         irc=1
         Go To 1 ! release memory and return
      End If
      Call GetMem('DiaErr','Free','Real',ip_Err,l_Err)
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
         Call Cho_PrtTim('Cholesky diagonal check',
     &                   TimSec(2,iSec),TimSec(1,iSec),
     &                   TimSec(4,iSec),TimSec(3,iSec),
     &                   1)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After diagonal check.')
#endif

C     Finalization.
C     =============

      iSec=8
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
         Write(LuPri,'(/,A)')
     &   '***** Starting Cholesky finalization *****'
         Call Cho_Flush(LuPri)
      End If
      Call Cho_TrcIdl_Final()
      Call Cho_PTS_Final(iWork(ip_NVT),l_NVT)
      If (iPrint .ge. Inf_Timing) Then
         Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
         Call Cho_PrtTim('Cholesky finalization',
     &                   TimSec(2,iSec),TimSec(1,iSec),
     &                   TimSec(4,iSec),TimSec(3,iSec),
     &                   1)
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After finalization')
#endif

C     Statistics.
C     ===========

      If (iPrint. ge. 1) Then
         iSec = 9
         If (iPrint .ge. Inf_Timing) Then
            Call Cho_Timer(TimSec(1,iSec),TimSec(3,iSec))
            Write(LuPri,'(/,A)')
     &      '***** Starting Cholesky statistics *****'
            Call Cho_Flush(LUPRI)
         End If
         Call Cho_PTS_Stat()
         Call Cho_GASync()
         If (iPrint .ge. Inf_Timing) Then
            Call Cho_Timer(TimSec(2,iSec),TimSec(4,iSec))
            Call Cho_PrtTim('Cholesky statistics',
     &                      TimSec(2,iSec),TimSec(1,iSec),
     &                      TimSec(4,iSec),TimSec(3,iSec),
     &                      1)
         End If
      End If
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem(SecNam//': After statistics')
#endif

C     Wrap it up and return.
C     ======================

      ! Close vector and restart files
      Call Cho_OpenVR(2,2)

      ! error termination point
    1 Continue
      ! check memory
      If (Abs(DumTst-Check(1)) .gt. DumTol) Then
         Write(LuPri,*) SecNam,': memory has been out of bounds [2]'
         irc=2
      End If
      Call GETMEM('KDIAG','FREE','REAL',KDIAG,iDum)
      Call mma_deallocate(Check)

      ! Print total timing
      If (iPrint.ge.Inf_Timing .and. irc.eq.0) Then
         Call Cho_Timer(tCPU1,tWall1)
         Call Cho_PrtTim('Cholesky Procedure',tCPU1,tCPU0,
     &                   tWall1,tWall0,1)
      End If

      Call Cho_Flush(LuPri)
#if defined (_DEBUGPRINT_)
      Call Cho_PrtMaxMem('End of '//SecNam)
#endif

      End

************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v.2.1.  *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2017 Varinia Bernales and Roland Lindh                 *
************************************************************************
      Subroutine Seward_DMET(ireturn,DMET_s,DMET_f,DMET_h,DMET_g,nBfn)
c.... please do not move this line further down (required for fortsplit)
************************************************************************
*                                                                      *
* (c) Copyright 1989, 1990, 1991, and 1992                             *
* Roland Lindh                                                         *
* Dept. of Theoretical Chemistry                                       *
* University of Lund, SWEDEN                                           *
*                                                                      *
* and                                                                  *
*                                                                      *
* IBM                                                                  *
* International Business Machine Corporation                           *
*                                                                      *
* No part of this code may be copied or redistributed without the      *
* written permission of the copyright owner.                           *
* The copyright owner does not take any responsibility for any         *
* errors in the code or documentation.                                 *
*                                                                      *
* All rights reserved.                                                 *
*                                                                      *
************************************************************************
************************************************************************
* In 1867, William Seward, for 2 cents per acre, purchased             *
* Alaska, a valueless wasteland of ice and snow.                       *
*                                                                      *
* In 1990, Roland Lindh and Ungsik Ryu worked on molecular             *
* integral evaluation, an exhausted scientific area with no            *
* room for innovation.                                                 *
*                                                                      *
* Bowen Liu                                                            *
* April, 1990                                                          *
************************************************************************
************************************************************************
*                                                                      *
*  Object: Driver for the one and two electron integral program        *
*          SEWARD. SEWARD computes integrals for cartesian and         *
*          spherical harmonic gaussian basis functions.                *
*                                                                      *
*                                                                      *
* Called from: None                                                    *
*                                                                      *
* Calling    : QEnter                                                  *
*              SetUp0                                                  *
*              DmpInf                                                  *
*              Input_Seward                                            *
*              SetUp                                                   *
*              Drv1El                                                  *
*              Drv2El                                                  *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
*          July '89 - May '90                                          *
*                                                                      *
* (c) Copyright 1989 IBM                                               *
* All rights reserved.                                                 *
* (c) Copyright 1990 IBM                                               *
* All rights reserved.                                                 *
*                                                                      *
*          Roland Lindh, Dept. of Theoretical Chemistry, University of *
*          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
*          screening, July 1991.                                       *
*                                                                      *
************************************************************************
      use Real_Spherical
      use Period
      use GeoList
      use MpmC
      Implicit Real*8 (A-H,O-Z)
      real*8 DMET_s(nBfn,nBfn), DMET_f(nBfn,nBfn)
      real*8 DMET_h(nBfn,nBfn), DMET_g(nBfn**4)
      External Integral_WrOut, Integral_WrOut2, Integral_RI_3
      Real*8, Dimension(:), Allocatable :: MemHide
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "lundio.fh"
#include "print.fh"
#include "gateway.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
      Integer iix(2), nChoV(8), GB
      Logical PrPrt_Save, Exist, DoRys, lOPTO
      Real*8  DiagErr(4), Dummy(2)
C-SVC: identify runfile with a fingerprint
      Character cDNA*256
      Logical IsBorn, Do_OneEl
      Integer IsGvMode
*                                                                      *
************************************************************************
*                                                                      *
      lOPTO = .False.
      nByte = iiLoc(iix(2)) - iiLoc(iix(1))
      Call CWTime(TCpu1,TWall1)
*
*     Prologue
*
      iRout=1
      Call qEnter('Seward')
      LuWr=6
      PrPrt_Save = .False. ! dummy initialize
*                                                                      *
************************************************************************
*                                                                      *
*     Figure out the run_mode
*
*     Seward can be run in two different modes!
*     GS_Mode: does the work of both Gateway and Seward
*     S_Mode:  only the work of Seward
*
*
*     Check if the run file is there
*
      Call f_Inquire('RUNFILE',Exist)
      If (Exist) Then
         Call Qpg_iScalar('Run_Mode',Exist)
         If (Exist) Then
*
*           The Run_mode of the runfile is either GS_Mode or G_Mode
*
            Call Get_iScalar('Run_Mode',Run_Mode)
*
*           If the Run_mode is that Gateway is in action then Seward
*           should be run in S_mode.
*
            If (Run_Mode.eq.G_Mode) Run_Mode=S_Mode
         Else
            Run_Mode=GS_Mode
         End If
      Else
*
*        Seward runs without Gateway
*
         Run_Mode=GS_Mode
         Call MkRun(iRC,0)
         Call Put_iScalar('Run_Mode',Run_Mode)
*
*     Determine and save the fingerprint of the runfile in a field with
*     label 'BirthCertificate' if it is empty.  This allows us to
*     uniquely identify the runfile and any later associated files.
*
         Call qpg_cArray('BirthCertificate',IsBorn,nDNA)
         If (.NOT.IsBorn) Then
           Call Get_Genome(cDNA,nDNA)
           Call Put_cArray('BirthCertificate',cDNA,nDNA)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Get the memory size available
*
      Call SetMem('Clear=Off')
*                                                                      *
************************************************************************
*                                                                      *
*     If Seward is run in S_mode most of the input is already on the
*     runfile. If Seward is run in GS_Mode it will handle the input and
*     runfile in the conventional way.
*
      If (Run_Mode.eq.S_Mode) Then
*
*        S_Mode
*
         Call Seward_Init()
         DoRys=.True.
         nDiff=0
         Call GetInf(Info,nInfo,DoRys,nDiff,1)
         Primitive_Pass=.True.
*                                                                      *
************************************************************************
*                                                                      *
      Else
*
*        GS_Mode
*                                                                      *
************************************************************************
*                                                                      *
*        Initialize common blocks
*
         Call Seward_Init()
         Call Funi_Init()
*                                                                      *
************************************************************************
*                                                                      *
*        Call GetMem to get pointer to first available core allocation.
*
         kB=2**10
         MB=kb*kB
         GB=kb*MB/8 ! adjust to real*8
         Call GetMem('Info','Max','Real',iDum,MaxM)
         nDInf=Max(MaxM/4,Min((9*MaxM)/10,GB))
         Call GetMem('Info','ALLO','REAL',Info,nDInf)
         Call FZero(Work(Info),nDInf)
         Info_Status=Active
         LctInf = Info
         nInfo = 0
*
      End If ! Run_Mode.eq.S_Mode

************************************************************************
*   columbus support: initialize additional items in Runfile
*   default: no mixed operation
      call Put_iScalar('Columbus',0)
      call Put_iScalar('colgradmode',0)
      dummy(1)=0.0d0
      dummy(2)=0.0d0
      call Put_dArray ('MR-CISD energy',dummy,2)
      Call NQGrid_Init()
*                                                                      *
************************************************************************
*                                                                      *
*     Spool the input
*
      LuSpool=21
      Call SpoolInp(LuSpool)
*     Read the input from input file
*
      Call RdCtl_DMET(Info,nInfo,LuSpool,lOPTO,Do_OneEl,
     &                  Work(Info),nDInf,nBfn)
      Call GvMode(IsGvMode)
      if(IsGvMode.gt.0) Onenly=.true.
*
      Call Close_LuSpool(LuSpool)
*
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
 199  Continue
*
*     Process the input.
*
      Call Input_Seward(lOPTO,Info,Work(Info),nDInf)
*
      If (Primitive_Pass) Then
         PrPrt_Save = PrPrt
         PrPrt=.False.
      Else
         PrPrt=PrPrt_Save
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Compute the Nuclear potential energy
*
*      If (.Not.Primitive_Pass) Then
*         Call Gen_RelPointers(-(Info-1))
*         Call DrvN0(Work(Info),nInfo)
*         Call Gen_RelPointers(Info-1)
*      End If
         Call DrvN0_DMET(ireturn)
*                                                                      *
************************************************************************
*                                                                      *
      If (Show) Then
*
*        Print out basis set information
*
         Write(6,*)
         Write(6,'(6X,A)')'Basis set specifications :'
         Write(6,'(6X,A,T30,8(2X,A))')
     &         'Symmetry species',     (lIrrep(i),i=0,nIrrep-1)
         Write(6,'(6X,A,T30,8I5)')'Basis functions',
     &                                 (nBas(i),i=0,nIrrep-1)
         Write(6,*)
*
      End If

*                                                                      *
************************************************************************
*                                                                      *
*     If only test case then clean up!
*
      If (Test) Go To 9999
*                                                                      *
************************************************************************
*                                                                      *
*     Write/update information on the run file.
*
      If (.Not.Primitive_Pass) Then
         Call Gen_RelPointers(-(Info-1))
         Call DmpInf(Work(Info),nInfo)
         Call basis2run(Work(Info),nInfo)
         Call Gen_RelPointers(Info-1)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     ONE-ELECTRON INTEGRAL SECTION
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _FDE_
      ! Embedding
      if (embPot.and..not.embPotInBasis) then
       Call embPotInit(.false.)
      end if
#endif

      Lu_One=2
      iOpt = 1
      iRC = -1
*
*     Generate primimitive integrals only if needed.
*
      If (Primitive_Pass.and.(DKroll.or.Nemo)) Then
         Call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
         Call OneBas('PRIM')
      Else
         Call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
      End If
      If (iRC.ne.0) Then
         Call WarningMessage(2,
     &                  ' *** Error in subroutine INPUT ***;'
     &                //'     Abend in subroutine OpnOne')
         Call Abend()
      End If
*
      If (Do_OneEl.and.
     &    (.Not.Primitive_Pass .or.
     &    (Primitive_Pass.and.(DKroll.or.NEMO)) ) )
     &   Call Drv1El_DMET(DMET_s,DMET_f,DMET_h,nBfn)
*
      iOpt = 0
      iRC = -1
      Call ClsOne(iRC,iOpt)
      If (iRC.ne.0) then
         Call WarningMessage(2,
     &              ' *** Error in SEWARD main ***;'
     &            //'  Abend in subroutine ClsOne')
         Call Abend()
      End If

#ifdef _FDE_
      ! Embedding
      if (embPot.and..not.embPotInBasis) Call embPotFreeMem
#endif

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     If a pass in which primitive integrals where computed do a second
*     pass.
*
      If (Primitive_Pass) Then
         Primitive_Pass=.False.
         Call Free_iSD()
         Go To 199
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Branch out if only one-electron integrals are to be computed!
*
      If (Onenly) Go To 9999
*                                                                      *
*     If ERIs/CD/RI already available, one may want not to redo it!
*
      If (Fake_ERIs) Then
         Call set_fake_ERIs()
         Go To 9999
      EndIf
************************************************************************
************************************************************************
*                                                                      *
*     TWO-ELECTRON INTEGRAL SECTION
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*
      Call mma_allocate(MemHide,Memhid)
*
      If ( iWRopt.eq.0 ) then
*
*------- Molcas format
*
         If (Cholesky) Then ! Cholesky decomposition
            Call Cho_MCA_Drv()
            Call Get_iArray('NumCho',nChoV,nIrrep)
            Write(6,'(6X,A,T30,8I5)')'Cholesky vectors',
     &               (nChoV(i),i=1,nIrrep)
            Write(6,*)
            Write(6,*)
         Else If (Do_RI) Then
            If (LocalDF) Then
               Call Drv2El_LocalDF()
            Else
               If (nPrint(iRout).ge.6) Then
                  Write (6,*)
                  Write (6,'(A)') 'Seward processing 2-center and '
     &                          //'3-center ERIs'
                  Write (6,*)
               End If
               Call Drv2El_3Center_RI(Integral_RI_3,Zero)
               Call Get_iArray('NumCho',nChoV,nIrrep)
               If (nPrint(iRout).ge.6) Then
                  Write(6,'(6X,A,T30,8I5)')'RI vectors',
     &                  (nChoV(i),i=1,nIrrep)
                  Write(6,*)
                  Write(6,*)
               End If
            End If
         Else
            iWrOpt_Save=iWrOpt
            iWrOpt=0
            Call Sort0
*
            Call Drv2El_DMET(Integral_WrOut2,Zero,DMET_g,nBfn)
*
            Call Sort1B
            Call Sort2
            Call Sort3(MaxDax)
*
            If (nPrint(iRout).ge. 6) Then
               Write (6,*)
               Write (6,'(A)')
     &           ' Integrals are written in MOLCAS2 format'
               Write (6,'(A,I10)')
     &           ' Total Number of integrals             '//
     &           '                = ',IntTot
               Write (6,'(A,I10)')
     &           ' Number of nonzero integrals passed to '//
     &           'packing routine = ',NotZer
               If ( iPack.ne.0 ) Then
                  Write (6,'(A)')
     &              ' No packing of integrals has been applied'
               Else
                  Write (6,'(A,G10.4)') ' Packing accuracy = ',
     &                                   PkAcc
                  Write (6,'(A,I10)')
     &             ' Highest disk address written',MaxDax
               End If
               If ( iSquar.eq.0 ) Then
                  Write (6,'(A,A)') ' Diagonal and subdiagonal, '
     &              //'symmetry allowed 2-el',
     &              ' integral blocks are stored on Disk'
               Else
                  Write (6,'(A,A)') ' All symmetry allowed 2-el '
     &              //'integral blocks are', ' stored on Disk'
               End If
            End If
            iWrOpt=iWrOpt_Save
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iWRopt.eq.1) Then
*
*------- Molecule format (Molcas 1.0)
*
         Lu_28=28
         Lu_28=isfreeunit(Lu_28)
         Call DaName_MF(Lu_28,'BASINT')
         iDisk=0
         lBuf=iiLoc(nUt)-idLoc(Buf)
         lBuf=(lBuf+nByte)/nByte
*
         Call Drv2El_DMET(Integral_WrOut,Zero)
*
         ip_Buf=ip_of_iWork(Buf)
         Call iDafile(Lu_28,1,iWork(ip_Buf),lBuf,iDisk)
         nUt=-1
         Call iDafile(Lu_28,1,iWork(ip_Buf),lBuf,iDisk)
         Write (6,*)
         Write (6,'(A)')' Integrals are written in MOLCAS1 format'
         Write (6,'(I10,A)') IntTot,' Integrals written on Disk'
*                                                                      *
************************************************************************
*                                                                      *
      Else
*
         Call WarningMessage(2,'Seward: Invalid value of iWRopt!')
         Call Abend()
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Dist) Then
         Write (6,*)
         Write (6,*)
     & ' Distribution of the Absolute Values of the Integrals'
         Write (6,*)
         Write (6,'(1x,10I8)')  (i,i=-20,-11)
         Write (6,'(1x,10I8)')  (NrInt(i),i=-20,-11)
         Write (6,*)
         Write (6,'(1x,10I8)')  (i,i=-10,-1)
         Write (6,'(1x,10I8)')  (NrInt(i),i=-10,-1)
         Write (6,*)
         Write (6,'(1x,10I8)')  (i,i=0,9)
         Write (6,'(1x,10I8)')  (NrInt(i),i=0,9)
         Write (6,*)
      End If
*
      Call mma_deallocate(MemHide)
*                                                                      *
************************************************************************
*                                                                      *
*     At the end of the calculation free all memory to check for
*     corruption of the memory.
*

 9999 Call ClsSew
      If (Allocated(AdCell)) Call mma_deallocate(AdCell)
      Call mma_deallocate(Coor_MPM)
      Call mma_deallocate(Chrg)
      Call mma_deallocate(Mass)
      Call mma_deallocate(Centr)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(4,TCpu2-TCpu1,TWall2-TWall1)
*                                                                      *
************************************************************************
*                                                                      *
*     Diagonal ERI check
*
      If (Cholesky .or. Do_RI) Then
         If (DiagCheck) Then
            write(6,*)' ==== Start Diagonal ERI check  ===='
            Call Cho_X_init(irc,ChFracMem)
            if (irc.ne.0) then
               Call WarningMessage(2,
     &                   ' Seward: Non-zero rc in Cho_X_init.')
               Call Abend()
            endif
            Call Cho_X_CheckDiag(irc,DiagErr)
            if (irc.ne.0) then
               Call WarningMessage(2,
     &                   ' Seward: Non-zero rc in Cho_X_CheckDiag.')
               Call Abend()
            endif
            Call Cho_X_Final(irc)
            if (irc.ne.0) then
               Call WarningMessage(2,
     &                   ' Seward: Non-zero rc in Cho_X_Final.')
               CALL Abend()
            endif
            write(6,*)
            write(6,*)' ====  End  Diagonal ERI check  ===='
         EndIf
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*     Automatic run of GuessOrb
*
*     If (Do_GuessOrb.and.Do_FckInt) Call GuessOrb(iReturn,.FALSE.)
*     If(IsGvMode.gt.0) then
*       Call DoGvMode(IsGvMode)
*     EndIf
*     If (.not.Prprt.and.Do_OneEl) Call Put_NucAttr()
*                                                                      *
************************************************************************
*                                                                      *
*
*     Epilogue
*
      Call qExit('Seward')
      If (nPrint(iRout).ge.6) Then
         Call qStat(' ')
         Call FastIO('STATUS')
      End If
*
*
      ireturn=_RC_ALL_IS_WELL_
      If (Test)  Then
         ireturn=_RC_EXIT_
      Else
         If (isGvMode.gt.0.or.lRP_Post)
     &       ireturn=_RC_INVOKED_OTHER_MODULE_
      End If
      Return
      End

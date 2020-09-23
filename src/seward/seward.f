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
* Copyright (C) 1989-1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      Subroutine Seward(ireturn)
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
*          Roland Lindh, Dept. of Theoretical Chemistry, University of *
*          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
*          screening, July 1991.                                       *
************************************************************************
      use Real_Spherical
      use Period
      use GeoList
      use MpmC
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: lIrrep
      use LundIO
      use Temporary_Parameters
      use Integral_parameters, only: iPack, iWROpt
      use DKH_Info, only: DKroll
      use Real_Info, only: PkAcc
      use RICD_Info, only: Do_RI, Cholesky, DiagCheck, LocalDF
      Implicit Real*8 (A-H,O-Z)
      External Integral_WrOut, Integral_WrOut2, Integral_RI_3
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "print.fh"
#include "gateway.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
      Integer iix(2), nChoV(8)
      Real*8 rrx(2)
      Logical PrPrt_Save, Exist, DoRys, lOPTO
      Real*8  DiagErr(4), Dummy(2)
C-SVC: identify runfile with a fingerprint
      Character cDNA*256
      Logical IsBorn, Do_OneEl
*                                                                      *
************************************************************************
*                                                                      *
C     Call Seward_Banner()
      lOPTO = .False.
      nByte = iiLoc(iix(2)) - iiLoc(iix(1))
      nByte_r = idloc(rrx(2))-idloc(rrx(1))
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
      Call Seward_Init()
      If (Run_Mode.eq.S_Mode) Then
*
*        S_Mode
*
         DoRys=.True.
         nDiff=0
         Call GetInf(DoRys,nDiff)
         Primitive_Pass=.True.
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        GS_Mode
*                                                                      *
************************************************************************
*                                                                      *
         Call Funi_Init()
         Call Basis_Info_Init()
         Call Center_Info_Init()
*
      End If ! Run_Mode.eq.S_Mode

*                                                                      *
************************************************************************
*                                                                      *
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
      Call RdCtl_Seward(LuSpool,lOPTO,Do_OneEl)
      If (Run_Mode.ne.S_Mode) Then
         Call Basis_Info_Dmp()
         Call Basis_Info_Free()
         Call Basis_Info_Get()
         Call Center_Info_Dmp()
         Call Center_Info_Free()
         Call Center_Info_Get()
      End If
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
      Call Input_Seward(lOPTO)
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
      If (.Not.Primitive_Pass) Call DrvN0()
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
         Call DmpInf()
         Call basis2run()
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
     &   Call Drv1El()
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
*
*     If ERIs/CD/RI already available, one may want not to redo it!
*
      If (Fake_ERIs) Then
         Call set_fake_ERIs()
         Go To 9999
      EndIf
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     TWO-ELECTRON INTEGRAL SECTION
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
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
            Call Drv2El(Integral_WrOut2,Zero)
*
            Call Sort1B
            Call Sort2
            Call Sort3(MaxDax)
*
            If (nPrint(iRout).ge. 6) Then
               Write (6,*)
               Write (6,'(A)')
     &           ' Integrals are written in MOLCAS2 format'
               If ( iPack.ne.0 ) Then
                  Write (6,'(A)')
     &              ' No packing of integrals has been applied'
               Else
                  Write (6,'(A,G11.4)') ' Packing accuracy =',
     &                                   PkAcc
                  Write (6,'(A,I10)')
     &             ' Highest disk address written',MaxDax
               End If
               Write (6,'(A,A)') ' Diagonal and subdiagonal, '
     &           //'symmetry allowed 2-el',
     &           ' integral blocks are stored on Disk'
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
         lBuf=idLoc(Buf%r_End)-idLoc(Buf%Buf(1))
         lBuf=(lBuf+nByte_r)/nByte_r - 1
*
         Call Drv2El(Integral_WrOut,Zero)
*
         Call dDafile(Lu_28,1,Buf%Buf,lBuf,iDisk)
         Buf%nUt=-1
         Call dDafile(Lu_28,1,Buf%Buf,lBuf,iDisk)
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
*     At the end of the calculation free all memory to check for
*     corruption of the memory.
*

 9999 Call ClsSew()
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
      If (Do_GuessOrb.and.Do_FckInt) Call GuessOrb(iReturn,.FALSE.)
      If (.not.Prprt.and.Do_OneEl) Call Put_NucAttr()
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
      If (Test)  Then
         ireturn=_RC_EXIT_EXPECTED_
      Else If (lRP_Post) Then
         ireturn=_RC_INVOKED_OTHER_MODULE_
      Else
         ireturn=_RC_ALL_IS_WELL_
      End If

      Return
      End

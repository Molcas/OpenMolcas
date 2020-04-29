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
* Copyright (C) 1989, Per Ake Malmqvist                                *
*               1989, Bjorn O. Roos                                    *
*               1991,1993, Markus P. Fuelscher                         *
*               1991,1993, Jeppe Olsen                                 *
*               1998, Roland Lindh                                     *
************************************************************************
      SUBROUTINE RASSCF(IRETURN)
************************************************************************
*                                                                      *
*           ######     #     #####   #####   #####  #######            *
*           #     #   # #   #     # #     # #     # #                  *
*           #     #  #   #  #       #       #       #                  *
*           ######  #     #  #####   #####  #       #####              *
*           #   #   #######       #       # #       #                  *
*           #    #  #     # #     # #     # #     # #                  *
*           #     # #     #  #####   #####   #####  #                  *
*                                                                      *
*                                                                      *
*                   A program for complete (CAS) and                   *
*                   restricted (RAS)SCF calculations                   *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and J. Olsen, P.Aa. Malmqvist and B.O. Roos       *
*     University of Lund, Sweden                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     MOLCAS version 1 by P.Aa. Malmqvist and B.O. Roos, 1989          *
*     MOLCAS version 2 by M.P. Fuelscher and J. Olsen  , 1991          *
*     MOLCAS version 3 by M.P. Fuelscher and J. Olsen  , 1993          *
*                                                                      *
*     Modified to process only unique symmetry blocks, R. Lindh,       *
*     March 1998.                                                      *
*                                                                      *
************************************************************************

#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_cfg
      use qcmaquis_interface_version
      use qcmaquis_interface_environment, only:
     &    finalize_dmrg, dump_dmrg_info
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use write_orbital_files, only : OrbFiles, putOrbFile

      use generic_CI, only: CI_solver_t
      use fciqmc, only: DoNECI, fciqmc_solver_t
      use CC_CI_mod, only: Do_CC_CI, CC_CI_solver_t
      use fcidump, only : make_fcidumps, transform, DumpOnly

      use orthonormalization, only : ON_scheme
      use print_RDMs_NECI_format, only: printRDMs_NECI

      Implicit Real*8 (A-H,O-Z)

#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.fh"
#include "input_ras.fh"
#include "rasscf.fh"
#include "rasrc.fh"
#include "general.fh"
#include "gas.fh"
#include "splitcas.fh"
#include "bk_approx.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='RASSCF  ')
#include "rctfld.fh"
#include "timers.fh"
#include "casvb.fh"
#include "rasscf_lucia.fh"
#include "lucia_ini.fh"
#include "csfbas.fh"
#include "gugx.fh"
#include "pamint.fh"
#include "davctl.fh"
#include "qnctl.fh"
#include "orthonormalize.fh"
#include "ciinfo.fh"
#ifdef _FDE_
#include "embpotdata.fh"
#endif
#include "raswfn.fh"

      Logical DSCF
      Logical lTemp, lOPTO
      Character*8 label
      Character*80 Line
      Character*1 CTHRE, CTHRSX, CTHRTE
      Logical DoQmat,DoActive, l_casdft
      Logical IfOpened
#ifdef _DMRG_
      Logical Do_ESPF
#endif

* --------- Cholesky stuff:
      Integer ALGO
      Logical DoCholesky
      Logical timings,DoLock,Deco
      Integer Nscreen
      COMMON /CHOTODO /DoActive,DoQmat,ipQmat
      COMMON /CHLCAS /DoCholesky,ALGO
      COMMON /CHOPAR/ ChFracMem
      COMMON /CHOTIME / timings
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen
* --------- End Cholesky stuff
      Logical Do_OFemb, KEonly, OFE_first
      COMMON  / OFembed_L / Do_OFemb,KEonly,OFE_first
      COMMON  / OFembed_I / ipFMaux, ip_NDSD, l_NDSD
      Character*8 EMILOOP
* --------- End Orbital-Free Embedding stuff
* --------- FCIDUMP stuff:
      real*8, allocatable :: orbital_E(:), folded_Fock(:)
* --------- End FCIDUMP stuff:
* --------- Procedure pointers for CI-solvers
        class(CI_solver_t), allocatable :: CI_solver
* --------- End Procedure pointers.

! actual_iter starts at 0, so iter 1A == 0, 1B == 1, 2 == 2, 3 == 3 and so on
      integer :: actual_iter

      Common /IDSXCI/ IDXCI(mxAct),IDXSX(mxAct)

      External Get_ProgName
      Character*100 ProgName, Get_ProgName
      Character*15 STLNE2
      External QEnter, QExit
      External RasScf_Init
      External Scan_Inp
      External Proc_Inp
#ifndef _DMRG_
      logical :: doDMRG = .false.
#else
      integer :: maxtrR
      real*8  :: maxtrW
#include "nevptp.fh"
#endif
      Dimension Dummy(1)

* Start the traceback utilities
*
      Call QENTER(ROUTINE)

* Set status line for monitor:
      Call StatusLine('RASSCF:',' Just started.')

* Set the return code(s)
      ITERM  = 0
      IRETURN=_RC_ALL_IS_WELL_

* Set some Cholesky stuff
      DoActive=.true.
      DoQmat=.false.
      lOPTO=.False.

* Set variable IfVB to check if this is a VB job.
      ProgName=Get_ProgName()
      IfVB=0
      If (ProgName(1:5).eq.'casvb') IfVB=2
* Default option switches and values, and initial data.
      EAV1=0.0d0
      Call RasScf_Init()
      Call Seward_Init()
* Open the one-olectron integral file:
       LuOne=77
       LuOne=isFreeUnit(LuOne)
       iRC=-1
       iOpt=0
       Call OpnOne(iRC,iOpt,'ONEINT',LuOne)
       If (iRC.ne.0) Then
        Write (6,*) 'Error when trying to open the one-electron'
        Write (6,*) 'integral file.'
        Call Quit(_RC_INTERNAL_ERROR_)
       End If
       Call StatusLine('RASSCF:',' Read-in ONEINT')
      If (IfVB.eq.2) go to 10

*
* Make a copy, upper-cased, left-adjusted, of the input between and including
* the '&RASSCF' and the 'End of input' markers, skipping all lines beginning
* with '*' or '!' or ' '  when left-adjusted, and replacing any rightmost
* substring beginning with '!' with blanks.
* That copy will be in file 'CleanInput', and its unit number is returned
* as LUInput in common (included file input_ras.fh) by the following call:
      Call cpinp(LUInput,iRc)
* If something wrong with input file:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       Call WarningMessage(2,'Input file is unusable.')
       Write(6,*)' RASSCF Error: Could not make a clean copy of'
       Write(6,*)' the input file. This is an unexpected bug.'
       IRETURN=_RC_INTERNAL_ERROR_
       GOTO 9990
      End If

* Scan the input file for keywords:
      Call Scan_Inp(iRc)
* If something wrong with input file:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       If (IPRLOC(1).GE.TERSE) Then
        Call WarningMessage(2,'Scanning input file failed.')
* Calling again, now with iRc indicating an error, will echo the keywords:
        Call Scan_Inp(iRc)
       End If
       IRETURN=_RC_INPUT_ERROR_
       GOTO 9990
      End If

* Local print level in this routine:
      IPRLEV=IPRLOC(1)
*
10    CONTINUE
* Open files
      Call OpnFls_RASSCF(DSCF,DoCholesky)

* Some preliminary input data:
      Call Rd1Int
      If ( .not.DSCF ) Call Rd2Int_RASSCF

* Printed program header:

* Process the input:
      Call StatusLine('RASSCF:',' Processing input')
      Call Proc_Inp(DSCF,Info,lOPTO,iRc)
* If something goes wrong in proc_inp:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       If (IPRLEV.ge.TERSE) Then
        Call WarningMessage(2,'Input processing failed.')
        Write(6,*)' RASSCF Error: Proc_Inp failed unexpectedly.'
        Write(6,*)' Check the output file for any previous messages'
        Write(6,*)' that can help explain the failure.'
        Write(6,*)' Here is a printing of the input file that'
        Write(6,*)' was processed:'
        Rewind(LUInput)
  15    Continue
        Read(LuInput,'(A80)',End=16,Err=16) Line
        Write(6,*) Line
        Go To 15
  16    Continue
       End If
       IRETURN=iRc
       GOTO 9990
      End If


* Local print level may have changed:
      IPRLEV=IPRLOC(1)


      Call InpPri(lOpto)

* Note that CI_solver subclasses can provide a final procedure
* (some people might call it destructor). Hence the deallocation and
* cleanup is automatically performed, when it goes out of scope.
      if (DoNECI) then
        allocate(fciqmc_solver_t :: CI_solver)
      else if (Do_CC_CI) then
        allocate(CC_CI_solver_t :: CI_solver)
      end if

      if (allocated(CI_solver)) call CI_solver%init()


*
* If this is not CASDFT make sure the DFT flag is unset
*
      If (KSDFT(1:3).eq.'SCF') Then
        Call Get_iScalar('System BitSwitch',iFlags)
        iFlags=iAnd(iFlags,Not(2**6))
        Call Put_iScalar('System BitSwitch',iFlags)
      End If

* If the ORBONLY option was chosen, then Proc_Inp just generated
*  orbitals from the JOBIPH file. Nothing more to do:
      IF(KeyORBO.or.(MAXIT.eq.0)) GOTO 9989
*                                                                 *
*******************************************************************
* Initialize global variable for mcpdft method                    *
       l_casdft = KSDFT(1:5).eq.'TLSDA'   .or.
     &            KSDFT(1:6).eq.'TLSDA5'  .or.
     &            KSDFT(1:5).eq.'TBLYP'   .or.
     &            KSDFT(1:6).eq.'TSSBSW'  .or.
     &            KSDFT(1:5).eq.'TSSBD'   .or.
     &            KSDFT(1:5).eq.'TS12G'   .or.
     &            KSDFT(1:4).eq.'TPBE'    .or.
     &            KSDFT(1:5).eq.'FTPBE'   .or.
     &            KSDFT(1:7).eq.'TREVPBE' .or.
     &            KSDFT(1:8).eq.'FTREVPBE'.or.
     &            KSDFT(1:6).eq.'FTLSDA'  .or.
     &            KSDFT(1:6).eq.'FTBLYP'
*******************************************************************
#ifdef _DMRG_
      if(l_casdft .and. doDMRG) domcpdftDMRG = .true.
#endif
*
* Allocate various matrices
*
      Call GetMem('FI','Allo','Real',LFI,NTOT1)
      Call GetMem('FA','Allo','Real',LFA,NTOT1)
      Call GetMem('D1I','Allo','Real',LD1I,NTOT2)
      Call GetMem('D1A','Allo','Real',LD1A,NTOT2)
      Call GetMem('D1tot','Allo','Real',LD1tot,NTOT1)
      Call GetMem('OCCN','Allo','Real',LOCCN,NTOT)
      Call GetMem('LCMO','Allo','Real',LCMO,NTOT2)
      Call GetMem('DIAF','Allo','Real',LDIAF,NTOT)
      lfi_cvb=lfi
      lfa_cvb=lfa
      ld1i_cvb=ld1i
      ld1a_cvb=ld1a
      ld1tot_cvb=ld1tot
      loccn_cvb=loccn
      lcmo_cvb=lcmo
      ldiaf_cvb=ldiaf
      Call FZero(Work(LFA),NTOT1)
      Call FZero(Work(LDIAF),NTOT)
*
      LTUVX=1
      LDMAT=1
      LDSPN=1
      LPMAT=1
      LPA  =1

      If (iCIRST.eq.1.and.DumpOnly) then
        write(6,*) 'ICIRST and DumpOnly flags are not compatible!'
        write(6,*) 'Choose only one.'
        Call QTrace
        Call Abend
      end if

      If(DumpOnly) then
        write(6,*) 'Dumping integrals.'
        write(6,*) 'Nothing else will be done.'
      end if

        If ( NAC.GT.0 ) then
         Call GetMem('TUVX','Allo','Real',LTUVX,NACPR2)
         Call FZero(Work(LTUVX),NACPR2)
         ltuvx_cvb=ltuvx
         if(.not.DumpOnly) then
           Call GetMem('DMAT','Allo','Real',LDMAT,NACPAR)
           Call GetMem('DSPN','Allo','Real',LDSPN,NACPAR)
           Call GetMem('PMAT','Allo','Real',LPMAT,NACPR2)
           Call GetMem('P2AS','Allo','Real',LPA,NACPR2)
           call dcopy_(NACPAR,[0.0d0],0,Work(LDMAT),1)
           call dcopy_(NACPAR,[0.0d0],0,Work(LDSPN),1)
           ldmat_cvb=ldmat
           ldspn_cvb=ldspn
           lpmat_cvb=lpmat
           lpa_cvb=lpa
         end if
        Else
         LTUVX = ip_Dummy
         ltuvx_cvb=ltuvx
         LDMAT = ip_Dummy
         LDSPN = ip_Dummy
         LPMAT = ip_Dummy
         LPA   = ip_Dummy
        End If
#ifdef _FDE_
      ! Embedding
      iDummyEmb=0
      Call Get_iScalar('embpot', iDummyEmb)
      if (iDummyEmb.eq.1) embPot=.true.
      if (embPot) then
       Call EmbPotRdRun
      end if
      if (embpot) then
       ! I have no idea why i need memory for x+4 entries
       ! and not just x...
       Call GetMem('Emb','ALLO','REAL',ipEmb,NTOT1+4)
       if (embPotInBasis) then
        ! If the potential is given in basis set representation it
        ! has not been calculated with a OneEl call and is just read
        ! from file here.
        iunit = isFreeUnit(1)
        call molcas_open(iunit, embPotPath)
        do iEmb=1, NTOT1
         read(iunit,*) Work(ipEmb+iEmb-1)
        end do
       else
        ! Read in the embedding potential one-electron integrals
        Label='embpot  '
        iRC=-1
        iOpt=0
        Call RdOne(iRC,iOpt,Label,1,Work(ipEmb),iSyLbl)
        If (iRC.ne.0) then
         Call WarningMessage(2,
     &                'Drv1El: Error reading ONEINT;'
     &              //'Label='//Label)
         Call Quit(_RC_IO_ERROR_READ_)
        End If
       end if
      end if
#endif

*
* Get start orbitals

* Initialize OCCN array, to prevent false alarms later from
* automated detection of using uninitialized variables:
      call dcopy_(NTot,[0.0D0],0,Work(lOCCN),1)

* PAM03: Note that removal of linear dependence may change the nr
* of secondary/deleted orbitals, affecting some of the global
* variables: NSSH(),NDEL(),NORB(),NTOT3, etc etc
      Call ReadVc(Work(LCMO),Work(lOCCN),
     & WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),ON_scheme)
      if (KeyORTH) then
! TODO(Oskar): Add fourth argument OCC
!   If the Occupation number is written properly as well.
        call putOrbFile(CMO=Work(lCMO : lCMO + nTot2 - 1),
     &                  orbital_E=work(lDIAF : lDIAF + nTot - 1),
     &                  iDoGAS=iDoGAS)
      end if
* Only now are such variables finally known.

      If ( IPRLEV.ge.DEBUG ) then
        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',Work(LDMAT),NAC)
        CALL TRIPRT('Averaged one-body spin density matrix DS, RASSCF',
     &              ' ',Work(LDSPN),NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',WORK(LPMAT),NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',WORK(LPA),NACPAR)
      END IF
*
* Allocate core space for dynamic storage of data
*
      CALL ALLOC
*
* Create job interphase on unit JOBIPH (FT15)
*
      if(ifvb.ne.2) then
        CALL CREIPH
        call cre_raswfn
      end if
      if(ifvb.eq.1)call casinfo2_cvb()

      Call Timing(Swatch,Swatch,Ebel_1,Swatch)

CGG03 Aug 03
      If(NAlter.gt.0) Call Alter_MO(Work(LCMO))

c At this point all is ready to potentially dump MO integrals... just do it if required.
      If(DumpOnly) goto 20
      if(ifvb.eq.2)goto 20

      if(dofcidump)then
        Write(LF,*)
        Write(LF,'(26X,A)')
     &  'Dumping integrals on file FCIDUMP - nothing else to be done'
        goto 20
      end if

************************************************************************
*
* Wave function section
*
************************************************************************

      Call StatusLine('RASSCF:',' Compute wave function.')
      If ( IPRLEV.GE.2 .AND..NOT.lOPTO) then
       Write(LF,*)
       Write(LF,'(6X,120A1)') ('*',i=1,120)
       Write(LF,'(6X,A,118X,A)') '*','*'
       Write(LF,'(6X,A,44X,A,45X,A)')
     &      '*','Wave function control section','*'
       Write(LF,'(6X,A,118X,A)') '*','*'
       Write(LF,'(6X,120A1)') ('*',i=1,120)
       Write(LF,*)
      End If

      If ( IPRLEV.GE.2 .AND..NOT.lOPTO) then
       IF (ICIONLY.eq.0) THEN
        Write(LF,*)
        if(doDMRG)then
          Write(LF,'(41X,A)')
     &         'DMRGSCF iterations: Energy and convergence statistics'
          Write(LF,'(41X,A)')
     &         '-----------------------------------------------------'
        else
          Write(LF,'(41X,A)')
     &         'RASSCF iterations: Energy and convergence statistics'
          Write(LF,'(41X,A)')
     &         '----------------------------------------------------'
        end if
        Write(LF,*)
       ELSE
        Write(LF,*)
        if(doDMRG)then
          Write(LF,'(41X,A)')
     &        'DMRGCI only, no orbital optimization will be done.'
          Write(LF,'(41X,A)')
     &        '--------------------------------------------------'
        else
          Write(LF,'(41X,A)')
     &        'CASCI only, no orbital optimization will be done.'
          Write(LF,'(41X,A)')
     &        '-------------------------------------------------'
        end if
        Write(LF,*)
       END IF
#ifdef _DMRG_
       if(doDMRG)then
          Write(LF,'(45x,a//,36x,a/,36x,a/,36x,a//,45x,a//,'//
     &               '36x,a/,36x,a/,36x,a//,36x,a/,36x,a,a/,36x,a//)')
     &        'Please cite for the QCMaquis-Molcas driver:',
     &        'Freitag L.; Keller S.; Knecht S.; Lindh R.; Ma Y.; ',
     &        'Stein C. J. and Reiher M., in preparation. (2018).',
     &        '------------------------------------------------'//
     &        '---------------',
     &        'Please cite for the QCMaquis DMRG software:',
     &        'S. Keller, M. Dolfi, M. Troyer, M. Reiher,',
     &        'J. Chem. Phys. 143, 244118 (2015)',
     &        '------------------------------------------------'//
     &        '---------------',
     &        trim(qcmaquis_interface_v),
     &        'git reference SHA        : ',
     &        trim(qcmaquis_interface_g(1:12)),
     &        '------------------------------------------------'//
     &        '---------------'
       end if
#endif
       IF(INOCALC.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' No calculation will be performed. Stopping in LUCIA'
       END IF
       IF(ISAVE_EXP.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' Information on the CI-vector will be written to ???.'
       END IF
       IF(IEXPAND.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' A shorter vector will be expanded in a longer'
       END IF
       If ( IPRLEV.LE.3) then
        If(DoSplitCAS) then
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI       SplitCAS       Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root      energy       change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        else if (DoBKAP) then
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI   RASSCF      CI    Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root  energy    energy  change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        else if (DoDMRG .and. ICIONLY == 0)then
         Write(LF,'(6X,A)')
     &         'Iter num  DMRG max tr DMRG  SX      DMRGSCF'//
     &         '       Energy    '//
     &         'max ROT   max BLB     max BLB  Level Ln srch  Step '//
     &         '  QN     CPU Time'
         Write(LF,'(6X,A)')
     &         '   sweeps/root weight/root iter     energy'//
     &         '        change    '//
     &         ' param    element      value   shift minimum  type '//
     &         'update   hh:mm:ss'
        else if (DoDMRG .and. ICIONLY /= 0)then

        else if( l_casdft ) then

        else
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI       RASSCF       Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root      energy       change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        end if
       End If
      End If
20    continue
*                                                                      *
************************************************************************
*                                                                      *
*     Start iterations
*                                                                      *
************************************************************************
*                                                                      *
      Rc_CI  = 0
      Rc_SX  = 0
      ECAS   = 0.0d0
      ROTMAX = 0.0d0
      ITER   = 0
      actual_iter = 0
      IFINAL = 0
      TMXTOT = 0.0D0
      Call GetMem('FOcc','ALLO','REAL',ipFocc,nTot1)
      ipfocc_cvb=ipfocc
*                                                                      *
************************************************************************
*                                                                      *
*     Entry point for second and successive iterations
*                                                                      *
************************************************************************
*                                                                      *
 1000 CONTINUE

      if( l_casdft ) then
        KSDFT_TEMP=KSDFT
        KSDFT='SCF'
        ExFac=1.0D0
      else
         KSDFT_TEMP=KSDFT
         ExFac=Get_ExFac(KSDFT)
      end if

      ITER=ITER+1
      Write(STLNE2,'(A12,I3)')' Iteration ',ITER
      Call StatusLine('RASSCF:',STLNE2)
      Call Timing(Swatch,Swatch,Certina_1,Swatch)
#ifdef _DMRG_
      ! Leon 27/11/2017: Skip the first CI iteration if we're using
      ! DMRGCI and CIOnly.It's enabled only for DMRGCI with QCMaquis
      ! now, (to exclude potential side effects)
      ! but consider extending it to other cases!
      call DecideOnESPF(Do_ESPF)
      !write(LF,*) ' |rasscf> DecideOnESPF == ',Do_ESPF
      If (( ITER.EQ.1 ).and.((.not.(DoDMRG.and.(ICIONLY.NE.0))).or.lRf
     &    .or.domcpdftDMRG.or.Do_ESPF))THEN
#else
      If ( ITER.EQ.1 ) THEN
#endif
************************************************************************
*     ^   First iteration
************************************************************************
*

* Print header to file containing informations on CI iterations.
*
        Write(IterFile,'(20A4)') ('****',i=1,20)
        Write(IterFile,'(15X,A)') 'RASSCF iteration: 1A'
*
        Start_Vectors=.True.
        lTemp = lRf
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
        Call Timing(Swatch,Swatch,Fortis_1,Swatch)

        If (.not.DoCholesky .or. ALGO.eq.1) Then
           Call GetMem('PUVX','Allo','Real',LPUVX,NFINT)
           Call FZero(Work(LPUVX),NFINT)
        EndIf

        Call Get_D1I_RASSCF(Work(LCMO),Work(lD1I))
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1I in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
             call wrtmat(Work(lD1I+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        End If

* Compute D1A from CMO coefficients and, if CIREstart, old DMAT.
        If (iCIRST.eq.1) Then

           Call GetMem('TmpDMAT','Allo','Real',ipTmpDMAT,NACPAR)
           call dcopy_(NACPAR,Work(LDMAT),1,Work(ipTmpDMAT),1)
           If (NASH(1).ne.NAC) Call DBLOCK(Work(ipTmpDMAT))
           Call Get_D1A_RASSCF(Work(LCMO),Work(ipTmpDMAT),WORK(LD1A))
           Call GetMem('TmpDMAT','Free','Real',ipTmpDMAT,NACPAR)

           DoActive = .true.

        Else

           lRf = .false.
           IF( .not.l_casdft )  then
             KSDFT='SCF'
             ExFac=1.0D0
           end IF
           Call dcopy_(NTOT2,[0.0D0],0,WORK(LD1A),1)

           DoActive = .false.

        End If

        DoQmat=.false.

        IPR=0
        IF(IPRLOC(2).EQ.4) IPR=5
        IF(IPRLOC(2).EQ.5) IPR=10

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' PUVX in rasscf bf first TRACTL2'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          call wrtmat(Work(LPUVX),1,nFint, 1, nFint)

          Write(LF,*)
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF bf TRACTL2 1'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        end if

*
* Transform two-electron integrals and compute the Fock matrices FI and FA
* FI and FA are output from TRACTL2...
        CALL TRACTL2(WORK(LCMO),WORK(LPUVX),WORK(LTUVX),WORK(LD1I),
     &               WORK(LFI),WORK(LD1A),WORK(LFA),IPR,lSquare,ExFac)

c         Write(6,*) ' TUVX after TRACTL2'
c         write(6,*) (WORK(LTUVX+ind),ind=0,NACPR2-1)
        IF (ITER.eq.1 .and. IfCRPR) Then
* Core shift applied to projection of WF with doubly occupied core
          Call MkCRVEC(Work(LCMO),Work(LCRVEC))
        END IF

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF af TRACTL2 1'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do

          Write(LF,*)
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          Write(LF,*) ' PUVX in rasscf af first TRACTL2'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          call wrtmat(Work(LPUVX),1,nFint, 1, nFint)

          Write(LF,*)
          Write(LF,*) ' ---------------------'

          Write(6,*) ' TUVX after TRACTL2'
          write(6,*) (WORK(LTUVX+ind),ind=0,NACPR2-1)
          Write(LF,*)
          Write(LF,*) ' ---------------------'
        end if

        If (.not.DoCholesky .or. ALGO.eq.1) Then
          Call GetMem('PUVX','Free','Real',LPUVX,NFINT)
        EndIf

        Call Timing(Swatch,Swatch,Fortis_2,Swatch)
        Fortis_2 = Fortis_2 - Fortis_1
        Fortis_3 = Fortis_3 + Fortis_2

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF bf first call to CICTL'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do
        End If

*
* Compute initial CI vectors and density matrices
*
        Call Timing(Swatch,Swatch,Zenith_1,Swatch)

        if (DumpOnly) then
          call mma_allocate(orbital_E, nTot)
          call mma_allocate(folded_Fock, nAcPar)
          call transform(iter,
     &                   CMO=work(LCMO : LCMO + nTot2 - 1),
     &                   DIAF=work(LDIAF : LDiaf + nTot - 1),
     &                   D1I_AO=work(lD1I : lD1I + nTot2 - 1),
     &                   D1A_AO=work(lD1A : lD1A + nTot2 - 1),
     &                   D1S_MO=work(lDSPN : lDSPN + nAcPar - 1),
     &                   F_IN=work(lFI : lFI + nTot1 - 1),
     &                   orbital_E=orbital_E,
     &                   folded_Fock=folded_Fock)
          call make_fcidumps('FCIDUMP', 'H5FCIDUMP',
     &      orbital_E, folded_Fock,
     &      TUVX=work(ltuvx : ltuvx + nAcPr2 - 1), core_energy=EMY)
          call mma_deallocate(orbital_E)
          call mma_deallocate(folded_Fock)
          write(6,*) "FCIDMP file generated. Here for serving you!"
          goto 2010
        end if

        if (allocated(CI_solver)) then
          call CI_solver%run(actual_iter=actual_iter,
     &                    CMO=work(LCMO : LCMO + nTot2 - 1),
     &                    DIAF=work(LDIAF : LDiaf + nTot - 1),
     &                    D1I_AO=work(lD1I : lD1I + nTot2 - 1),
     &                    D1A_AO=work(lD1A : lD1A + nTot2 - 1),
     &                    TUVX=work(ltuvx : ltuvx + nAcPr2 - 1),
     &                    F_IN=work(lFI : lFI + nTot1 - 1),
     &                    D1S_MO=work(lDSPN : lDSPN + nAcPar - 1),
     &                    DMAT=work(lDMAT : lDMAT + nAcPar - 1),
     &                    PSMAT=work(lpmat : lPMat + nAcpr2 - 1),
     &                    PAMAT=work(lpa : lpa + nAcPr2 - 1))

#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
        else If(DoBlockDMRG) then
          CALL DMRGCTL(WORK(LCMO),
     &                 WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &                 WORK(LFI),WORK(LD1I),WORK(LD1A),
     &                 WORK(LTUVX),IFINAL,0)
#endif
        else
          CALL CICTL(WORK(LCMO),
     &               WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &               WORK(LFI),WORK(LD1I),WORK(LD1A),
     &               WORK(LTUVX),IFINAL)

          if(dofcidump)then
           write(LF,*) " FCIDUMP file generated. This is the end..."
           goto 9990
          end if
#ifdef _FDE_
          !Thomas Dresselhaus
          if (embpot) then
            !Eemb=DDot_(NACPAR,Work(ipEmb),1,Work(LDMAT),1)
!           Eemb=embPotEne(Work(LD1I), Work(LD1A), Work(ipEmb),
!    &                   Work(LCMO), nBasFunc, nFrozenOrbs, .true.)
            Eemb=embPotEneMODensities(Work(LD1I), Work(LD1A),
     &            Work(ipEmb), nBas, nTot2, nFrozenOrbs, nSym)
            Write(LF,*) "Energy from embedding potential with the"
            Write(LF,*) "initial CI vectors: ", Eemb
          end if
        !!!!!!!!!!!!!!!!!!!
#endif
* PAM 2015: Additional output line.
          If ( IPRLEV.ge.USUAL .and. .not. doDMRG) Then
             write(6,'(a,i4)')' Nr of preliminary CI iterations:',ITERCI
          End If
        end if

c.. dongxia testing jobiph
c.. upt to here, jobiph are all zeros at iadr15(2)
* If CASVB job, go directly to return.
        invec_cvb=invec
        If (DSCF) NewFock=1
        If (IfVB.eq.2) GoTo 9990

        EAV=0.0d0
        If(DoSplitCAS) then
          EAV = ENER(lRootSplit,ITER)
        Else
          DO KROOT=1,NROOTS
            EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
          END DO
        End If

        Call Get_D1A_RASSCF(WORK(LCMO),WORK(LDMAT),WORK(LD1A))

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF af Get_D1A_RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        end if

        Call Timing(Swatch,Swatch,Zenith_2,Swatch)
        Zenith_2 = Zenith_2 - Zenith_1
        Zenith_3 = Zenith_3 + Zenith_2
        lRf = lTemp

        IF( .not.l_casdft ) then
          KSDFT=KSDFT_TEMP
          ExFac=Get_ExFac(KSDFT)
        end IF

*     v GLM for MC-PDFT
      If(KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM'.or.l_casdft) Then
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF bf call NATORB_RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do
        End If

        Call GetMem('CMON','Allo','Real',LCMON,NTOT2)
        Call GetMem('OCCX','Allo','Real',LOCCX,NTOT)
        noscr1=max(nacpar,no2m)
        Call GetMem('NOscr1','Allo','Real',lNOscr1,noscr1)
        Call GetMem('NOscr2','Allo','Real',lNOscr2,NO2M)
        Call GetMem('SMAT','Allo','Real',LSMAT,NTOT1)
        CALL NATORB_RASSCF(Work(LCMO),
     &            Work(lNOscr1),Work(lNOscr2),Work(LSMAT),
     &            Work(LCMON),Work(LOCCX))
        Call dCopy_(NTOT2,WORK(LCMON),1,WORK(LCMO),1)
        Call Put_CMO(WORK(LCMO),ntot2)
        Call GetMem('NOscr1','Free','Real',lNOscr1,NOscr1)
        Call GetMem('NOscr2','Free','Real',lNOscr2,NO2M)
        Call GetMem('SMAT','Free','Real',LSMAT,NTOT1)
        Call GetMem('OCCX','Free','Real',LOCCX,NTOT)
        Call GetMem('CMON','Free','Real',LCMON,NTOT2)
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF af call NATORB_RASSCF & bf 2 CICTL'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
           iBas = nBas(iSym)
           if(iBas.ne.0) then
             write(6,*) 'Sym =', iSym
             do i= 1,iBas
               write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
             end do
             iOff = iOff + (iBas*iBas)
           end if
         End Do
        End If
      End If
*     ^ GLM End If for MC-PDFT

        EAV=0.0d0
        DO KROOT=1,NROOTS
         EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
        END DO
        If ( IPRLEV.ge.DEBUG ) then
          IF( l_casdft ) then
            write(6,*) 'EAV value in RASSCF after first call to CICTL:'
            write(6,*) EAV
          END if
        End if
      END IF
************************************************************************
*     ^ End First iteration
************************************************************************
*
* Print header to file containing informations on CI iterations.
*
      actual_iter = actual_iter + 1
      Write(IterFile,*)
      Write(IterFile,'(20A4)') ('****',i=1,20)
      IF (Iter .Eq. 1) Then
         Write(IterFile,'(15X,A)') 'RASSCF iteration: 1B'
      Else
         Write(IterFile,'(15X,A,I3)') 'RASSCF iteration: ',Iter
      End If
*
      If ( IPRLEV.ge.DEBUG .and. l_casdft) then
        write(6,*) ('*',i=1,70)
        write(6,*) 'we are done withe first standard CAS-CI iteration  '
        write(6,*) 'CI coeffs are known and mantained fix in next stage'
        write(6,*) 'We are now going to remove exchange from FI and FA '
        write(6,*) 'in TRACTL2 --> TRA_CTL2 --> TRADRV --> FTWO        '
        write(6,*) 'the ExFac is going to be set to 0.0d0 as DT asked! '
        write(6,*) 'FI and FA are going to change... '
        write(6,*) 'Check with previous printout to see differences.\  '
        write(6,*) ('*',i=1,70)
      End if


      IF( l_casdft ) then
        KSDFT=KSDFT_TEMP
        ExFac=0.0d0
      end IF


      IF(ICIONLY.NE.0) IFINAL=1
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
      Call Timing(Swatch,Swatch,Fortis_1,Swatch)
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call GetMem('PUVX','Allo','Real',LPUVX,NFINT)
         Call FZero(Work(LPUVX),NFINT)
      EndIf
      Call Get_D1I_RASSCF(Work(LCMO),Work(lD1I))

      DoActive = .true.

      If (DoCholesky.and.ALGO.eq.2) Then
         DoQmat=.true. ! to be used in the subsequent SX-section
         NTav=0
         do iSym=1,nSym
            NTav = NTav + nBas(iSym)*nAsh(iSym)
         end do
         Call GetMem('Q-mat','Allo','Real',ipQmat,NTav)
         Call Fzero(Work(ipQmat),NTav)
      EndIf

       IPR=0
       IF(IPRLOC(2).EQ.4) IPR=5
       IF(IPRLOC(2).EQ.5) IPR=10
       If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in RASSCF bf TRACTL2 2'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End Do
       end if
       CALL TRACTL2(WORK(LCMO),WORK(LPUVX),WORK(LTUVX),WORK(LD1I),
     &              WORK(LFI),WORK(LD1A),WORK(LFA),IPR,lSquare,ExFac)

      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in RASSCF af TRACTL2 2'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End Do
      end if

      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call GetMem('PUVX','Free','Real',LPUVX,NFINT)
      EndIf

      Call Timing(Swatch,Swatch,Fortis_2,Swatch)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2

*
* Compute the CI vectors and density matrices
*
      IF (.not. l_casdft) THEN !the following is skipped in CASDFT-GLM

        If(KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Then
          Call Put_CMO(WORK(LCMO),ntot2)
        End If

        Call Timing(Swatch,Swatch,Zenith_1,Swatch)
        if (allocated(CI_solver)) then
          call CI_solver%run(actual_iter=actual_iter,
     &                    CMO=work(LCMO : LCMO + nTot2 - 1),
     &                    DIAF=work(LDIAF : LDiaf + nTot - 1),
     &                    D1I_AO=work(lD1I : lD1I + nTot2 - 1),
     &                    D1A_AO=work(lD1A : lD1A + nTot2 - 1),
     &                    TUVX=work(ltuvx : ltuvx + nAcPr2 - 1),
     &                    F_IN=work(lFI : lFI + nTot1 - 1),
     &                    D1S_MO=work(lDSPN : lDSPN + nAcPar - 1),
     &                    DMAT=work(lDMAT : lDMAT + nAcPar - 1),
     &                    PSMAT=work(lpmat : lPMat + nAcpr2 - 1),
     &                    PAMAT=work(lpa : lpa + nAcPr2 - 1))
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
        else If(DoBlockDMRG) Then
            CALL DMRGCTL(WORK(LCMO),
     &             WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &             WORK(LFI),WORK(LD1I),WORK(LD1A),
     &             WORK(LTUVX),IFINAL,1)
#endif
        else
          CALL CICTL(WORK(LCMO),
     &               WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &               WORK(LFI),WORK(LD1I),WORK(LD1A),
     &               WORK(LTUVX),IFINAL)
        end if

c      call triprt('twxy',' ',WORK(LTUVX),nAc*(nAc+1)/2)
c      call triprt('P-mat 2',' ',WORK(LPMAT),nAc*(nAc+1)/2)

        EAV=0.0d0
        If (DoSplitCAS) Then
          EAV = ENER(lRootSplit,ITER)
        Else
         DO KROOT=1,NROOTS
           EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
         END DO
        End If

        IF (IPRLEV.ge.DEBUG) THEN
          write(6,*) 'EAV value in RASSCF after second call to CICTL:'
          write(6,*) EAV

          write(6,*) 'Printing matrices in RASSCF'
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do

          Write(LF,*)
          Write(LF,*) ' D1I in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(Work(lD1I+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
          write(6,*)
          write(6,*) 'Total Charge :', Tot_Charge

          Call GetMem('Fcore','Allo','Real',iTmp1,nTot1)
          iComp  =  1
          iSyLbl =  1
          iRc    = -1
          iOpt   =  6
          Call RdOne(iRc,iOpt,'OneHam',iComp,Work(iTmp1),iSyLbl)
          If ( iRc.ne.0 ) then
           Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
           Write(LF,*) 'Label = ',Label
           Write(LF,*) 'iRc = ',iRc
           Call QTrace
           Call Abend
          End if

          Write(LF,*)
          Write(LF,*) ' OneHam in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=0
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Work(iTmp1+iOff),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do

          Call GetMem('Fcore','Free','Real',iTmp1,nTot1)
          Call Get_dScalar('PotNuc',potNuc)

          write(6,*)
          write(6,*) 'PotNuc :', PotNuc

          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(Work(lD1A+ioff-1),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        End if
      else
        CALL GETMEM('CASDFT_Fock','ALLO','REAL',LFOCK,NACPAR)
* To fix the DS bug... I forgot to transform it to the AO basis... Agrrrrhhh!
        if(iSpin.eq.1) then
          If ( IPRLEV.ge.DEBUG ) then
            write(6,*) 'running a singlet. LDSPN set to zero!'
          end if
          Call dcopy_(NACPAR,[0.0d0],0,Work(LDSPN),1)
        end if
        CALL GETMEM('TmpDS_DFT' ,'Allo','REAL',ipTmpDS_DFT ,NACPAR)
        CALL GETMEM('TmpD1S_DFT','Allo','REAL',ipTmpD1S_DFT,NTOT2)
        Call dcopy_(NACPAR,Work(LDSPN),1,Work(ipTmpDS_DFT),1)
        IF ( NASH(1).NE.NAC ) then
          CALL DBLOCK(Work(ipTmpDS_DFT))
        end if
        Call Get_D1A_RASSCF(Work(LCMO),Work(ipTmpDS_DFT),
     &                      Work(ipTmpD1S_DFT))
        CALL GETMEM('TmpDS_DFT' ,'Free','REAL',ipTmpDS_DFT ,NACPAR)
        call CASDFT_terms(WORK(LCMO),WORK(LFOCK),WORK(LFI),WORK(LD1I),
     &                    WORK(LD1A),Work(ipTmpD1S_DFT))
        CALL GETMEM('TmpD1S_DFT','Free','REAL',ipTmpD1S_DFT,NTOT2)
        CALL GETMEM('CASDFT_Fock','FREE','REAL',LFOCK,NACPAR)
* to fix complains from garble option on borr machines... initialize THMAX to zero.
        THMAX = 0.0d0
      end if

c        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
c     &              ' ',Work(LDMAT),NAC)
c        CALL TRIPRT('Averaged one-body spin density matrix DS, RASSCF',
c     &              ' ',Work(LDSPN),NAC)
c        CALL TRIPRT('Averaged two-body density matrix, P',
c     &              ' ',WORK(LPMAT),NACPAR)
c        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
c     &              ' ',WORK(LPA),NACPAR)

      Call Get_D1A_RASSCF(WORK(LCMO),WORK(LDMAT),WORK(LD1A))
       If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' D1A in AO basis in RASSCF bf SXCTL'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(Work(LD1A+ioff-1),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
        End Do
       End If
      Call Timing(Swatch,Swatch,Zenith_2,Swatch)
      Zenith_2 = Zenith_2 - Zenith_1
      Zenith_3 = Zenith_3 + Zenith_2

*
c      Call rasscf_xml(Iter)
      Call rasscf_mcontrol(Iter)
*
* SX-section
*
      Call Timing(Swatch,Swatch,Gucci_1,Swatch)

       If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*) ' In RASSCF bf SXCTL'
        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',Work(LDMAT),NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',WORK(LPMAT),NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',WORK(LPA),NACPAR)
       end if
      CALL SXCTL(WORK(LCMO),WORK(LOCCN),
     &           WORK(LDMAT),WORK(LPMAT),WORK(LPA),
     &           WORK(LFI),WORK(LFA),WORK(LD1A),THMAX,IFINAL)
      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*)
       Write(LF,*) ' FI+FA in RASSCF after SXCTL'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=0
       Do iSym = 1,nSym
        iBas = nBas(iSym)
        Call TriPrt(' ',' ',Work(LFA+iOff),iBas)
        iOff = iOff + (iBas*iBas+iBas)/2
       End Do
      End If

cGLM   write(6,*) 'ECAS in RASSCF after call to SXCTL', ECAS
      If (DoCholesky.and.ALGO.eq.2) Then
         Call GetMem('Q-mat','Free','Real',ipQmat,NTav)
      EndIf

      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' D1A in AO basis in RASSCF af SXCTL'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(Work(LD1A+ioff-1),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
        End Do
      end if
      Call Get_D1A_RASSCF(WORK(LCMO),WORK(LDMAT),WORK(LD1A))
      CASDFT_Funct=0.0d0
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') then
        Call Get_dScalar('CASDFT energy',CASDFT_Funct)
cGLM        write(6,*) 'CASDFT energy :', CASDFT_Funct

      end IF
      DE=(ECAS+CASDFT_Funct)-ECAS1
      ECAS1=ECAS+CASDFT_Funct

      IF(NAC.EQ.0)   EAV=ECAS

      If (KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM') Then
        IF (nConf.EQ.1) Then
          EAV = ECAS
        Else
          EAV = EAV - VIA_DFT - HALFQ
        End If
      End IF

      IF (ITER.EQ.1) DE=0.0d0
      Call Timing(Swatch,Swatch,Gucci_2,Swatch)
      Gucci_2 = Gucci_2 - Gucci_1
      Gucci_3 = Gucci_3 + Gucci_2
      TMXTOT=MAX(TMXTOT,THMAX)
*
* Save energies and convergence parameters
*
      CONV(1,ITER)=ECAS
      CONV(2,ITER)=ESX
      CONV(3,ITER)=CMAX
      CONV(4,ITER)=DE
      CONV(5,ITER)=CBLBM
      CONV(6,ITER)=ROTMAX
      IAD15 = IADR15(6)
      CALL DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)
      CALL DDAFILE(JOBIPH,1,CONV,6*mxIter,IAD15)
#ifdef _HDF5_
      call mh5_put_attr (wfn_iter, Iter)
      call mh5_put_dset_array_real(wfn_energy, ENER(1,Iter))
#endif
*
* Print output of energies and convergence parameters
*
      Call Timing(Swatch,Swatch,Certina_2,Swatch)
      Certina_3 = Certina_2 - Certina_1
* Character indicating unconvergence/convergence criterion fulfilled:
      CTHRE=' '
      CTHRSX=' '
      CTHRTE=' '
      IF(ABS(DE).GT.THRE) CTHRE='*'
      IF(ABS(CBLBM).GT.THRSX) CTHRSX='*'
      IF(ABS(ROTMAX).GT.THRTE) CTHRTE='*'
      IROT=0
      IF(NROOTS.EQ.1) IROT=IROOT(1)
      IF ( IPRLEV.GE.2 .and. IPRLEV.LE.3 .AND..NOT.lOPTO) THEN
*----------------------------------
* Shift total energies (BOR 070411)
        if(iter.eq.1 .and. ICIONLY == 0)then
         kau=int(ECAS/1000.d0)
         EVAC=1000.d0*DBLE(kau)
         if(kau.ne.0) THEN
          write(6,'(6x,A,f23.2,A)')
     &       'Total energies have been shifted. Add ', EVAC,' au'
         end if
        endif
*----------------------------------
        ihh=int(Certina_3/3600)
        imm=int(Certina_3-ihh*3600)/60
        iss=int(Certina_3-ihh*3600-imm*60)
        if (DoSplitCAS) then
         Write(LF,'(6X,I3,I4,I5,I5,F15.8,E12.2,A1,E10.2,A1,2I4,I2,'//
     &            'E10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)')
     &        ITER,iterSplit,
     &        ITERSX,IROT,EAV,
     &        DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
        else if(DoBKAP) then
      Write(LF,'(3X,I3,I4,I2,I2,F15.8,F15.8,E12.2,A1,E10.2,A1,2I4,I2,'//
     &            'E10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)')
     &        ITER,ITERCI,
     &        ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,EAV,DE,CTHRE,
     &        ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss

        else
          if(doDMRG.and.KeyCION)then ! If DMRG only -- yma
            Write(LF,'(/6X,a,F19.10)')
     &      'DMRGCI energy              :',EAV+CASDFT_Funct
            Write(LF,'(6X,a,I5,A1,I2.2,A1,I2.2/)')
     &      'Total time spent (hh:mm:ss):        ',ihh,':',imm,':',iss
          else
            if(doDMRG)then
#ifdef _DMRG_
              maxtrW = 0.0d0
              maxtrR = -1
              ITERCI = MAXVAL(dmrg_energy%num_sweeps)
              IROT   = MAXLOC(dmrg_energy%num_sweeps,nroots)
              maxtrW = MAXVAL(dmrg_energy%max_truncW)
              maxtrR = MAXLOC(dmrg_energy%max_truncW,nroots)
#endif
         Write(LF,'(6X,I3,I3,I4,E12.2,I4,I5,F15.8,E12.2,A1,E9.2,A1,'//
     &   '2I4,I2,E10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I7,A1,I2.2,A1,I2.2)')
     &        ITER,ITERCI,IROT,maxtrW,maxtrR,
     &        ITERSX,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,
     &        ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
            else
            Write(LF,'(6X,I3,I4,I5,I5,F15.8,E12.2,A1,E10.2,A1,2I4,I2,'//
     &          'E10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)')
     &          ITER,ITERCI,
     &          ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,
     &          ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &          SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
            end if
          end if
        end if
      ELSE IF ( IPRLEV.GE.4) THEN
        Write(LF,'(6X,A,I4)') 'Energy statistics and convergence'//
     &                                          ' in iteration',ITER
        Write(LF,'(6X,A,10X,F26.6)') 'Average of CI energies',EAV-EVAC
        Write(LF,'(6X,A,F26.6)') 'Complete active space SCF energy',
     &                                                       ECAS-EVAC
        Write(LF,'(6X,A,17X,F26.6)') 'Super-CI energy',ESX
        Write(LF,'(6X,A,12X,F26.6,3X,A1)') 'RASSCF energy change',
     &                                                        DE,CTHRE
        Write(LF,'(6X,A,I1,A1,2I3,F19.6)') 'Maximum BLB matrix'
     &                  //' element(sym=',ISYMBB,')',IBLBM,JBLBM,CBLBM
        Write(LF,'(6X,A,10X,F26.6,3X,A1)') 'Max rotation parameter',
     &                                               ROTMAX,CTHRSX
        Write(LF,'(6X,A,14X,F26.6,3X,A1)') 'Max rotation angle',
     &                                               THMAX,CTHRTE
        Write(LF,'(6X,A,3X,F26.6)') 'Max change in MO coefficients',CMAX
      END IF
#ifdef _FDE_
      ! Embedding
      if (embpot) then
       Eemb=embPotEneMODensities(Work(LD1I), Work(LD1A), Work(ipEmb),
     &                nBas, nTot2, nFrozenOrbs, nSym)
       Write(LF,*)"E from embedding potential (<Psi|v_emb|Psi>): ",Eemb
      end if
#endif
      Call XFlush(6)
cGLM some additional printout for MC-PDFT

      If( l_casdft ) then
             CASDFT_E = ECAS-EVAC+CASDFT_Funct
             call print_mcpdft(CASDFT_E)
      end if
*
* Compare RASSCF energy and average CI energy
*
        DIFFE = ABS((ECAS-EAV)/ECAS)
      if(.not.DoSplitCAS) then
        if(DoBKAP) then
         if(iter.eq.1) then
          IF (DIFFE.GT.1.D-10 .AND. NROOTS.EQ.1) THEN
            Write(LF,'(6X,120A1)') ('=',i=1,120)
            Call WarningMessage(2,'Rasscf and CI energies will differ.')
            Write(LF,*)'This is the price you pay by the diagonal '
     &   //'approximation over the BB block in the SplitCAS method.'
            Write(LF,*)'The RASSCF energy might also diverge!'
            Write(LF,'(80A1)') ('#',i=1,80)
          END IF
         end if
        else if( l_casdft ) then
          Write(LF,'(6X,80A)') ('=',i=1,80)
          Write(LF,'(10X,A)') 'This is a POST-SCF correction using a '
     & //'modified  Hamiltonian.'
          write(LF,'(10X,A)') 'The RASSCF energy has been corrected and'
     & //' it will differ from'
          write(LF,'(10X,A)') 'the preceding CI energy.'
          Write(LF,'(6X,80A)') ('=',i=1,80)
        else
          IF(doDMRG)then

#ifdef _DMRG_DEBUG_
            write(lf,*) "DMRG-SCF energy    ",ECAS
            write(lf,*) "DMRG sweeped energy",EAV
#endif

            if(KeyCION)then
            else
              IF (DIFFE.GT.1.D-6 .AND. NROOTS.EQ.1) THEN
                Write(LF,'(6X,120A1)') ('=',i=1,120)
              Call WarningMessage(2,'DMRGSCF and DMRG energies differ.')
                Write(LF,'(6X,A,I11)')    'iteration           ',ITER
                Write(LF,'(6X,A,F22.10)') 'DMRGSCF energy      ',ECAS
                Write(LF,'(6X,A,F22.10)') 'DMRG energy         ',EAV
                Write(LF,'(6X,A,F22.10)') 'relative difference ',DIFFE
                Write(LF,*)'About this difference:'
                Write(LF,*)'1) If possible, consider a larger M value'
               Write(LF,*)'2) Severe convergence problems. Maybe active'
                Write(LF,*)'   space is unsuitable for this system?'
                Write(LF,'(6X,120A1)') ('=',i=1,120)
                IF (DIFFE.GT.5.D-04 .AND. NROOTS.EQ.1) THEN
                  Write(LF,*)
                  Write(LF,*)"Warning : "
                  Write(LF,*)" Relative difference is near unacceptable"
                  Write(LF,*)" If possible, consider a larger M value"
                  Write(LF,*)
                end if
              end if
            end if
          else
            IF (DIFFE.GT.1.D-10 .AND. NROOTS.EQ.1) THEN
              Write(LF,'(6X,120A1)') ('=',i=1,120)
              Call WarningMessage(2,'Rasscf and CI energies differ.')
              Write(LF,'(6X,A,I11)')    'iteration           ',ITER
              Write(LF,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
              Write(LF,'(6X,A,F22.10)') 'CI energy           ',EAV
              Write(LF,'(6X,A,F22.10)') 'relative difference ',DIFFE
              Write(LF,*)'Severe convergence problems. Maybe the active'
              Write(LF,*)'   space is unsuitable for this system?'
              Write(LF,'(6X,120A1)') ('=',i=1,120)
              IF(DIFFE.GT.1.D-04.AND.NROOTS.EQ.1.AND. .not.l_casdft)THEN
                Write(LF,*)
                Write(LF,'(6X,A)') 'The program has to stop !!!'
                Write(LF,*)
                Write(LF,'(6X,120A1)') ('=',i=1,120)
                Write(LF,*)
                ITERM=99
                GOTO 2000
              END IF
            end if
          END IF
        end if
      else
*          Write(LF,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
        IF (DIFFE.GT.5.0D-03) THEN
          Write(LF,'(6X,120A1)') ('*',i=1,120)
          Write(LF,'(6X,A)') 'The Split-RASSCF and Split-CI '//
     &                       'energies differ !!!'
*          Write(LF,'(6X,A,I11)')    'iteration           ',ITER
          Write(LF,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
          Write(LF,'(6X,A,F22.10)') 'Split-CI energy        ',EAV
          Write(LF,'(6X,A,F22.10)') 'Relative difference    ', DIFFE
          Write(LF,*)'     Smaller is the difference more realiable is',
     &               ' the result.'
          Write(LF,*)'     To make the difference smaller try',
     &               ' to select a bigger AA block or use firstOrder ',
     &               ' keyword.'
          Write(LF,'(6X,120A1)') ('*',i=1,120)
        END IF
      end if

*
* Convergence check:
* check is done on largest BLB matrix
* element (CBLBM), on difference in
* average energy, DE and on maximum value of non-
* diagonal rotation matrix element.
*
************************************************************************
************************************************************************
* IF CIONLY calculation the convergence is skipped and goes to line 2000
************************************************************************
************************************************************************

      IF (IFINAL.EQ.1) GOTO 2000
      IF (DE.GT.1.0D0) THEN
        Call StatusLine('RASSCF:','No convergence.')
        Write(LF,*)
        Write(LF,'(6X,120A1)') ('=',i=1,120)
        Call WarningMessage(2,'Rasscf energy diverges.')
        Write(LF,'(6X,A,I11)')    'iteration           ',ITER
        Write(LF,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
        Write(LF,'(6X,A,F22.10)') 'energy difference   ',DE
        Write(LF,'(6X,120A1)') ('=',i=1,120)
        Write(LF,*)
        Write(LF,'(6X,A)') '!!! The program was forced to stop !!!'
        Write(LF,*)
        Write(LF,'(6X,120A1)') ('=',i=1,120)
        Write(LF,*)
        ITERM=99
        GOTO 2000
      ENDIF
      IF(ITER.LT.MAXIT) THEN
        Write(STLNE2,'(A12,I3)')' Iteration ',ITER
        Call StatusLine('RASSCF converged:',STLNE2)
        IF(ABS(DE).GT.THRE) GO TO 1000
        IF(ABS(CBLBM).GT.THRSX) GO TO 1000
        IF(ABS(ROTMAX).GT.THRTE) GO TO 1000
        IF(ITER.LE.3.AND.ICIONLY.EQ.0) GO TO 1000   ! 3->0 checking
        IF(IPRLEV.ge.TERSE) Write(LF,'(6X,A,I3,A)')
     &        'Convergence after',ITER,' iterations'
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF after convergence printout'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If
        IFINAL=1
        Call Add_Info('RASSCF_ITER',[DBLE(ITER)],1,8)
*    Call Add_Info('RASSCF_THMX',TMXTOT,1,5)
        GOTO 1000
      ELSE
        IF(IPRLEV.ge.TERSE) Write(LF,'(6X,A,I3,A)')
     &        'No convergence after',ITER,' iterations'
        Write(STLNE2,'(A12,I3)')' Iteration ',ITER
        Call StatusLine('RASSCF max iter:',STLNE2)
        IFINAL=1
        ITERM=16
        GOTO 1000
      ENDIF
*                                                                      *
************************************************************************
*                                                                      *
*     Compute Final CI vectors
*                                                                      *
************************************************************************
*                                                                      *
 2000 IFINAL=2
      ICICH=0
************************************************************************
******************           Closing up MC-PDFT      *******************
************************************************************************

c Clean-close as much as you can the CASDFT stuff...
      if( l_casdft ) goto 2010

** IPT2 = 1 for OUTO, CANOnical keyword...
      IF(IPT2.EQ.1) THEN
        IAD=IADR15(9)
        CALL DDAFILE(JOBIPH,2,WORK(LCMO),NTOT2,IAD)
      ELSE
        IAD=IADR15(2)
        CALL DDAFILE(JOBIPH,2,WORK(LCMO),NTOT2,IAD)
      ENDIF
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF after DDAFILE'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (Work(LCMO+ioff-1+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If
      IF (NROOTS.GT.1) THEN
       Call StatusLine('RASSCF:','Compute final CI vectors')
      ELSE
       Call StatusLine('RASSCF:','Compute final CI vector')
      END IF
*
* Transform two-electron integrals
*
      Call Timing(Swatch,Swatch,Fortis_1,Swatch)
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call GetMem('PVX1','Allo','Real',LPUVX,NFINT)
         Call FZero(Work(LPUVX),NFINT)
      EndIf

      Call Get_D1I_RASSCF(Work(LCMO),Work(lD1I))

       DoQmat=.false.

       IPR=0
       IF(IPRLOC(2).EQ.4) IPR=5
       IF(IPRLOC(2).EQ.5) IPR=10
       CALL TRACTL2(WORK(LCMO),WORK(LPUVX),WORK(LTUVX),WORK(LD1I),
     &              WORK(LFI),WORK(LD1A),WORK(LFA),IPR,lSquare,ExFac)
*
       If (.not.DoCholesky .or. ALGO.eq.1) Then
          Call GetMem('PVX1','Free','Real',LPUVX,NFINT)
       EndIf

      Call Timing(Swatch,Swatch,Fortis_2,Swatch)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2
*
* CI-section (to obtain final wave function)
* 1st and 2nd order density matrix in MO basis
* for the gradients.
*
*
* Print header to file containing informations on CI iterations.
*
       Write(IterFile,*)
       Write(IterFile,'(20A4)') ('****',i=1,20)
       Write(IterFile,'(15X,A)') 'RASSCF iteration: Final'
*
      Call Timing(Swatch,Swatch,Zenith_1,Swatch)

      if (allocated(CI_solver)) then
          call CI_solver%run(actual_iter=actual_iter,
     &                    CMO=work(LCMO : LCMO + nTot2 - 1),
     &                    DIAF=work(LDIAF : LDiaf + nTot - 1),
     &                    D1I_AO=work(lD1I : lD1I + nTot2 - 1),
     &                    D1A_AO=work(lD1A : lD1A + nTot2 - 1),
     &                    TUVX=work(ltuvx : ltuvx + nAcPr2 - 1),
     &                    F_IN=work(lFI : lFI + nTot1 - 1),
     &                    D1S_MO=work(lDSPN : lDSPN + nAcPar - 1),
     &                    DMAT=work(lDMAT : lDMAT + nAcPar - 1),
     &                    PSMAT=work(lpmat : lPMat + nAcpr2 - 1),
     &                    PAMAT=work(lpa : lpa + nAcPr2 - 1))
#if defined _ENABLE_BLOCK_DMRG_ || defined _ENABLE_CHEMPS2_DMRG_
      else If(DoBlockDMRG) Then
        CALL DMRGCTL(WORK(LCMO),
     &           WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &           WORK(LFI),WORK(LD1I),WORK(LD1A),
     &           WORK(LTUVX),IFINAL,1)
#endif
#ifdef _DMRG_
! Leon 27/11/2017: Skip the final CI iteration if we're using DMRGCI
! and CIOnly. It's enabled only for DMRGCI with QCMaquis now
! (to exclude potential side effects)
! but consider extending it to other cases!
      else if (doDMRG .and. (ICIONLY.NE.0)) then
        continue
#endif
      else
        CALL CICTL(WORK(LCMO),
     &           WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA),
     &           WORK(LFI),WORK(LD1I),WORK(LD1A),
     &           WORK(LTUVX),IFINAL)
      end if

      EAV=0.0d0
      If(DoSplitCAS) then
        EAV = ENER(lRootSplit,ITER)
      else
       DO KROOT=1,NROOTS
         EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
       END DO
      end if
      IF(NAC.EQ.0) EAV=ECAS
      IF(NCRVEC.gt.0) then
* Core shift has been used
        Call GetMem('CRVEC','Free','Real',LCRVEC,NCRVEC)
        Call GetMem('CRPROJ','Free','Real',LCRPROJ,NCRPROJ)
      END IF
      Call Timing(Swatch,Swatch,Zenith_2,Swatch)
      Zenith_2 = Zenith_2 - Zenith_1
      Zenith_3 = Zenith_3 + Zenith_2
      Call Timing(Swatch,Swatch,Ebel_2,Swatch)
*
* Calculation of natural orbitals. These orbitals are stored on
* JOBIPH in IADR15(12), followed by the occupation numbers.
* Usage: Only for one-electron properties
* Note: in an average calculation natural orbitals are obtained
* for all roots. Stored sequentially at IADR15(12) on JOBPIH.
*PAM2009 This is done in NATORB. Note that it places the resulting
* natural orbitals and occ nos into JOBIPH. Thus the arguments nr
* 5 and 6 are simply those for the last root that NATORB processed,
* and they should not be reused for anything after return from NATORB.
* NATORB args: Arg1 is current CMO coeffs, used in CI;
*  all the rest should be regarded as scratch.
*
      Call StatusLine('RASSCF:','Compute natural orbitals')
      IPR=0
      IF(IPRLOC(6).EQ.4) IPR=5
      IF(IPRLOC(6).EQ.5) IPR=10
      Call GetMem('CMON','Allo','Real',LCMON,NTOT2)
      Call GetMem('OCCX','Allo','Real',LOCCX,NTOT)
      Call GetMem('NOscr1','Allo','Real',lNOscr1,MAX(NACPAR,NO2M))
      Call GetMem('NOscr2','Allo','Real',lNOscr2,NO2M)
      Call GetMem('SMAT','Allo','Real',LSMAT,NTOT1)
*PAM2009 NATORB args: Arg1 is current CMO coeffs, used in CI;
*  all the rest should be regarded as scratch.
      CALL NATORB_RASSCF(Work(LCMO),
     &            Work(lNOscr1),Work(lNOscr2),Work(LSMAT),
     &            Work(LCMON),Work(LOCCX))
      Call GetMem('NOscr1','Free','Real',lNOscr1,MAX(NACPAR,NO2M))
      Call GetMem('NOscr2','Free','Real',lNOscr2,NO2M)
*PAM2009 Deallocate CMON, OCCX.
      Call GetMem('CMON','Free','Real',LCMON,NTOT2)
      Call GetMem('OCCX','Free','Real',LOCCX,NTOT)

*
* Compute transition density matrices
      If (KeyTDM) Then
#ifdef _HDF5_
         Call GetMem('TMP','ALLO','REAL',iTmp,NConf)
         Call GetMem('LVEC','ALLO','REAL',iVecL,NConf)
         Call GetMem('RVEC','ALLO','REAL',iVecR,NConf)
         Call GetMem('KCNF','ALLO','INTE',ivkcnf,NACTEL)
         Call GetMem('Dtmp','ALLO','REAL',LW6,NAC*NAC)
         Call GetMem('SDtmp','ALLO','REAL',LW7,NAC*NAC)
         jDisk=IADR15(4)
         Call DDafile(JOBIPH,2,Work(iTmp),nConf,jDisk)
         Do jRoot=2,lRoots
*           Read and reorder the left CI vector
            Call DDafile(JOBIPH,2,Work(iTmp),nConf,jDisk)
            Call Reord2(NAC,NACTEL,LSYM,1,
     &                  iWork(KICONF(1)),iWork(KCFTP),
     &                  Work(iTmp),Work(iVecL),iWork(ivkcnf))
            C_Pointer=iVecL
            kDisk=IADR15(4)
            Do kRoot=1,jRoot-1
*              Read and reorder the right CI vector
               Call DDafile(JOBIPH,2,Work(iTmp),nConf,kDisk)
               Call Reord2(NAC,NACTEL,LSYM,1,
     &                     iWork(KICONF(1)),iWork(KCFTP),
     &                     Work(iTmp),Work(iVecR),iWork(ivkcnf))
*              Compute TDM and store in h5 file
               Call Lucia_Util('Densi',iVecR,iDummy,Dummy)
               idx=(jRoot-2)*(jRoot-1)/2+kRoot
               Call mh5_put_dset_array_real(wfn_transdens, Work(LW6),
     &              [NAC,NAC,1], [0,0,idx-1])
               If (iSpin.gt.1)
     &         Call mh5_put_dset_array_real(wfn_transsdens, Work(LW7),
     &              [NAC,NAC,1], [0,0,idx-1])
            End Do
         End Do
         Call GetMem('TMP','FREE','REAL',iTmp,NConf)
         Call GetMem('LVEC','FREE','REAL',iVecL,NConf)
         Call GetMem('RVEC','FREE','REAL',iVecR,NConf)
         Call GetMem('KCNF','FREE','INTE',ivkcnf,NACTEL)
         Call GetMem('Dtmp','FREE','REAL',LW6,NAC*NAC)
         Call GetMem('SDtmp','FREE','REAL',LW7,NAC*NAC)
#else
         Call WarningMessage(1,'HDF5 support disabled, '//
     &                         'TDM keyword ignored.')
#endif
      End If
*
*****************************************************************
* Export all information relevant to geometry optimizations.
* Save also the reaction field operator.
      Call Timing(Swatch,Swatch,Oris_1,Swatch)
      If (iRlxRoot.eq.0) iRlxRoot=iRoot(1)
*
* Replace average occ Fock with occ Fock for state iRlxRoot
* and densities with the densities for state iRLXRoot
c      write(6,*) 'I am in RASSCF before call to PutRlx!'
      If ( ITERM.ne.99 ) Then
         Call GetMem('Dens','ALLO','REAL',ipDens,nTot1)
         Call PutRlx(Work(LDMAT),Work(LDSPN),WORK(LPMAT),
     &            Work(ipDens),Work(LCMO))
         Call Export1(IFINAL,WORK(LCMO),WORK(LDMAT),WORK(LPMAT),
     &             Work(ipDens),Work(ipFocc))
         Call GetMem('Dens','FREE','REAL',ipDens,nTot1)
      End If
      Call Timing(Swatch,Swatch,Oris_2,Swatch)
      Oris_2 = Oris_2 - Oris_1
*****************************************************************
*
      EMY=EMY+CASDFT_Funct
*                                                                      *
************************************************************************
*                                                                      *
*     Output section
*                                                                      *
************************************************************************
*                                                                      *

      Call StatusLine('RASSCF:','Printing results')
      IF (IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
        Write(LF,*)
        Write(LF,'(6X,120A1)') ('*',i=1,120)
        Write(LF,'(6X,A,118X,A)') '*','*'
        Write(LF,'(6X,A,52X,A,53X,A)')
     &        '*','Final results','*'
        Write(LF,'(6X,A,118X,A)') '*','*'
        Write(LF,'(6X,120A1)') ('*',i=1,120)
        Write(LF,*)
      END IF
#ifdef _FDE_
      ! Embedding
      if (embpot) then
       Eemb=embPotEneMODensities(Work(LD1I), Work(LD1A), Work(ipEmb),
     &                nBas, nTot2, nFrozenOrbs, nSym)
       Write(LF,*) "Final energy from embedding potential: ", Eemb
       ! Write out ESP on grid if requested
       if (embWriteEsp) then
        Call Get_iScalar('Unique atoms',nNuc)
        Call embPotOutputMODensities(nNuc,nSym,LD1I,LD1A,nBas,nTot2)
       end if
      end if
#endif

      If ( ITERM.ne.99 ) THEN
       If (.not.DoSplitCAS) then
        CALL OUTCTL(WORK(LCMO),WORK(LOCCN),WORK(LSMAT),lOPTO)
       else
        CALL OUTCTLSplit(WORK(LCMO),WORK(LOCCN),WORK(LSMAT),lOPTO)
       end if
      End If

      Call GetMem('SMAT','Free','Real',LSMAT,NTOT1)

*
* Write information for MOLDEN
*
c  i_root=0 gives state- and spin-averaged natural orbitals
c  i_root>0 gives natural spin orbitals for that root
      mroots=nroots
      if ((iSpin.eq.1).and.(nroots.eq.1)) mroots=0
      Do i_root=0,mroots
        Call Interf(i_root,FDIAG,1,0)
      End Do

* Create output orbital files:
      Call OrbFiles(JOBIPH,IPRLEV)
*
************************************************************************
************ Priniting final RDMs in NECI format       *****************
************************************************************************
      If ( IPRLEV.ge.DEBUG ) then
       Call printRDMs_NECI(Work(LDMAT),NAC,Work(LPMAT),Work(LPA),NACPAR)
      End If
************************************************************************
******************           Closing up RASSCF       *******************
************************************************************************

2010   continue

c deallocating TUVX memory...
      IF(NAC.GT.0) THEN
         Call GetMem('TUVX','Free','Real',LTUVX,NACPR2)
         if(.not.DumpOnly) Then
            Call GetMem('P2AS','Free','Real',LPA,NACPR2)
            Call GetMem('PMAT','Free','Real',LPMAT,NACPR2)
            Call GetMem('DMAT','Free','Real',LDMAT,NACPAR)
            Call GetMem('DSPN','Free','Real',LDSPN,NACPAR)
         endif
      END IF
*
* release SEWARD
      Call ClsSew
* ClsSew is needed for releasing memory used by integral_util, rys... which is allocated when MC-PDFT run is performed.

*---  Finalize Cholesky information if initialized
      if (DoCholesky)then
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
            Write(LF,*)'RASSCF: Cho_X_Final fails with return code ',irc
            Write(LF,*)' Try to recover. Calculation continues.'
         endif
         If (Do_OFemb) Then
            Call GetMem('FMaux','Free','Real',ipFMaux,nTot1)
            Call OFE_print(EAV)
         EndIf
      endif

!      do i=1,NTOT2
!        write(*,*) "A,I",work(LD1A+i-1),work(LD1I+i-1)
!      end do

*  Release  some memory allocations
      Call GetMem('DIAF','Free','Real',LDIAF,NTOT)
      Call GetMem('FOCC','FREE','REAL',ipFocc,idum)
      Call GetMem('FI','Free','Real',LFI,NTOT1)
      Call GetMem('FA','Free','Real',LFA,NTOT1)
      Call GetMem('D1I','Free','Real',LD1I,NTOT2)
      Call GetMem('D1A','Free','Real',LD1A,NTOT2)
      Call GetMem('D1tot','Free','Real',lD1tot,NTOT1)
      Call GetMem('OCCN','Free','Real',LOCCN,NTOT)
      Call GetMem('LCMO','Free','Real',LCMO,NTOT2)
      If (iClean.eq.1) Call Free_iWork(ipCleanMask)

*
* Skip Lucia stuff if NECI or BLOCK-DMRG is on
      If (.not. any([allocated(CI_solver), DumpOnly,
     &              doDMRG, doBlockDMRG])) then
        Call Lucia_Util('CLOSE',iDummy,iDummy,Dummy)
      else if (allocated(CI_solver)) then
        call CI_solver%cleanup()
        deallocate(CI_solver)
      end if


      Call StatusLine('RASSCF:','Finished.')
      If (IPRLEV.GE.2) Write(LF,*)
      if(ifvb.eq.1) call make_close_rvb
cvv call to grid is moved up, in order to call clssew safely..
c      If (iCIonly.eq.0) Then
c        Call Grid_driver(-1,'RASSCF','RASORB',iR)
c      End If

      Call Timing(Swatch,Swatch,Ebel_3,Swatch)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim
       Call FastIO('STATUS')
      END IF
      Call ClsFls_RASSCF()

*
      Rc_RAS = ITERM
      Rc_RAS = Max(RC_RAS,Rc_CI)
      Rc_RAS = Max(RC_RAS,Rc_SX)
      If (Rc_Ras.eq.0) then
         ireturn=_RC_ALL_IS_WELL_
      Else If (Rc_Ras.eq.16) then
         ireturn=_RC_NOT_CONVERGED_
      Else
         Call WarningMessage(2,'Something is wrong: Did CI fail?')
         ireturn=_RC_GENERAL_ERROR_
      End If
*
      If (Do_OFemb) Then
         Call GetEnvF('EMIL_InLoop',EMILOOP)
         If (EMILOOP.eq.' ') EMILOOP='0'
         If (EMILOOP(1:1).ne.'0') Then
            If (iReturn.ne._RC_ALL_IS_WELL_) Then
               Call WarningMessage(1,'RASSCF: non-zero return code.')
            EndIf
            iReturn=_RC_CONTINUE_LOOP_
            Call Check_FThaw(iReturn)
         EndIf
      EndIf

      if (.not. (iDoGas .or. doDMRG .or. doBlockDMRG
     &          .or. allocated(CI_solver) .or. DumpOnly)) then
        Call MKGUGA_FREE
      end if

!Leon: The velociraptor comes! xkcd.com/292/
 9989 Continue
* DMRG: Save results for other use
! ==========================================================
      if(doDMRG)then
#ifdef _DMRG_
        !Leon: Generate 4-RDM evaluation templates for NEVPT2

        if (DoNEVPT2Prep) then
          do i=1,NROOTS
            if (DoEvaluateRDM) then
              Write (6,'(a,i4)') 'Evaluating 4-RDM for state ', i
            else
        Write (6,'(a,i4)') 'Writing 4-RDM QCMaquis template'//
     &   ' for state ', i
            end if
            call dmrg_interface_ctl(task='w 4rdmin',state=i-1,
     &       msproj=MPSCompressM,rdm4=DoEvaluateRDM)
      !> MPSCompressM is passed to dmrg_interface_ctl via msproj,
      !> if it is > 0, then MPS compression is triggered there
          end do
          ! Generate 3-TDM templates
          do i=1,NROOTS
            do j=i+1,NROOTS
              if (DoEvaluateRDM) then
            Write (6,'(a,i4,i4)') 'Evaluating 3-TDM for states ', i,j
            else
        Write (6,'(a,i4,i4)') 'Writing 3-TDM QCMaquis'//
     &   ' template for state ', i,j
            end if

              call dmrg_interface_ctl(task='w 3tdmin',state=i-1,
     &         stateL=j-1,rdm3=DoEvaluateRDM)
            end do
          end do
        end if
!       call dump_dmrg_info()
        call finalize_dmrg()
#endif
      end if
! ==========================================================
! Exit

 9990 Continue
C Close the one-electron integral file:
      iRC=-1
      iOpt=0
      Call ClsOne(iRC,iOpt)
      If (iRC.ne.0) Then
         Write (6,*) 'Error when trying to close the one-electron'
         Write (6,*) 'integral file.'
         Call Quit(_RC_INTERNAL_ERROR_)
      End If

      If (IfVB.ne.2) Then
        DO I=10,99
          INQUIRE(UNIT=I,OPENED=IfOpened)
          IF (IfOpened.and.I.ne.19) CLOSE (I)
        END DO
        Close(LUInput)
      End If
      Call qExit(ROUTINE)

      return

      end subroutine rasscf

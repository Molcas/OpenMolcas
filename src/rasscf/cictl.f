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
* Copyright (C) Bjorn O. Roos                                          *
*               Per Ake Malmqvist                                      *
*               1991, Jeppe Olsen                                      *
*               1991,1996, Markus P. Fuelscher                         *
************************************************************************
*  CICtl
*
*> @brief
*>   CI Control
*> @author B. O. Roos
*> @author P. &Aring;. Malmqvist
*> @modified_by P. &Aring;. Malmqvist
*>
*> @details
*> Depends on \p IFINAL, which is set in ::RASSCF. If \p IFINAL = ``0``, repeated
*> calculations with orbital optimization before each call. If \p IFINAL = ``1``,
*> there has been no orbital optimization, or the calculation is
*> converged. \p IFINAL = ``2`` means this is a final CI calculation, using the
*> final orbitals. For meaning of global variables \c NTOT1, \c NTOT2, \c NACPAR
*> and \c NACPR2, see src/Include/general.fh and src/Include/rasscf.fh.
*>
*> @param[in]     CMO    MO coefficients
*> @param[out]    D      Average 1-dens matrix
*> @param[out]    DS     Average spin 1-dens matrix
*> @param[out]    P      Average symm. 2-dens matrix
*> @param[out]    PA     Average antisymm. 2-dens matrix
*> @param[out]    FI     Fock matrix from inactive density
*> @param[in,out] D1I    Inactive 1-dens matrix
*> @param[in,out] D1A    Active 1-dens matrix
*> @param[in]     TUVX   Active 2-el integrals
*> @param[in]     IFINAL Calculation status switch
************************************************************************
      Subroutine CICtl(CMO,D,DS,P,PA,FI,D1I,D1A,TUVX,IFINAL)
* ****************************************************************
* history:                                                       *
* updated to use determinant based CI-procedures                 *
* J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
* updated for MOLCAS version 3                                   *
* J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
* updated for integral direct and reaction field calculations    *
* M.P. Fuelscher, University of Lund, Sweden, 1996               *
* ****************************************************************
#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_main, only: file_name_generator
#endif

      Implicit Real* 8 (A-H,O-Z)

      Dimension CMO(*),D(*),DS(*),P(*),PA(*),FI(*),D1I(*),D1A(*),
     &          TUVX(*)
      Logical Exist,Do_ESPF,l_casdft
*JB   variables for state rotation on final states
      Logical do_rotate
*JB   end of variables for state rotation calculation

#include "rasdim.fh"
#include "rasscf.fh"
#include "splitcas.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='CICTL   ')
#include "csfbas.fh"
#include "gugx.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "rctfld.fh"
#include "timers.fh"
#include "casvb.fh"
#include "wadr.fh"
#include "rasscf_lucia.fh"
#include "pamint.fh"
#include "input_ras.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#  include "raswfn.fh"
      real*8, allocatable :: density_square(:,:)
#endif
      Common /IDSXCI/ IDXCI(mxAct),IDXSX(mxAct)

#ifndef _DMRG_
      logical :: doDMRG = .false.
#else
      character(len=2300) :: maquis_name_states
      character(len=2300) :: maquis_name_results
      logical             :: rfh5DMRG
      logical             :: twordm_qcm
#endif
      Dimension rdum(1)

*PAM05      SymProd(i,j)=1+iEor(i-1,j-1)
      Call qEnter('CICTL')
C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

* PAM 2017-05-23 Modify TUVX by adding a shift vector to TUVX, which has
* the effect of adding a scalar times a projector for doubly-occupied
* core states.
          IF(IfCRPR) Then
*      write(6,*)' CICTL calling MKPROJ.'
*      call xflush(6)
            CALL MKPROJ(Work(LCRVEC),CMO,TUVX)
*      write(6,*)' CICTL back from MKPROJ.'
*      call xflush(6)
          END IF


* set up flag 'IFCAS' for GAS option, which is set up in gugatcl originally.
* IFCAS = 0: This is a CAS calculation
* IFCAS = 1: This is a RAS calculation
*
      if(iDoGas.or.ifcas.gt.2) call setsxci
      If ( IPRLEV.gt.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Enter CI section, CICTL routine'
        Write(LF,*) ' ================'
        Write(LF,*)
        Write(LF,*) ' iteration count =',ITER
      End If
      if(ifinal.ne.0) PamGen1=.True.

!      do i=1,NTOT2  ! yma
!        write(*,*)"ifinal CMO",ifinal,i,CMO(i)
!      end do

*
* SOME DIRTY SETUPS
*
      S=0.5D0*DBLE(ISPIN-1)
*
* COMPUTE ONE ELECTRON INTEGRALS IN MO BASIS
* AND ADD CORE INTERACTION
*
* LW1: FOCK MATRIX IN MO-BASIS
* LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION
*
      CALL GETMEM('CICTL1','ALLO','REAL',LW1,NACPAR)
      IF (IPRLEV.GE.DEBUG) THEN
        Write(LF,*) ' WORK SPACE VARIABLES IN SUBR. CICTL: '
        Write(LF,*) ' SGFCIN ',LW1
      END IF
      Call DecideOnESPF(Do_ESPF)

*                                                                      *
************************************************************************
* Global variable for MCPDFT functionals                               *
      l_casdft = KSDFT(1:5).eq.'TLSDA'   .or.
     &           KSDFT(1:6).eq.'TLSDA5'  .or.
     &           KSDFT(1:5).eq.'TBLYP'   .or.
     &           KSDFT(1:6).eq.'TSSBSW'  .or.
     &           KSDFT(1:5).eq.'TSSBD'  .or.
     &           KSDFT(1:5).eq.'TS12G'  .or.
     &           KSDFT(1:4).eq.'TPBE'    .or.
     &           KSDFT(1:5).eq.'FTPBE'   .or.
     &           KSDFT(1:7).eq.'TREVPBE' .or.
     &           KSDFT(1:8).eq.'FTREVPBE'.or.
     &           KSDFT(1:6).eq.'FTLSDA'  .or.
     &           KSDFT(1:6).eq.'FTBLYP'
*                                                                      *
************************************************************************
*                                                                      *

      If ( lRf .or. KSDFT.ne.'SCF' .or. Do_ESPF) THEN
*
* In case of a reaction field in combination with an average CAS
* select the potential of the appropriate state.
*
        jDisk = IADR15(3)
        Do i=1,IPCMROOT-1
          Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
        End Do
*
        CALL GETMEM('D1A_FULL','ALLO','REAL',LRCT_F,NTOT2)
        CALL GETMEM('D1S_FULL','ALLO','REAL',LRCT_FS,NTOT2)
        If (IFinal.eq.0) Then
*
* Use normal MOs
*
           CALL GETMEM('D1A_RCT','ALLO','REAL',LRCT,NACPAR)
           CALL GETMEM('P2MO','ALLO','REAL',ipP2MO,NACPR2)
*
* Get the total density in MOs
*
           Call DDafile(JOBIPH,2,Work(LRCT),NACPAR,jDisk)
           Call Put_D1MO(Work(LRCT),NACPAR)  ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LRCT))
* Transform to AOs
           Call Get_D1A_RASSCF(CMO,WORK(LRCT),WORK(LRCT_F))
*
* Get the spin density in MOs
*
           IF (NACTEL.EQ.0) THEN
             CALL DCOPY_(NTOT2,[0.0D0],0,WORK(LRCT_FS),1)
           ELSE
             CALL GETMEM('D1S_RCT','ALLO','REAL',LRCT_S,NACPAR)
             Call DDafile(JOBIPH,2,Work(LRCT_S),NACPAR,jDisk)
             IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LRCT_S))
* Transform to AOs
             Call Get_D1A_RASSCF(CMO,WORK(LRCT_S),WORK(LRCT_FS))
             CALL GETMEM('D1S_RCT','FREE','REAL',LRCT_S,NACPAR)
           END IF
*
* Get the 2-particle density in MO
*
           Call DDafile(JOBIPH,2,Work(ipP2MO),NACPR2,jDisk)
           Call Put_P2MO(Work(ipP2MO),NACPR2) ! Put it on the RUNFILE
*
           CALL SGFCIN(CMO,WORK(LW1),FI,D1I,Work(LRCT_F),Work(LRCT_FS))
*
           CALL GETMEM('P2MO','FREE','REAL',ipP2MO,NACPR2)
           CALL GETMEM('D1A_RCT','FREE','REAL',LRCT,NACPAR)
*
        Else
*
* Here the pseudo-natural orbitals are in CMO and we need to
* get the D1A of the selected state in this basis.
*
*
* Compute the density of the particular state
*
           Call GetMem('CIVEC','ALLO','REAL',LW4,NCONF)
           If (NACTEL.EQ.0) THEN
             Work(LW4)=1.0D0
           Else
             if(.not.(doDMRG))then
!               write(*,*)"run the load back CI vector part" ! yma
               iDisk = IADR15(4)
               Do jRoot = 1,IPCMROOT
                   iOpt=0
                   If (jRoot.eq.IPCMROOT) iOpt=2
* load back one CI vector at the time
                   Call DDafile(JOBIPH,iOpt,Work(LW4),nConf,iDisk)
                  !call DVcPrt('BLUBB-start CI PCM',' ',Work(LW4),nConf)
               End Do
             end if
           End If

* compute density matrices

           Call GetMem('Dtmp ','ALLO','REAL',LW6,NACPAR)
           Call GetMem('DStmp','ALLO','REAL',LW7,NACPAR)
           Call GetMem('Ptmp ','ALLO','REAL',LW8,NACPR2)
           If ( NAC.ge.1 ) Then

              If (NACTEL.eq.0) THEN
                 call dcopy_(NACPAR,[0.0D0],0,WORK(LW6),1)
                 call dcopy_(NACPAR,[0.0D0],0,WORK(LW7),1)
                 call dcopy_(NACPR2,[0.0D0],0,WORK(LW8),1)
              Else

                if(doDMRG)then
#ifdef _DMRG_
                call dmrg_interface_ctl(
     &                                  task  = 'imp rdmX',
     &                                  x1    = work(lw6:lw6+NACPAR-1),
     &                                  x2    = work(lw8:lw8+NACPR2-1),
     &                                  ndim  = nacpar,
     &                                  mdim  = nacpr2,
     &                                  state = ipcmroot,
     &                                  rdm1  = .true.,
     &                                  rdm2  = (.not.KeyCION)
     &                                 )

                !> import 1p-spin density
                call dcopy_(NACPAR,[0.0D0],0,work(lw7),1)
                call dmrg_interface_ctl(
     &                                  task  = 'imp spdX',
     &                                  x1    = work(lw7:lw7+NACPAR-1),
     &                                  ndim  = nacpar,
     &                                  state = ipcmroot,
     &                                  rdm1  = .true.
     &                                 )
#endif
               else
                 Call GetMem('PAtmp','ALLO','REAL',LW9,NACPR2)
                 Call GetMem('Pscr','ALLO','REAL',LW10,NACPR2)
                 C_Pointer = Lw4
                 CALL Lucia_Util('Densi',0,iDummy,rdum)
                 If (IFCAS.GT.2 .OR. iDoGAS) Then
                   Call CISX(IDXSX,Work(LW6),Work(LW7),Work(LW8),
     &                     Work(LW9),Work(LW10))
                 End If
                 Call GetMem('Pscr','FREE','REAL',LW10,NACPR2)
                 Call GetMem('PAtmp','FREE','REAL',LW9,NACPR2)
               end if ! doDMRG/doBLOK or CI

             End If
           Else
              call dcopy_(NACPAR,[0.0D0],0,WORK(LW6),1)
              call dcopy_(NACPAR,[0.0D0],0,WORK(LW7),1)
              call dcopy_(NACPR2,[0.0D0],0,WORK(LW8),1)
           End If
* Modify the symmetric 2-particle density if only partial
* "exact exchange" is included.
c          n_Det=2
c          n_unpaired_elec=(iSpin-1)
c          n_paired_elec=nActEl-n_unpaired_elec
c          If(n_unpaired_elec+n_paired_elec/2.eq.nac) n_Det=1
*                  write(6,*) 'n_Det =', n_Det
           If (ExFac.ne.1.0D0.AND.(.not.l_casdft))
     &    Call Mod_P2(Work(LW8),NACPR2,
     &                                   Work(LW6),NACPAR,
     &                                   Work(LW7),ExFac,n_Det)
*
           Call Put_P2MO(Work(LW8),NACPR2) ! Put it on the RUNFILE
*
           Call GetMem('Ptmp ','FREE','REAL',LW8,NACPR2)
*
           Call Put_D1MO(Work(LW6),NACPAR) ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LW6))
           Call Get_D1A_RASSCF(CMO,Work(LW6),Work(LRCT_F))
*
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LW7))
           Call Get_D1A_RASSCF(CMO,Work(LW7),Work(LRCT_FS))
*
!           do i=1,NACPAR  !yma
!              write(*,*)"i-rdms1 1",i,Work(LW6+i-1)
!           end do

           Call GetMem('DStmp','FREE','REAL',LW7,NACPAR)
           Call GetMem('Dtmp ','FREE','REAL',LW6,NACPAR)
           Call GetMem('CIVEC','FREE','REAL',LW4,NCONF)
*
           Call SGFCIN(CMO,Work(LW1),FI,D1I,Work(LRCT_F),Work(LRCT_FS))
*
        End If
        CALL GETMEM('D1S_FULL','FREE','REAL',LRCT_FS,NTOT2)
        CALL GETMEM('D1A_RCT','FREE','REAL',LRCT_F,NTOT1)
*
      ELSE
*********************************************************************************************
* Normal case: no reaction field, no ESP to be added to the Fock matrix or ksdft /= SCF
*********************************************************************************************
         IF (IPRLEV.GE.DEBUG) THEN
           CALL TRIPRT('DS input  ',' ',DS,NAC)
           CALL TRIPRT('D input  ',' ',D,NAC)
         END IF
        CALL GETMEM('TmpDS' ,'Allo','REAL',ipTmpDS ,NACPAR)
        CALL GETMEM('TmpD1S','Allo','REAL',ipTmpD1S,NTOT2 )
        call dcopy_(NACPAR,DS,1,Work(ipTmpDS),1)
        IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(ipTmpDS))
        Call Get_D1A_RASSCF(CMO,Work(ipTmpDS),Work(ipTmpD1S))

        CALL GETMEM('TmpDS' ,'Free','REAL',ipTmpDS ,NACPAR)
*
        CALL SGFCIN(CMO,WORK(LW1),FI,D1I,D1A,Work(ipTmpD1S))
        CALL GETMEM('TmpD1S','Free','REAL',ipTmpD1S,NTOT2 )
*
      END IF

!     stknecht: the core fock matrix (pointer lw1)
!               total core energy (emy) as output of sgfcin
!     task: write the FCIDUMP file
      if(doDMRG)then
#ifdef _DMRG_
        call dmrg_interface_ctl(
     &                          task   = 'fci dump',
     &                          x1     = work(lw1),
     &                          x2     = tuvx,
     &                          energy = emy
     &                         )
        if(dofcidump)then
          CALL GETMEM('CICTL1','FREE','REAL',LW1,NACPAR)
          goto 9000
        end if
#endif
      end if
*
      lw1_cvb=lw1
      If (IfVB.eq.2) GoTo 9000

*
* C
* DAVIDSON DIAGONALIZATION
* C
*
C     kh0_pointer is used in Lucia to retrieve H0 from Molcas.
      kh0_pointer = lw1
      if(IfVB.eq.1)then
        call cvbmn_rvb(max(ifinal,1))
      else
        If (KSDFT(1:3).ne.'SCF'
     &      .and.DFTFOCK(1:4).eq.'DIFF'.and.nac.ne.0) Then
          nTmpPUVX=nFint
          Call GetMem('TmpPUVX','Allo','Real',ipTmpPUVX,nTmpPUVX)
          Call GetMem('TmpTUVX','Allo','Real',ipTmpTUVX,NACPR2)
          Call dCopy_(NACPR2,[0.0d0],0,Work(ipTmpTUVX),1)
          Call Get_dArray('DFT_TwoEl',Work(ipTmpPUVX),nTmpPUVX)
          Call Get_TUVX(Work(ipTmpPUVX),Work(ipTmpTUVX))
          Call DaXpY_(NACPR2,1.0d0,TUVX,1,Work(ipTmpTUVX),1)
         if (DoSplitCAS) then  ! (GLMJ)
           Call SplitCtl(Work(LW1),Work(ipTmpTUVX),IFINAL,iErrSplit)
           Call GetMem('TmpTUVX','Free','Real',ipTmpTUVX,NACPR2)
           Call GetMem('TmpPUVX','Free','Real',ipTmpPUVX,nTmpPUVX)
           if (iErrSplit.eq.1) then
            write(LF,*) ('*',i=1,120)
            write(LF,*)'WARNING!!!'
            write(LF,*) 'SplitCAS iterations don''t converge.'
            write(LF,*) 'The program will continue'
            write(LF,*) 'Hopefully your calculation will converge',
     &                 'next iteration!'
            write(LF,*) ('*',i=1,120)
           end if
           if (iErrSplit.eq.2) then
            write(LF,*) ('*',i=1,120)
            write(LF,*)'WARNING!!!'
            write(LF,*) 'SplitCAS iterations don''t converge.'
          write(LF,*)'Try to increase MxIterSplit or SplitCAS threshold'
            write(LF,*) 'The program will STOP!!!'
            write(LF,*) ('*',i=1,120)
            call xQuit(96)
           end if
         end if
         If (.not.DoSplitCAS) then
           Call   DavCtl(Work(LW1),Work(ipTmpTUVX),IFINAL)
           Call GetMem('TmpTUVX','Free','Real',ipTmpTUVX,NACPR2)
           Call GetMem('TmpPUVX','Free','Real',ipTmpPUVX,nTmpPUVX)
         end if
        Else
         if (DoSplitCAS) then !(GLMJ)
           Call SplitCtl(Work(LW1),TUVX,IFINAL,iErrSplit)
           if (iErrSplit.eq.1) then
            write(LF,*) ('*',i=1,120)
            write(LF,*)'WARNING!!!'
            write(LF,*) 'SplitCAS iterations don''t converge.'
            write(LF,*) 'The program will continue'
            write(LF,*) 'Hopefully your calculation will converge',
     &                 'next iteration!'
            write(LF,*) ('*',i=1,120)
           end if
           if (iErrSplit.eq.2) then
            write(LF,*) ('*',i=1,120)
            write(LF,*)'WARNING!!!'
            write(LF,*) 'SplitCAS iterations don''t converge.'
            write(LF,*) 'MxIterSplit', MxIterSplit
            write(LF,*) 'SplitCAS ThreShold', ThrSplit
          write(LF,*)'Try to increase MxIterSplit or SplitCAS threshold'
            write(LF,*) 'The program will STOP'
            write(LF,*) ('*',i=1,120)
            call xQuit(96)
           end if
         end if
         if(.not.DoSplitCAS) then
! Call DMRG staff in Molcas - yingjin
           if(doDMRG)then
#ifdef _DMRG_
                                                  twordm_qcm = .true.
             if(KeyCION .and. .not. domcpdftDMRG) twordm_qcm = .false.
             call dmrg_interface_ctl(
     &                               task = 'run DMRG',
     &                               Key_CION = .not.twordm_qcm,
     &                               IterSCF  = Iter
     &                              )
! Keep the root energies
             Do jRoot = 1,lRoots
                ENER(jRoot,ITER)=dmrg_energy%dmrg_state_specific(jroot)
             End Do
#endif
           else
! Normal Davidson algorithm
             Call DavCtl(Work(LW1),TUVX,IFINAL)
           end if
         end if
        End If
      endif

*
* C
* CALCULATE DENSITY MATRICES
* SAVE DENSITY MATRICES ON FILE
* COMPUTE AVERAGE DENSITY MATRICES
* C
*
* LW6: ONE-BODY DENSITY
* LW7: ONE-BODY SPIN DENSITY
* LW8: SYMMETRIC TWO-BODY DENSITY
* LW9: ANTISYMMETRIC TWO-BODY DENSITY
*
      Call Timing(Rado_1,Swatch,Swatch,Swatch)
      Call dCopy_(NACPAR,[0.0D0],0,D,1)
      Call dCopy_(NACPAR,[0.0D0],0,DS,1)
      Call dCopy_(NACPR2,[0.0D0],0,P,1)
      Call dCopy_(NACPR2,[0.0D0],0,PA,1)
      CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
      CALL GETMEM('Dtmp ','ALLO','REAL',LW6,NACPAR)
      CALL GETMEM('DStmp','ALLO','REAL',LW7,NACPAR)
      CALL GETMEM('Ptmp ','ALLO','REAL',LW8,NACPR2)
      CALL GETMEM('PAtmp','ALLO','REAL',LW9,NACPR2)
      CALL GETMEM('Pscr','ALLO','REAL',LW10,NACPR2)
#ifdef _HDF5_
      call mma_allocate(density_square, nac, nac)
#endif
      iDisk = IADR15(4)
      jDisk = IADR15(3)
      IF (.not.DoSplitCAS) THEN
*JB   Instead of RASSCF/RASCI energy, print out energy for rotated
*JB   states
       do_rotate=.False.
       If (ifinal.eq.2) Then
        If(IRotPsi==1) Then
         CALL f_inquire('ROT_VEC',Do_Rotate)
        End If
        If(Do_Rotate) Then
         CALL RotState()
        Else
         If(IRotPsi==1) Then
          write(LF,'(6X,A,A)')'Do_Rotate.txt is not found. ',
     &   'MCSCF states will not be rotated'
         End If
        End If
*JB    End of condition 'Do_Rotate' to initialize rotated states
       End If
*JB    End If for ifinal=2
       Do jRoot = 1,lRoots
* load back one CI vector at the time
*JB      If do_rotate=.true., then we read CI vectors from Work(LRCIVec)
*JB      Otherwise we read if from JOBIPH
         Call DDafile(JOBIPH,2,Work(LW4),nConf,iDisk)
         IF (IPRLEV.GE.DEBUG) THEN
          call DVcPrt('CI-Vec in CICTL',' ',Work(LW4),nConf )
          Write(LF,*) ' WORK SPACE VARIABLES IN SUBR. CICTL: '
          Write(LF,'(1x,A,5I10)') 'DENSI',LW4,LW6,LW7,LW8,LW9
         END IF
* compute density matrices

         If ( NAC.ge.1 ) Then
           C_Pointer = Lw4
           if(.not.(doDMRG))
     &       CALL Lucia_Util('Densi',0,iDummy,rdum)
           IF ( IPRLEV.GE.INSANE  ) THEN
             write(6,*) 'At root number =', jroot
             CALL TRIPRT('D after lucia  ',' ',Work(LW6),NAC)
             CALL TRIPRT('DS after lucia  ',' ',Work(LW7),NAC)
             CALL TRIPRT('P after lucia',' ',Work(LW8),NACPAR)
             CALL TRIPRT('PA after lucia',' ',Work(LW9),NACPAR)
           END IF
         EndIf
         IF (.not.doDMRG .and. (IFCAS.GT.2 .OR. iDoGAS))
     &   CALL CISX(IDXSX,Work(LW6),Work(LW7),
     &               Work(LW8),Work(LW9),Work(LW10))
! 1,2-RDMs importing from DMRG calculation -- Stefan/Yingjin
         if(doDMRG)then
#ifdef _DMRG_
           call dmrg_interface_ctl(
     &                             task  = 'imp rdmX',
     &                             x1    = work(lw6:lw6+NACPAR-1),
     &                             x2    = work(lw8:lw8+NACPR2-1),
     &                             ndim  = nacpar,
     &                             mdim  = nacpr2,
     &                             state = jroot,
     &                             rdm1  = .true.,
     &                             rdm2  = twordm_qcm
     &                            )

           !> import 1p-spin density
           call dcopy_(NACPAR,[0.0D0],0,work(lw7),1)
           call dmrg_interface_ctl(
     &                             task  = 'imp spdX',
     &                             x1    = work(lw7:lw7+NACPAR-1),
     &                             ndim  = nacpar,
     &                             state = jroot,
     &                             rdm1  = .true.
     &                            )
#ifdef _DMRG_DEBUG_
           write(6,*)"==============================================="
           write(6,*)"  Set all elems in anti-symmetric 2-RDM to zero"
           write(6,*)"==============================================="
#endif
           call dcopy_(NACPR2,[0.0D0],0,work(lw9),1)

           IF ( IPRLEV.GE.INSANE  ) THEN
             CALL TRIPRT('D after  DMRG',' ',Work(LW6),NAC)
             CALL TRIPRT('DS after DMRG',' ',Work(LW7),NAC)
             CALL TRIPRT('P after  DMRG',' ',Work(LW8),NACPAR)
             CALL TRIPRT('PA after DMRG',' ',Work(LW9),NACPAR)
           END IF
#endif
         end if
* Modify the symmetric 2-particle density if only partial
* "exact exchange" is included.
c        n_Det=2
c        n_unpaired_elec=(iSpin-1)
c        n_paired_elec=nActEl-n_unpaired_elec
c        If(n_unpaired_elec+n_paired_elec/2.eq.nac) n_Det=1
c
c           Write(LF,*) ' iSpin=',iSpin
c           Write(LF,*) ' n_unpaired_elec',n_unpaired_elec
c           Write(LF,*) ' n_paired_elec',  n_paired_elec
c           Write(LF,*) ' n_unpaired_elec+n_paired_elec/2',
c     &                  n_unpaired_elec+n_paired_elec/2
c           Write(LF,*) ' n_Det=',n_Det
c
c

*          write(6,*) 'second call to Mod_P2'

           If (ExFac.ne.1.0D0.AND.(.not.l_casdft))
     &                     Call Mod_P2(Work(LW8),NACPR2,
     &                                 Work(LW6),NACPAR,
     &                                 Work(LW7),ExFac,n_Det)

* update average density matrices
         Scal = 0.0d0
         Do kRoot = 1,nRoots
           If ( iRoot(kRoot).eq.jRoot ) then
             Scal = Weight(kRoot)
           End If
         End Do
         Call daXpY_(NACPAR,Scal,Work(LW6),1,D,1)
         Call daXpY_(NACPAR,Scal,Work(LW7),1,DS,1)
         Call daXpY_(NACPR2,Scal,Work(LW8),1,P,1)
cGLM Put the D1MO and the P2MO values in RUNFILE
*
         Call Put_D1MO(Work(LW6),NACPAR) ! Put it on the RUNFILE
         Call Put_P2MO(Work(LW8),NACPR2) ! Put it on the RUNFILE
         Call daXpY_(NACPR2,Scal,Work(LW9),1,PA,1)
* save density matrices on disk
         Call DDafile(JOBIPH,1,Work(LW6),NACPAR,jDisk)
         Call DDafile(JOBIPH,1,Work(LW7),NACPAR,jDisk)
         Call DDafile(JOBIPH,1,Work(LW8),NACPR2,jDisk)
         Call DDafile(JOBIPH,1,Work(LW9),NACPR2,jDisk)
CSVC: store a single column instead of the whole array (which is for each root!)
C and for now don't bother with 2-electron active density matrices
#ifdef _HDF5_
         call square(work(lw6),density_square,1,nac,nac)
         call mh5_put_dset_array_real(wfn_dens, density_square,
     $           [nac,nac,1], [0,0,jRoot-1])
         call square(work(lw7),density_square,1,nac,nac)
         call mh5_put_dset_array_real(wfn_spindens, density_square,
     $           [nac,nac,1], [0,0,jRoot-1])
#endif
       End Do

      ELSE  ! SplitCAS run
        Call DDafile(JOBIPH,2,Work(LW4),nConf,iDisk)
        IF (IPRLEV.GE.DEBUG) then
          call DVcPrt('CI-Vec in CICTL SplitCAS sect',' ',
     &              Work(LW4),nConf)
        end if
* compute density matrices
        If ( NAC.ge.1 ) Then
           C_Pointer = Lw4
           CALL Lucia_Util('Densi',0,iDummy,rdum)
           IF ( IPRLEV.GE.INSANE  ) THEN
             CALL TRIPRT('D after lucia',' ',Work(LW6),NAC)
             CALL TRIPRT('DS after lucia',' ',Work(LW7),NAC)
             CALL TRIPRT('P after lucia',' ',Work(LW8),NACPAR)
             CALL TRIPRT('PA after lucia',' ',Work(LW9),NACPAR)
           END IF
        EndIf
        IF (IDoGAS.or.ifcas.gt.2) CALL CISX(IDXSX,Work(LW6),Work(LW7),
     &              Work(LW8),Work(LW9),Work(LW10))
           If (ExFac.ne.1.0D0.AND.(.not.l_casdft))
     &                      Call Mod_P2(Work(LW8),NACPR2,
     &                                Work(LW6),NACPAR,
     &                                Work(LW7),ExFac,n_Det)
        Scal = 1.0d0
        call daxpy_(NACPAR,Scal,Work(LW6),1,D,1)
        call daxpy_(NACPAR,Scal,Work(LW7),1,DS,1)
        call daxpy_(NACPR2,Scal,Work(LW8),1,P,1)
        call daxpy_(NACPR2,Scal,Work(LW9),1,PA,1)
* save density matrices on disk
        Call DDafile(JOBIPH,1,Work(LW6),NACPAR,jDisk)
        Call DDafile(JOBIPH,1,Work(LW7),NACPAR,jDisk)
        Call DDafile(JOBIPH,1,Work(LW8),NACPR2,jDisk)
        Call DDafile(JOBIPH,1,Work(LW9),NACPR2,jDisk)
      END IF

#ifdef _HDF5_
      call mma_deallocate(density_square)
#endif
      CALL GETMEM('Pscr','FREE','REAL',LW10,NACPR2)
      CALL GETMEM('PAtmp','FREE','REAL',LW9,NACPAR)
      CALL GETMEM('Ptmp ','FREE','REAL',LW8,NACPAR)
      CALL GETMEM('DStmp','FREE','REAL',LW7,NACPR2)
      CALL GETMEM('Dtmp ','FREE','REAL',LW6,NACPR2)
*
* print matrices
      IF ( IPRLEV.GE.INSANE  ) THEN
        CALL TRIPRT('Averaged one-body density matrix, D',
     &              ' ',D,NAC)
        CALL TRIPRT('Averaged one-body spin density matrix, DS',
     &              ' ',DS,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',P,NACPAR)
        CALL TRIPRT('Averaged antisymmetric two-body density matrix,PA',
     &              ' ',PA,NACPAR)
      END IF
      Call Put_D1MO(D,NACPAR) ! Put it on the RUNFILE
c
      IF ( NASH(1).NE.NAC ) CALL DBLOCK(D)
      Call Timing(Rado_2,Swatch,Swatch,Swatch)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2
*
* C
* IF FINAL ITERATION REORDER THE WAVEFUNCTION ACCORDING TO
* THE SPLIT GRAPH GUGA CONVENTIONS AND PRINT IT.
* C
*
* LW11: Temporary copy of a CI vector
*
      IF (IFINAL.EQ.2 .AND. NAC.GT.0 ) THEN
       IF (IPRLEV.ge.USUAL) THEN
        Write(LF,*)
        Write(LF,'(6X,120(1H*))')
        Write(LF,'(54X,A)') 'Wave function printout:'
        Write(LF,'(23X,A)') 'occupation of active orbitals, and '//
     &                      'spin coupling of open shells '//
     &                      '(u,d: Spin up or down)'
        Write(LF,'(6X,120(1H*))')
        Write(LF,*)
        Write(6,'(6x,A)') 'Note: transformation to natural orbitals'
        Write(6,'(6x,A)')
     &     'has been made, which may change the order of the CSFs.'
       END IF
       Call GetMem('PrSel','Allo','Inte',LW12,nConf)
       Call iCopy(nConf,[0],0,iWork(LW12),1)
       Call GetMem('CIVtmp','Allo','Real',LW11,nConf)
       iDisk = IADR15(4)

       IF (.Not.DoSplitCAS) THEN
        Do i = 1,lRoots
          jDisk=iDisk
* load back one CI vector at the time
*          Call DDafile(JOBIPH,2,Work(LW4),nConf,iDisk)
           Call DDafile(JOBIPH,2,Work(LW4),nConf,iDisk)
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL last cycle',' ',
     &        Work(LW4),nConf)
           Write(LF,*) ' WORK SPACE VARIABLES IN SUBR. CICTL: '
           Write(LF,'(1x,A,2I10)') 'REORD',LW4,LW11
          END IF
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
         if(.not.iDoGas)then
          Call Reord2(NAC,NACTEL,LSYM,0,
     &                iWork(KICONF(1)),iWork(KCFTP),
     &                Work(LW4),Work(LW11),iWork(ivkcnf))
c        end if
c         call getmem('kcnf','free','inte',ivkcnf,nactel)

* save reorder CI vector on disk
c         if(.not.iDoGas)then
          Call DDafile(JOBIPH,1,Work(LW11),nConf,jDisk)
#ifdef _HDF5_
          call mh5_put_dset_array_real
     $            (wfn_cicoef,Work(LW11),[nconf,1],[0,i-1])

#endif
c         else
c         call DDafile(JOBIPH,1,Work(LW4),nConf,jDisk)
c         end if
* printout of the wave function
           if(doDMRG)then
             ! If DMRG, the SRCAS can give the CI-coefficients
           else
             IF (IPRLEV.GE.USUAL) THEN
              Write(LF,*)
              Write(LF,'(6X,A,F6.2,A,I3)')
     &                'printout of CI-coefficients larger than',
     &                 PRWTHR,' for root',i
              Write(LF,'(6X,A,F15.6)')
     &                 'energy=',ENER(I,ITER)
               CALL SGPRWF(iWork(LW12),IWORK(LNOCSF),IWORK(LIOCSF),
     &                  IWORK(LNOW),IWORK(LIOW),WORK(LW11))
             End If
           end if
         else ! for iDoGas
          Write(LF,'(1x,a)') 'WARNING: true GAS, JOBIPH not compatible!'
c.. save CI vector on disk
          Call DDafile(JOBIPH,1,Work(LW4),nconf,jDisk)
CSVC: store CI as a column array of the on-disk CI (which is for all roots!)
#ifdef _HDF5_
          call mh5_put_dset_array_real
     $            (wfn_cicoef,Work(LW4),[nconf,1],[0,i-1])
#endif
C.. printout of the wave function
          IF (IPRLEV.GE.USUAL) THEN
            Write(LF,*)
            Write(LF,'(6X,A,F6.2,A,I3)')
     c                'printout of CI-coefficients larger than',
     c                 prwthr,' for root', i
            Write(LF,'(6X,A,F15.6)')
     c                'energy=',ener(i,iter)
          call gasprwf(iwork(lw12),nac,nactel,lsym,iwork(kiconf(1)),
     c                 iwork(kcftp),work(lw4),iwork(ivkcnf))
          End If
         end if
          call getmem('kcnf','free','inte',ivkcnf,nactel)
*          END IF
        End Do

        ELSE !RUN SPLITCAS

          jDisk=iDisk
* load back one CI vector at the time
          Call DDafile(JOBIPH,2,Work(LW4),nConf,iDisk)
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL SplitCAS last cycle',' ',
     &        Work(LW4),nConf)
          END IF
          IF (IPRLEV.GE.DEBUG) THEN
           Write(LF,*) ' WORK SPACE VARIABLES IN SUBR. CICTL: '
           Write(LF,'(1x,A,2I10)') 'REORD',LW4,LW11
          END IF
* reorder it according to the split graph GUGA conventions
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
          Call Reord2(NAC,NACTEL,LSYM,0,
     &                iWork(KICONF(1)),iWork(KCFTP),
     &                Work(LW4),Work(LW11),iWork(ivkcnf))
          call getmem('kcnf','free','inte',ivkcnf,nactel)
* save reorder CI vector on disk
          Call DDafile(JOBIPH,1,Work(LW11),nConf,jDisk)
#ifdef _HDF5_
          call mh5_put_dset_array_real
     $            (wfn_cicoef,Work(LW11),[nconf,1],[0,i-1])
#endif
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL after Reord',' ',
     &        Work(LW11),nConf)
          END IF
* printout of the wave function
          IF (IPRLEV.GE.USUAL) THEN
            Write(LF,*)
            Write(LF,'(6X,A,F6.2,A,I3)')
     &                'printout of CI-coefficients larger than',
     &                 PRWTHR,' for root',lRootSplit
            Write(LF,'(6X,A,F15.6)')
     &                'Split-energy=',ENER(lRootSplit,ITER)
            CALL SGPRWF(iWork(LW12),IWORK(LNOCSF),IWORK(LIOCSF),
     &                IWORK(LNOW),IWORK(LIOW),WORK(LW11))
          END IF
        END IF

        Call GetMem('PrSel','Free','Inte',LW12,nConf)
        Call GetMem('CIVtmp','Free','Real',LW11,nConf)
      ENDIF
#ifdef _DMRG_
          call mh5_put_dset_array_str
     &         (wfn_dmrg_checkpoint,dmrg_file%qcmaquis_checkpoint_file)
#endif

      CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)
      CALL GETMEM('CICTL1','FREE','REAL',LW1,NACPAR)

 9000 Continue
*
*     For RF calculations make sure that the we are following the
*     correct root.
*
*     In the current implementation the overlap between the CI vectors
*     of different macro iterations is used. This criterion stricktly
*     only hold if the orbitals are not changed in between the
*     interations, to make sure that this approximately holds the
*     comparision is only to be considered to be valid if the rotmax
*     parameter is below an empirical threshold. In the future the
*     procedure for indentifying root flipping has to be made more
*     robust.
*
CSVC: if CISElect is used, roots should be traced regardless of orbital
C     rotation, so in that case ignore the automatic tracing and follow
C     the relative CISE root given in the input by the 'CIRF' keyword.
*
      If (lRF.and.KeyCISE.and.KeyCIRF) Then
        JPCMROOT=IPCMROOT
        IPCMROOT=IROOT(ICIRFROOT)
        Call Put_iScalar("RF CASSCF root",IPCMROOT)
        If (JPCMROOT.ne.IPCMROOT) Then
          Write (6,'(1X,A,I3,A,I3)') 'RF Root has flipped from ',
     &                 JPCMROOT, ' to ',IPCMROOT
        End If
      Else If (lRF) Then
        Call Qpg_iScalar('RF CASSCF root',Exist)
        If (.NOT.Exist) Then
*
*          We are here since we are using the default values.
*
           Call Put_iScalar("RF CASSCF root",IPCMROOT)
           Call Put_iScalar("RF0CASSCF root",IPCMROOT)
        End If
*
        mconf = 0
        Call Allocate_Work(ipRF,nConf)
        Call Qpg_dArray("RF CASSCF Vector",Exist,mConf)

        !> check whether the rf target h5 file exists (needed at this
        !point for numerical gradient calculations)
#ifdef _DMRG_
        if(doDMRG.and.exist)then
          inquire(file="rf.results_state.h5", exist=rfh5DMRG)
          if(.not.rfh5DMRG)then
            maquis_name_states  = ""
            maquis_name_results = ""
            call file_name_generator(IPCMROOT-1,"checkpoint_state.",
     &                               17,".h5",3,maquis_name_states)
            call file_name_generator(IPCMROOT-1,"results_state.",
     &                               14,".h5",3,maquis_name_results)

          !> copy current target wave function to local wave function
            call systemf(
     & "cp -f "//trim(maquis_name_results)//" rf.results_state.h5 && "//
     & "rm -rf rf.checkpoint_state.h5 && "//
     & "cp -r "//trim(maquis_name_states)//" rf.checkpoint_state.h5",
     &                 iErr)
          end if
        end if
#endif

        If (Exist
     &      .and. mConf .eq. nConf
     &      .and. iFinal.ne.2
     &      .and. (ABS(RotMax).lt.1.0D-3 .or. KeyCISE)
     &     ) Then

           if(doDMRG)then
#ifdef _DMRG_
             overlap = 0.0d0
             call dmrg_interface_ctl(
     &                               task   = 'overlapR',
     &                               energy = overlap,
     &                               stateL = -1
     &                              )
             rnorm = sqrt(overlap)
#endif
           else
             Call Get_dArray("RF CASSCF Vector",Work(ipRF),nConf)
             rNorm=Sqrt(DDot_(nConf,Work(ipRF),1,Work(ipRF),1))
           end if
*          Write (6,*) 'rNorm=',rNorm
           JPCMROOT=IPCMROOT
           If (rNorm.gt.1.0D-10) Then
              Call Allocate_Work(ipTemp,nConf)
              rMax=0.0D0
              jDisk = IADR15(4)
              Do i = 1, lRoots
                 if(doDMRG)then
#ifdef _DMRG_
                   overlap = 0.0d0
                   call dmrg_interface_ctl(
     &                                     task   = 'overlapR',
     &                                     energy = overlap,
     &                                     stateL = i
     &                                    )
                   qmax = abs(overlap)
#endif
                 else
                   Call DDafile(JOBIPH,2,Work(ipTemp),nConf,jDisk)
                   qMax=Abs(DDot_(nConf,Work(ipTemp),1,Work(ipRF),1))
                 end if
*                Write (6,*) 'qMax=',qMax
                 If (qMax.gt.rMax .and.
     &               qMax.gt.0.5D0) Then
                    rMax = qMax
                    JPCMROOT=i
                 End If
              End Do
              Call Free_Work(ipTemp)
           End If
        Else
           JPCMROOT=IPCMROOT
        End If
*
        If (JPCMROOT.ne.IPCMROOT) Then
           Write (6,*) ' RF Root has flipped from ',IPCMROOT, ' to ',
     &                                              JPCMROOT
           IPCMROOT=JPCMROOT
           Call Put_iScalar("RF CASSCF root",IPCMROOT)
        End If
*
        if(doDMRG)then
#ifdef _DMRG_
          maquis_name_states  = ""
          maquis_name_results = ""
          call file_name_generator(IPCMROOT-1,"checkpoint_state.",
     &                             17,".h5",3,maquis_name_states)
          call file_name_generator(IPCMROOT-1,"results_state.",
     &                             14,".h5",3,maquis_name_results)

!         print *, 'my string: ',
!    & "cp -f "//trim(maquis_name_results)//" rf.results_state.h5 && "//
!    & "rm -rf rf.checkpoint_state.h5 && "//
!    & "cp -r "//trim(maquis_name_states)//" rf.checkpoint_state.h5"

          !> copy current target wave function to local wave function
          call system(
     & "cp -f "//trim(maquis_name_results)//" rf.results_state.h5 && "//
     & "rm -rf rf.checkpoint_state.h5 && "//
     & "cp -r "//trim(maquis_name_states)//" rf.checkpoint_state.h5",
     &               iErr)
#endif
        else
          jDisk = IADR15(4)
          Do i=1,IPCMROOT-1
            Call DDafile(JOBIPH,0,rdum,nConf,jDisk)
          End Do
          Call DDafile(JOBIPH,2,Work(ipRF),nConf,jDisk)
        end if

        Call Put_dArray("RF CASSCF Vector",Work(ipRF),nConf)
        Call Free_Work(ipRF)
      End If

      Call qExit('CICTL')
      Return
      End

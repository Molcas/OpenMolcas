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
*> @param         FA
*> @param[in,out] D1I    Inactive 1-dens matrix
*> @param[in,out] D1A    Active 1-dens matrix
*> @param[in]     TUVX   Active 2-el integrals
*> @param[in]     IFINAL Calculation status switch
************************************************************************
      Subroutine CICtl(CMO,D,DS,P,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
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
      use qcmaquis_interface
      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines,
     &     only: fiedlerorder_length, file_name_generator,
     &           qcmaquis_interface_fcidump
#endif
#ifdef _HDF5_
      use mh5, only: mh5_put_dset
#endif
      use csfbas, only: CONF
      use glbbas, only: CFTP
      use casvb_global, only: ifvb
      use CMS, only: iCMSOpt,CMSGiveOpt
      use rctfld_module
      use rasscf_lucia, only: PAtmp, Pscr, CIVEC, PTmp, DStmp, Dtmp
#ifdef _DMRG_
      use rasscf_lucia, only: RF1, RF2
#endif
      use Lucia_Interface, only: Lucia_Util
      use wadr, only: FMO
      use gugx, only: IFCAS,  NOCSF,  IOCSF, NOW1, IOW1
      use sxci, only: IDXSX

      Implicit Real* 8 (A-H,O-Z)

      Dimension CMO(*),D(*),DS(*),P(*),PA(*),FI(*),FA(*),D1I(*),D1A(*),
     &          TUVX(*)
      Logical Exist,Do_ESPF
*JB   variables for state rotation on final states
      Logical do_rotate
#ifdef _DMRG_
      ! function defined in misc_util/pcm_on.f
      Logical, external :: PCM_On
#endif
      ! Filename used to write GronOR vecdet files (tps/cdg 20210430)
      character(len=128) :: filename
      integer, external :: IsFreeUnit

#include "rasdim.fh"
#include "rasscf.fh"
#include "splitcas.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Character*16 ROUTINE
      Parameter (ROUTINE='CICTL   ')
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "timers.fh"
#include "pamint.fh"
#include "input_ras.fh"
#include "stdalloc.fh"
#ifdef _HDF5_
#include "raswfn.fh"
      real*8, allocatable :: density_square(:,:)
#endif

#ifdef _DMRG_
      character(len=2300) :: maquis_name_states
      character(len=2300) :: maquis_name_results
      logical             :: rfh5DMRG
      ! check if we calculate entanglement and spin density
      ! we do it only in the last iteration
      logical             :: doEntanglement
      character(len=:), allocatable :: fiedler_order_str

      ! arrays for 1- and 2-RDMs and spin-1-RDMs, size: nrdm x nroots
      real*8, allocatable :: d1all(:,:), d2all(:,:), spd1all(:,:)
c #include "nevptp.fh"
#endif
      Dimension rdum(1)
      Real*8, Allocatable:: CIV(:)
      Integer, Allocatable:: PrSel(:)

*PAM05      SymProd(i,j)=1+iEor(i-1,j-1)
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
* FMO: FOCK MATRIX IN MO-BASIS
* LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION
*
      CALL mma_allocate(FMO,NACPAR,Label='FMO')
      Call DecideOnESPF(Do_ESPF)

! initialize RDM arrays for QCMaquis
#ifdef _DMRG_
      if (doDMRG) then
        ! Calculate entanglement and spin densities
        ! only in the last iteration
        ! i.e. only for IFINAL.eq.2 if we do DMRG-SCF
        ! otherwise always
        doEntanglement = merge(.true.,IFINAL.eq.2,KeyCION)

        call mma_allocate(d1all, NACPAR, lRoots)
        d1all = 0.0d0
        if (twordm_qcm) then
          call mma_allocate(d2all, NACPR2, lRoots)
          d2all = 0.0d0
        end if
        ! Allocate spin density only for the last iteration
        if(doEntanglement) then
          call mma_allocate(spd1all, NACPAR, lRoots)
          spd1all = 0.0d0
        end if
      end if
#endif

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
           Call Put_dArray('D1mo',Work(LRCT),NACPAR)  ! Put on RUNFILE
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
           Call Put_dArray('P2mo',Work(ipP2MO),NACPR2) ! Put on RUNFILE
*
           CALL SGFCIN(CMO,FMO,FI,D1I,Work(LRCT_F),Work(LRCT_FS))
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
           Call mma_allocate(CIVEC,NCONF,Label='CIVEC')
           If (NACTEL.EQ.0) THEN
             CIVEC(1)=1.0D0
           Else
             if(.not.(doDMRG))then
!               write(*,*)"run the load back CI vector part" ! yma
               iDisk = IADR15(4)
               Do jRoot = 1,IPCMROOT
                   iOpt=0
                   If (jRoot.eq.IPCMROOT) iOpt=2
* load back one CI vector at the time
                   Call DDafile(JOBIPH,iOpt,CIVEC,nConf,iDisk)
                  !call DVcPrt('BLUBB-start CI PCM',' ',CIVEC,nConf)
               End Do
             end if
           End If

* compute density matrices

           Call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
           Call mma_allocate(DStmp,NAC**2,Label='DStmp')
           Call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
           If ( NAC.ge.1 ) Then

              If (NACTEL.eq.0) THEN
                 Dtmp(:)=0.0D0
                 DStmp(:)=0.0D0
                 Ptmp(:)=0.0D0
              Else

                if(doDMRG)then
#ifdef _DMRG_
                ! copy the DMs from d1rf/d2rf for ipcmroot
                call dcopy_(NACPAR,rf1,1,Dtmp,1)
                if (twordm_qcm) then
                  call dcopy_(NACPR2,rf2,1,Ptmp,1)
                end if

        ! Import RDMs from QCMaquis that we've got from the last optimization
        ! Here we should import one-particle spin density.
        ! However, the spin density has been temporarily disabled here:
        ! For performance reasons, it is calculated
        ! only once in the last iteration of DMRG-SCF optimisation.
        ! If you need it at every iteration for some reason
        ! please change this code accordingly
                DStmp(:)=0.0D0
#endif
               else
                 Call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
                 Call mma_allocate(Pscr,NACPR2,Label='Pscr')
                 CALL Lucia_Util('Densi',
     &                           CI_Vector=CIVEC(:))
                 If (IFCAS.GT.2 .OR. iDoGAS) Then
                   Call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
                 End If
                 Call mma_deallocate(Pscr)
                 Call mma_deallocate(PAtmp)
               end if ! doDMRG/doBLOK or CI

             End If
           Else
              Dtmp(:)=0.0D0
              DStmp(:)=0.0D0
              Ptmp(:)=0.0D0
           End If
* Modify the symmetric 2-particle density if only partial
* "exact exchange" is included.
c          n_Det=2
c          n_unpaired_elec=(iSpin-1)
c          n_paired_elec=nActEl-n_unpaired_elec
c          If(n_unpaired_elec+n_paired_elec/2.eq.nac) n_Det=1
*                  write(6,*) 'n_Det =', n_Det
           If (ExFac.ne.1.0D0.AND.(.not.l_casdft))
     &    Call Mod_P2(Ptmp,NACPR2,
     &                                   Dtmp,NACPAR,
     &                                   DStmp,ExFac,n_Det)
*
           Call Put_dArray('P2mo',Ptmp,NACPR2) ! Put on RUNFILE
*
           Call mma_deallocate(Ptmp)
*
           Call Put_dArray('D1mo',Dtmp,NACPAR) ! Put on RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Dtmp)
           Call Get_D1A_RASSCF(CMO,Dtmp,Work(LRCT_F))
*
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(DStmp)
           Call Get_D1A_RASSCF(CMO,DStmp,Work(LRCT_FS))
*
!           do i=1,NACPAR  !yma
!              write(*,*)"i-rdms1 1",i,Dtmp(i)
!           end do

           Call mma_deallocate(DStmp)
           Call mma_deallocate(Dtmp)
           Call mma_deallocate(CIVEC)
*
           Call SGFCIN(CMO,FMO,FI,D1I,Work(LRCT_F),Work(LRCT_FS))
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
        CALL SGFCIN(CMO,FMO,FI,D1I,D1A,Work(ipTmpD1S))
        CALL GETMEM('TmpD1S','Free','REAL',ipTmpD1S,NTOT2 )
*
      END IF

      if(doDMRG)then
#ifdef _DMRG_
        ! update integrals for QCMaquis

        call qcmaquis_interface_update_integrals(FMO,tuvx,emy)

        !!! Fiedler order/CI-DEAS run
        if (dmrg_warmup%dofiedler.or.dmrg_warmup%docideas) then
        ! allocate string where Fiedler ordering will be returned with the correct length
          ilen = fiedlerorder_length(qcmaquis_param%L)
          allocate(character(len=ilen) :: fiedler_order_str)
          fiedler_order_str(:) = ' '


        ! if HF guess is present, use it (required for CI-DEAS)
        ! If not, error handling should be done in QCMaquis
          if(sum(dmrg_orbital_space%initial_occ) > 0)then
              call qcmaquis_interface_run_starting_guess(nRoots,
     &     dmrg_warmup%dofiedler,
     &     dmrg_warmup%docideas,
     &     fiedler_order_str,
           ! pass the HF occupation as 1D array to QCMaquis
     &     reshape(dmrg_orbital_space%initial_occ,
     &          (/1, sum(nash)*nroots/)))
          else
            call qcmaquis_interface_run_starting_guess(nRoots,
     &       dmrg_warmup%dofiedler,
     &       dmrg_warmup%docideas,
     &       fiedler_order_str)
          endif

          if (dmrg_warmup%dofiedler)
     &      call qcmaquis_interface_set_param("orbital_order",
     &      fiedler_order_str)
            write (6,*) "Fiedler orbital ordering: "//fiedler_order_str

          dmrg_warmup%dofiedler = .false.
          dmrg_warmup%docideas = .false.
          if(allocated(fiedler_order_str)) deallocate(fiedler_order_str)
        end if
        if(dofcidump)then
          ! Produce a FCIDUMP file
          ! TODO:
          ! We already have the fcidump module in rasscf, which is called elsewhere
          ! so ensure the compatibility of the FCIDUMP files produced by this module
          ! and remove the code below
          call qcmaquis_interface_fcidump(FMO,tuvx,emy)
          Call mma_deallocate(FMO)
          goto 9000
        end if
#endif
      end if
*
      If (IfVB.eq.2) GoTo 9000

*
* C
* DAVIDSON DIAGONALIZATION
* C
*
      if(IfVB.eq.1)then
        call cvbmn_rvb(max(ifinal,1))
      else
         if (DoSplitCAS) then !(GLMJ)
           Call SplitCtl(FMO,TUVX,IFINAL,iErrSplit)
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
           if(doDMRG)then
#ifdef _DMRG_
             ! Get also spin density at the last iteration
             ! Please help me call it more cleanly than with these if clauses
             ! and different optional arguments
             if (doEntanglement) then
               if (twordm_qcm) then
                 call qcmaquis_interface_run_dmrg(nstates=lroots,
     &             d1=d1all, d2=d2all, spd=spd1all,
     &             entanglement=doEntanglement)
               else
                 call qcmaquis_interface_run_dmrg(nstates=lroots,
     &             d1=d1all, spd=spd1all, entanglement=doEntanglement)
               end if
             else
               if (twordm_qcm) then
                 call qcmaquis_interface_run_dmrg(nstates=lroots,
     &             d1=d1all, d2=d2all, entanglement=doEntanglement)
               else
                 call qcmaquis_interface_run_dmrg(nstates=lroots,
     &               d1=d1all, entanglement=doEntanglement)
               end if
             end if

             ! For PCM calculations: copy RDMs for the PCM root
             if (PCM_On()) then
               call dcopy_(NACPAR,d1all(:,ipcmroot),1,rf1,1)
               if (twordm_qcm) then
                 call dcopy_(NACPR2,d2all(:,ipcmroot),1,rf2,1)
               end if
             end if
! Keep the root energies
             Do jRoot = 1,lRoots
                ENER(jRoot,ITER)=dmrg_energy%dmrg_state_specific(jroot)
             End Do
! The new QCMaquis interface requires that the density matrices are calculated immediately after the DMRG run
! So we either need to keep them all in memory, or move the saving routines up here.
! The 2nd option requires code refactoring, so for now we keep them all in memory.
#endif
           else
! Normal Davidson algorithm
             Call DavCtl(FMO,TUVX,IFINAL)
           end if
         end if
      endif

*
* C
* CALCULATE DENSITY MATRICES
* SAVE DENSITY MATRICES ON FILE
* COMPUTE AVERAGE DENSITY MATRICES
* C
*
* Dtmp: ONE-BODY DENSITY
* DStmp: ONE-BODY SPIN DENSITY
* Ptmp: SYMMETRIC TWO-BODY DENSITY
* PAtmp: ANTISYMMETRIC TWO-BODY DENSITY
*
      Call Timing(Rado_1,dum1,dum2,dum3)
      Call dCopy_(NACPAR,[0.0D0],0,D,1)
      Call dCopy_(NACPAR,[0.0D0],0,DS,1)
      Call dCopy_(NACPR2,[0.0D0],0,P,1)
      Call dCopy_(NACPR2,[0.0D0],0,PA,1)
      CALL mma_allocate(CIVEC,NCONF,Label='CIVEC')
      CALL mma_allocate(Dtmp,NAC**2,Label='Dtmp')
      CALL mma_allocate(DStmp,NAC**2,Label='DStmp')
      CALL mma_allocate(Ptmp,NACPR2,Label='Ptmp')
      CALL mma_allocate(PAtmp,NACPR2,Label='PAtmp')
      CALL mma_allocate(Pscr,NACPR2,Label='Pscr')
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
        IF(IXMSP.eq.1) THEN
         CALL XMSRot(CMO,FI,FA)
        End If
        IF(ICMSP.eq.1) THEN
         If(trim(CMSStartMat).eq.'XMS') Then
          CALL XMSRot(CMO,FI,FA)
         End If
         If(.not.CMSGiveOpt) Then
          if(lRoots.eq.2) iCMSOpt=2
          if(lRoots.ge.3) iCMSOpt=1
         End If
         If(iCMSOpt.eq.1) Then
          CALL CMSOpt(TUVX)
         Else If (iCMSOpt.eq.2) Then
          CALL CMSRot(TUVX)
         End If
        END IF
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
         Call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
         IF (IPRLEV.GE.DEBUG) THEN
          call DVcPrt('CI-Vec in CICTL',' ',CIVEC,nConf )
         END IF
* compute density matrices

         If ( NAC.ge.1 ) Then
           if(.not.(doDMRG)) CALL Lucia_Util('Densi',
     &                                       CI_Vector=CIVEC(:))
           IF ( IPRLEV.GE.INSANE  ) THEN
             write(6,*) 'At root number =', jroot
             CALL TRIPRT('D after lucia  ',' ',Dtmp,NAC)
             CALL TRIPRT('DS after lucia  ',' ',DStmp,NAC)
             CALL TRIPRT('P after lucia',' ',Ptmp,NACPAR)
             CALL TRIPRT('PA after lucia',' ',PAtmp,NACPAR)
           END IF
         EndIf
         IF (.not.doDMRG .and. (IFCAS.GT.2 .OR. iDoGAS))
     &   CALL CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
! 1,2-RDMs importing from DMRG calculation -- Stefan/Yingjin
         if(doDMRG)then
#ifdef _DMRG_
          ! for QCMaquis, just copy the RDMs
          ! actually, copying is not needed! TODO
          call dcopy_(NACPAR,d1all(:,jroot),1,Dtmp,1)
          if (twordm_qcm) then
            call dcopy_(NACPR2,d2all(:,jroot),1,Ptmp,1)
          end if

           !> import 1p-spin density
           ! disable spin density if not in the last iteration
           if (doEntanglement) then
             call dcopy_(NACPAR,spd1all(:,jroot),1,DStmp,1)
           else
             DStmp(:)=0.0D0
           end if

           ! disable antisymmetric 2-RDM
           PAtmp(:)=0.0D0

           IF ( IPRLEV.GE.INSANE  ) THEN
             CALL TRIPRT('D after  DMRG',' ',Dtmp,NAC)
             CALL TRIPRT('DS after DMRG',' ',DStmp,NAC)
             CALL TRIPRT('P after  DMRG',' ',Ptmp,NACPAR)
             CALL TRIPRT('PA after DMRG',' ',PAtmp,NACPAR)
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
     &                     Call Mod_P2(Ptmp,NACPR2,
     &                                 Dtmp,NACPAR,
     &                                 DStmp,ExFac,n_Det)

* update average density matrices
         Scal = 0.0d0
         Do kRoot = 1,nRoots
           If ( iRoot(kRoot).eq.jRoot ) then
             Scal = Weight(kRoot)
           End If
         End Do
         Call daXpY_(NACPAR,Scal,Dtmp,1,D,1)
         Call daXpY_(NACPAR,Scal,DStmp,1,DS,1)
         Call daXpY_(NACPR2,Scal,Ptmp,1,P,1)
cGLM Put the D1MO and the P2MO values in RUNFILE
*
         Call Put_dArray('D1mo',Dtmp,NACPAR) ! Put on RUNFILE
         Call Put_dArray('P2mo',Ptmp,NACPR2) ! Put on RUNFILE
         Call daXpY_(NACPR2,Scal,PAtmp,1,PA,1)
* save density matrices on disk
         Call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
         Call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
         Call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
         Call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
CSVC: store a single column instead of the whole array (which is for each root!)
C and for now don't bother with 2-electron active density matrices
#ifdef _HDF5_
         call square(Dtmp,density_square,1,nac,nac)
         call mh5_put_dset(wfn_dens, density_square,
     $           [nac,nac,1], [0,0,jRoot-1])
         call square(DStmp,density_square,1,nac,nac)
         call mh5_put_dset(wfn_spindens, density_square,
     $           [nac,nac,1], [0,0,jRoot-1])
#endif
       End Do

      ELSE  ! SplitCAS run
        Call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
        IF (IPRLEV.GE.DEBUG) then
          call DVcPrt('CI-Vec in CICTL SplitCAS sect',' ',
     &              CIVEC,nConf)
        end if
* compute density matrices
        If ( NAC.ge.1 ) Then
           CALL Lucia_Util('Densi',
     &                     CI_Vector=CIVEC(:))
           IF ( IPRLEV.GE.INSANE  ) THEN
             CALL TRIPRT('D after lucia',' ',Dtmp,NAC)
             CALL TRIPRT('DS after lucia',' ',DStmp,NAC)
             CALL TRIPRT('P after lucia',' ',Ptmp,NACPAR)
             CALL TRIPRT('PA after lucia',' ',PAtmp,NACPAR)
           END IF
        EndIf
        IF (IDoGAS.or.ifcas.gt.2) CALL CISX(IDXSX,Dtmp,DStmp,
     &              Ptmp,PAtmp,Pscr)
           If (ExFac.ne.1.0D0.AND.(.not.l_casdft))
     &                      Call Mod_P2(Ptmp,NACPR2,
     &                                Dtmp,NACPAR,
     &                                DStmp,ExFac,n_Det)
        Scal = 1.0d0
        call daxpy_(NACPAR,Scal,Dtmp,1,D,1)
        call daxpy_(NACPAR,Scal,DStmp,1,DS,1)
        call daxpy_(NACPR2,Scal,Ptmp,1,P,1)
        call daxpy_(NACPR2,Scal,PAtmp,1,PA,1)
* save density matrices on disk
        Call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
        Call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
      END IF

#ifdef _HDF5_
      call mma_deallocate(density_square)
#endif
      Call mma_deallocate(Pscr)
      Call mma_deallocate(PAtmp)
      Call mma_deallocate(Ptmp)
      Call mma_deallocate(DStmp)
      Call mma_deallocate(Dtmp)
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
      Call Put_dArray('D1mo',D,NACPAR) ! Put on RUNFILE
c
      IF ( NASH(1).NE.NAC ) CALL DBLOCK(D)
      Call Timing(Rado_2,dum1,dum2,dum3)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2
*
* C
* IF FINAL ITERATION REORDER THE WAVEFUNCTION ACCORDING TO
* THE SPLIT GRAPH GUGA CONVENTIONS AND PRINT IT.
* C
*
* CIV: Temporary copy of a CI vector
*
      IF (IFINAL.EQ.2 .AND. NAC.GT.0 ) THEN
       IF (IPRLEV.ge.USUAL) THEN
        Write(LF,*)
        Write(LF,'(6X,120("*"))')
        Write(LF,'(54X,A)') 'Wave function printout:'
        Write(LF,'(23X,A)') 'occupation of active orbitals, and '//
     &                      'spin coupling of open shells '//
     &                      '(u,d: Spin up or down)'
        Write(LF,'(6X,120("*"))')
        Write(LF,*)
        Write(6,'(6x,A)') 'Note: transformation to natural orbitals'
        Write(6,'(6x,A)')
     &     'has been made, which may change the order of the CSFs.'
       END IF
       Call mma_allocate(PrSel,nConf,Label='PrSel')
       PrSel(:)=0
       Call mma_allocate(CIV,nConf,Label='CIV')
       iDisk = IADR15(4)

       if (.not.doDMRG) then
       IF (.Not.DoSplitCAS) THEN
        Do i = 1,lRoots
          jDisk=iDisk
* load back one CI vector at the time
           Call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL last cycle',' ',
     &        CIVEC,nConf)
          END IF
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
         if(.not.iDoGas)then
          Call Reord2(NAC,NACTEL,STSYM,0,
     &                CONF,CFTP,
     &                CIVEC,CIV,iWork(ivkcnf))
c        end if
c         call getmem('kcnf','free','inte',ivkcnf,nactel)

* save reorder CI vector on disk
c         if(.not.iDoGas)then
          Call DDafile(JOBIPH,1,CIV,nConf,jDisk)
#ifdef _HDF5_
          call mh5_put_dset(wfn_cicoef,CIV(1:nConf),[nconf,1],[0,i-1])

#endif
c         else
c         call DDafile(JOBIPH,1,CIVEC,nConf,jDisk)
c         end if
* printout of the wave function
            IF (IPRLEV.GE.USUAL) THEN
              Write(LF,*)
              Write(LF,'(6X,A,F6.2,A,I3)')
     &                  'printout of CI-coefficients larger than',
     &                   PRWTHR,' for root',i
              Write(LF,'(6X,A,F15.6)')
     &             'energy=',ENER(I,ITER)
              If (KeyPRSD) Then
!     Define filename to write GronOR vecdet files (tps/cdg 20210430)
                write(filename,'(a7,i1)') 'VECDET.',i
!     filename = 'VECDET.'//merge(str(i), 'x', i.lt.999)
                LuVecDet=39
                LuVecDet=IsFreeUnit(LuVecDet)
                call Molcas_open(LuVecDet,filename)
                write(LuVecDet,'(8i4)') nish
              End If
              CALL SGPRWF(PrSel,NOCSF,IOCSF,NOW1,IOW1,CIV)
!     Close GronOR vecdet file (tps/cdg 20210430)
              If (KeyPRSD) close(LuVecDet)
            End If
         else ! for iDoGas
          Write(LF,'(1x,a)') 'WARNING: true GAS, JOBIPH not compatible!'
c.. save CI vector on disk
          Call DDafile(JOBIPH,1,CIVEC,nconf,jDisk)
CSVC: store CI as a column array of the on-disk CI (which is for all roots!)
#ifdef _HDF5_
          call mh5_put_dset(wfn_cicoef,CIVEC(1:nconf),
     &                      [nconf,1],[0,i-1])
#endif
C.. printout of the wave function
          IF (IPRLEV.GE.USUAL) THEN
            Write(LF,*)
            Write(LF,'(6X,A,F6.2,A,I3)')
     c                'printout of CI-coefficients larger than',
     c                 prwthr,' for root', i
            Write(LF,'(6X,A,F15.6)')
     c                'energy=',ener(i,iter)

            call gasprwf(PrSel,nac,nactel,stsym,conf,
     c           cftp,CIVEC,iwork(ivkcnf))
          End If
         end if
          call getmem('kcnf','free','inte',ivkcnf,nactel)
*          END IF
        End Do

        ELSE !RUN SPLITCAS

          jDisk=iDisk
* load back one CI vector at the time
          Call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL SplitCAS last cycle',' ',
     &        CIVEC,nConf)
          END IF
* reorder it according to the split graph GUGA conventions
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
          Call Reord2(NAC,NACTEL,STSYM,0,
     &                CONF,CFTP,
     &                CIVEC,CIV,iWork(ivkcnf))
          call getmem('kcnf','free','inte',ivkcnf,nactel)
* save reorder CI vector on disk
          Call DDafile(JOBIPH,1,CIV,nConf,jDisk)
#ifdef _HDF5_
          call mh5_put_dset(wfn_cicoef,CIV(1:nConf),[nconf,1],[0,i-1])
#endif
          IF (IPRLEV.GE.DEBUG) THEN
           call DVcPrt('CI-Vec in CICTL after Reord',' ',CIV,nConf)
          END IF
* printout of the wave function
          IF (IPRLEV.GE.USUAL) THEN
            Write(LF,*)
            Write(LF,'(6X,A,F6.2,A,I3)')
     &                'printout of CI-coefficients larger than',
     &                 PRWTHR,' for root',lRootSplit
            Write(LF,'(6X,A,F15.6)')
     &           'Split-energy=',ENER(lRootSplit,ITER)
!     Open GronOR vecdet file (tps/cdg 20210430)
            write(filename,'(a7,i1)') 'VECDET.',i
!     filename = 'VECDET.'//merge(str(i),'x',i.lt.999)
            LuVecDet=39
            LuVecDet=IsFreeUnit(LuVecDet)
            call Molcas_open(LuVecDet,filename)
            write(LuVecDet,'(8i4)') nish
            CALL SGPRWF(PrSel,NOCSF,IOCSF,NOW1,IOW1,CIV)
!     Close GronOR vecdet file (tps/cdg 20210430)
            close(LuVecDet)
          END IF
        END IF
        endif

        Call mma_deallocate(PrSel)
        Call mma_deallocate(CIV)
      ENDIF

#ifdef _DMRG_
      if(doDMRG)then
        call mh5_put_dset
     &         (wfn_dmrg_checkpoint,dmrg_file%qcmaquis_checkpoint_file)
        call mma_deallocate(d1all)
        if(twordm_qcm) call mma_deallocate(d2all)
        if(doEntanglement) then
          if(allocated(spd1all)) call mma_deallocate(spd1all)
        end if
      end if
#endif

      Call mma_deallocate(CIVEC)
      Call mma_deallocate(FMO)

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
          call f_inquire('rf.results_state.h5', rfh5DMRG)
          if(.not.rfh5DMRG)then
            maquis_name_states  = ""
            maquis_name_results = ""
            call file_name_generator(IPCMROOT-1,"checkpoint_state.",
     &                               ".h5",maquis_name_states)
            call file_name_generator(IPCMROOT-1,"results_state.",
     &                               ".h5",maquis_name_results)

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

           rNorm = 1.0d0
           ! Shouldn't the overlap in this case be always 1?
           ! For DMRG it seems it is...
           ! But just to make sure we calculate it anyway
           ! in case of non-DMRG calculation
           if (.not.doDMRG) then
             Call Get_dArray("RF CASSCF Vector",Work(ipRF),nConf)
             rNorm=Sqrt(DDot_(nConf,Work(ipRF),1,Work(ipRF),1))
           end if
*          Write (6,*) 'rNorm=',rNorm
           JPCMROOT=IPCMROOT
           If (rNorm.gt.1.0D-10) Then
              Call Allocate_Work(ipTemp,nConf)
              rMax=0.0D0
              qMax=0.0d0
              jDisk = IADR15(4)
              Do i = 1, lRoots
                 if(doDMRG)then
#ifdef _DMRG_
                   qmax = abs(qcmaquis_interface_get_overlap(i))
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
     &                             ".h5",maquis_name_states)
          call file_name_generator(IPCMROOT-1,"results_state.",
     &                             ".h5",maquis_name_results)

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

      Return
      End

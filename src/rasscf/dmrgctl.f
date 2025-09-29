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
* Copyright (C) Naoki Nakatani                                         *
************************************************************************
*  DMRGCtl
*
*> @brief
*>   DMRG Control
*> @author N. Nakatani
*>
*> @details
*> Many functionality is same as ::CICtl (copy most part of ::CICtl to here)
*> To use DMRG as CAS-CI solver.
*> \p IRst = ``0`` is used for the first and the last (DMRG-) CI calculations
*> to fully optimize the DMRG wavefunction.
*> \p IRst = ``1`` enables to restart DMRG calculation from last iteration
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
*> @param[in]     IRst   DMRG restart status switch
************************************************************************

#include "compiler_features.h"

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      Subroutine DMRGCtl(CMO,D,DS,P,PA,FI,D1I,D1A,TUVX,IFINAL,IRst)

      use wadr, only: FMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use rctfld_module, only: lRF
      Use casvb_global, Only: ifvb
      use timers, only: TimeDens
      use lucia_data, only: PAtmp, Pscr, Ptmp, DStmp, Dtmp
      use gas_data, only: iDoGAS
      use Constants, only: Zero
      use rasscf_global, only: KSDFT, ExFac, iPCMRoot, ITER, lRoots,
     &                         n_Det, NAC, NACPAR, NACPR2, nFint,
     &                         nRoots, S, iAdr15, iRoot, Weight,
     &                         DFTFOCK
      use printlevel, only: DEBUG,INSANE
      use output_ras, only: LF,IPRLOC
      use general_data, only: ISPIN,jobiph,nactel,ntot2,nash

      Implicit None
      Integer iFinal, IRst
      Real*8 CMO(*),D(*),DS(*),P(*),PA(*),FI(*),D1I(*),D1A(*),
     &          TUVX(*)
c     Logical Exist
      Logical Do_ESPF
      Real*8 rdum(1)
      Real*8, Allocatable:: RCT_F(:), RCT_FS(:), RCT(:), RCT_S(:),
     &                      P2MO(:), TmpDS(:), TmpD1S(:),
     &                      TmpPUVX(:), TmpTUVX(:)

      Character(LEN=16), Parameter:: ROUTINE='DMRGCTL '
C Local print level (if any)
      Integer iPrLev, i, jDisk, jRoot, kRoot, NACT4, nTmpPUVX
      Real*8 dum1, dum2, dum3, Scal, Time(2)

      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

* set up flag 'IFCAS' for GAS option, which is set up in gugatcl originally.
* IFCAS = 0: This is a CAS calculation
* IFCAS = 1: This is a RAS calculation - This might cause an error in DMRG-CASSCF.
*
      if(iDoGas) call setsxci()
      If(IPRLEV.gt.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Enter DMRG section'
        Write(LF,*) ' =================='
        Write(LF,*)
        Write(LF,*) ' iteration count =',ITER
      End If
*
* SOME DIRTY SETUPS
*
      S=0.5D0*DBLE(ISPIN-1)
*
* COMPUTE ONE ELECTRON INTEGRALS IN MO BASIS
* AND ADD CORE INTERACTION
*
* FMO FOCK MATRIX IN MO-BASIS
* LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION
*
      CALL mma_allocate(FMO,NACPAR,Label='FMO')
      Call DecideOnESPF(Do_ESPF)
      If ( lRF .or. KSDFT.ne.'SCF' .or. Do_ESPF) THEN
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
        CALL mma_allocate(RCT_F,NTOT2,Label='RCT_F')
        CALL mma_allocate(RCT_FS,NTOT2,Label='RCT_FS')
        If (IFINAL.eq.0) Then
*
* Use normal MOs
*
           CALL mma_allocate(RCT,NACPAR,Label='RCT')
           CALL mma_allocate(P2MO,NACPR2,Label='P2MO')
*
* Get the total density in MOs
*
           Call DDafile(JOBIPH,2,RCT,NACPAR,jDisk)
           Call Put_dArray('D1mo',RCT,NACPAR)  ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(RCT)
* Transform to AOs
           Call Get_D1A_RASSCF(CMO,RCT,RCT_F)
*
* Get the spin density in MOs
*
           IF (NACTEL.EQ.0) THEN
             CALL DCOPY_(NTOT2,[0.0D0],0,RCT_FS,1)
           ELSE
             CALL mma_allocate(RCT_S,NACPAR,Label='RCT_S')
             Call DDafile(JOBIPH,2,RCT_S,NACPAR,jDisk)
             IF ( NASH(1).NE.NAC ) CALL DBLOCK(RCT_S)
* Transform to AOs
             Call Get_D1A_RASSCF(CMO,RCT_S,RCT_FS)
             CALL mma_deallocate(RCT_S)
           END IF
*
* Get the 2-particle density in MO
*
           Call DDafile(JOBIPH,2,P2MO,NACPR2,jDisk)
           Call Put_dArray('P2mo',P2MO,NACPR2) ! Put it on the RUNFILE
*
           CALL SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)
*
           CALL mma_deallocate(P2MO)
           CALL mma_deallocate(RCT)
*
        Else
*
* Here the pseudo-natural orbitals are in CMO and we need to
* get the D1A of the selected state in this basis.
*
*
* Compute the density of the particular state
*

           Call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
           Call mma_allocate(DStmp,NAC**2,Label='DStmp')
           Call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
           If ( NAC.ge.1 ) Then
              If (NACTEL.eq.0) THEN
                 Dtmp(:)=0.0D0
                 DStmp(:)=0.0D0
                 Ptmp(:)=0.0D0
              Else
* load back 1- and 2-RDMs from previous DMRG run
                 NACT4=NAC**4
                 Call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
                 CALL mma_allocate(Pscr,NACT4,Label='Pscr')
#ifdef _ENABLE_BLOCK_DMRG_
                 CALL block_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                   Ptmp,PAtmp,Pscr)
#elif _ENABLE_CHEMPS2_DMRG_
                 CALL chemps2_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                     Ptmp,PAtmp,Pscr)
#elif _ENABLE_DICE_SHCI_
                 CALL dice_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                  Ptmp,PAtmp,Pscr)
#endif

* NN.14 NOTE: IFCAS must be 0 for DMRG-CASSCF
c                If (IFCAS.GT.2) Call CISX(IDXSX,Dtmp,
c    &                                     DStmp,Ptmp,PAtmp,Pscr)
                 Call mma_deallocate(Pscr)
                 Call mma_deallocate(PAtmp)
              EndIf
*
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
           If (ExFac.ne.1.0D0) Call Mod_P2(Ptmp,NACPR2,
     &                                   Dtmp,NACPAR,
     &                                   DStmp,ExFac,n_Det)
*
           Call Put_dArray('P2mo',Ptmp,NACPR2) ! Put it on the RUNFILE
*
           Call mma_deallocate(Ptmp)
*
           Call Put_dArray('D1mo',Dtmp,NACPAR) ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Dtmp)
           Call Get_D1A_RASSCF(CMO,Dtmp,RCT_F)
*
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(DStmp)
           Call Get_D1A_RASSCF(CMO,DStmp,RCT_FS)
*
           Call mma_deallocate(Dtmp)
           Call mma_deallocate(DStmp)
*
           Call SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)
*
        End If
        CALL mma_deallocate(RCT_FS)
        CALL mma_deallocate(RCT_F)
*
      ELSE
*
* Normal case
*
*
*
*
        CALL mma_allocate(TmpDS ,NACPAR,Label='TmpDS')
        CALL mma_allocate(TmpD1S,NTOT2,Label='TmpD1S')
        call dcopy_(NACPAR,DS,1,TmpDS,1)
        IF ( NASH(1).NE.NAC ) CALL DBLOCK(TmpDS)
        Call Get_D1A_RASSCF(CMO,TmpDS,TmpD1S)
        CALL mma_deallocate(TmpDS)
*
        CALL SGFCIN(CMO,FMO,FI,D1I,D1A,TmpD1S)
        CALL mma_deallocate(TmpD1S)
*
      END IF
*
      If (IfVB.eq.2) GoTo 9000
*
* SOLVE DMRG WAVEFUNCTION
*
      if(IfVB.eq.1)then
* NN.14 FIXME: I'm not sure whether this option should work?
        call cvbmn_rvb(max(ifinal,1))
      else
        If (KSDFT(1:3).ne.'SCF'.
     &      and.DFTFOCK(1:4).eq.'DIFF'.and.nac.ne.0) Then
          nTmpPUVX=nFint
          Call mma_allocate(TmpPUVX,nTmpPUVX,Label='TmpPUVX')
          Call mma_allocate(TmpTUVX,NACPR2,Label='TmpTUVX')
          TmpTUVX(:)=0.0D0
          Call Get_dArray('DFT_TwoEl',TmpPUVX,nTmpPUVX)
          Call Get_TUVX(TmpPUVX,TmpTUVX)
          Call DaXpY_(NACPR2,1.0d0,TUVX,1,TmpTUVX,1)
#ifdef _ENABLE_BLOCK_DMRG_
          Call BlockCtl(FMO,TmpTUVX,IFINAL,IRst)
#elif _ENABLE_CHEMPS2_DMRG_
          Call Chemps2Ctl(FMO,TmpTUVX,IFINAL,IRst)
#elif _ENABLE_DICE_SHCI_
          Call DiceCtl(FMO,TmpTUVX,IFINAL,IRst)
#endif

          Call mma_deallocate(TmpTUVX)
          Call mma_deallocate(TmpPUVX)
        Else
#ifdef _ENABLE_BLOCK_DMRG_
          Call BlockCtl(FMO,TUVX,IFINAL,IRst)
#elif _ENABLE_CHEMPS2_DMRG_
          Call Chemps2Ctl(FMO,TUVX,IFINAL,IRst)
#elif _ENABLE_DICE_SHCI_
          Call DiceCtl(FMO,TUVX,IFINAL,IRst)
#endif
        End If
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
      Call Timing(Time(1),dum1,dum2,dum3)
      Call dCopy_(NACPAR,[Zero],0,D,1)
      Call dCopy_(NACPAR,[Zero],0,DS,1)
      Call dCopy_(NACPR2,[Zero],0,P,1)
      Call dCopy_(NACPR2,[Zero],0,PA,1)
      CALL mma_allocate(Dtmp,NAC**2,Label='Dtmp')
      CALL mma_allocate(DStmp,NAC**2,Label='DStmp')
      CALL mma_allocate(Ptmp,NACPR2,Label='Ptmp')
      CALL mma_allocate(PAtmp,NACPR2,Label='PAtmp')
      jDisk = IADR15(3)

      Do jRoot = 1,lRoots
* load density matrices from DMRG run
        If ( NAC.ge.1 ) Then
          NACT4=NAC**4
          CALL mma_allocate(Pscr,NACT4,Label='Pscr')
#ifdef _ENABLE_BLOCK_DMRG_
          CALL block_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#elif _ENABLE_CHEMPS2_DMRG_
          CALL chemps2_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#elif _ENABLE_DICE_SHCI_
          CALL dice_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#endif
          Call mma_deallocate(Pscr)
        EndIf
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
        If (ExFac.ne.1.0D0) Call Mod_P2(Ptmp,NACPR2,
     &                                Dtmp,NACPAR,
     &                                DStmp,ExFac,n_Det)

* update average density matrices
        Scal = 0.0d0
        Do kRoot = 1,nRoots
          If ( iRoot(kRoot).eq.jRoot ) then
            Scal = Weight(kRoot)
          End If
        End Do
        call daxpy_(NACPAR,Scal,Dtmp,1,D,1)
        call daxpy_(NACPAR,Scal,DStmp,1,DS,1)
        call daxpy_(NACPR2,Scal,Ptmp,1,P,1)
        call daxpy_(NACPR2,Scal,PAtmp,1,PA,1)
* save density matrices on disk
        Call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
        Call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
      End Do

      Call mma_deallocate(PAtmp)
      Call mma_deallocate(Ptmp)
      Call mma_deallocate(DStmp)
      Call mma_deallocate(Dtmp)
*
* C
* PREPARE DENSITY MATRICES AS USED BY THE SUPER CI SECTION
* C
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
      IF ( NASH(1).NE.NAC ) CALL DBLOCK(D)
      Call Timing(Time(2),dum1,dum2,dum3)
      TimeDens = TimeDens + Time(2) - Time(1)
*
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
* ======================================================================
* NN.14 FIXME:
*     The overlap b/w the DMRG wave of different macro iterations
*     hasn't yet been implemented. Eventually this must be fixed to
*     chose the correct root for RF calculation with DMRG-CASSCF.
*     For the time, just skip following.
* ======================================================================
*
c     If (lRF.and.KeyCISE.and.KeyCIRF) Then
c       JPCMROOT=IPCMROOT
c       IPCMROOT=IROOT(ICIRFROOT)
c       Call Put_iScalar("RF CASSCF root",IPCMROOT)
c       If (JPCMROOT.ne.IPCMROOT) Then
c         Write (6,'(1X,A,I3,A,I3)') 'RF Root has flipped from ',
c    &                 JPCMROOT, ' to ',IPCMROOT
c       End If
c     Else If (lRF) Then
c       Call Qpg_iScalar('RF CASSCF root',Exist)
c       If (.NOT.Exist) Then
*
*          We are here since we are using the default values.
*
c          Call Put_iScalar("RF CASSCF root",IPCMROOT)
c          Call Put_iScalar("RF0CASSCF root",IPCMROOT)
c       End If
*
c       Call mma_allocate(RF,nConf)
c       Call Qpg_dArray("RF CASSCF Vector",Exist,mConf)
*       Write (6,*) 'Exist=',Exist
c       If (Exist
c    &      .and. mConf .eq. nConf
c    &      .and. iFinal.ne.2
c    &      .and. (ABS(RotMax).lt.1.0D-3 .or. KeyCISE)
c    &     ) Then
c          Call Get_dArray("RF CASSCF Vector",RF,nConf)
c          rNorm=Sqrt(DDot_(nConf,RF,1,RF,1))
*          Write (6,*) 'rNorm=',rNorm
c          JPCMROOT=IPCMROOT
c          If (rNorm.gt.1.0D-10) Then
c             Call mma_allocate(Temp,nConf,Label='Temp')
c             rMax=0.0D0
c             jDisk = IADR15(4)
c             Do i = 1, lRoots
c                Call DDafile(JOBIPH,2,Temp,nConf,jDisk)
c                qMax=Abs(DDot_(nConf,Temp,1,RF,1))
*                Write (6,*) 'qMax=',qMax
c                If (qMax.gt.rMax .and.
c    &               qMax.gt.0.5D0) Then
c                   rMax = qMax
c                   JPCMROOT=i
c                End If
c             End Do
c             Call mma_deallocate(Temp)
c          End If
c       Else
c          JPCMROOT=IPCMROOT
c       End If
*
c       If (JPCMROOT.ne.IPCMROOT) Then
c          Write (6,*) ' RF Root has flipped from ',IPCMROOT, ' to ',
c    &                                              JPCMROOT
c          IPCMROOT=JPCMROOT
c          Call Put_iScalar("RF CASSCF root",IPCMROOT)
c       End If
*
c       jDisk = IADR15(4)
c       Do i=1,IPCMROOT-1
c         Call DDafile(JOBIPH,0,rdum,nConf,jDisk)
c       End Do
c       Call DDafile(JOBIPH,2,RF,nConf,jDisk)
c       Call Put_dArray("RF CASSCF Vector",RF,nConf)
c       Call mma_deallocate(RF)
c     End If
*
      End Subroutine DMRGCtl

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#     include "macros.fh"
      subroutine empty_DMRGCtl()
      end subroutine empty_DMRGCtl

#endif

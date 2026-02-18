!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Bjorn O. Roos                                          *
!               Per Ake Malmqvist                                      *
!               1991, Jeppe Olsen                                      *
!               1991,1996, Markus P. Fuelscher                         *
!***********************************************************************
!  CICtl
!
!> @brief
!>   CI Control
!> @author B. O. Roos
!> @author P. &Aring;. Malmqvist
!> @modified_by P. &Aring;. Malmqvist
!>
!> @details
!> Depends on \p IFINAL, which is set in ::RASSCF. If \p IFINAL = ``0``, repeated
!> calculations with orbital optimization before each call. If \p IFINAL = ``1``,
!> there has been no orbital optimization, or the calculation is
!> converged. \p IFINAL = ``2`` means this is a final CI calculation, using the
!> final orbitals.
!>
!> @param[in]     CMO    MO coefficients
!> @param[out]    D      Average 1-dens matrix
!> @param[out]    DS     Average spin 1-dens matrix
!> @param[out]    P      Average symm. 2-dens matrix
!> @param[out]    PA     Average antisymm. 2-dens matrix
!> @param[out]    FI     Fock matrix from inactive density
!> @param         FA
!> @param[in,out] D1I    Inactive 1-dens matrix
!> @param[in,out] D1A    Active 1-dens matrix
!> @param[in]     TUVX   Active 2-el integrals
!> @param[in]     IFINAL Calculation status switch
!***********************************************************************

subroutine CICtl(CMO,D,DS,P,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
! ****************************************************************
! history:                                                       *
! updated to use determinant based CI-procedures                 *
! J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
! updated for MOLCAS version 3                                   *
! J. Olsen and M.P. Fuelscher, University of Lund, Sweden, 1991  *
! updated for integral direct and reaction field calculations    *
! M.P. Fuelscher, University of Lund, Sweden, 1996               *
! ****************************************************************

#ifdef _DMRG_
use qcmaquis_interface, only: qcmaquis_interface_get_overlap, qcmaquis_interface_run_dmrg, qcmaquis_interface_run_starting_guess, &
                              qcmaquis_interface_set_param, qcmaquis_interface_update_integrals
use qcmaquis_interface_cfg, only: dmrg_energy, dmrg_file, dmrg_orbital_space, dmrg_warmup, qcmaquis_param
use qcmaquis_interface_utility_routines, only: fiedlerorder_length, file_name_generator, qcmaquis_interface_fcidump
use lucia_data, only: RF1, RF2
use RASWfn, only: wfn_dmrg_checkpoint
use rasscf_global, only: DOFCIDump, Emy, TwoRDM_qcm
#endif
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use RASWfn, only: wfn_cicoef, wfn_dens, wfn_spindens
#endif
use csfbas, only: CONF
use casvb_global, only: ifvb
use CMS, only: CMSGiveOpt, iCMSOpt
use rctfld_module, only: lRF
use lucia_data, only: CFTP, DStmp, Dtmp, PAtmp, Pscr, PTmp
use Lucia_Interface, only: Lucia_Util
use wadr, only: FMO
use gugx, only: CIS, SGS
use sxci, only: IDXSX
use gas_data, only: iDoGAS
use input_ras, only: Key
use timers, only: TimeDens
use rasscf_global, only: CMSStartMat, DoDMRG, Ener, ExFac, IADR15, iCIRFRoot, ICMSP, IFCRPR, iPCMRoot, iRoot, iRotPsi, ITER, &
                         IXMSP, KSDFT, l_casdft, lroots, n_Det, NAC, NACPAR, NACPR2, nRoots, PrwThr, RotMax, S, Weight
use SplitCas_Data, only: DoSPlitCas, lRootSplit, MxIterSplit, ThrSplit
use PrintLevel, only: DEBUG, INSANE, USUAL
use output_ras, only: IPRLOC
use general_data, only: CRVec, ISPIN, JOBIPH, NACTEL, NASH, NCONF, NISH, NTOT2, STSYM
use DWSol, only: DWSolv
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*), D(*), DS(*), P(*), PA(*), FI(*), FA(*), D1I(*), D1A(*), TUVX(*)
integer(kind=iwp) :: iFinal
integer(kind=iwp) :: i, iDisk, iErrSplit, iOpt, iPrLev, jDisk, jPCMRoot, jRoot, kRoot, LuVecDet, mconf
real(kind=wp) :: dum1, dum2, dum3, qMax, rdum(1), rMax, rNorm, Scal, Time(2)
logical(kind=iwp) :: Do_ESPF, do_rotate, Exists, Skip
character(len=128) :: filename
real(kind=wp), allocatable :: CIV(:), CIVec(:), P2MO(:), RCT(:), RCT_F(:), RCT_FS(:), RCT_S(:), RF(:), Temp(:), TmpD1S(:), TmpDS(:)
integer, allocatable :: kCnf(:)
#ifdef _HDF5_
real(kind=wp), allocatable :: density_square(:,:)
#endif
#ifdef _DMRG_
integer(kind=iwp) :: iErr, ilen
logical(kind=iwp) :: doEntanglement, rfh5DMRG
character(len=2300) :: maquis_name_results, maquis_name_states
character(len=:), allocatable :: fiedler_order_str
real(kind=wp), allocatable :: d1all(:,:), d2all(:,:), spd1all(:,:)
logical(kind=iwp), external :: PCM_On
#endif
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: DDot_
#include "warnings.h"

! Local print level (if any)
IPRLEV = IPRLOC(3)
if (IPRLEV >= DEBUG) write(u6,*) ' Entering CICTL'

! PAM 2017-05-23 Modify TUVX by adding a shift vector to TUVX, which has
! the effect of adding a scalar times a projector for doubly-occupied
! core states.
if (IfCRPR) call MKPROJ(CRVEC,CMO,TUVX)

! set up flag 'IFRAS' for GAS option, which is set up in gugatcl originally.
! IFRAS = 0: This is a CAS calculation
! IFRAS = 1: This is a RAS calculation

if (iDoGas .or. (SGS%IFRAS > 2)) call setsxci()
if (IPRLEV > DEBUG) then
  write(u6,*)
  write(u6,*) ' Enter CI section, CICTL routine'
  write(u6,*) ' ================'
  write(u6,*)
  write(u6,*) ' iteration count =',ITER
end if

!do i=1,NTOT2 ! yma
!  write(u6,*) 'ifinal CMO',ifinal,i,CMO(i)
!end do

! SOME DIRTY SETUPS

S = Half*real(ISPIN-1,kind=wp)

! COMPUTE ONE ELECTRON INTEGRALS IN MO BASIS AND ADD CORE INTERACTION

! FMO: FOCK MATRIX IN MO-BASIS
! LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION

call mma_allocate(FMO,NACPAR,Label='FMO')
call DecideOnESPF(Do_ESPF)

! initialize RDM arrays for QCMaquis
#ifdef _DMRG_
if (doDMRG) then
  ! Calculate entanglement and spin densities
  ! only in the last iteration
  ! i.e. only for IFINAL == 2 if we do DMRG-SCF
  ! otherwise always
  doEntanglement = merge(.true.,IFINAL == 2,Key('CION'))

  call mma_allocate(d1all,NACPAR,lRoots)
  d1all(:) = Zero
  if (twordm_qcm) then
    call mma_allocate(d2all,NACPR2,lRoots)
    d2all(:) = Zero
  end if
  ! Allocate spin density only for the last iteration
  if (doEntanglement) then
    call mma_allocate(spd1all,NACPAR,lRoots)
    spd1all(:) = Zero
  end if
end if
#endif

if ((lRf .or. (KSDFT /= 'SCF') .or. Do_ESPF) .and. IPCMROOT > 0) then

  ! In case of a reaction field in combination with an average CAS
  ! select the potential of the appropriate state.

  jDisk = IADR15(3)
  do i=1,IPCMROOT-1
    call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
    call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
  end do

  call mma_allocate(RCT_F,NTOT2,Label='RCT_F')
  call mma_allocate(RCT_FS,NTOT2,Label='RCT_FS')

  if ((IPCMROOT > 0) .and. (DWSolv%DWZeta /= Zero)) then
    call DWDens_RASSCF(CMO,D1A,RCT_FS,IFINAL)
    call SGFCIN(CMO,FMO,FI,D1I,D1A,RCT_FS)
  else if (IFinal == 0) then

    ! Use normal MOs

    call mma_allocate(RCT,NACPAR,Label='RCT')
    call mma_allocate(P2MO,NACPR2,Label='P2MO')

    ! Get the total density in MOs

    call DDafile(JOBIPH,2,RCT,NACPAR,jDisk)
    call Put_dArray('D1mo',RCT,NACPAR)  ! Put on RUNFILE
    if (NASH(1) /= NAC) call DBLOCK(RCT)
    ! Transform to AOs
    call Get_D1A_RASSCF(CMO,RCT,RCT_F)

    ! Get the spin density in MOs

    if (NACTEL == 0) then
      RCT_FS(:) = Zero
    else
      call mma_allocate(RCT_S,NACPAR,Label='RCT_S')
      call DDafile(JOBIPH,2,RCT_S,NACPAR,jDisk)
      if (NASH(1) /= NAC) call DBLOCK(RCT_S)
      ! Transform to AOs
      call Get_D1A_RASSCF(CMO,RCT_S,RCT_FS)
      call mma_deallocate(RCT_S)
    end if

    ! Get the 2-particle density in MO

    call DDafile(JOBIPH,2,P2MO,NACPR2,jDisk)
    call Put_dArray('P2mo',P2MO,NACPR2) ! Put on RUNFILE

    call SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)

    call mma_deallocate(P2MO)
    call mma_deallocate(RCT)

  else

    ! Here the pseudo-natural orbitals are in CMO and we need to
    ! get the D1A of the selected state in this basis.

    ! Compute the density of the particular state

    call mma_allocate(CIVEC,NCONF,Label='CIVEC')
    if (NACTEL == 0) then
      CIVEC(1) = One
    else
      if (.not. doDMRG) then
        !write(u6,*) 'run the load back CI vector part' ! yma
        iDisk = IADR15(4)
        do jRoot=1,IPCMROOT
          iOpt = 0
          if (jRoot == IPCMROOT) iOpt = 2
          ! load back one CI vector at the time
          call DDafile(JOBIPH,iOpt,CIVEC,nConf,iDisk)
          !call DVcPrt('BLUBB-start CI PCM',' ',CIVEC,nConf)
        end do
      end if
    end if

    ! compute density matrices

    call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
    call mma_allocate(DStmp,NAC**2,Label='DStmp')
    call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
    if (NAC >= 1) then

      if (NACTEL == 0) then
        Dtmp(:) = Zero
        DStmp(:) = Zero
        Ptmp(:) = Zero
      else

        if (doDMRG) then
#         ifdef _DMRG_
          ! copy the DMs from d1rf/d2rf for ipcmroot
          Dtmp(1:NACPAR) = rf1(1:NACPAR)
          if (twordm_qcm) Ptmp(1:NACPR2) = rf2(1:NACPR2)

          ! Import RDMs from QCMaquis that we've got from the last optimization
          ! Here we should import one-particle spin density.
          ! However, the spin density has been temporarily disabled here:
          ! For performance reasons, it is calculated
          ! only once in the last iteration of DMRG-SCF optimisation.
          ! If you need it at every iteration for some reason
          ! please change this code accordingly
          DStmp(:) = Zero
#         endif
        else
          call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
          call mma_allocate(Pscr,NACPR2,Label='Pscr')
          call Lucia_Util('Densi',CI_Vector=CIVEC)
          if ((SGS%IFRAS > 2) .or. iDoGAS) call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
          call mma_deallocate(Pscr)
          call mma_deallocate(PAtmp)
        end if ! doDMRG/doBLOK or CI

      end if
    else
      Dtmp(:) = Zero
      DStmp(:) = Zero
      Ptmp(:) = Zero
    end if
    ! Modify the symmetric 2-particle density if only partial
    ! "exact exchange" is included.
    !n_Det = 2
    !n_unpaired_elec = iSpin-1
    !n_paired_elec = nActEl-n_unpaired_elec
    !if (n_unpaired_elec+n_paired_elec/2 == nac) n_Det = 1
    !write(u6,*) 'n_Det =',n_Det
    if ((ExFac /= One) .and. (.not. l_casdft)) call Mod_P2(Ptmp,NACPR2,Dtmp,NACPAR,DStmp,ExFac,n_Det)

    call Put_dArray('P2mo',Ptmp,NACPR2) ! Put on RUNFILE

    call mma_deallocate(Ptmp)

    call Put_dArray('D1mo',Dtmp,NACPAR) ! Put on RUNFILE
    if (NASH(1) /= NAC) call DBLOCK(Dtmp)
    call Get_D1A_RASSCF(CMO,Dtmp,RCT_F)

    if (NASH(1) /= NAC) call DBLOCK(DStmp)
    call Get_D1A_RASSCF(CMO,DStmp,RCT_FS)

    !do i=1,NACPAR ! yma
    !  write(u6,*) 'i-rdms1 1',i,Dtmp(i)
    !end do

    call mma_deallocate(DStmp)
    call mma_deallocate(Dtmp)
    call mma_deallocate(CIVEC)

    call SGFCIN(CMO,FMO,FI,D1I,RCT_F,RCT_FS)

  end if
  call mma_deallocate(RCT_FS)
  call mma_deallocate(RCT_F)

else
  !**************************************************************************************
  ! Normal case: no reaction field, no ESP to be added to the Fock matrix or ksdft /= SCF
  !**************************************************************************************
  if (IPRLEV >= DEBUG) then
    call TRIPRT('DS input  ',' ',DS,NAC)
    call TRIPRT('D input  ',' ',D,NAC)
  end if
  call mma_allocate(TmpDS,NACPAR,Label='TmpDS')
  call mma_allocate(TmpD1S,NTOT2,Label='TmpD1S')
  TmpDS(:) = DS(1:NACPAR)
  if (NASH(1) /= NAC) call DBLOCK(TmpDS)
  call Get_D1A_RASSCF(CMO,TmpDS,TmpD1S)

  call mma_deallocate(TmpDS)

  call SGFCIN(CMO,FMO,FI,D1I,D1A,TmpD1S)
  call mma_deallocate(TmpD1S)

end if

Skip = .false.
#ifdef _DMRG_
if (doDMRG) then
  ! update integrals for QCMaquis

  call qcmaquis_interface_update_integrals(FMO,tuvx,emy)

  !!! Fiedler order/CI-DEAS run
  if (dmrg_warmup%dofiedler .or. dmrg_warmup%docideas) then
    ! allocate string where Fiedler ordering will be returned with the correct length
    ilen = fiedlerorder_length(qcmaquis_param%L)
    allocate(character(len=ilen) :: fiedler_order_str)
    fiedler_order_str(:) = ' '

    ! if HF guess is present, use it (required for CI-DEAS)
    ! If not, error handling should be done in QCMaquis
    ! pass the HF occupation as 1D array to QCMaquis
    if (sum(dmrg_orbital_space%initial_occ) > 0) then
      call qcmaquis_interface_run_starting_guess(nRoots,dmrg_warmup%dofiedler,dmrg_warmup%docideas,fiedler_order_str, &
                                                 reshape(dmrg_orbital_space%initial_occ,[1,sum(nash)*nroots]))
    else
      call qcmaquis_interface_run_starting_guess(nRoots,dmrg_warmup%dofiedler,dmrg_warmup%docideas,fiedler_order_str)
    end if

    if (dmrg_warmup%dofiedler) call qcmaquis_interface_set_param('orbital_order',fiedler_order_str)
    write(u6,*) 'Fiedler orbital ordering: '//fiedler_order_str

    dmrg_warmup%dofiedler = .false.
    dmrg_warmup%docideas = .false.
    if (allocated(fiedler_order_str)) deallocate(fiedler_order_str)
  end if
  if (dofcidump) then
    ! Produce a FCIDUMP file
    ! TODO:
    ! We already have the fcidump module in rasscf, which is called elsewhere
    ! so ensure the compatibility of the FCIDUMP files produced by this module
    ! and remove the code below
    call qcmaquis_interface_fcidump(FMO,tuvx,emy)
    call mma_deallocate(FMO)
    Skip = .true.
  end if
end if
#endif

if ((.not. Skip) .and. (IfVB /= 2)) then

  ! DAVIDSON DIAGONALIZATION

  if (IfVB == 1) then
    call cvbmn_rvb(max(ifinal,1))
  else
    if (DoSplitCAS) then !(GLMJ)
      call SplitCtl(FMO,TUVX,IFINAL,iErrSplit)
      if (iErrSplit == 1) then
        write(u6,*) repeat('*',120)
        write(u6,*) 'WARNING!!!'
        write(u6,*) 'SplitCAS iterations don''t converge.'
        write(u6,*) 'The program will continue'
        write(u6,*) 'Hopefully your calculation will converge next iteration!'
        write(u6,*) repeat('*',120)
      end if
      if (iErrSplit == 2) then
        write(u6,*) repeat('*',120)
        write(u6,*) 'WARNING!!!'
        write(u6,*) 'SplitCAS iterations don''t converge.'
        write(u6,*) 'MxIterSplit',MxIterSplit
        write(u6,*) 'SplitCAS ThreShold',ThrSplit
        write(u6,*) 'Try to increase MxIterSplit or SplitCAS threshold'
        write(u6,*) 'The program will STOP'
        write(u6,*) repeat('*',120)
        call xQuit(_RC_NOT_CONVERGED_)
      end if
    end if
    if (.not. DoSplitCAS) then
      if (doDMRG) then
#       ifdef _DMRG_
        ! Get also spin density at the last iteration
        ! Please help me call it more cleanly than with these if clauses
        ! and different optional arguments
        if (doEntanglement) then
          if (twordm_qcm) then
            call qcmaquis_interface_run_dmrg(nstates=lroots,d1=d1all,d2=d2all,spd=spd1all,entanglement=doEntanglement)
          else
            call qcmaquis_interface_run_dmrg(nstates=lroots,d1=d1all,spd=spd1all,entanglement=doEntanglement)
          end if
        else
          if (twordm_qcm) then
            call qcmaquis_interface_run_dmrg(nstates=lroots,d1=d1all,d2=d2all,entanglement=doEntanglement)
          else
            call qcmaquis_interface_run_dmrg(nstates=lroots,d1=d1all,entanglement=doEntanglement)
          end if
        end if

        ! For PCM calculations: copy RDMs for the PCM root
        if (PCM_On()) then
          rf1(1:NACPAR) = d1all(:,ipcmroot)
          if (twordm_qcm) rf2(1:NACPR2) = d2all(:,ipcmroot)
        end if
        ! Keep the root energies
        do jRoot=1,lRoots
          ENER(jRoot,ITER) = dmrg_energy%dmrg_state_specific(jroot)
        end do
        ! The new QCMaquis interface requires that the density matrices are calculated immediately after the DMRG run
        ! So we either need to keep them all in memory, or move the saving routines up here.
        ! The 2nd option requires code refactoring, so for now we keep them all in memory.
#       endif
      else
        ! Normal Davidson algorithm
        call DavCtl(FMO,TUVX,IFINAL)
      end if
    end if
  end if

  ! CALCULATE DENSITY MATRICES
  ! SAVE DENSITY MATRICES ON FILE
  ! COMPUTE AVERAGE DENSITY MATRICES

  ! Dtmp: ONE-BODY DENSITY
  ! DStmp: ONE-BODY SPIN DENSITY
  ! Ptmp: SYMMETRIC TWO-BODY DENSITY
  ! PAtmp: ANTISYMMETRIC TWO-BODY DENSITY

  call Timing(Time(1),dum1,dum2,dum3)
  D(1:NACPAR) = Zero
  DS(1:NACPAR) = Zero
  P(1:NACPR2) = Zero
  PA(1:NACPR2) = Zero
  call mma_allocate(CIVEC,NCONF,Label='CIVEC')
  call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
  call mma_allocate(DStmp,NAC**2,Label='DStmp')
  call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
  call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
  call mma_allocate(Pscr,NACPR2,Label='Pscr')
# ifdef _HDF5_
  call mma_allocate(density_square,nac,nac)
# endif

  !if (DWSCF%do_DW) call DWSol_wgt(1,ENER(:,ITER),weight)
  iDisk = IADR15(4)
  jDisk = IADR15(3)
  if (.not. DoSplitCAS) then
    !JB Instead of RASSCF/RASCI energy, print out energy for rotated states
    do_rotate = .false.
    if (ifinal == 2) then
      if (IXMSP == 1) call XMSRot(CMO,FI,FA)
      if (ICMSP == 1) then
        if (trim(CMSStartMat) == 'XMS') call XMSRot(CMO,FI,FA)
        if (.not. CMSGiveOpt) then
          if (lRoots == 2) iCMSOpt = 2
          if (lRoots >= 3) iCMSOpt = 1
        end if
        if (iCMSOpt == 1) then
          call CMSOpt(TUVX)
        else if (iCMSOpt == 2) then
          call CMSRot(TUVX)
        end if
      end if
      if (IRotPsi == 1) call f_inquire('ROT_VEC',Do_Rotate)
      if (Do_Rotate) then
        call RotState()
      else
        if (IRotPsi == 1) write(u6,'(6X,A)') 'Do_Rotate.txt is not found. MCSCF states will not be rotated'
      end if
      !JB End of condition 'Do_Rotate' to initialize rotated states
    end if
    !JB End If for ifinal=2
    do jRoot=1,lRoots
      ! load back one CI vector at the time
      !JB If do_rotate=.true., then we read CI vectors from CIVec
      !JB Otherwise we read if from JOBIPH
      call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
      if (IPRLEV >= DEBUG) call DVcPrt('CI-Vec in CICTL',' ',CIVEC,nConf)
      ! compute density matrices

      if (NAC >= 1) then
        if (.not. doDMRG) call Lucia_Util('Densi',CI_Vector=CIVEC)
        if (IPRLEV >= INSANE) then
          write(u6,*) 'At root number =',jroot
          call TRIPRT('D after lucia  ',' ',Dtmp,NAC)
          call TRIPRT('DS after lucia  ',' ',DStmp,NAC)
          call TRIPRT('P after lucia',' ',Ptmp,NACPAR)
          call TRIPRT('PA after lucia',' ',PAtmp,NACPAR)
        end if
      end if
      if ((.not. doDMRG) .and. ((SGS%IFRAS > 2) .or. iDoGAS)) call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
      ! 1,2-RDMs importing from DMRG calculation -- Stefan/Yingjin
      if (doDMRG) then
#       ifdef _DMRG_
        ! for QCMaquis, just copy the RDMs
        ! actually, copying is not needed! TODO
        Dtmp(1:NACPAR) = d1all(:,jroot)
        if (twordm_qcm) Ptmp(1:NACPR2) = d2all(:,jroot)

        !> import 1p-spin density
        ! disable spin density if not in the last iteration
        if (doEntanglement) then
          DStmp(1:NACPAR) = spd1all(:,jroot)
        else
          DStmp(:) = Zero
        end if

        ! disable antisymmetric 2-RDM
        PAtmp(:) = Zero

        if (IPRLEV >= INSANE) then
          call TRIPRT('D after  DMRG',' ',Dtmp,NAC)
          call TRIPRT('DS after DMRG',' ',DStmp,NAC)
          call TRIPRT('P after  DMRG',' ',Ptmp,NACPAR)
          call TRIPRT('PA after DMRG',' ',PAtmp,NACPAR)
        end if
#     endif
      end if
      ! Modify the symmetric 2-particle density if only partial "exact exchange" is included.
      !n_Det = 2
      !n_unpaired_elec = iSpin-1
      !n_paired_elec = nActEl-n_unpaired_elec
      !if (n_unpaired_elec+n_paired_elec/2 == nac) n_Det = 1
      !  write(u6,*) ' iSpin=',iSpin
      !  write(u6,*) ' n_unpaired_elec',n_unpaired_elec
      !  write(u6,*) ' n_paired_elec', n_paired_elec
      !  write(u6,*) ' n_unpaired_elec+n_paired_elec/2',n_unpaired_elec+n_paired_elec/2
      !  write(u6,*) ' n_Det=',n_Det
      !end if

      !write(u6,*) 'second call to Mod_P2'

      if ((ExFac /= One) .and. (.not. l_casdft)) call Mod_P2(Ptmp,NACPR2,Dtmp,NACPAR,DStmp,ExFac,n_Det)

      ! update average density matrices
      Scal = Zero
      do kRoot=nRoots,1,-1
        if (iRoot(kRoot) == jRoot) then
          Scal = Weight(kRoot)
          exit
        end if
      end do
      D(1:NACPAR) = D(1:NACPAR)+Scal*Dtmp(1:NACPAR)
      DS(1:NACPAR) = DS(1:NACPAR)+Scal*DStmp(1:NACPAR)
      P(1:NACPR2) = P(1:NACPR2)+Scal*Ptmp(1:NACPR2)
      PA(1:NACPR2) = PA(1:NACPR2)+Scal*PAtmp(1:NACPR2)
      !GLM Put the D1MO and the P2MO values in RUNFILE

      call Put_dArray('D1mo',Dtmp,NACPAR) ! Put on RUNFILE
      call Put_dArray('P2mo',Ptmp,NACPR2) ! Put on RUNFILE
      ! save density matrices on disk
      call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
      call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
      call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
      call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
#     ifdef _HDF5_
      !SVC: store a single column instead of the whole array (which is for each root!)
      ! and for now don't bother with 2-electron active density matrices
      call square(Dtmp,density_square,1,nac,nac)
      call mh5_put_dset(wfn_dens,density_square,[nac,nac,1],[0,0,jRoot-1])
      call square(DStmp,density_square,1,nac,nac)
      call mh5_put_dset(wfn_spindens,density_square,[nac,nac,1],[0,0,jRoot-1])
#     endif
    end do

  else  ! SplitCAS run
    call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
    if (IPRLEV >= DEBUG) call DVcPrt('CI-Vec in CICTL SplitCAS sect',' ',CIVEC,nConf)
    ! compute density matrices
    if (NAC >= 1) then
      call Lucia_Util('Densi',CI_Vector=CIVEC)
      if (IPRLEV >= INSANE) then
        call TRIPRT('D after lucia',' ',Dtmp,NAC)
        call TRIPRT('DS after lucia',' ',DStmp,NAC)
        call TRIPRT('P after lucia',' ',Ptmp,NACPAR)
        call TRIPRT('PA after lucia',' ',PAtmp,NACPAR)
      end if
    end if
    if (IDoGAS .or. (SGS%IFRAS > 2)) call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
    if ((ExFac /= One) .and. (.not. l_casdft)) call Mod_P2(Ptmp,NACPR2,Dtmp,NACPAR,DStmp,ExFac,n_Det)
    D(1:NACPAR) = D(1:NACPAR)+Dtmp(1:NACPAR)
    DS(1:NACPAR) = DS(1:NACPAR)+DStmp(1:NACPAR)
    P(1:NACPR2) = P(1:NACPR2)+Ptmp(1:NACPR2)
    PA(1:NACPR2) = PA(1:NACPR2)+PAtmp(1:NACPR2)
    ! save density matrices on disk
    call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
    call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
    call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
    call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
  end if

# ifdef _HDF5_
  call mma_deallocate(density_square)
# endif
  call mma_deallocate(Pscr)
  call mma_deallocate(PAtmp)
  call mma_deallocate(Ptmp)
  call mma_deallocate(DStmp)
  call mma_deallocate(Dtmp)

  ! print matrices
  if (IPRLEV >= INSANE) then
    call TRIPRT('Averaged one-body density matrix, D',' ',D,NAC)
    call TRIPRT('Averaged one-body spin density matrix, DS',' ',DS,NAC)
    call TRIPRT('Averaged two-body density matrix, P',' ',P,NACPAR)
    call TRIPRT('Averaged antisymmetric two-body density matrix,PA',' ',PA,NACPAR)
  end if
  call Put_dArray('D1mo',D,NACPAR) ! Put on RUNFILE
  if (lRf .and. (IPCMROOT <= 0)) call Put_dArray('P2mo',P,NACPR2) ! Put on RUNFILE

  if (NASH(1) /= NAC) call DBLOCK(D)
  call Timing(Time(2),dum1,dum2,dum3)
  TimeDens = TimeDens+Time(2)-Time(1)

  ! IF FINAL ITERATION REORDER THE WAVEFUNCTION ACCORDING TO
  ! THE SPLIT GRAPH GUGA CONVENTIONS AND PRINT IT.

  ! CIV: Temporary copy of a CI vector

  if ((IFINAL == 2) .and. (NAC > 0)) then
    if (IPRLEV >= USUAL) then
      write(u6,*)
      write(u6,'(6X,A)') repeat('*',120)
      write(u6,'(54X,A)') 'Wave function printout:'
      write(u6,'(23X,A)') 'occupation of active orbitals, and spin coupling of open shells (u,d: Spin up or down)'
      write(u6,'(6X,A)') repeat('*',120)
      write(u6,*)
      write(u6,'(6x,A)') 'Note: transformation to natural orbitals'
      write(u6,'(6x,A)') 'has been made, which may change the order of the CSFs.'
    end if
    call mma_allocate(CIV,nConf,Label='CIV')
    iDisk = IADR15(4)

    if (.not. doDMRG) then
      if (.not. DoSplitCAS) then
        do i=1,lRoots
          jDisk = iDisk
          ! load back one CI vector at the time
          call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
          if (IPRLEV >= DEBUG) call DVcPrt('CI-Vec in CICTL last cycle',' ',CIVEC,nConf)
          call mma_allocate(kcnf,nactel,Label='kCnf')
          if (.not. iDoGas) then
            call Reord2(NAC,NACTEL,STSYM,0,CONF,CFTP,CIVEC,CIV,kcnf)
            !end if
            !call mma_deallocate(kcnf)

            ! save reorder CI vector on disk
            !if (.not. iDoGas) then
            call DDafile(JOBIPH,1,CIV,nConf,jDisk)
#           ifdef _HDF5_
            call mh5_put_dset(wfn_cicoef,CIV(1:nConf),[nconf,1],[0,i-1])
#           endif
            !else
            !  call DDafile(JOBIPH,1,CIVEC,nConf,jDisk)
            !end if
            ! printout of the wave function
            if (IPRLEV >= USUAL) then
              write(u6,*)
              write(u6,'(6X,A,F6.2,A,I3)') 'printout of CI-coefficients larger than',PRWTHR,' for root',i
              write(u6,'(6X,A,F15.6)') 'energy=',ENER(I,ITER)
              if (Key('PRSD')) then
                ! Define filename to write GronOR vecdet files (tps/cdg 20210430)
                write(filename,'(a7,i0)') 'VECDET.',i
                !filename = 'VECDET.'//merge(str(i),'x',i <= 999)
                LuVecDet = IsFreeUnit(39)
                call Molcas_open(LuVecDet,filename)
                write(LuVecDet,'(8i4)') nish
              end if
              call SGPRWF(SGS,CIS,STSYM,PRWTHR,iSpin,CIV,nConf,Key('PRSD'),LUVECDET)
              ! Close GronOR vecdet file (tps/cdg 20210430)
              if (Key('PRSD')) close(LuVecDet)
            end if
          else ! for iDoGas
            write(u6,'(1x,a)') 'WARNING: true GAS, JOBIPH not compatible!'
            !.. save CI vector on disk
            call DDafile(JOBIPH,1,CIVEC,nconf,jDisk)
#           ifdef _HDF5_
            !SVC: store CI as a column array of the on-disk CI (which is for all roots!)
            call mh5_put_dset(wfn_cicoef,CIVEC(1:nconf),[nconf,1],[0,i-1])
#           endif
            !.. printout of the wave function
            if (IPRLEV >= USUAL) then
              write(u6,*)
              write(u6,'(6X,A,F6.2,A,I3)') 'printout of CI-coefficients larger than',prwthr,' for root',i
              write(u6,'(6X,A,F15.6)') 'energy=',ener(i,iter)

              call gasprwf(nac,nactel,stsym,conf,cftp,CIVEC,kcnf)
            end if
          end if
          call mma_deallocate(kcnf)
          !end if
        end do

      else !RUN SPLITCAS

        jDisk = iDisk
        ! load back one CI vector at the time
        call DDafile(JOBIPH,2,CIVEC,nConf,iDisk)
        if (IPRLEV >= DEBUG) call DVcPrt('CI-Vec in CICTL SplitCAS last cycle',' ',CIVEC,nConf)
        ! reorder it according to the split graph GUGA conventions
        call mma_allocate(kcnf,nactel,Label='kcnf')
        call Reord2(NAC,NACTEL,STSYM,0,CONF,CFTP,CIVEC,CIV,kcnf)
        call mma_deallocate(kcnf)
        ! save reorder CI vector on disk
        call DDafile(JOBIPH,1,CIV,nConf,jDisk)
#       ifdef _HDF5_
        call mh5_put_dset(wfn_cicoef,CIV(1:nConf),[nconf,1],[0,i-1])
#       endif
        if (IPRLEV >= DEBUG) call DVcPrt('CI-Vec in CICTL after Reord',' ',CIV,nConf)
        ! printout of the wave function
        if (IPRLEV >= USUAL) then
          write(u6,*)
          write(u6,'(6X,A,F6.2,A,I3)') 'printout of CI-coefficients larger than',PRWTHR,' for root',lRootSplit
          write(u6,'(6X,A,F15.6)') 'Split-energy=',ENER(lRootSplit,ITER)
          ! Open GronOR vecdet file (tps/cdg 20210430)
          write(filename,'(a7,i0)') 'VECDET.',i
          !filename = 'VECDET.'//merge(str(i),'x',i <= 999)
          LuVecDet = 39
          LuVecDet = IsFreeUnit(LuVecDet)
          call Molcas_open(LuVecDet,filename)
          write(LuVecDet,'(8i4)') nish
          call SGPRWF(SGS,CIS,STSYM,PRWTHR,iSpin,CIV,nConf,Key('PRSD'),LUVECDET)
          ! Close GronOR vecdet file (tps/cdg 20210430)
          close(LuVecDet)
        end if
      end if
    end if

    call mma_deallocate(CIV)
  end if

# ifdef _DMRG_
  if (doDMRG) then
    call mh5_put_dset(wfn_dmrg_checkpoint,dmrg_file%qcmaquis_checkpoint_file)
    call mma_deallocate(d1all)
    if (twordm_qcm) call mma_deallocate(d2all)
    if (doEntanglement) call mma_deallocate(spd1all,safe='*')
  end if
# endif

  call mma_deallocate(CIVEC)
  call mma_deallocate(FMO)

end if

! For RF calculations make sure that the we are following the
! correct root.
!
! In the current implementation the overlap between the CI vectors
! of different macro iterations is used. This criterion stricktly
! only hold if the orbitals are not changed in between the
! interations, to make sure that this approximately holds the
! comparision is only to be considered to be valid if the rotmax
! parameter is below an empirical threshold. In the future the
! procedure for indentifying root flipping has to be made more
! robust.

!SVC: if CISElect is used, roots should be traced regardless of orbital
!     rotation, so in that case ignore the automatic tracing and follow
!     the relative CISE root given in the input by the 'CIRF' keyword.

if (lRF .and. Key('CISE') .and. Key('CIRF') .and. (IPCMROOT > 0)) then
  JPCMROOT = IPCMROOT
  IPCMROOT = IROOT(ICIRFROOT)
  call Put_iScalar('RF CASSCF root',IPCMROOT)
  if (JPCMROOT /= IPCMROOT) write(u6,'(1X,A,I3,A,I3)') 'RF Root has flipped from ',JPCMROOT,' to ',IPCMROOT
else if (lRF .and. (IPCMROOT > 0)) then
  call Qpg_iScalar('RF CASSCF root',Exists)
  if (.not. Exists) then

    ! We are here since we are using the default values.

    call Put_iScalar('RF CASSCF root',IPCMROOT)
    call Put_iScalar('RF0CASSCF root',IPCMROOT)
  end if

  mconf = 0
  call mma_allocate(RF,nConf,Label='RF')
  call Qpg_dArray('RF CASSCF Vector',Exists,mConf)

  !> check whether the rf target h5 file exists (needed at this
  ! point for numerical gradient calculations)
# ifdef _DMRG_
  if (doDMRG .and. Exists) then
    call f_inquire('rf.results_state.h5',rfh5DMRG)
    if (.not. rfh5DMRG) then
      maquis_name_states = ''
      maquis_name_results = ''
      call file_name_generator(IPCMROOT-1,'checkpoint_state.','.h5',maquis_name_states)
      call file_name_generator(IPCMROOT-1,'results_state.','.h5',maquis_name_results)

      !> copy current target wave function to local wave function
      call systemf('cp -f '//trim(maquis_name_results)//' rf.results_state.h5 && rm -rf rf.checkpoint_state.h5 && cp -r '// &
                   trim(maquis_name_states)//' rf.checkpoint_state.h5',iErr)
    end if
  end if
# endif

  if (Exists .and. (mConf == nConf) .and. (iFinal /= 2) .and. ((abs(RotMax) < 1.0e-3_wp) .or. Key('CISE'))) then

    rNorm = One
    ! Shouldn't the overlap in this case be always 1?
    ! For DMRG it seems it is...
    ! But just to make sure we calculate it anyway
    ! in case of non-DMRG calculation
    if (.not. doDMRG) then
      call Get_dArray('RF CASSCF Vector',RF,nConf)
      rNorm = sqrt(DDot_(nConf,RF,1,RF,1))
    end if
    !write(u6,*) 'rNorm=',rNorm
    JPCMROOT = IPCMROOT
    if (rNorm > 1.0e-10_wp) then
      call mma_allocate(Temp,nConf,Label='Temp')
      rMax = Zero
      qMax = Zero
      jDisk = IADR15(4)
      do i=1,lRoots
        if (doDMRG) then
#         ifdef _DMRG_
          qmax = abs(qcmaquis_interface_get_overlap(i))
#         endif
        else
          call DDafile(JOBIPH,2,Temp,nConf,jDisk)
          qMax = abs(DDot_(nConf,Temp,1,RF,1))
        end if
        !write(u6,*) 'qMax=',qMax
        if ((qMax > rMax) .and. (qMax > Half)) then
          rMax = qMax
          JPCMROOT = i
        end if
      end do
      call mma_deallocate(Temp)
    end if
  else
    JPCMROOT = IPCMROOT
  end if

  if (JPCMROOT /= IPCMROOT) then
    write(u6,*) ' RF Root has flipped from ',IPCMROOT,' to ',JPCMROOT
    IPCMROOT = JPCMROOT
    call Put_iScalar('RF CASSCF root',IPCMROOT)
  end if

  if (doDMRG) then
#   ifdef _DMRG_
    maquis_name_states = ''
    maquis_name_results = ''
    call file_name_generator(IPCMROOT-1,'checkpoint_state.','.h5',maquis_name_states)
    call file_name_generator(IPCMROOT-1,'results_state.','.h5',maquis_name_results)

    !> copy current target wave function to local wave function
    call systemf('cp -f '//trim(maquis_name_results)//' rf.results_state.h5 && rm -rf rf.checkpoint_state.h5 && cp -r '// &
                 trim(maquis_name_states)//' rf.checkpoint_state.h5',iErr)
#   endif
  else
    jDisk = IADR15(4)
    do i=1,IPCMROOT-1
      call DDafile(JOBIPH,0,rdum,nConf,jDisk)
    end do
    call DDafile(JOBIPH,2,RF,nConf,jDisk)
  end if

  call Put_dArray('RF CASSCF Vector',RF,nConf)
  call mma_deallocate(RF)
end if

end subroutine CICtl

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
! Copyright (C) 2015, Kamal Sharkas                                    *
!               2019, Thomas J. Duignan                                *
!               2021, Rulin Feng                                       *
!               2026, Dong Q. Le                                       *
!***********************************************************************
!
! Note: The hyperfine code is based on the analogous
! pre-existing G-tensor functionality

module hfcop

#ifdef _HDF5_
use mh5, only: mh5_put_dset
use RASSIWfn, only: wfn_h_hfc_rms
#endif
use Molcas, only: LenIn
use spin_data, only: free_spin_data, get_first_nonzero_GNUC, GNUC_by_nucspin, GNUC_NUCSPIN_by_nucmass, init_spin_data, &
                     NUCSPIN_by_gnuc
use Cntrl, only: AngMom_idx, ASD_idx, Atens_Req, AutoSelect_GFac, DEGEN_ETHR, GNuc, GNuc_set, HypF_rms_Req, HypoIso, LCSTATES, &
                 LPRPR, MLTPLT, NATens_Calc, NAtoms, NCOUP, NMass_set, NPNMR_Calc, NPROP, NSpin_set, NSTATE, NTP, NucMass,     &
                 NucSpin, pNMR_req, PSO_idx, TMAXP, TMINP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Twelve, Half, cZero, auTocm, auToHz, auTokJ, c_in_au, gElectron, kBoltzmann, &
                     proton_mass_in_au
use Definitions, only: iwp, wp, u6

implicit none
private

! RASSI runtime variables
! NSS, ESO

! EPR Calc. :: A tensor, principal and isotropic values
! Atens_fac, prin_vals, signs_resolved

! PNMR Calc :: Shielding tensor, Chemical Shift
! Temp_in_K, pBoltz, Z_HFC_int_oper, Curie_ChemShift, LinRes_ChemShift, LR_tens, C_tens, Z_HFC_over_dE, dE_inv, degen_group,
! n_uniq_ener

! State information and Wigner-Eckart theorem
! MAPST, ETHR_in_cm, CGO_mat, CGx_mat, CGy_mat

! Hamiltonian
! h_FC, h_SD, h_FCSD, h_PSO, h_TOT, USO, h_Zeeman, h_rms_nuc, h_hfc_rms

! Iterators and indices
! iACalc, ipNMR_Calc, degen_start_idx, degen_end_idx, xyz, contrib_lab

! Isotopes and electronic structure
! e_spin, LAtNumb, LAtomLbl, LStability

! Routing variables
! do_calc, do_EPR, do_pNMR

integer(kind=iwp) :: iACalc, ipNMR_Calc, n_uniq_ener, NSS
real(kind=wp) :: Atens_fac, e_spin, ETHR_in_cm
logical(kind=iwp) :: do_calc, do_EPR, do_pNMR
integer(kind=iwp), allocatable :: degen_end_idx(:), degen_group(:), degen_start_idx(:), LAtNumb(:), MAPST(:)
real(kind=wp), allocatable :: C_tens(:,:,:), CGo_mat(:,:), CGx_mat(:,:), CGy_mat(:,:), C_shifts(:,:,:), dE_inv(:,:), &
                              ESO(:), h_hfc_rms(:,:), h_rms_nuc(:,:), LR_shifts(:,:,:), LR_tens(:,:,:), pBoltz(:,:), &
                              prin_vals(:,:,:), Temp_in_K(:), Z_HFC_int_oper(:,:,:,:), Z_HFC_over_dE(:,:,:,:)
complex(kind=wp), allocatable :: h_FC(:,:,:), h_FCSD(:,:,:), h_PSO(:,:,:), h_SD(:,:,:), h_TOT(:,:,:), h_Zeeman(:,:,:), USO(:,:)
logical(kind=iwp), allocatable :: signs_resolved(:)
character(len=LenIn), allocatable :: LAtomLbl(:)
character, allocatable :: LStability(:)

! Conversion factors
real(kind=wp), parameter :: alpha2 = One/(c_in_au*c_in_au), au2J = auTokJ*1.0e3_wp, beta_e = One/(Two*c_in_au), &
                            beta_n = beta_e/proton_mass_in_au, con_to_MHz = -gElectron*beta_e*beta_n*auToHz*1.0e-6_wp, &
                            kBoltzman_in_cm = kBoltzmann*auTocm/au2J, to_ppm = 1.0e6_wp*auTocm*alpha2, TwoThird = Two/Three
character(len=*), parameter :: contrib_lab(5) = [character(len=7) :: 'FC', 'SD','FCSD','PSO','TOTAL']
character, parameter :: xyz(3) = ['x','y','z']

public :: Hyperfine_Oper

contains

subroutine Hyperfine_Oper(PROP,USOR,USOI,JBNUM)

  real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP), USOR(:,:), USOI(:,:)
  integer(kind=iwp), intent(in) :: JBNUM(NSTATE)
  integer(kind=iwp) :: iAtom

  ! Setup reused variables------------------------------
  call setup_hfc_calc(JBNUM,USOR,USOI)
  if (allocated(pNMR_req)) call setup_pNMR_calc(PROP)
  !----------------------------------------------------

  ! do_calc :: calculate Hyperfine hamiltonian
  ! do_EPR  :: calculate A tensor, principal values for EPR spectroscopy
  ! do_pNMR :: calculate paramagnetic NMR chemical shift [Curie and Linear Response]
  ! Calculation begins----------------------------------
  iACalc = 0
  ipNMR_Calc = 0
  do iAtom=1,NAtoms
    call route_calc(iAtom)
    if (do_calc) call calc_h_HFC(iAtom,PROP)
    if (HypF_rms_Req) call update_h_HFC_RMS(iAtom,h_TOT)
  end do

  ! Printing final results------------------------------
  if (HypF_rms_Req) call save_h_rms()
  if (allocated(Atens_Req)) call print_EPR_summary()
  if (allocated(pNMR_req)) call print_pNMR_summary()

  call cleanup_hfcop()

end subroutine Hyperfine_Oper

subroutine setup_hfc_calc(JBNUM,USOR,USOI)

  use wigner_util, only: dclebs

  integer(kind=iwp), intent(in) :: JBNUM(NSTATE)
  real(kind=wp), intent(in) :: USOR(:,:), USOI(:,:)
  integer(kind=iwp) :: ISS, ISTATE, JOB, JSS, MPLET, MSPROJ
  real(kind=wp) :: CGm, CGp, FACT, MPLET1, MPLET2, MSPROJ1, MSPROJ2, S1, S2, SM1, SM2
  integer(kind=iwp), allocatable :: MAPMS(:), MAPSP(:)
  real(kind=wp), allocatable :: rtemp(:)

  ! Form transformation matrix USO with complex numbers
  NSS = size(USOR,1)
  call mma_allocate(USO,NSS,NSS,Label='USO')
  USO(:,:) = cmplx(USOR(:,:),USOI(:,:),kind=wp)

  call mma_allocate(h_FC,3,NSS,NSS,Label='h_FC')
  call mma_allocate(h_SD,3,NSS,NSS,Label='h_SD')
  call mma_allocate(h_FCSD,3,NSS,NSS,Label='h_FCSD')
  call mma_allocate(h_PSO,3,NSS,NSS,Label='h_PSO')
  call mma_allocate(h_TOT,3,NSS,NSS,Label='h_TOT')

  if (HypF_rms_Req) then
    call mma_allocate(h_hfc_rms,NSS,NSS,Label='h_hfc_rms')
    call mma_allocate(h_rms_nuc,NSS,NSS,Label='h_rms_nuc')
    ! IMPORTANT: h_hfc_rms must be initialized to zero, as it is updated via addition within the loop.
    !            h_rms_nuc will be re-assigned later, therefore it does not need to be initialized
    h_hfc_rms(:,:) = Zero
  end if

  if (allocated(Atens_Req)) then
    call mma_allocate(prin_vals,NATens_Calc,5,3,Label='prin_vals')
    call mma_allocate(signs_resolved,NATens_Calc,Label='signs_resolved')
    signs_resolved(:) = .false.
  end if

  ! GET: AtomLbl
  call mma_allocate(LAtomLbl,LenIn*nAtoms,Label='LAtNumb')
  call Get_cArray('Unique Atom Names',LAtomLbl,LenIn*nAtoms)

  ! GET: AtNumb(Z)
  call mma_allocate(rtemp,NAtoms,Label='rtemp_hfcop')
  call mma_allocate(LAtNumb,NAtoms,Label='LAtNumb')
  call Get_dArray('Nuclear charge',rtemp,nAtoms)
  LAtNumb(:) = nint(rtemp(:))
  call mma_deallocate(rtemp)

  ! GET: Energy of SO states (ESO)
  call mma_allocate(ESO,NSS,Label='ESO')
  call get_dArray('ESO_SINGLE',ESO,NSS)

  ! MAP: from spin-free and to spin states:
  call mma_allocate(MAPST,NSS,Label='MAPST')
  call mma_allocate(MAPSP,NSS,Label='MAPSP')
  call mma_allocate(MAPMS,NSS,Label='MAPMS')

  ISS = 0
  do ISTATE=1,NSTATE
    JOB = JBNUM(ISTATE)
    MPLET = MLTPLT(JOB)
    do MSPROJ=-MPLET+1,MPLET-1,2
      ISS = ISS+1
      MAPST(ISS) = ISTATE
      MAPSP(ISS) = MPLET
      MAPMS(ISS) = MSPROJ
    end do
  end do

  ! CALC: Clebsch-Gordan coefficients. Store CGo, CGx, CGy to matrices
  call mma_allocate(CGo_mat,NSS,NSS,Label='CGo')
  call mma_allocate(CGx_mat,NSS,NSS,Label='CGx')
  call mma_allocate(CGy_mat,NSS,NSS,Label='CGy')

  ! Wigner-Eckart theorem
  do ISS=1,NSS

    MPLET1 = MAPSP(ISS)
    MSPROJ1 = MAPMS(ISS)
    S1 = Half*real(MPLET1-1,kind=wp)
    SM1 = Half*real(MSPROJ1,kind=wp)
    do JSS=1,NSS

      MPLET2 = MAPSP(JSS)
      MSPROJ2 = MAPMS(JSS)
      S2 = Half*real(MPLET2-1,kind=wp)
      SM2 = Half*real(MSPROJ2,kind=wp)
      FACT = One/sqrt(real(MPLET1,kind=wp))
      if (MPLET1 == MPLET2-2) FACT = -FACT
      CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
      CGo_mat(ISS,JSS) = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
      CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
      CGx_mat(ISS,JSS) = sqrt(Half)*(CGm-CGp)
      CGy_mat(ISS,JSS) = sqrt(Half)*(CGm+CGp)
    end do
  end do
  ! MAPSP and MAPMS are not needed anymore, but MAPST is required to compute h_HFC
  call mma_deallocate(MAPSP)
  call mma_deallocate(MAPMS)

  ! Process spin data (Nuclear spin + g-Factor)
  if (HypF_rms_Req .or. allocated(Atens_Req)) then
    call proc_spin_data()
    call print_isotope_info()
  end if

  ! GET: Coupled states used to calculate A_tensors and/or pNMR tensors
  ETHR_in_cm = DEGEN_ETHR*auTocm
  if (allocated(Atens_Req) .or. allocated(pNMR_req)) call proc_coupl_states()

end subroutine setup_hfc_calc

subroutine setup_pNMR_calc(PROP)

  real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
  integer(kind=iwp) :: ISS, iT, JSS, lmb
  real(kind=wp) :: dlt_E, dlt_T, Zstat
  real(kind=wp), parameter :: dltE_cutoff = 1.0e-6_wp, pBoltz_cutoff = 1.0e-100_wp

  call mma_allocate(C_shifts,NPNMR_Calc,4,NTP,Label='Curie_ChemShift')
  call mma_allocate(LR_shifts,NPNMR_Calc,4,NTP,Label='LinRes_ChemShift')

  ! Initialize temperature grid
  call mma_allocate(Temp_in_K,NTP,Label='Temp_in_K')
  dlt_T = (TMAXP-TMINP)/(real(NTP-1,kind=wp))
  if (TMINP == Zero) then
    write(u6,*) 'WARNING: TMINP is set to zero. Adjusting TMINP to 0.1K to avoid numerical issues.'
    Temp_in_K(1) = 0.1_wp
  else
    Temp_in_K(1) = TMINP
  end if
  do iT=2,NTP-1
    Temp_in_K(iT) = TMINP+dlt_T*real(iT-1,kind=wp)
  end do
  Temp_in_K(NTP) = TMAXP

  ! Initialize Boltzmann factors
  call mma_allocate(pBoltz,NTP,NSS,Label='p_Boltz')
  do iT=1,NTP
    pBoltz(iT,:) = exp(-ESO(:)/kBoltzman_in_cm/Temp_in_K(iT))
    ! Truncate :: Prevent overflow/underflow when excited states have very high energies
    where(pBoltz(iT,:) < pBoltz_cutoff) pBoltz(iT,:) = Zero
    Zstat = sum(pBoltz(iT,:))
    pBoltz(iT,:) = pBoltz(iT,:)/Zstat
  end do

  ! Calc 1/dE
  call mma_allocate(dE_inv,NSS,NSS,Label='dE_inv')
  dE_inv = Zero
  do ISS=1,NSS
    do JSS=ISS,NSS
      dlt_E = ESO(ISS)-ESO(JSS)
      ! Truncate :: Prevent dividing by small numbers
      if (abs(dlt_E) > dltE_cutoff) dE_inv(ISS,JSS) = -One/dlt_E
    end do
  end do

  call mma_allocate(LR_tens,NTP,3,3,Label='LR_tens')
  call mma_allocate(C_tens,NTP,3,3,Label='C_tens')
  call mma_allocate(Z_HFC_over_dE,3,3,NSS,NSS,Label='Z_HFC_over_dE')
  call mma_allocate(Z_HFC_int_oper,3,3,NSS,NSS,Label='Z_HFC_int_oper')

  ! GROUPING DEGENERATE STATES--------------------------------------
  ! This is only needed for pNMR calculations (to include degenerate excited states).
  ! For EPR (A_tensor), we only take degenerate ground states.
  ! Therefore, it does not appear in setup_hfc_calc subroutine
  call mma_allocate(degen_group,NSS,Label='degen_group')
  n_uniq_ener = size(degen_start_idx)
  do lmb=1,n_uniq_ener
    degen_group(degen_start_idx(lmb):degen_end_idx(lmb)) = lmb
  end do

  call calc_h_Zeeman(PROP)

end subroutine setup_pNMR_calc

subroutine calc_h_HFC(iAtom,PROP)

  integer(kind=iwp), intent(in) :: iAtom
  real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
  integer(kind=iwp) :: idx(6), ISS, iState, JSS, jState
  real(kind=wp) :: A_tens(3,3,5)
  real(kind=wp), allocatable :: ASD(:,:,:)

  idx(:) = ASD_idx(iAtom,:)
  call mma_allocate(ASD,6,NSS,NSS,Label='ASD')
  do ISS=1,NSS
    iState = MAPST(ISS)
    do JSS=ISS,NSS
      jState = MAPST(JSS)
      ASD(:,ISS,JSS) = PROP(iState,jState,idx(:))
      ASD(:,JSS,ISS) = ASD(:,ISS,JSS)
    end do
  end do

  if(do_EPR .or. do_pNMR) then
    write(u6,*) ""
    write(u6,'(3X,A30)') repeat("=",30)
    if (do_EPR .and. do_pNMR)       write(u6,'(3X,A24,A6)') "HFC & pNMR Calc. for :: ", LAtomLbl(iAtom)
    if (do_EPR .and. .not. do_pNMR) write(u6,'(6X,A17,A6)') "HFC Calc. for :: ", LAtomLbl(iAtom)
    if (.not. do_EPR .and. do_pNMR) write(u6,'(6X,A18,A6)') "pNMR Calc. for :: ", LAtomLbl(iAtom)
    write(u6,'(3X,A30)') repeat("=",30)
    write(u6,*) ""
    write(u6,*) ""
  endif


! CALCULATE HAMILTONIAN
  call calc_h_FC(ASD(6,:,:))
  call calc_h_SD(ASD)
  call mma_deallocate(ASD)
  h_FCSD(:,:,:) = h_FC(:,:,:) + h_SD(:,:,:)
  call calc_h_PSO(iAtom, PROP)
  h_TOT(:,:,:) = h_FCSD(:,:,:) + h_PSO(:,:,:)


! TRANSFORM TO SPIN-ORIBT BASIS HAMILTONIAN
  call to_cmpl_SO_states(h_FC)
  call to_cmpl_SO_states(h_SD)
  call to_cmpl_SO_states(h_FCSD)
  call to_cmpl_SO_states(h_PSO)
  call to_cmpl_SO_states(h_TOT)


! PRINT SPIN-ORIBT BASIS HAMILTONIAN
  if(LPRPR) then
    call save_h_hfc(h_FC,   iAtom,  1)
    call save_h_hfc(h_SD,   iAtom,  2)
    call save_h_hfc(h_FCSD, iAtom,  3)
    call save_h_hfc(h_PSO,  iAtom,  4)
    call save_h_hfc(h_TOT,  iAtom,  5)
  endif


! CALCULATE A_TENSOR (EPR)
  if(do_EPR) then
    ! 1. Total first
    call calc_A_tens(A_tens(:,:,5),h_TOT)

    ! 2. Contributions decomposition
    call calc_A_tens(A_tens(:,:,1),h_FC)
    call calc_A_tens(A_tens(:,:,2),h_SD)
    call calc_A_tens(A_tens(:,:,3),h_FCSD)
    call calc_A_tens(A_tens(:,:,4),h_PSO)

    call calc_prin_val(iAtom,A_tens)
  endif


! CALCULATE PNMR_TENSOR
  if(do_pNMR) then
    ! Total
    call calc_pNMR_Tensor(iAtom,h_TOT,5)

    ! Contributions decomposition
    call calc_pNMR_Tensor(iAtom,h_FC,1)
    call calc_pNMR_Tensor(iAtom,h_SD,2)
    ! Skip do_pNMR for FCSD (just a sum)
    call calc_pNMR_Tensor(iAtom,h_PSO,4)
  endif

end subroutine calc_h_HFC

subroutine calc_h_FC(ASD_zz)

  real(kind=wp), intent(in) :: ASD_zz(NSS,NSS)

  h_FC(1,:,:) = cmplx(CGx_mat(:,:)*ASD_zz(:,:),Zero,kind=wp)
  h_FC(2,:,:) = cmplx(Zero,CGy_mat(:,:)*ASD_zz(:,:),kind=wp)
  h_FC(3,:,:) = cmplx(CGo_mat(:,:)*ASD_zz(:,:),Zero,kind=wp)
  h_FC(:,:,:) = TwoThird*h_FC(:,:,:)

end subroutine calc_h_FC

subroutine calc_h_SD(ASD)

  real(kind=wp), intent(out) :: ASD(6,NSS,NSS)

  ASD(:,:,:) = -ASD(:,:,:)
  ASD(6,:,:) = -ASD(1,:,:)-ASD(4,:,:)

  h_SD(1,:,:) = cmplx(CGx_mat(:,:)*ASD(1,:,:)+CGo_mat(:,:)*ASD(3,:,:),CGy_mat(:,:)*ASD(2,:,:),kind=wp)
  h_SD(2,:,:) = cmplx(CGx_mat(:,:)*ASD(2,:,:)+CGo_mat(:,:)*ASD(5,:,:),CGy_mat(:,:)*ASD(4,:,:),kind=wp)
  h_SD(3,:,:) = cmplx(CGx_mat(:,:)*ASD(3,:,:)+CGo_mat(:,:)*ASD(6,:,:),CGy_mat(:,:)*ASD(5,:,:),kind=wp)

end subroutine calc_h_SD

subroutine proc_nuc_spin(iAtom,icase)

  integer(kind=iwp), intent(in) :: iAtom, icase

  write(u6,'(7X,A18,A6)') 'Process for atom: ',LAtomLbl(iAtom)
  select case (icase)
    case (1)
      write(u6,'(11X,A38)') 'USE: Isotopic mass [SEWARD or GATEWAY]'
      call GNUC_NUCSPIN_by_nucmass(LAtNumb(iAtom),NucMass(iAtom),GNuc(iAtom),NucSpin(iAtom),LStability(iAtom))
    case (2)
      write(u6,'(11X,A31,I0)') 'USE: NMASs [RASSI input]   A = ',NucMass(iAtom)
      call GNUC_NUCSPIN_by_nucmass(LAtNumb(iAtom),NucMass(iAtom),GNuc(iAtom),NucSpin(iAtom),LStability(iAtom))
    case (3)
      write(u6,'(11X,A31,F3.1)') 'USE: NSPIn [RASSI input]   I = ',NucSpin(iAtom)
      call GNUC_by_nucspin(LAtNumb(iAtom),NucMass(iAtom),GNuc(iAtom),NucSpin(iAtom),LStability(iAtom))
    case (4)
      write(u6,'(11X,A37,F12.8)') 'USE: GNUC [RASSI input]   g-factor = ',GNuc(iAtom)
      call NUCSPIN_by_gnuc(LAtNumb(iAtom),NucMass(iAtom),GNuc(iAtom),NucSpin(iAtom),LStability(iAtom))
    case (5)
      write(u6,'(11X,A37)') 'USE: Most abundance non-zero g-factor'
      call get_first_nonzero_GNUC(LAtNumb(iAtom),NucMass(iAtom),GNuc(iAtom),NucSpin(iAtom),LStability(iAtom))
  end select

end subroutine proc_nuc_spin

subroutine proc_spin_data()
  !PURPOSE: To find g-factor and NucSpin for given atoms in EZSpin database.

  ! Case 1: Use SEWARD, GATEWAY mass  --> get g-factor and NucSpin
  ! Case 2: Use Nuclear mass     [NMASs, RASSI] --> get g-factor and NucSpin
  ! Case 3: Use Nuclear spin     [NSPIn, RASSI] --> get g-factor
  ! Case 4: Use Nuclear g-factor [GNUC, RASSI] --> get NucSpin
  ! Case 5: Use Most abundant non-zero g-factor (useful for EPR-HFCC)

  ! Note: NucSpin and mass number (integer) has higher priority [IF-ELSE] than g-factor (the last case).

  integer(kind=iwp) :: iAtom, icase
  logical(kind=iwp) :: use_seward_mass
  real(kind=wp), allocatable :: Weights(:)

  call init_spin_data()

  if (.not. allocated(NucSpin)) then
    call mma_allocate(NucSpin,NAtoms,'NucSpin')
    NucSpin(:) = -100.0_wp
  end if

  if (.not. allocated(GNuc)) then
    call mma_allocate(GNuc,NAtoms,'gNuc')
    GNuc(:) = -100.0_wp
  end if

  call mma_allocate(Weights,NAtoms,Label='Weights_hfcop')
  call Get_dArray('Weights',Weights,NAtoms)
  if (NMass_set) then
    ! CASE 1: RASSI and SEWARD might be inconsistent because of user-input NMASs
    do iAtom=1,nAtoms
      if (NucMass(iAtom) == -100.0_wp) NucMass(iAtom) = nint(Weights(iAtom),kind=iwp)
    end do
  else
    ! CASE 2: RASSI and SEWARD mass numbers is consistent at first.
    !         However, users might define wrong spin for unreal isotope
    !         --> Mass number might be incosistent later --> Warnings instead of termination
    call mma_allocate(NucMass,NAtoms,'NucMass_hfcop')
    NucMass(:) = nint(Weights(:),kind=iwp)
  end if

  ! Stability will be re-assigned after processing spin data
  call mma_allocate(LStability,NAtoms,'Stability')
  LStability(:) = '?'

  ! DEFAULT SETTINGS (if users do NOT specify anything : NMASs, NSPIn or GNUC in RASSI input)
  !------------------
  !       1. HFCOperator --> CASE 1 use SEWARD mass
  !                      REASON: to be consistent with molecular dynamics simulations where mass affects the forces on nuclei.
  !       2. HFCAtoms    --> CASE 5 use most abundant non-zero g-factor
  !                      REASON: Users want to compute EPR parameters without looking up the isotopic information.
  !                              This is convienient because the code also print out the conversion factor.
  use_seward_mass = .false.
  ! Setting DEFAULT case based on user input
  if (HypF_rms_Req) then
    if (.not.(AutoSelect_GFac .or. NMass_set .or. NSpin_set .or. GNuc_set)) use_seward_mass = .true.
  else if (allocated(Atens_Req)) then
    if (.not.(AutoSelect_GFac .or. NMass_set .or. NSpin_set .or. GNuc_set)) AutoSelect_GFac = .true.
  end if

  if (use_seward_mass) icase = 1
  if (NMass_set) icase = 2
  if (NSpin_set) icase = 3
  if (GNuc_set) icase = 4
  if (AutoSelect_GFac) icase = 5

  ! HYPOTHEICAL ISOTOPE-------------------------------------------------------------
  if (.not. allocated(HypoIso)) then
    call mma_allocate(HypoIso,NAtoms,'HypoIso_hfcop')
    HypoIso(:) = .false.
  end if
  write(u6,*)
  if (HypF_rms_Req) then
    do iAtom=1,NAtoms
      if (.not. HypoIso(iAtom)) call proc_nuc_spin(iAtom,icase)
      if (NucMass(iAtom) /= nint(Weights(iAtom))) then
        write(u6,'(11X,A28,I3,A24,I3)') 'Warning: Mass number RASSI= ',NucMass(iAtom),' does NOT match SEWARD= ', &
                                        nint(Weights(iAtom))
        write(u6,*) ''
      end if
    end do
  else
    do iAtom=1,NAtoms
      if (.not. HypoIso(iAtom) .and. Atens_Req(iAtom)) call proc_nuc_spin(iAtom,icase)
      if (NucMass(iAtom) /= nint(Weights(iAtom))) then
        write(u6,'(11X,A28,I3,A24,I3)') 'Warning: Mass number RASSI= ',NucMass(iAtom),' does NOT match SEWARD= ', &
                                        nint(Weights(iAtom))
        write(u6,*) ''
      end if
    end do
  end if

  call mma_deallocate(Weights)

  call free_spin_data()

end subroutine proc_spin_data

subroutine route_calc(iAtom)
  ! PURPOSE: Route the calculation when user requests HFC_RMS, A-tensor, pNMR_tensor togehter or separately for different atoms.
  !          Also, it updates iACalc and iPNMR_Calc
  !------------------------------------------------
  !           |HFCOperator | HFCAToms | NMRAToms  |
  !-----------------------------------|-----------|
  ! do_calc   |     TRUE   |    TRUE  |   TRUE    |
  ! do_EPR    |     ___    |    TRUE  |   ___     |
  ! do_pNMR   |     ___    |    ___   |   TRUE    |
  !------------------------------------------------
  !   --- : to be determined by this subroutine.

  ! Note: iACalc and ipNMR_Calc are global variables controlled by route_calc for long-passing.
  !       Do not use iACalc, ipNMR_Calc for loops.

  integer(kind=iwp), intent(in) :: iAtom

  do_calc = .false.
  do_EPR = .false.
  do_pNMR = .false.

  if (allocated(Atens_Req)) then
    do_EPR = Atens_Req(iAtom)
    if (do_EPR) iACalc = iACalc+1
  end if

  if (allocated(pNMR_req)) then
    do_pNMR = pNMR_req(iAtom)
    if (do_pNMR) ipNMR_Calc = ipNMR_Calc+1
  end if

  if (HypF_rms_Req) then
    do_calc = .true.
  else
    do_calc = do_EPR .or. do_pNMR
  end if

end subroutine route_calc

subroutine print_isotope_info()

  integer(kind=iwp) :: iAtom

  write(u6,*)
  write(u6,*)
  write(u6,'(7X,A21)') 'Isotope specification'
  if (HypF_rms_Req) then
    ! Print full table (Spin G-factor)
    write(u6,'(7X,A45)') repeat('-',45)
    write(u6,'(7X,A45)') 'Atom     A   Spin       g-Factor    Stability'
    write(u6,'(7X,A45)') repeat('-',45)
    do iAtom=1,NAtoms
      write(u6,'(7X,A6,1X,I3,4x,F3.1,3x,F12.8,12X,A1)') LAtomLbl(iAtom),NucMass(iAtom),NucSpin(iAtom),GNuc(iAtom),LStability(iAtom)
    end do
    write(u6,'(7X,A45)') repeat('-',45)

  else
    ! Otherwise, only print g-factor. Spin is not needed for EPR calc.
    write(u6,'(7X,A38)') repeat('-',38)
    write(u6,'(7X,A38)') 'Atom     A       g-Factor    Stability'
    write(u6,'(7X,A38)') repeat('-',38)
    do iAtom=1,NAtoms
      if (Atens_Req(iAtom)) write(u6,'(7X,A6,1X,I3,3x,F12.8,12X,A1)') LAtomLbl(iAtom),NucMass(iAtom),GNuc(iAtom),LStability(iAtom)
    end do
    write(u6,'(7X,A38)') repeat('-',38)
  end if
  write(u6,'(7X,A53)') '[ ] = stable, [*] = radioactive, [?] = no information'
  write(u6,*)

  call mma_deallocate(LStability)
  call mma_deallocate(NucMass)

end subroutine print_isotope_info

subroutine print_pNMR_summary()

  integer(kind=iwp) :: iAtom, iT, iContr
  real(kind=wp) :: total_shifts(4)

  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(3X,A96)') repeat('=',96)
  write(u6,'(36X,A28)') 'Summary pNMR Chemical Shifts'
  write(u6,'(3X,A96)') repeat('=',96)
  write(u6,*)
  write(u6,*)

  ipNMR_Calc = 0
  do iAtom=1,nAtoms

    if (pNMR_req(iAtom)) then
      ipNMR_Calc = ipNMR_Calc + 1
      write(u6,*) ""
      write(u6,'(3X,A10,A6)') '>>> ATOM: ', adjustl(LAtomLbl(iAtom))
      write(u6,*) ""

      do iT = 1, NTP

        total_shifts(:) = LR_shifts(ipNMR_Calc,:,iT) + C_shifts(ipNMR_Calc,:,iT)

        write(u6,'(12X,A17,F6.1)')   'Temperature (K): ', Temp_in_K(iT)
        write(u6,'(12X,A75)') repeat('=',75)


        write(u6,'(30X,A12,5X,A12,2(2X,A12))')                                   &
        '   FC+SD+PSO', "          FC", "          SD", "         PSO"
        write(u6,'(30X,A12,5X,A12,2(2X,A12))') (repeat('-',12), iContr=1,4)
        ! Linear Response
        write(u6,'(12x,A17,1x,F12.2,5x,F12.2,2(2x,F12.2))') "LinRes (ppm)     ", LR_shifts(ipNMR_Calc,4,iT), &
                                (LR_shifts(ipNMR_Calc,iContr,iT), iContr=1,3)
        ! Curie
        write(u6,'(12x,A17,1x,F12.2,5x,F12.2,2(2x,F12.2))') "Curie (ppm)      ", C_shifts(ipNMR_Calc,4,iT),  &
                                    (C_shifts(ipNMR_Calc,iContr,iT), iContr=1,3)
        ! Total
        write(u6,'(12x,A17,1x,F12.2,5x,F12.2,2(2x,F12.2))') "Total pNMR shifts", total_shifts(4), &
                                    (total_shifts(iContr), iContr=1,3)
        write(u6,*) ""
        write(u6,*) ""
      enddo
    endif
  enddo
end subroutine print_pNMR_summary

subroutine print_EPR_summary()

  integer(kind=iwp) :: iAtom, iAxis, iCheckVal, iContr
  real(kind=wp) :: Aiso_tot, conv
  character(len=12) :: string_val
  real(kind=wp), allocatable :: checkfile_vals(:)
  character(len=*), parameter :: undef_res = '-/+/-/+/-', unt = '(MHz)'

  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(3X,A96)') repeat('=',96)
  write(u6,'(36X,A28)') 'Summary HFC Principal Values'
  write(u6,'(3X,A96)') repeat('=',96)
  write(u6,*)
  write(u6,*)

  !--> Checkfile
  !        (1) Checkfile should be generated before converting to MHz.
  !            The reason is nuclear gfactor might be updated every year,
  !            which might generate numerical noise in test cases (if spin_data is updated)
  !        (2) Sign determination is experimental. Unsigned principal values are better.
  !        (3) The diagonalization step in A_tensor is sometimes sensitive to compilers and geometry [also symmetry breaking].
  !            Test case should be independent of axis choices.
  !
  ! ----------> checkfile_vals = (ABS(A_xx) + ABS(A_yy) + ABS(A_zz))/3     [note. this not prvl or isotropic vals]
  call mma_allocate(checkfile_vals,NATens_Calc*5,Label='EPR_checkfile')
  iACalc = 0
  iCheckVal = 0
  do iAtom=1,nAtoms
    if (Atens_Req(iAtom)) then
      iACalc = iACalc+1
      do iContr=1,5
        iCheckVal = iCheckVal+1
        checkfile_vals(iCheckVal) = sum(abs(prin_vals(iACalc,iContr,1:3)))/Three
      end do
    end if
  end do
  call Add_Info('EPR_HFCOP_CHKVAL',checkfile_vals,NATens_Calc*5,2)
  call mma_deallocate(checkfile_vals)

  !--> Unit conversion to MHz
  iACalc = 0
  do iAtom=1,nAtoms
    if (Atens_Req(iAtom)) then
      iACalc = iACalc+1
      conv = con_to_MHz*GNuc(iAtom)
      prin_vals(iACalc,:,:) = conv*prin_vals(iACalc,:,:)
      call assign_hfc_prvl_signs()
    end if
  end do

  !--> Print Electronic State information
  write(u6,'(3X,A29,F6.2,2X,A27,I0)') '>>> Electronic Pseudospin :: ',e_spin,'Number of Coupling States: ',NCOUP
  write(u6,*)
  write(u6,*)

  !--> Print summary table of principal values in MHz unit
  iACalc = 0
  do iAtom=1,nAtoms
    if (Atens_Req(iAtom)) then
      iACalc = iACalc+1
      write(u6,*)
      write(u6,'(3X,A10,A6,A19,F12.8)') '>>> ATOM: ',adjustl(LAtomLbl(iAtom)),'Nuclear g-factor = ',GNuc(iAtom)
      if (.not. HypoIso(iAtom)) write(u6,'(19X,A19,2X,F3.1)') 'Nuclear Spin     = ',NucSpin(iAtom)
      write(u6,*)
      write(u6,'(3X,A4,2X,A7,A5,5X,A7,A5,3(2X,A7,A5))') 'Comp',' TOTAL ',unt,'    FC ',unt,'    SD ',unt,'  FCSD ',unt,'   PSO ',unt
      write(u6,'(3X,A4,2X,A12,5X,A12,3(2X,A12))') repeat('-',4),(repeat('-',12),iContr=1,5)

      do iAxis=1,3
        write(u6,'(3X,A2,A1,A1,2x,F12.2,5x,F12.2,3(2x,F12.2))') 'A_',xyz(iAxis),xyz(iAxis),prin_vals(iACalc,5,iAxis), &
                                                                (prin_vals(iACalc,iContr,iAxis),iContr=1,4)
      end do

      !--> Print isotropic values
      if (signs_resolved(iACalc)) then
        Aiso_tot = sum(prin_vals(iACalc,5,:))/Three
        write(u6,'(3X,A4,2X,A12,5X,A12,3(2X,A12))') repeat('-',4),(repeat('-',12),iContr=1,5)
        write(u6,'(3X,A5,1x,F12.2,5x,F12.2,3(2x,F12.2))') 'A_iso',abs(Aiso_tot), &
                                                          (abs(sum(prin_vals(iACalc,iContr,:))/Three),iContr=1,4)
        write(u6,*)
        write(u6,'(3X,A33,F13.3,1X,A5)') '>>>>  Isotropic HFCCs (Total)  = ',abs(Aiso_TOT),unt
        write(u6,*)
      else
        write(u6,'(3X,A5,1x,A12,5x,A12,3(2x,A12))') 'A_iso',undef_res,undef_res,undef_res,undef_res,undef_res
        write(u6,*)
        write(u6,'(17x,A57)') 'NOTE: Signs of principal values cannot be determined.'
      end if
    end if
  end do

  ! FUTURE IMPROVEMENT: Accurate sign determination needs spin density at nucleus. There will be an update for this feature.

  write(u6,*)
  write(u6,*)
  write(u6,*) '   ------------'
  write(u6,*) '     ( Note )'
  write(u6,*) '   ------------'
  write(u6,'(A93)') '       1. All isotropic hyperfine coupling constants (HFCCs) are reported as absolute values.'
  write(u6,'(A88)') '       2. If a different isotope is required, you can convert the reported values to MHz'
  write(u6,'(A65)') '          without restarting by using the provided formula below:'
  write(u6,*)

  write(string_val,'(F12.4)')-gElectron*beta_e*beta_n*1.0e9_wp
  write(u6,'(A52,A7,A4)') '                A(MHz) = A(au) * nuclear-g-factor * ',adjustl(string_val(6:12)),'E-09'
  write(string_val,'(F12.4)') con_to_MHz
  write(u6,'(A60,A7)') '             or A(MHz) = sqrt(eigvval) * nuclear-g-factor * ',adjustl(string_val(6:12))
  write(u6,*)
  write(u6,*)

end subroutine print_EPR_summary

subroutine to_cmpl_SO_states(h)

  complex(kind=wp), intent(out) :: h(3,NSS,NSS)
  integer(kind=iwp) :: u
  complex(kind=wp), allocatable :: tmp_matr(:,:)

  call mma_allocate(tmp_matr,NSS,NSS)
  tmp_matr = cmplx(Zero,Zero,kind=wp)
  do u=1,3
    call zgemm_('n','n',NSS,NSS,NSS,cmplx(1.0_wp,0.0_wp,kind=wp),h(u,:,:),NSS,USO,NSS,cmplx(0.0_wp,0.0_wp,kind=wp),tmp_matr,NSS)
    call zgemm_('c','n',NSS,NSS,NSS,cmplx(1.0_wp,0.0_wp,kind=wp),USO,NSS,tmp_matr,NSS,cmplx(0.0_wp,0.0_wp,kind=wp),h(u,:,:),NSS)
  end do
  call mma_deallocate(tmp_matr)

end subroutine to_cmpl_SO_states

subroutine proc_coupl_states()

  integer(kind=iwp) :: ISS
  real(kind=wp) :: degen, min_energy

  ! PRIORITY of KEYWORDS
  ! NCOU  >  COUPL  > DETH

  ! Note (1): Users might run both HFC and NMR. The results should be consistent, therefore,
  ! both of calculations should couple states by energy threshold. In the otherhand, NCOUP or COUP
  ! keywords give users the freedom to choose states for EPR calculations.

  call get_degen_states(1)

  ! User specifies NCOU without COUP ---> First NCOUP lowest-energy states will be coupled
  if (NCOUP > 0_iwp .and. .not.(allocated(LCSTATES))) then
    call mma_allocate(LCSTATES,NCOUP)
    do ISS=1,NCOUP
      LCSTATES(ISS) = ISS
    end do
  end if

  ! SHIFTING ENERGY by NCOUP/COUP-------------------------------------------
  ! then re-determine the coupled states based on energy threshold.
  if (allocated(LCSTATES)) then
    min_energy = ESO(1)
    do ISS=1,size(LCSTATES)
      ESO(LCSTATES(ISS)) = min_energy
    end do
    call get_degen_states(2)
    ! The assignment NCOUP here is different from lines below in else block.
    ! This line is "after" allocated(LCSTATES) to make consistency.
    ! NCOUP in else block is "before" allocating
    NCOUP = degen_end_idx(1)-degen_start_idx(1)+1_iwp
  else
    ! NO SHIFTING ENERGY----------------------------------------------------
    ! If NCOUP or COUP are not specified, coupling states for the A tensor
    ! are derived from GROUND STATES via get_degen_states(1), controlled by DETH.

    ! Get NCOUP from previous call get_degen_states(1)
    NCOUP = degen_end_idx(1)-degen_start_idx(1)+1_iwp
    call mma_allocate(LCSTATES,NCOUP,'LCStates')
    do ISS=degen_start_idx(1),degen_end_idx(1)
      LCSTATES(ISS) = ISS
    end do
  end if

  if (allocated(Atens_Req)) then
    ! Print states used to calculate A_HFC
    write(u6,*)
    write(u6,*)
    write(u6,'(7X,A25)') ' A TENSOR Coupling States'
    write(u6,'(7X,A25)') '-------------------------'
    write(u6,'(7X,A25)') 'SO-State           ENERGY'
    write(u6,'(7X,A25)') '-------------------------'
    do ISS=1,NCOUP
      write(u6,'(7X,I4,6x,F15.6)') LCSTATES(ISS),ESO(LCSTATES(ISS))
    end do
    write(u6,'(7X,A25)') '-------------------------'
    write(u6,*)
    write(u6,*)

    ! PSEUDOSPIN APPROACH
    !--------------------
    ! Purpose: Calculate Atens_fac used to calculate Atensor
    degen = real(NCOUP,kind=wp)
    e_spin = (degen-One)/Two
    ! REF: 10.1021/acs.jctc.0c01005    [eq. 15]
    ! Plug S = e_spin = (degen-1)/2 into eq. 15, one gets:
    ! Atens_fac = 6/[4*S(S+1)(2S+1)] = 12/(degen^3-degen)
    Atens_fac = Twelve/(degen**3-degen)
  end if

end subroutine proc_coupl_states

subroutine get_degen_states(opt)

  integer(kind=iwp), intent(in) :: opt
  integer(kind=iwp) :: degeneracy, first_state, ISS, last_state, lmb
  real(kind=wp) :: ener,prev_ener

  ! Count number of degenerate groups
  prev_ener = ESO(1)
  n_uniq_ener = 1
  do ISS=2,NSS
    if (abs(ESO(ISS)-prev_ener) >= ETHR_in_cm) then
      n_uniq_ener = n_uniq_ener+1
      prev_ener = ESO(ISS)
    end if
  end do

  ! get_degen_states() can be called multiple times, after shifting energy by NCOU or COUPLedstates.
  if (allocated(degen_start_idx)) call mma_deallocate(degen_start_idx)
  if (allocated(degen_end_idx)) call mma_deallocate(degen_end_idx)

  call mma_allocate(degen_start_idx,n_uniq_ener,Label='beg_idx')
  call mma_allocate(degen_end_idx,n_uniq_ener,Label='end_idx')

  ! Assign begin/end indices for each degenerate group
  degen_start_idx(1) = 1
  lmb = 1
  prev_ener = ESO(1)

  do ISS=2,NSS
    if (abs(ESO(ISS)-prev_ener) > ETHR_in_cm) then
      degen_end_idx(lmb) = ISS-1
      lmb = lmb+1
      degen_start_idx(lmb) = ISS
      prev_ener = ESO(ISS)
    end if
  end do
  degen_end_idx(n_uniq_ener) = NSS

  write(u6,*)
  write(u6,*)
  select case (opt)
    case (1)
      write(u6,'(7x,A34,E9.2,A3)') 'RASSI DEGENERATE SO-STATES within ',DEGEN_ETHR,' au'
    case (2)
      write(u6,'(12x,A35)') '<< SHIFTED by NCOU/COUP keywords >>'
  end select
  write(u6,'(7X,A45)') repeat('-',45)
  write(u6,'(7X,A45)') 'Range [SO-State]  Degen.       Energy (cm^-1)'
  write(u6,'(7X,A45)') repeat('-',45)
  do lmb=1,n_uniq_ener
    first_state = degen_start_idx(lmb)
    last_state = degen_end_idx(lmb)
    degeneracy = last_state-first_state+1_iwp
    ener = ESO(degen_start_idx(lmb))
    write(u6,'(10X,I3,A2,I3,8X,I2,9X,F15.6)') first_state,' - ',last_state,degeneracy,ener
  end do
  write(u6,'(7X,A45)') repeat('-',45)

end subroutine get_degen_states

subroutine calc_A_tens(A_tens,h_HFC)

  real(kind=wp), intent(out) :: A_tens(3,3)
  complex(kind=wp), intent(in) :: h_HFC(3,NSS,NSS)
  integer(kind=iwp) :: u, w
  complex(kind=wp), allocatable :: h_temp(:,:,:)

  ! VARIBLE DESCRIPTION:
  ! u,w    :: are cartesian components {x,y,z}
  ! h_temp :: slice of h_HFC contains only coupling states [degenerate ground states]

  !------> UNOPTIMIZED VERSION (for reference, debug only)
  ! REF: 10.1021/acs.jctc.0c01005    [eq. 14-15]
  ! h_HFC = Atens_fac * sum(<i| h_u |j><j| h_w |i>)
  ! We need two kinds of loop (1) u,w over {x,y,z}
  !                           (2) i,j over {coupled states}
  !
  ! do i=1,degen
  !   do j=1,degen
  !     do u=1,3
  !       do w=1,3
  !         A_tens(u,w)=A_tens(u,w)+real(h_HFC(u,i,j)*h_HFC(w,j,i),kind=wp)
  !       enddo
  !     enddo
  !   enddo
  ! enddo

  call mma_allocate(h_temp,3,NCOUP,NCOUP,Label='h_temp_Atens')
  h_temp(:,1:NCoup,1:NCoup) = h_HFC(:,LCSTATES,LCSTATES)

  !------> VECTORIZED VERSION
  do u=1,3
    do w=1,3
      A_tens(u,w) = sum(real(h_temp(u,:,:)*transpose(h_temp(w,:,:)),kind=wp))
    end do
  end do

  call mma_deallocate(h_temp)
  ! Atens_fac is determined by subroutine get_Atens_coupled_states()
  ! PSEUDOSPIN APPROACH :  REF: 10.1021/acs.jctc.0c01005
  A_tens(:,:) = Atens_fac*A_tens(:,:)

end subroutine calc_A_tens

subroutine assign_hfc_prvl_signs()

  integer(kind=iwp) :: iAxis
  logical(kind=iwp) :: is_determined(3)

  !----------------------------------------------
  ! STEP 1: determine FC  +/- SD  =  +/- FCSD
  is_determined(:) = .false.
  do iAxis=1,3
    call assign_abc_signs(prin_vals(iACalc,1,iAxis),prin_vals(iACalc,2,iAxis),prin_vals(iACalc,3,iAxis),is_determined(iAxis))
  end do

  if (all(is_determined)) then
    !----------------------------------------------
    ! STEP 2: determine FCSD  +/- PSO  =  +/- TOTAL
    is_determined(:) = .false.
    do iAxis=1,3
      call assign_abc_signs(prin_vals(iACalc,3,iAxis),prin_vals(iACalc,4,iAxis),prin_vals(iACalc,5,iAxis),is_determined(iAxis))
    end do
    if (all(is_determined)) signs_resolved(iACalc) = .true.
  end if

  if (.not. signs_resolved(iACalc)) prin_vals(iACalc,:,:) = abs(prin_vals(iACalc,:,:))

end subroutine assign_hfc_prvl_signs

subroutine assign_abc_signs(a,b,c,is_determined)
  ! PURPOSE: Determine signs of this equatiion:
  !           a  +/- |b| = +/- |c|
  ! where sign(a) is fixed (always positive or negative). OUTPUT: sign(b), sign(c)
  ! There are two steps:
  !           STEP 1: Compare  MAX(|a|,|b|) and |c|
  !           STEP 2: Compare with the tolerance [5% of the abs maximum]

  real(kind=wp), intent(in) :: a
  real(kind=wp), intent(out) :: b, c
  logical(kind=iwp), intent(out) :: is_determined
  real(kind=wp) :: abs_a, abs_b, abs_c, max_left_absval
  real(kind=wp), parameter :: tol = 0.05_wp

  ! Tolerance for determining signs
  ! 5% of the maximum absolute value among a, b, c
  abs_a = abs(a)
  abs_b = abs(b)
  abs_c = abs(c)
  max_left_absval = max(abs_a,abs_b)
  is_determined = .false.

  if (abs_c > max_left_absval) then
    ! CASE 1: abs_c increases
    !------------------------
    ! a, b should have the SAME sign --> 2 cases: [+] + [+] = [+]
    !                                          or [-] + [-] = [-]
    if (abs((abs_a+abs_b)/abs_c-One) <= tol) then
      is_determined = .true.
      ! a = [+] --> b = [+]  , c = [+]
      ! a = [-] --> b = [-]  , c = [-]
      if (a < Zero) then
        b = -abs_b
        c = -abs_c
      end if
    end if

  else if (abs_c < max_left_absval) then
    ! CASE 2: abs_c decreases
    !------------------------
    ! a, b should have the DIFFERENT signs --> 2 cases: [+] + [-]
    !                                                or [-] + [+]

    if (abs(abs(abs_a-abs_b)/abs_c-One) <= tol) then
      is_determined = .true.

      ! b needs to have the OPOSITE sign of a. Sign of c depends on a & b
      ! a = [+] --> b = [-]  , c = [+] if |a| > |b| and vice versa
      ! a < [-] --> b = [+]  , c = [-] if |a| > |b| and vice versa

      if (a >= Zero) then
        b = -abs_b
        if (abs_a < abs_b) c = -abs_c
      else if (a < Zero) then
        if (abs_a > abs_b) c = -abs_c
      end if
    end if
  end if

end subroutine assign_abc_signs

subroutine transf_prin_axes(A, X, a_sm)

  real(kind=wp), intent(in) :: A(3,3,5), X(3,3)
  real(kind=wp),intent(out) :: a_sm(3,3,5)
  integer(kind=iwp) :: iContr
  real(kind=wp) :: tmpmat(3,3)

  do iContr=1,5
    tmpmat(:,:) = Zero
    call dgemm_('n','n',3,3,3,1.0_wp,A(:,:,iContr),3,X,3,0.0_wp,tmpmat,3)
    call dgemm_('t','n',3,3,3,1.0_wp,X,3,tmpmat,3,0.0_wp,a_sm(:,:,iContr),3)
  enddo

end subroutine transf_prin_axes

subroutine calc_prin_val(iAtom,A_tens)

  integer(kind=iwp), intent(in) :: iAtom
  real(kind=wp), intent(in) :: A_tens(3,3,5)
  real(kind=wp) :: a_small(3,3,5)
  integer(kind=iwp) :: iAxis, IERR, iContr
  real(kind=wp) :: EVI(3), EVR(3), fnorm_diag, fnorm_off_diag, prvl, tmpmat(3,3), X(3,3)
  real(kind=wp), parameter :: to_au = -gElectron*beta_e*beta_n
  real(kind=wp), external :: dnrm2_

  ! VARIABLE DESCRIPTION:
  !   fnorm_diag:       Frobenius norm of diag(sm_a)
  !   fnorm_off_diag:   Frobenius norm of off-diag(sm_a)
  ! iContr              : Contribution: FC = 1
  !                                     SD = 2
  !                                   FCSD = 3
  !                                    PSO = 4
  !                                  Total = 5
  tmpmat(:,:) = A_tens(:,:,5)
  X(:,:) = Zero
  EVR(:) = Zero
  EVI(:) = Zero
  call XEIGEN(1,3,3,tmpmat,EVR,EVI,X,IERR)

  call transf_prin_axes(A_tens, X, a_small)

  write(u6,'(4X,A14)') "PRINCIPAL AXES"
  write(u6,'(4X,A38)') repeat('-',38)
  write(u6,'(12x,3(A1,12x))') xyz(1:3)
  do iAxis=1,3
    write(u6,'(4X,3(ES12.3,1X))') X(iAxis,1:3)
  end do
  write(u6,*)


  ! Print A-tensor, a-matrix, and principal values for each contribution
  do iContr=1,5

    write(u6,'(3X,A82)') repeat('-',82)
    write(u6,'(3X,A10,A6,A22,A7)') '>>> ATOM: ',LAtomLbl(iAtom),'HYPERFINE COUPLING :: ',contrib_lab(iContr)
    write(u6,'(3X,A82)') repeat('-',82)
    write(u6,'(14X,A17,32X,A8)') 'A-tensor (A=aa^T)','a-matrix'
    write(u6,'(12x,3(A1,12x),3x,3(A1,12x))') xyz(1:3),xyz(1:3)
    do iAxis=1,3
      write(u6,'(3X,A1,3(1x,ES12.3),3x,3(1x,ES12.3))') xyz(iAxis),A_tens(iAxis,1:3,iContr),a_small(iAxis,1:3,iContr)
    end do

    ! CHECK DIAG/OFF-DIAG NORMS OF a-matrix
    ! After transformation, the a-matrix should be diagonal. If not, it indicates numerical instability.
    tmpmat(:,:) = a_small(:,:,iContr)
    ! Calculate diagonal elements norm
    fnorm_diag = sqrt(tmpmat(1,1)**2+tmpmat(2,2)**2+tmpmat(3,3)**2)
    ! Calculate off-diagonal elements norm
    tmpmat(1,1) = Zero
    tmpmat(2,2) = Zero
    tmpmat(3,3) = Zero
    fnorm_off_diag = dnrm2_(9,tmpmat,1)

    if (fnorm_off_diag/fnorm_diag > 0.05_wp) then
      call WarningMessage(1,'Relative Frobenius diag/off-diag norm > 5%')
    end if

    do iAxis=1,3
      EVR(iAxis) = a_small(iAxis,iAxis,iContr)
    end do

    if (any(EVR(:) < Zero)) then
      call WarningMessage(2,'Negative eigenvalues found. Cannot take square root.')
      call AbEnd()
    end if

    ! principal values (sqrt of diagonal elements)
    prin_vals(iACalc,iContr,:) = sqrt(EVR)

    ! Print absolute principal values without signs
    write(u6,*) ''
    write(u6,'(20X,A31,A7)') 'ABS. PRINCIPAL VALUES [+/-] :: ',contrib_lab(iContr)
    write(u6,'(12X,A52)') repeat('.',52)
    write(u6,'(22X,A10,12X,A4,11X,A5)') 'sqrt(a_ii)','(au)','(MHz)'
    write(u6,'(12X,A52)') repeat('.',52)
    do iAxis=1,3
      prvl = prin_vals(iACalc,iContr,iAxis)
      write(u6,'(12X,A1,A1,5X,E13.6,3x,E13.6,3x,E13.6)') xyz(iAxis),xyz(iAxis),prvl,prvl*to_au,prvl*con_to_MHz*GNuc(iAtom)
    end do
    write(u6,*) ''
    write(u6,*) ''
  end do

end subroutine calc_prin_val

subroutine update_h_HFC_RMS(iAtom,h_HFC)

  integer(kind=iwp), intent(in) :: iAtom
  complex(kind=wp), intent(in) :: h_HFC(3,NSS,NSS)
  real(kind=wp) :: fac, I_sq, NSpin_I

  ! VERSION 1: PARTLY VECTORIZED
  ! h_rms_nuc =  abs(h_HFC(1,:,:))**2 + abs(h_HFC(2,:,:))**2 + abs(h_HFC(3,:,:))**2

  ! VERSION 2: LEGACY COMPILERS
  ! h_rms_nuc(:,:) = Zero
  ! do i = 1, 3
  !   h_rms_nuc = h_rms_nuc + real(h_HFC(i,:,:))**2 + aimag(h_HFC(i,:,:))**2
  ! end do

  ! OPTIMIZED VERSION :: FULLY VECTORIZED
  h_rms_nuc(:,:) = sum(abs(h_HFC(:,:,:))**2,dim=1)

  ! Scale with factor || I_sq ||^2
  NSpin_I = NucSpin(iAtom)
  I_sq = NSpin_I*(NSpin_I+One)*(Two*NSpin_I+One)/Three
  fac = I_sq*GNuc(iAtom)**2
  h_rms_nuc(:,:) = fac*h_rms_nuc(:,:)

  ! Update h_RMS for this nuclei
  h_hfc_rms(:,:) = h_hfc_rms(:,:)+h_rms_nuc(:,:)

  if (iAtom == NAtoms) then
    ! Note: con_to_MHz = gElectron*beta_e*beta_n*auToHz*1e-6_wp
    fac = -gElectron*beta_e*beta_n
    h_hfc_rms(:,:) = fac*sqrt(h_hfc_rms(:,:))
  end if
end subroutine update_h_HFC_RMS

subroutine calc_h_PSO(iAtom,PROP)

  integer(kind=iwp), intent(in) :: iAtom
  real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
  integer(kind=iwp) :: u
  real(kind=wp), allocatable :: Im_h_PSO(:,:,:)

  call mma_allocate(Im_h_PSO,3,NSS,NSS,Label='Im_h_PSO')
  Im_h_PSO(:,:,:) = Zero
  do u=1,3
    call SMMAT(PROP,Im_h_PSO(u,:,:),NSS,PSO_idx(iAtom,u),u)
  end do
  h_PSO(:,:,:) = cmplx(0.0_wp,Im_h_PSO(:,:,:),kind=wp)
  call mma_deallocate(Im_h_PSO)

end subroutine calc_h_PSO

subroutine save_h_hfc(h_HFC,iAtom, iContr)

  complex(kind=wp), intent(in) :: h_HFC(3,NSS,NSS)
  integer(kind=iwp), intent(in) :: iAtom, iContr
  integer(kind=iwp) :: LU, JSS, ISS, u, istatus
  integer(kind=iwp), External:: IsFreeUnit
  character(len=35) :: file_name
  logical(kind=iwp) :: is_error


  do u = 1, 3
    file_name = 'h_' // trim(contrib_lab(iContr)) // '_' // trim(LAtomLbl(iAtom)) // '_' // xyz(u) // '.txt'
    Lu = IsFreeUnit(88)
    istatus = 100
    call molcas_open_ext2(Lu,file_name,'SEQUENTIAL','FORMATTED',istatus,.false.,1,'REPLACE',is_error)
    write(Lu,*) "NSS = ",NSS
    write(Lu,*) "#NROW NCOL REAL IMAG"
    do JSS=1,NSS
      do ISS=1,NSS
      write(Lu,'(I6,1X,I6,1X,ES25.16,1X,ES25.16)') ISS,JSS, real(h_HFC(u,ISS,JSS)), aimag(h_HFC(u,ISS,JSS))
      end do
    end do
  CLOSE(Lu)
  enddo
end subroutine

subroutine save_h_rms()

  integer(kind=iwp) :: ISS, istatus, JSTA, LU
  logical(kind=iwp) :: is_error
  integer(kind=iwp), external :: IsFreeUnit

# ifdef _HDF5_
  call mh5_put_dset(wfn_h_hfc_rms,h_hfc_rms)
# endif

  Lu = IsFreeUnit(88)
  istatus = 100
  call molcas_open_ext2(Lu,'h_RMS.txt','SEQUENTIAL','FORMATTED',istatus,.false.,1,'REPLACE',is_error)
  write(Lu,*) 'NSS= ',NSS
  write(Lu,*) '#NROW NCOL REAL'
  do JSTA=1,NSS
    do ISS=1,NSS
      write(Lu,'(I6,1X,I6,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ',h_hfc_rms(ISS,JSTA)
    end do
  end do
  close(Lu)

end subroutine save_h_rms

subroutine calc_h_Zeeman(PROP)

  real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
  real(kind=wp), allocatable :: Angmom(:,:,:), Im_h_Zeeman(:,:,:), Re_h_Zeeman(:,:,:)

  call mma_allocate(Re_h_Zeeman,3,NSS,NSS,Label='Re_h_Zeeman')
  call mma_allocate(Im_h_Zeeman,3,NSS,NSS,Label='Im_h_Zeeman')
  call mma_allocate(Angmom,3,NSS,NSS,Label='Angmom')
  call mma_allocate(h_Zeeman,3,NSS,NSS,Label='h_Zeeman')

  h_Zeeman(:,:,:) = cZero
  Re_h_Zeeman(:,:,:) = Zero
  Im_h_Zeeman(:,:,:) = Zero

  Angmom(:,:,:) = Zero

  call SMMAT(PROP,Angmom(1,:,:),NSS,AngMom_idx(1),0)
  call SMMAT(PROP,Angmom(2,:,:),NSS,AngMom_idx(2),0)
  call SMMAT(PROP,Angmom(3,:,:),NSS,AngMom_idx(3),0)

  call SMMAT(PROP,Re_h_Zeeman(1,:,:),NSS,0,1)
  call SMMAT(PROP,Im_h_Zeeman(2,:,:),NSS,0,2)
  call SMMAT(PROP,Re_h_Zeeman(3,:,:),NSS,0,3)

  Re_h_Zeeman(1,:,:) = -gElectron*Re_h_Zeeman(1,:,:)
  Im_h_Zeeman(2,:,:) = -gElectron*Im_h_Zeeman(2,:,:)
  Re_h_Zeeman(3,:,:) = -gElectron*Re_h_Zeeman(3,:,:)
  Im_h_Zeeman(:,:,:) = Im_h_Zeeman(:,:,:)+Angmom(:,:,:)

  h_Zeeman(:,:,:) = cmplx(Re_h_Zeeman(:,:,:),Im_h_Zeeman(:,:,:),kind=wp)
  call to_cmpl_SO_states(h_Zeeman)

  h_Zeeman(:,:,:) = h_Zeeman(:,:,:)/Two

  call mma_deallocate(Re_h_Zeeman)
  call mma_deallocate(Im_h_Zeeman)
  call mma_deallocate(Angmom)
  ! h_Zeeman needs to deallocated at the end [after pNMR calc. finishes]

end subroutine calc_h_Zeeman

subroutine calc_pNMR_Tensor(iAtom,h_HFC,iContr)

  integer(kind=iwp), intent(in) :: iAtom, iContr
  complex(kind=wp), intent(in) :: h_HFC(3,NSS,NSS)

  integer(kind=iwp) :: ISS, iT, lmb, lmb_a, lmb_ap, u, w, ipNMR_contr
  real(kind=wp) :: fac

  ! VARIABLES DESCRIPTION
  ! REF: DOI: 10.1021/acs.jctc.6b00462   [eq. 1]
  ! lmb                 : "lambda" is the group of all degenerate states.
  ! lmb_a & lmb_ap      : states between lmb_a - lmb_ap are degenerate
  ! Z_HFC_int_oper      : <i |h_Zeeman| j> <j |h_HFC| i>
  ! Z_HFC_over_dE       : Z_HFC_int_oper / delta_E(ij)
  ! LR_tens             : Linear Response tensor [first term,  eq 1]
  ! C_tens              : Curie tensor           [second term, eq 1]
  !               iContr (HFC)   ipNMR_Contr
  !            FC    1              1
  !            SD    2              2
  !            FCSD  3              Skip!
  !            PSO   4              3
  !            TOTAL 5              4
  !
  ! We skip FCSD, therefore ipNMR_Contr needs to be updated (to prevent out of bounds)
  ipNMR_contr = iContr
  if (iContr >= 4) ipNMR_contr = ipNMR_contr - 1

  ! Calculate temperature-indepedent terms (numerator) in ! REF: DOI: 10.1021/acs.jctc.6b00462 Eq. 1
  do u=1,3
    do w=1,3
      Z_HFC_int_oper(u,w,:,:) = real(h_Zeeman(u,:,:)*conjg(h_HFC(w,:,:)),kind=wp)
      Z_HFC_over_dE(u,w,:,:) = Z_HFC_int_oper(u,w,:,:)*dE_inv(:,:)
    end do
  end do

  ! Linear response, Curie pNMR tensor are temperature-dependent
  LR_tens(:,:,:) = Zero
  C_tens(:,:,:) = Zero
  do iT=1,NTP
    do ISS=1,NSS
      ! lambda (lmb) is a group of degenereate states
      ! All degenerate states have the SAME lambda
      lmb = degen_group(ISS)
      lmb_a = degen_start_idx(lmb)
      lmb_ap = degen_end_idx(lmb)

      ! Vectorized Curie tensor is trivial, loop between lmb_a and lmb_ap [degenerate states]
      C_tens(iT,:,:) = C_tens(iT,:,:)+pBoltz(iT,ISS)*sum(Z_HFC_int_oper(:,:,ISS,lmb_a:lmb_ap),dim=3)

      ! Non-degenerate states will contribute to Linear Response term
      ! Two if-logic has been used to avoid empty array (lmb_a == 1 means 1:lmb_a-1 == 1:0)
      if (lmb_a > 1) LR_tens(iT,:,:) = LR_tens(iT,:,:)+pBoltz(iT,ISS)*sum(Z_HFC_over_dE(:,:,ISS,1:lmb_a-1),dim=3)
      if (lmb_ap < NSS) LR_tens(iT,:,:) = LR_tens(iT,:,:)+pBoltz(iT,ISS)*sum(Z_HFC_over_dE(:,:,ISS,lmb_ap+1:NSS),dim=3)

    end do
    fac = One/(kBoltzman_in_cm*Temp_in_K(iT))
    C_tens(iT,:,:) = fac*C_tens(iT,:,:)
  end do !Temperature

  ! Unit conversion to ppm
  C_tens(:,:,:) = to_ppm*C_tens(:,:,:)
  fac = to_ppm*Two
  LR_tens(:,:,:) = fac*LR_tens(:,:,:)

  write(u6,'(3X,A91)') repeat('-',91)
  write(u6,'(3X,A10,A6,A9,A7)') '>>> ATOM: ', LAtomLbl(iAtom), ' pNMR :: ', contrib_lab(iContr)
  write(u6,'(3X,A91)') repeat('-',91)
  do iT=1,NTP
    call print_pNMR_tens(Temp_in_K(iT),LR_tens(iT,:,:), C_tens(iT,:,:))
  end do

  ! Calculate chemical shift by average diagonal elements
  ! Store chemical shift to print summary later
  ! Ref: 10.1016/bs.arcc.2015.09.006
  C_shifts(ipNMR_Calc,ipNMR_contr,:)  = (C_tens(:,1,1)+C_tens(:,2,2)+C_tens(:,3,3))/Three
  LR_shifts(ipNMR_Calc,ipNMR_contr,:) = (LR_tens(:,1,1)+LR_tens(:,2,2)+LR_tens(:,3,3))/Three

end subroutine calc_pNMR_Tensor

subroutine print_pNMR_tens(temp, LinRes,Curie)

  real(kind=wp), intent(in) :: temp, LinRes(3,3), Curie(3,3)
  real(kind=wp) :: lr_shift, curie_shift
  integer(kind=iwp) :: u

  lr_shift = (LinRes(1,1) + LinRes(2,2) + LinRes(3,3))/Three
  curie_shift = (Curie(1,1) + Curie(2,2) + Curie(3,3))/Three

  write(u6,'(3X,A8,13X,A15,34X,A5)') 'Temp (K)', 'Linear Response', 'Curie'
  write(u6,'(2X,F6.1,3(2x,ES12.3),2X,3(2x,ES12.3))') temp, LinRes(1,1:3), Curie(1,1:3)
  do u = 2, 3
    write(u6,'(8X,3(2x,ES12.3),2X,3(2x,ES12.3))') LinRes(u,1:3), Curie(u,1:3)
  end do
  write(u6,*) ''
  write(u6,'(12X,A15,ES12.3,17X,A15,ES12.3)') 'LinRes Shift : ', lr_shift, 'Curie Shift : ', curie_shift
  write(u6,*) ''
  write(u6,*) ''

end subroutine print_pNMR_tens

subroutine cleanup_hfcop()

  call mma_deallocate(MAPST)
  call mma_deallocate(CGx_mat)
  call mma_deallocate(CGy_mat)
  call mma_deallocate(CGo_mat)
  call mma_deallocate(LAtNumb)
  call mma_deallocate(LAtomLbl)
  call mma_deallocate(USO)
  call mma_deallocate(ESO)
  call mma_deallocate(h_FC)
  call mma_deallocate(h_SD)
  call mma_deallocate(h_FCSD)
  call mma_deallocate(h_PSO)
  call mma_deallocate(h_TOT)

  if (HypF_rms_Req) then
    call mma_deallocate(h_hfc_rms)
    call mma_deallocate(h_rms_nuc)
  end if

  if (allocated(pNMR_req)) then
    call mma_deallocate(C_shifts)
    call mma_deallocate(LR_shifts)
    call mma_deallocate(pNMR_req)
    call mma_deallocate(pBoltz)
    call mma_deallocate(Temp_in_K)

    call mma_deallocate(Z_HFC_int_oper)
    call mma_deallocate(Z_HFC_over_dE)
    call mma_deallocate(AngMom_idx)
    call mma_deallocate(h_Zeeman)

    call mma_deallocate(dE_inv)
    call mma_deallocate(LR_tens)
    call mma_deallocate(C_tens)
    call mma_deallocate(degen_group)
  end if

  call mma_deallocate(HypoIso,safe='*')

  call mma_deallocate(NucSpin,safe='*')
  call mma_deallocate(GNuc,safe='*')
  call mma_deallocate(LCSTATES,safe='*')

  call mma_deallocate(ASD_idx,safe='*')
  call mma_deallocate(PSO_idx,safe='*')

  call mma_deallocate(degen_start_idx,safe='*')
  call mma_deallocate(degen_end_idx,safe='*')

  if (allocated(Atens_Req)) then
    call mma_deallocate(Atens_Req)
    call mma_deallocate(prin_vals)
    call mma_deallocate(signs_resolved)
  end if

end subroutine cleanup_hfcop

end module hfcop

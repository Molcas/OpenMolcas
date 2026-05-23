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
!
module hyperfine
#ifdef _HDF5_
  use mh5,         only: mh5_put_dset
  use RASSIWfn,    only: wfn_h_hfc_rms
#endif
  use Molcas,      only: LenIn
  use Definitions, only: iwp, wp
  use Constants,   only: Zero, Half, One, Two, Three, Four, Six, Twelve, auTocm, c_in_au,      &
                         gElectron, auToHz, kBoltzmann, auTokJ
  use stdalloc,    only: mma_allocate, mma_deallocate
  use hfc_data
  use Cntrl,       only: NSTATE, NPROP, PNAME,NAtoms, MULTIP, ICOMP, MLTPLT, ASD_idx, PSO_idx, &
                         NPNMR_Calc, pNMR_req, HypF_rms_Req, NATens_Calc, Atens_Req,           &
                         NucMass, NMass_set, NucSpin, NSpin_set, NucGFac, NGFac_set, Hypo_Iso, &
                         AutoSelect_GFac, LCSTATES, NCOUP, NTP, TMINP, TMAXP, DEGEN_ETHR

  use hfc_logical, only: MAG_X2C

  implicit none

  ! RASSI runtime variables
  integer(kind=iwp)               :: NSS
  real(kind=wp), allocatable      :: ESO(:)


  ! Conversion factors
  real(kind=wp), parameter :: e_proton_mass_ratio = 1836.152673426d0,                          &
                              beta_e = One/(Two*c_in_au), TwoThird = Two/Three,                &
                              beta_n = beta_e/e_proton_mass_ratio,                             &
                              con_to_MHz = -gElectron*beta_e*beta_n*auToHz*1.0e-6_wp,          &
                              alpha2 = One/(c_in_au * c_in_au),                                &
                              to_ppm = 1.0e6_wp * auTocm * alpha2,                             &
                              au2J = auTokJ*1.0e3_wp, kBoltzman_in_cm = kBoltzmann*auTocm/au2J


  ! EPR Calc. :: A tensor, principal and isotropic values
  real(kind=wp)                   :: Atens_fac
  real(kind=wp), allocatable      :: prin_vals(:,:,:)
  logical(kind=iwp), allocatable  :: signs_resolved(:)


  ! PNMR Calc :: Shielding tensor, Chemical Shift
  real(kind=wp),allocatable       :: Temp_in_K(:), pBoltz(:,:), Z_HFC_int_oper(:,:,:,:),       &
                                     Curie_ChemShift(:,:,:), LinRes_ChemShift(:,:,:),          &
                                     LR_tens(:,:,:), C_tens(:,:,:),  Z_HFC_over_dE(:,:,:,:),   &
                                     dE_inv(:,:)
  integer(kind=iwp), allocatable  ::  degen_group(:)
  integer(kind=iwp)               ::  n_uniq_ener


  ! State information and Wigner-Eckart theorem
  integer(kind=iwp),allocatable   :: MAPST(:)
  real(kind=wp) :: ETHR_in_cm
  real(kind=wp),allocatable       :: CGo_mat(:,:), CGx_mat(:,:), CGy_mat(:,:)


  ! Hamiltonian
  complex(kind=wp),allocatable    :: h_FC(:,:,:) , h_SD(:,:,:), h_FCSD(:,:,:), h_PSO(:,:,:),   &
                                     h_TOT(:,:,:), USO(:,:), h_Zeeman(:,:,:)
  real(kind=wp),allocatable       :: h_rms_nuc(:,:), h_hfc_rms(:,:)


  ! Iterators and indices
  integer(kind=iwp)               :: JOB, IERR, nData, iACalc, ipNMR_Calc
  integer(kind=iwp), allocatable  :: degen_start_idx(:), degen_end_idx(:)
  character(LEN=1)                :: xyz(3) = ['x', 'y', 'z']
  character(len=7)                :: contrib_lab(5) = [character(len=7) :: '[TOTAL]', '[FC]',  &
                                    '[SD]', '[FCSD]', '[PSO]']


  ! Isotopes and electronic structure
  real(kind=wp)                      ::  e_spin
  integer(kind=iwp), allocatable     ::  LAtNumb(:)
  character(len=LenIn), allocatable  ::  LAtomLbl(:)

  ! Routing variables
  logical(kind=iwp)                  :: do_calc, do_EPR, do_pNMR


  public :: HFCOP

!----------------------------------------------------------------------
  Contains



  subroutine HFCOP(PROP,USOR,USOI,JBNUM)
    integer(kind=iwp)               :: iAtom
    integer(kind=iwp), intent(in)   ::  JBNUM(NSTATE)
    real(kind=wp), intent(in)       ::  PROP(NSTATE,NSTATE,NPROP), USOR(:,:), USOI(:,:)


! Setup reused variables------------------------------
    call setup_hfc_calc(JBNUM,USOR,USOI)
    if(allocated(pNMR_req)) call setup_pNMR_calc(PROP)
!-----------------------------------------------------

! do_calc :: calculate Hyperfine hamiltonian
! do_EPR  :: calculate A tensor, principal values for EPR spectroscopy
! do_pNMR :: calculate paramagnetic NMR chemical shift [Curie and Linear Response]
! Calculation begins----------------------------------
    iACalc = 0
    ipNMR_Calc = 0
    do iAtom = 1, NAtoms
      call route_calc(iAtom)
      if(do_calc)       call calc_h_HFC(iAtom,PROP)
      if(HypF_rms_Req)  call update_h_HFC_RMS(iAtom,h_TOT)
    enddo
!-----------------------------------------------------

! Printing final results------------------------------
    if (HypF_rms_Req) call save_h_rms()
    if(allocated(Atens_Req)) call print_EPR_summary()
    if(allocated(pNMR_req))  call print_pNMR_summary()
!-----------------------------------------------------

    call cleanup_hfcop()

    end subroutine HFCOP
!======================================================================


!======================================================================
  subroutine setup_hfc_calc(JBNUM, USOR, USOI)
    integer(kind=iwp), intent(in)   :: JBNUM(NSTATE)
    real(kind=wp)                   :: rtemp(NAtoms)
    real(kind=wp)                   :: DCLEBS
    real(kind=wp),intent(in)        :: USOR(:,:), USOI(:,:)

    integer(kind=iwp),allocatable   :: MAPSP(:), MAPMS(:)
    integer(kind=iwp)               :: MSPROJ,  MPLET, ISS, JSS, ISTATE, JSTATE
    real(kind=wp)                   :: S1, S2, SM1, SM2,FACT, MPLET1,MSPROJ1, MPLET2, MSPROJ2,   &
                                       IMLTPL, CGm, Co, CGp, Cx, Cy

    NSS = size(USOR,1)
    call mma_allocate(USO, NSS, NSS, Label='USO')
    USO(:,:) = cmplx(USOR(:,:),USOI(:,:),kind=wp)


    call mma_allocate(h_FC,   3, NSS, NSS, Label='h_FC')
    call mma_allocate(h_SD,   3, NSS, NSS, Label='h_SD')
    call mma_allocate(h_FCSD, 3, NSS, NSS, Label='h_FCSD')
    call mma_allocate(h_PSO,  3, NSS, NSS, Label='h_PSO')
    call mma_allocate(h_TOT,  3, NSS, NSS, Label='h_TOT')

    if(HypF_rms_Req) then
      call mma_allocate(h_hfc_rms, NSS, NSS, Label='h_hfc_rms')
      call mma_allocate(h_rms_nuc, NSS, NSS, Label='h_rms_nuc')
      ! IMPORTANT: h_hfc_rms must be initialized to zero, as it is updated via addition within the loop.
      !            h_rms_nuc will be re-assigned later, therefore it does not need to be initialized
      h_hfc_rms(:,:) = Zero
    endif

    if(allocated(Atens_Req)) then
      call mma_allocate(prin_vals, NATens_Calc, 5, 3, Label="prin_vals")
      call mma_allocate(signs_resolved, NATens_Calc, Label="signs_resolved")
      signs_resolved(:) = .false.
    endif

! GET: AtomLbl
    call mma_allocate(LAtomLbl,LenIn*nAtoms, Label='LAtNumb')
    call Get_cArray('Unique Atom Names',LAtomLbl,LenIn*nAtoms)

! GET: AtNumb(Z)
    call mma_allocate(LAtNumb,NAtoms, Label='LAtNumb')
    call Get_dArray('Nuclear charge', rtemp, nAtoms)
    LAtNumb(:) = int(rtemp(:))

! GET: Eenergy of SO states (ESO)
    call mma_allocate(ESO,NSS, Label='ESO')
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


! CALC: Clebsch-Gordan coefficients
  call mma_allocate(CGo_mat,NSS,NSS,Label='CGo')
  call mma_allocate(CGx_mat,NSS,NSS,Label='CGx')
  call mma_allocate(CGy_mat,NSS,NSS,Label='CGy')


! Wigner-Eckart theorem
  do ISS=1,NSS
    ISTATE = MAPST(ISS)
    MPLET1 = MAPSP(ISS)
    MSPROJ1 = MAPMS(ISS)
    S1 = Half*real(MPLET1-1,kind=wp)
    SM1 = Half*real(MSPROJ1,kind=wp)
    do JSS=1,NSS
      JSTATE = MAPST(JSS)
      MPLET2 = MAPSP(JSS)
      MSPROJ2 = MAPMS(JSS)
      S2 = Half*real(MPLET2-1,kind=wp)
      SM2 = Half*real(MSPROJ2,kind=wp)
      FACT = One/sqrt(real(MPLET1,kind=wp))
      if (MPLET1 == MPLET2-2) FACT = -FACT
      CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
      CGo_mat(ISS,JSS)=FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
      CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
      CGx_mat(ISS,JSS) = SQRT(Half)*(CGm-CGp)
      CGy_mat(ISS,JSS) = SQRT(Half)*(CGm+CGp)
    end do
  end do
! MAPSP and MAPMS are not needed anymore, but MAPST is required to compute h_HFC
  call mma_deallocate(MAPSP)
  call mma_deallocate(MAPMS)

! Process spin data (Nuclear spin + g-Factor)
  if (HypF_rms_Req .or. allocated(Atens_Req)) then
    call proc_spin_data()
    call print_isotope_info()
  endif

! GET: Coupled states used to calculate A_tensor
    ETHR_in_cm = DEGEN_ETHR * auTocm
    if(allocated(Atens_Req) .or. allocated(pNMR_req)) call proc_coupl_states()

! These variables are not needed
    if(allocated(NucMass)) call mma_deallocate(NucMass)
  end subroutine
!======================================================================


!======================================================================
  subroutine setup_pNMR_calc(PROP)
    real(kind=wp), intent(in)   ::  PROP(NSTATE,NSTATE,NPROP)
    real(kind=wp)               :: dlt_T, Zstat, dlt_E
    real(kind=wp),parameter     :: pBoltz_cutoff = 10.0_wp **(-100), dltE_cutoff = 10.0_wp ** (-6)

    integer(kind=iwp)           :: iT, ISS, JSS, lmb

    call mma_allocate(Curie_ChemShift, NPNMR_Calc, 4, NTP, Label="Curie_ChemShift")
    call mma_allocate(LinRes_ChemShift, NPNMR_Calc, 4, NTP, Label="LinRes_ChemShift")

! Initialized temperature grid
    call mma_allocate(Temp_in_K,NTP, Label="Temp_in_K")
    dlt_T=(TMAXP-TMINP)/(real(NTP-1 ,kind=wp))
    if(TMINP == Zero) then
      write(6,*) "WARNING: TMINP is set to zero. Adjusting TMINP to 0.1K to avoid numerical issues."
      Temp_in_K(1)=0.1_wp
    else
      Temp_in_K(1)=TMINP
    endif
    do iT=2,NTP-1
      Temp_in_K(iT)=TMINP+dlt_T*real(iT-1,kind=wp)
    enddo
    Temp_in_K(NTP)=TMAXP

! Initialized Boltzmann factors
    call mma_allocate(pBoltz,NTP,NSS, Label="p_Boltz")
    do iT=1,NTP
      pBoltz(iT,:) = exp(-ESO(:)/kBoltzman_in_cm/Temp_in_K(iT))
      ! Truncate :: Prevent overflow/underflow when excited states have very high energies
      where (pBoltz(iT,:) < pBoltz_cutoff)
        pBoltz(iT,:) = Zero
      end where
      Zstat        = sum(pBoltz(iT,:))
      pBoltz(iT,:) = pBoltz(iT,:)/Zstat
    enddo


! Calc 1/dE
    call mma_allocate(dE_inv,NSS,NSS,Label="dE_inv")
    dE_inv = Zero
    do ISS = 1, NSS
      do JSS = ISS, NSS
        dlt_E = ESO(ISS) - ESO(JSS)
        ! Truncate :: Prevent dividing by small numbers
        if(ABS(dlt_E) > dltE_cutoff) then
          dE_inv(ISS,JSS) = One/dlt_E
          dE_inv(ISS,JSS) = -dE_inv(ISS,JSS)
        end if
      enddo
    enddo

    call mma_allocate(LR_tens, NTP, 3, 3, Label="LR_tens")
    call mma_allocate(C_tens, NTP, 3, 3, Label="C_tens")
    call mma_allocate(Z_HFC_over_dE,3,3 ,NSS,NSS, Label="Z_HFC_over_dE" )
    call mma_allocate(Z_HFC_int_oper,3,3,NSS,NSS, Label="Z_HFC_int_oper")


    call mma_allocate(degen_group,NSS, Label="degen_group")
    n_uniq_ener = size(degen_start_idx)
    do lmb = 1, n_uniq_ener
      degen_group(degen_start_idx(lmb):degen_end_idx(lmb)) = lmb
    end do

    call calc_h_Zeeman(PROP)

  end subroutine

!======================================================================


!======================================================================
  subroutine calc_h_HFC(iAtom, PROP)
    integer(kind=iwp), intent(in) :: iAtom
    real(kind=wp), intent(in)     :: PROP(NSTATE,NSTATE,NPROP)
    real(kind=wp) ,allocatable    :: ASD(:,:,:)
    integer(kind=iwp)             :: idx(6)
    real(kind=wp)                 :: A_tens(3,3)

    integer(kind=iwp)             :: ISS, JSS, iState, jState

    idx(:) = ASD_idx(iAtom,:)
    call mma_allocate(ASD,6,NSS,NSS,Label="ASD")
    do ISS = 1, NSS
      iState = MAPST(ISS)
      do JSS = 1,NSS
        jState = MAPST(JSS)
        ASD(:,ISS,JSS) = PROP(iState,jState,idx(:))
        ASD(:,JSS,ISS) = ASD(:,ISS,JSS)
      enddo
    enddo

    if(do_EPR .or. do_pNMR) then
      write(6,*) ""
      write(6,*) ""
      write(6,*) ""
      write(6,*) ""
      write(6,'(3X,A30)') repeat("=",30)
      if (do_EPR .and. do_pNMR)       write(6,'(3X,A24,A6)') "HFC & pNMR Calc. for :: ", LAtomLbl(iAtom)
      if (do_EPR .and. .not. do_pNMR) write(6,'(6X,A17,A6)') "HFC Calc. for :: ", LAtomLbl(iAtom)
      if (.not. do_EPR .and. do_pNMR) write(6,'(6X,A18,A6)') "pNMR Calc. for :: ", LAtomLbl(iAtom)
      write(6,'(3X,A30)') repeat("=",30)
    endif

! NOTE: While it is possible to make the code concise by enclosing if(do_EPR) or do_pNMR
!    for different contributions FC, SD,... it has been made by designing that
!    the code always prints partial output A_tensor or principal values
!    in case of abnormal termination. Most users need A_tens and prinval more
!    than hamiltonian.
!-----------------------------------
    call calc_h_FC(ASD(6,:,:))
    if(do_EPR) then
      call calc_A_tens(A_tens,h_FC)
      call calc_prin_val(iAtom,A_tens,1)
    endif
    if(do_pNMR) call calc_pNMR_Tensor(iAtom,h_FC,1)
!-----------------------------------
    call calc_h_SD(ASD)
    if(do_EPR) then
      call calc_A_tens(A_tens,h_SD)
      call calc_prin_val(iAtom,A_tens,2)
    endif
    if(do_pNMR) call calc_pNMR_Tensor(iAtom,h_SD,2)
!-----------------------------------
    h_FCSD(:,:,:) = h_FC(:,:,:) + h_SD(:,:,:)
    if(do_EPR) then
      call calc_A_tens(A_tens,h_FCSD)
      call calc_prin_val(iAtom,A_tens,3)
    endif
    ! Skip do_pNMR for FCSD (just a sum)
!-----------------------------------
    call calc_h_PSO(iAtom, PROP)
    if(do_EPR) then
      call calc_A_tens(A_tens,h_PSO)
      call calc_prin_val(iAtom,A_tens,4)
    endif
    if(do_pNMR) call calc_pNMR_Tensor(iAtom,h_PSO,3)
!-----------------------------------
    h_TOT(:,:,:) = h_FCSD(:,:,:) + h_PSO(:,:,:)
    if(do_EPR) then
      call calc_A_tens(A_tens,h_TOT)
      call calc_prin_val(iAtom,A_tens,5)
    endif
    if(do_pNMR) call calc_pNMR_Tensor(iAtom,h_TOT,4)

    call mma_deallocate(ASD)

  end subroutine
!======================================================================


!======================================================================
  subroutine calc_h_FC(ASD_zz)
    real(kind=wp), intent(in)  ::  ASD_zz(NSS,NSS)

    h_FC(1,:,:)=cmplx(CGx_mat(:,:)*ASD_zz(:,:),0.0_wp,kind=wp)
    h_FC(2,:,:)=cmplx(0.0_wp,CGy_mat(:,:)*ASD_zz(:,:),kind=wp)
    h_FC(3,:,:)=cmplx(CGo_mat(:,:)*ASD_zz(:,:),0.0_wp,kind=wp)
    h_FC(:,:,:) = TwoThird * h_FC(:,:,:)

    call to_cmpl_SO_states(h_FC)

  end subroutine
!======================================================================


!======================================================================
  subroutine calc_h_SD(ASD)
    real(kind=wp),intent(out)   :: ASD(6,NSS,NSS)

    ASD(:,:,:)  = -ASD(:,:,:)
    ASD(6,:,:)  = -ASD(1,:,:) - ASD(4,:,:)

    h_SD(1,:,:) = cmplx(CGx_mat(:,:) * ASD(1,:,:) + CGo_mat(:,:) * ASD(3,:,:) , CGy_mat(:,:) * ASD(2,:,:), kind=wp)
    h_SD(2,:,:) = cmplx(CGx_mat(:,:) * ASD(2,:,:) + CGo_mat(:,:) * ASD(5,:,:) , CGy_mat(:,:) * ASD(4,:,:), kind=wp)
    h_SD(3,:,:) = cmplx(CGx_mat(:,:) * ASD(3,:,:) + CGo_mat(:,:) * ASD(6,:,:) , CGy_mat(:,:) * ASD(5,:,:), kind=wp)

    call to_cmpl_SO_states(h_SD)

  end subroutine
!======================================================================


!======================================================================
  subroutine proc_nuc_spin_data(iAtom,case,seward_mass)
    integer(kind=iwp), intent(in)             :: iAtom, case
    integer(kind=iwp), intent(in), optional   :: seward_mass

    write(6,'(7X,A18,A6)') 'Process for atom: ', LAtomLbl(iAtom)
    select case (case)
    case(1)
      write(6,'(11X,A38)') 'USE: Isotopic mass [SEWARD or GATEWAY]'
      call gfac_spin_by_mass(LAtNumb(iAtom), seward_mass, NucGFac(iAtom),NucSpin(iAtom))
    case(2)
      write(6,'(11X,A31,I0)') 'USE: NMASs [RASSI input]   A = ', NucMass(iAtom)
      call gfac_spin_by_mass(LAtNumb(iAtom), NucMass(iAtom), NucGFac(iAtom), NucSpin(iAtom))
    case(3)
      write(6,'(11X,A31,F6.2)') 'USE: NSPIn [RASSI input]   I = ', NucSpin(iAtom)
      NucGFac(iAtom) = gfactor_by_NucSpin(LAtNumb(iAtom), NucSpin(iAtom))
    case(4)
      write(6,'(11X,A37,F12.8)') 'USE: GNUC [RASSI input]   g-factor = ', NucGFac(iAtom)
      NucSpin(iAtom) = NucSpin_by_gfactor(LAtNumb(iAtom), NucGFac(iAtom))
    case(5)
      write(6,'(11X,A28,F12.8)') 'USE: Most abundance non-zero'
      call get_first_nonzero_gfactor(LAtNumb(iAtom),NucGFac(iAtom),NucSpin(iAtom))
    end select

  end subroutine
!======================================================================

!======================================================================
  subroutine proc_spin_data()
    real(kind=wp)           :: Weights(NAtoms)
    integer(kind=iwp)       :: case
    logical(kind=iwp)       :: use_seward_mass, not_set_by_users
    integer(kind=iwp)       :: iAtom
!PURPOSE: To find g-factor and NucSpin for given atoms in EZSpin database.

! Case 1: Use SEWARD, GATEWAY mass  --> get g-factor and NucSpin
! Case 2: Use ISOMAss [RASSI] --> get g-factor and NucSpin
! Case 3: Use NucSpin [RASSI] --> get g-factor
! Case 4: Use GNUC    [RASSI] --> get NucSpin
! Case 5: Use Most abundant non-zero g-factor (useful for EPR-HFCC)

! Note: NucSpin and mass number (integer) has higher priority [IF-ELSE] than g-factor (the last case).

    call init_isotope_data()

    if (.not. allocated(NucSpin))     call mma_allocate(NucSpin,NAtoms,"NucSpin")
    if (.not. allocated(NucGFac))     call mma_allocate(NucGFac,NAtoms,"gNuc")

    call Get_dArray('Weights', Weights, NAtoms)


! DEFAULT SETTINGS (if users do NOT specify anything : ISOMass, NucSpin or GNUC in RASSI input)
!------------------
!       1. HFCOperator --> CASE 1 use SEWARD mass to get g-factor and NucSpin
!                      REASON: to be consistent with molecular dynamics simulations where mass affects the forces on nuclei.
!       2. HFCAtoms    --> CASE 5 use most abundant non-zero g-factor to get g-factor and NucSpin
!                      REASON: Users want to compute EPR parameters without looking up the isotopic information.
!                              This is convienient because the code also print out the conversion factor.
    use_seward_mass = .false.
    ! Setting DEFAULT case based on user input
    if (HypF_rms_Req) then
      if(.not. (AutoSelect_GFac .or. NMass_set .or. NSpin_set .or. NGFac_set)) use_seward_mass = .true.
    else if (allocated(Atens_Req)) then
      if(.not. (AutoSelect_GFac .or. NMass_set .or. NSpin_set .or. NGFac_set)) AutoSelect_GFac   = .true.
    endif

    if (use_seward_mass)  case = 1
    if (NMass_set)        case = 2
    if (NSpin_set)        case = 3
    if (NGFac_set)        case = 4
    if (AutoSelect_GFac)  case = 5

!   If users specify both NucSpin and GNUC, they should use "HISO" keyword
!   to indicate hypothetical isotopes rather than finding real isotopes in EZSpin data.
!   Hypo_Iso will be turned on automatically if both NucSpin and GNUC are set by users.
    if (NSpin_set .and. NGFac_set) Hypo_Iso = .true.

! When HFC_RMS hamiltonian is requested, we need to process the nuclear data for ALL atoms.
    if (HypF_rms_Req) then
      do iAtom = 1, NAtoms
        not_set_by_users = (NSpin_set  .and. NucSpin(iAtom)  == -100.0_wp)  .or. &
                            (NGFac_set .and. NucGFac(iAtom)  == -100.0_wp)  .or. &
                            (NMass_set .and. NucMass(iAtom)  == -100_iwp)

        if(not_set_by_users) then
          call proc_nuc_spin_data(iAtom, 1, int(Weights(iAtom)))
        else
          if(.not. Hypo_Iso) call proc_nuc_spin_data(iAtom, case, int(Weights(iAtom)))
        endif
      enddo
! Only process atoms for which A tensors (EPR parameters) are requested.
    else
      do iAtom = 1, NAtoms
        if (Atens_Req(iAtom)) then
          if (.not. Hypo_Iso) call proc_nuc_spin_data(iAtom, case, int(Weights(iAtom)))
        end if
      enddo
    endif

    call free_isotope_data()
  end subroutine
!======================================================================


!======================================================================
  subroutine route_calc(iAtom)
    integer(kind=iwp), intent(in) :: iAtom
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

    do_calc  = .false.
    do_EPR   = .false.
    do_pNMR  = .false.

    if (allocated(Atens_Req)) then
      do_EPR = Atens_Req(iAtom)
      if(do_EPR) iACalc = iACalc + 1
    end if

    if (allocated(pNMR_req)) then
      do_pNMR = pNMR_req(iAtom)
      if(do_pNMR) ipNMR_Calc = ipNMR_Calc + 1
    end if

    if (HypF_rms_Req) then
      do_calc = .true.
    else
      do_calc = do_EPR .or. do_pNMR
    endif

  end subroutine route_calc
!======================================================================


!======================================================================
  subroutine print_isotope_info()
    integer(kind=iwp)  :: iAtom
    write(6,*)
    write(6,*)

    if(HypF_rms_Req) then

      ! Print full table (Spin G-factor)
      write(6,'(7X,A32)') repeat('-',32)
      write(6,'(7X,A32)')'Atom     Spin           g-Factor'
      write(6,'(7X,A32)') repeat('-',32)
      do iAtom=1,NAtoms
        write(6,'(7X,A6,1x,F6.2,7x,F12.8)') LAtomLbl(iAtom), NucSpin(iAtom), NucGFac(iAtom)
      enddo
      write(6,'(7X,A32)') repeat('-',32)

    else
    ! Otherwise, only print g-factor. Spin is not needed for EPR calc.
      write(6,'(7X,A22)') repeat('-',22)
      write(6,'(7X,A22)')'Atom          g-Factor'
      write(6,'(7X,A22)') repeat('-',22)
      do iAtom=1,NAtoms
        if(Atens_Req(iAtom)) write(6,'(7X,A6,4x,F12.8)') LAtomLbl(iAtom), NucGFac(iAtom)
      enddo
      write(6,'(7X,A22)') repeat('-',22)
    endif
  end subroutine
!======================================================================


!======================================================================
subroutine print_pNMR_summary()
  integer(kind=iwp)     :: iAtom, iT, iContr

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,'(3X,A96)') REPEAT('=',96)
  write(6,'(36X,A28)') 'Summary pNMR Chemical Shifts'
  write(6,'(3X,A96)') REPEAT('=',96)
  write(6,*)
  write(6,*)

  ipNMR_Calc = 0
  do iAtom = 1, nAtoms

    if (pNMR_req(iAtom)) then
      ipNMR_Calc = ipNMR_Calc + 1
      write(6,*) ""
      write(6,'(3X,A10,A6)') '>>> ATOM: ', adjustl(LAtomLbl(iAtom))
      write(6,*) ""

      write(6,'(3X,A66)') repeat('-',66)
      write(6,'(3X,A10,A6,A20)') '>>> ATOM: ', adjustl(LAtomLbl(iAtom)), 'CURIE CHEMICAL SHIFT'
      write(6,'(3X,A66)') repeat('-',66)
      write(6,'(3X,A7,2X,A12,5X,A12,2(2X,A12))') 'Temp(K)',                                  &
      ' Curie (ppm)', "        [FC]", "        [SD]",  "       [PSO]"
      write(6,'(3X,A7,2X,A12,5X,A12,2(2X,A12))') repeat('-',7),(repeat('-',12), iContr=1,4)
      do iT = 1, NTP
        write(6,'(3X,F7.1,2x,F12.2,5x,F12.2,2(2x,F12.2))') Temp_in_K(iT) ,                        &
        Curie_ChemShift(ipNMR_Calc,4,iT), (Curie_ChemShift(ipNMR_Calc,iContr,iT), iContr=1,3)
      enddo


      write(6,'(3X,A66)') repeat('-',66)
      write(6,'(3X,A10,A6,A30)') '>>> ATOM: ',adjustl(LAtomLbl(iAtom)), 'LINEAR RESPONSE CHEMICAL SHIFT'
      write(6,'(3X,A66)') repeat('-',66)
      write(6,'(3X,A7,2X,A12,5X,A12,2(2X,A12))') 'Temp(K)',                                  &
      'LinRes (ppm)', "        [FC]", "        [SD]", "       [PSO]"
      write(6,'(3X,A7,2X,A12,5X,A12,2(2X,A12))') repeat('-',7),(repeat('-',12), iContr=1,4)

      do iT = 1, NTP
        write(6,'(3X,F7.1,2x,F12.2,5x,F12.2,2(2x,F12.2))') Temp_in_K(iT) ,                        &
        LinRes_ChemShift(ipNMR_Calc,4,iT), (LinRes_ChemShift(ipNMR_Calc,iContr,iT), iContr=1,3)
      enddo
    endif
  enddo
end subroutine
!======================================================================


!======================================================================
subroutine print_EPR_summary()
  real(kind=wp)             :: conv, Aiso_tot
  character(len=5)          :: unit ="(MHz)"
  character(len=12)         :: undef_res
  character(len=16)         :: string_val

  integer(kind=iwp)         :: iAtom, iContr, iAxis

  write(6,*)
  write(6,*)
  write(6,*)
  write(6,'(3X,A96)') REPEAT('=',96)
  write(6,'(36X,A28)') 'Summary HFC Principal Values'
  write(6,'(3X,A96)') REPEAT('=',96)
  write(6,*)
  write(6,*)

  call assign_hfc_prvl_signs()

!--> Unit conversion to MHz
!--------------------------
  iACalc = 0
  do iAtom = 1, nAtoms
    if (Atens_Req(iAtom)) then
      iACalc = iACalc + 1
      conv = con_to_MHz * NucGFac(iAtom)
      prin_vals(iACalc,:,:) = conv * prin_vals(iACalc,:,:)
    end if
  end do
!--------------------------


!--> Print summary table in MHz unit
  write(6,'(3X,A29,F6.2,2X,A27,I0)') '>>> Electronic Pseudospin :: ',e_spin, "Number of Coupling States: ",NCOUP
  write(6,*)
  write(6,*)
  iACalc = 0
  do iAtom = 1, nAtoms
    if (Atens_Req(iAtom)) then
      iACalc = iACalc + 1

!--> Print principal values
      write(6,*) ""
      write(6,'(3X,A10,A6,A12,F12.4,A12,F3.1)') '>>> ATOM: ',adjustl(LAtomLbl(iAtom)), &
              ' g-factor = ',NucGFac(iAtom),'  NucSpin = ',NucSpin(iAtom)
      write(6,*) ""
      write(6,'(3X,A4,2X,A7,A5,5X,A7,A5,3(2X,A7,A5))') 'Comp',                         &
      ' TOTAL ',unit,'    FC ',unit,'    SD ',unit,'  FCSD ',unit,'   PSO ',unit
      write(6,'(3X,A4,2X,A12,5X,A12,3(2X,A12))') repeat('-',4),(repeat('-',12), iContr=1,5)

      do iAxis = 1, 3
        write(6,'(3X,A2,A1,A1,2x,F12.2,5x,F12.2,3(2x,F12.2))') 'A_',xyz(iAxis),xyz(iAxis), &
        prin_vals(iACalc,5,iAxis) , (prin_vals(iACalc,iContr,iAxis), iContr=1,4)
      enddo

!--> Print isotropic values
      if (signs_resolved(iACalc)) then
        Aiso_tot = sum(prin_vals(iACalc,5,:))/Three
        write(6,'(3X,A4,2X,A12,5X,A12,3(2X,A12))') repeat('-',4),(repeat('-',12), iContr=1,5)
        write(6,'(3X,A5,1x,F12.2,5x,F12.2,3(2x,F12.2))') 'A_iso', &
        abs(Aiso_tot) , (abs(sum(prin_vals(iACalc,iContr,:))/Three), iContr = 1, 4)
        write(6,*) ""
        write(6,'(3X,A33,F13.3,1X,A5)') ">>>>  Isotropic HFCCs (Total)  = ", abs(Aiso_TOT) , unit
        write(6,*) ""
      else
        undef_res = "-/+/-/+/-"
        write(6,'(3X,A5,1x,A12,5x,A12,3(2x,A12))') 'A_iso', undef_res, undef_res, undef_res, undef_res, undef_res
        write(6,*) ""
        write(6,'(17x,A57)') 'NOTE: Signs of principal values cannot be determined.'
      endif
    endif
  enddo

  ! FUTURE IMPROVEMENT: Accurate sign determination needs spin density at nucleus. There will be an update for this feature.

  write(6,*)
  write(6,*)
  write(6,*) "   ------------"
  write(6,*) "     ( Note )"
  write(6,*) "   ------------"
  write(6,'(A93)')  "       1. All isotropic hyperfine coupling constants (HFCCs) are reported as absolute values."
  write(6,'(A88)')  "       2. If a different isotope is required, you can convert the reported values to MHz"
  write(6,'(A65)')  "          without restarting by using the provided formula below:"
  write(6,*)
  write(string_val,'(F12.4)') con_to_MHz
  write(6,'(A52,E11.4)') "                A(MHz) = A(au) * nuclear-g-factor * ", -gElectron * beta_e * beta_n
  write(6,'(A60,A16)') "             or A(MHz) = sqrt(eigvval) * nuclear-g-factor * ", adjustl(string_val)
  write(6,*)
  write(6,*)

end subroutine
!======================================================================


!======================================================================
subroutine to_cmpl_SO_states(h)
  complex(kind=wp), intent(out) :: h(3,NSS,NSS)
  complex(kind=wp),allocatable  :: tmp_matr(:,:)
  integer(kind=iwp)             :: u

  call mma_allocate(tmp_matr,NSS,NSS)

  tmp_matr=cmplx(Zero,Zero,kind=wp)
  do u = 1, 3
    call zgemm_('n','n',NSS,NSS,NSS,cmplx(1.0_wp,0.0_wp,kind=wp),h(u,:,:),NSS,USO,NSS,cmplx(0.0_wp,0.0_wp,kind=wp),tmp_matr,NSS)
    call zgemm_('c','n',NSS,NSS,NSS,cmplx(1.0_wp,0.0_wp,kind=wp),USO,NSS,tmp_matr,NSS,cmplx(0.0_wp,0.0_wp,kind=wp),h(u,:,:),NSS)
  enddo

  call mma_deallocate(tmp_matr)

end subroutine
!======================================================================


!======================================================================
  subroutine proc_coupl_states()
    real(kind=wp)       :: degen, min_energy
    integer(kind=iwp)   :: ISS
! PRIORITY of KEYWORDS
! NCOUP  >  COUPL  > EPRA

! Note (1): Users might run both HFC and NMR. The results should be consistent, therefore,
! both of calculations should couple states by energy threshold. In the otherhand, NCOUPLED or COUPLEDstates
! keywords give users the freedom to choose states for EPR calculations.

    call get_degen_states(1)

    ! NCOUP > COUPL :: NCOUP has higher priority than CPOUP
    if(NCOUP > 0_iwp) then
      if(allocated(LCSTATES)) then
        write(6,'(7X, A86)') "WARNING: Both NCOUP/COUP are set. --> NCOUP will be used and LCSTATES will be ignored."
        call mma_deallocate(LCSTATES)
      end if
      call mma_allocate(LCSTATES,NCOUP,'LCStates')
      do ISS = 1, NCOUP
        LCSTATES(ISS) = ISS
      enddo
    endif

    ! SHIFT ENERGY by NCOUPLED or COUPLEDstates if specified,
    ! then re-determine the coupled states based on energy threshold.
    if(allocated(LCSTATES)) then
      min_energy = ESO(1)
      do ISS = 1, size(LCSTATES)
        ESO(LCSTATES(ISS)) = min_energy
      enddo
      call get_degen_states(2)
      ! The assignment NCOUP here is different from lines below in else block.
      ! This line is "after" allocated(LCSTATES) to make consistency.
      ! NCOUP in else block is "before" allocating
      NCOUP = degen_end_idx(1) - degen_start_idx(1) + 1_iwp
    else
      ! If NCOUP or COUP are not specified, coupling states for the A tensor
      ! are derived from GROUND STATES via get_coupled_degen_states, controlled by ETHR.

      ! Get NCOUP from previous call get_degen_states(1)
      NCOUP = degen_end_idx(1) - degen_start_idx(1) + 1_iwp
      call mma_allocate(LCSTATES,NCOUP,"LCStates")
      do ISS = degen_start_idx(1), degen_end_idx(1)
        LCSTATES(ISS) = ISS
      end do
    endif


    if(allocated(Atens_Req)) then
      ! Print states used to calculate A_HFC
      write(6,*)
      write(6,*)
      write(6,'(7X,A25)') " A TENSOR Coupling States"
      write(6,'(7X,A25)') '-------------------------'
      write(6,'(7X,A25)') 'SO-State           ENERGY'
      write(6,'(7X,A25)') '-------------------------'
      do ISS = 1, NCOUP
        write(6,'(7X,I4,6x,F15.6)') LCSTATES(ISS), ESO(LCSTATES(ISS))
      enddo
      write(6,'(7X,A25)') '-------------------------'
      write(6,*)
      write(6,*)

      ! PSEUDOSPIN APPROACH
      !--------------------
      ! Purpose: Calculate Atens_fac used to calculate Atensor
      degen = real(NCOUP,kind=wp)
      e_spin=(degen - One)/Two
      ! REF: 10.1021/acs.jctc.0c01005    [eq. 15]
      ! Plug S = e_spin = (degen-1)/2 into eq. 15, one gets:
      ! Atens_fac = 6/[4*S(S+1)(2S+1)] = 12/(degen^3-degen)
      Atens_fac = Twelve/(degen**3 - degen)
    endif

  end subroutine
!======================================================================


!======================================================================
subroutine get_degen_states(opt)
  integer(kind=iwp) :: lmb, opt, degeneracy, first_state, last_state
  real(kind=wp)     :: prev_ener, ener
  integer(kind=iwp) :: ISS

! Similar to get_Atens_coupled_states() but for pNMR shift calculation
! Count number of degenerate groups
  prev_ener=ESO(1)
  n_uniq_ener=1
  do ISS = 2, NSS
    if (abs(ESO(ISS) - prev_ener) >= ETHR_in_cm) then
      n_uniq_ener = n_uniq_ener + 1
      prev_ener = ESO(ISS)
    endif
  end do

! get_coupled_degen_states() can be called multiple times, after shifting energy by NCOU or COUPLedstates.
  if (allocated(degen_start_idx)) call mma_deallocate(degen_start_idx)
  if (allocated(degen_end_idx)) call mma_deallocate(degen_end_idx)

  call mma_allocate(degen_start_idx,n_uniq_ener,Label='beg_idx')
  call mma_allocate(degen_end_idx,n_uniq_ener,Label='end_idx')

  ! Assign begin/end indices for each degenerate group
  degen_start_idx(1)=1
  lmb = 1
  prev_ener = ESO(1)

  do ISS = 2, NSS
    if(abs(ESO(ISS) - prev_ener) > ETHR_in_cm) then
      degen_end_idx(lmb) = ISS-1
      lmb = lmb+1
      degen_start_idx(lmb) = ISS
      prev_ener = ESO(ISS)
    endif
  end do
  degen_end_idx(n_uniq_ener) = NSS

  write(6,*)
  write(6,*)
  if (opt == 1) write(6,'(7x,A34,E9.2,A3)') "RASSI DEGENERATE SO-STATES within ", DEGEN_ETHR," au"
  if (opt == 2) write(6,'(12x,A35)') "<< SHIFTED by NCOU/COUP keywords >>"
  write(6,'(7X,A45)') repeat('-',45)
  write(6,'(7X,A45)') 'Range [SO-State]  Degen.       Energy (cm^-1)'
  write(6,'(7X,A45)') repeat('-',45)
  do lmb = 1, n_uniq_ener
    first_state = degen_start_idx(lmb)
    last_state  = degen_end_idx(lmb)
    degeneracy  = last_state - first_state + 1_iwp
    ener        = ESO(degen_start_idx(lmb))
    write(6,'(10X,I3,A2,I3,8X,I2,9X,F15.6)') first_state, " - ", last_state, degeneracy, ener
  enddo
  write(6,'(7X,A45)') repeat('-',45)
end subroutine
 !======================================================================


!======================================================================
subroutine calc_A_tens(A_tens,h_HFC)
  ! INPUT VARIABLE
  complex(kind=wp), intent(in)  :: h_HFC(3,NSS,NSS)

  ! OUTPUT VARIABLE
  real(kind=wp),intent(out)     :: A_tens(3,3)

  ! LOCAL VARIABLE
  integer(kind=iwp)             :: u, w
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

  call mma_allocate(h_temp,3,NCOUP,NCOUP,Label="h_temp_Atens")
  h_temp(:,1:NCoup,1:NCoup) = h_HFC(:,LCSTATES,LCSTATES)

  !------> VECTORIZED VERSION
  do u=1,3
    do w=1,3
      A_tens(u,w) = sum(real(h_temp(u,:,:)*transpose(h_temp(w,:,:)), kind=wp))
    enddo
  enddo

  call mma_deallocate(h_temp)
  ! Atens_fac is determined by subroutine get_Atens_coupled_states()
  A_tens(:,:) = Atens_fac * A_tens(:,:)

end subroutine
!======================================================================


!======================================================================
  subroutine assign_hfc_prvl_signs()
    logical(kind=iwp)     :: is_determined(3)
    integer(kind=iwp)     :: iAxis

    do iACalc = 1, NATens_Calc

      !----------------------------------------------
      ! STEP 1: determine FC  +/- SD  =  +/- FCSD
      is_determined(:) = .false.
      do iAxis = 1, 3
        ! call assign_abc_signs(FC_prvl(iACalc,i), SD_prvl(iACalc,i),FCSD_prvl(iACalc,i),is_determined(i))
        call assign_abc_signs(prin_vals(iACalc,1,iAxis), prin_vals(iACalc,2,iAxis),prin_vals(iACalc,3,iAxis),is_determined(iAxis))
      enddo
      if (all(is_determined)) then

      !----------------------------------------------
      ! STEP 2: determine FCSD  +/- PSO  =  +/- TOTAL
        is_determined(:) = .false.
        do iAxis = 1, 3
          ! call assign_abc_signs(FCSD_prvl(iACalc,i),PSO_prvl(iACalc,i),Tot_prvl(iACalc,i),is_determined(i))
          call assign_abc_signs(prin_vals(iACalc,3,iAxis), prin_vals(iACalc,4,iAxis),prin_vals(iACalc,5,iAxis),is_determined(iAxis))
        enddo
      endif
      if(all(is_determined)) signs_resolved(iACalc) = .true.
    enddo

  end subroutine assign_hfc_prvl_signs
!======================================================================


!======================================================================
subroutine assign_abc_signs(a,b,c,is_determined)
  real(kind=wp), intent(out)    :: a, b, c
  logical,intent(out)           :: is_determined
  real(kind=wp)                 :: tol, max_left_absval, abs_a,abs_b,abs_c

! PURPOSE: Determine signs of this equatiion:
!           a  +/- |b| = +/- |c|
! where sign(a) is fixed (always positive or negative). OUTPUT: sign(b), sign(c)
! There are two steps:
!           STEP 1: Compare  MAX(|a|,|b|) and |c|
!           STEP 2: Compare with the tolerance [5% of the abs maximum]


! Tolerance for determining signs
! 5% of the maximum absolute value among a, b, c
  tol = 0.05_wp
  abs_a = abs(a)
  abs_b = abs(b)
  abs_c = abs(c)
  max_left_absval = maxval([abs_a,abs_b])


  if( abs_c > max_left_absval ) then
    ! CASE 1: abs_c increases
    !------------------------
    ! a, b should have the SAME sign --> 2 cases: [+] + [+] = [+]
    !                                          or [-] + [-] = [-]
    if (abs((abs_a + abs_b)/abs_c - One ) <= tol) then
      is_determined = .true.
    ! a = [+] --> b = [+]  , c = [+]
    ! a = [-] --> b = [-]  , c = [-]
      if (a < Zero) then
        b = -abs_b
        c = -abs_c
      endif
    endif

  else if ( abs_c < max_left_absval ) then
    ! CASE 2: abs_c decreases
    !------------------------
    ! a, b should have the DIFFERENT signs --> 2 cases: [+] + [-]
    !                                                or [-] + [+]

    if (abs(abs(abs_a - abs_b)/abs_c - One ) <= tol) THEN
      is_determined = .true.

      ! b needs to have the OPOSITE sign of a. Sign of c depends on a & b
      ! a = [+] --> b = [-]  , c = [+] if |a| > |b| and vice versa
      ! a < [-] --> b = [+]  , c = [-] if |a| > |b| and vice versa

      if (a >= Zero) then
        b = - abs_b
        if (abs_a < abs_b) c = - abs_c
      else ! a < Zero
        if (abs_a > abs_b) c = - abs_c
      endif
    endif
  endif

  if(.not. is_determined) then
    a = abs_a
    b = abs_b
    c = abs_c
  endif

end subroutine assign_abc_signs
!======================================================================


!======================================================================
  subroutine calc_prin_val(iAtom,A_tens,contrib)
    real(kind=wp), intent(in)    :: A_tens(3,3)
    integer(kind=iwp),intent(in) :: contrib
    real(kind=wp)                :: X(3,3),EVR(3), EVI(3), prvl, fnorm_diag, fnorm_off_diag, tmpmat(3,3), dnrm2_
    real(kind=wp),parameter      :: to_au = -gElectron * beta_e * beta_n

    integer(kind=iwp)            :: iAxis
    integer(kind=iwp),intent(in) :: iAtom

! VARIABLE DESCRIPTION:
!   fnorm_diag:       Frobenius norm of diag(sm_a)
!   fnorm_off_diag:   Frobenius norm of off-diag(sm_a)

    tmpmat(:,:) = A_tens(:,:)
    X(:,:) = Zero
    EVR(:) = Zero
    EVI(:) = Zero
    call XEIGEN(1,3,3,tmpmat,EVR,EVI,X,IERR)

    ! Calcate diagonal elements norm
    fnorm_diag = sqrt(tmpmat(1,1)**2 + tmpmat(2,2)**2 + tmpmat(3,3)**2)
    ! Calcate off-diagonal elements norm
    tmpmat(1,1) = Zero
    tmpmat(2,2) = Zero
    tmpmat(3,3) = Zero
    fnorm_off_diag = dnrm2_(9, tmpmat, 1)

    if (fnorm_off_diag/fnorm_diag > 0.05_wp) then
      write(6,*) ''
      write(6,*) ''
      write(6,'(3X,A48,F6.4,A5)') 'WARNING: Relative Frobenius diag/off-diag norm: ', 100.0_wp * fnorm_off_diag/fnorm_diag, " > 5%"
    end if

    ! principal values (sqrt of diagonal elements)
    prin_vals(iACalc,contrib,:) = SQRT(EVR(:))

    ! Print A-tensor, the principal axes and eigenvalues
    write(6,*) ''
    write(6,*) ''
    write(6,'(3X,A96)') REPEAT('-',96)
    write(6,'(3X,A10,A6,A22,A7)') '>>> ATOM: ', LAtomLbl(iAtom), 'HYPERFINE COUPLING :: ', contrib_lab(contrib)
    write(6,'(3X,A96)') REPEAT('-',96)
    write(6,'(14X,A17,29X,A14,14X,A11)') "A-tensor (A=aa^T)","principal axes", "eigenvalues"
    write(6,'(12x,3(A1,12x),3x,3(A1,12x))') xyz(1:3),xyz(1:3)
    do iAxis = 1, 3
      write(6,'(3X,A1,3(1x,ES12.3),3x,3(1x,ES12.3),2x,ES12.3)') xyz(iAxis), A_tens(iAxis,1:3), X(iAxis,1:3), EVR(iAxis)
    end do

    ! Print absolute principal values without signs
    write(6,*) ''
    write(6,'(20X,A31,A7)') 'ABS. PRINCIPAL VALUES [+/-] :: ', contrib_lab(contrib)
    write(6,'(12X,A52)') repeat('.',52)
    write(6,'(20X,A12,12X,A4,11X,A5)')  "sqrt(eigval)", "(au)", "(MHz)"
    write(6,'(12X,A52)') repeat('.',52)
    do iAxis = 1, 3
      prvl = prin_vals(iACalc,contrib,iAxis)
      write(6,'(12X,A2,A1,A1,3X,E13.6,3x,E13.6,3x,E13.6)') 'A_', xyz(iAxis),xyz(iAxis), &
      prvl,  prvl*to_au,  prvl*con_to_MHz*NucGFac(iAtom)
    end do

  end subroutine calc_prin_val
!======================================================================


!======================================================================
subroutine update_h_HFC_RMS(iAtom, h_HFC)
  integer(kind=iwp),intent(in)  :: iAtom
  complex(kind=wp), intent(in)  :: h_HFC(3,NSS,NSS)
  real(kind=wp)                 :: I_sq, NSpin_I, fac

  ! VERSION 1: PARTLY VECTORIZED
  ! h_rms_nuc =  abs(h_HFC(1,:,:))**2 + abs(h_HFC(2,:,:))**2 + abs(h_HFC(3,:,:))**2

  ! VERSION 2: LEGACY COMPILERS
  ! h_rms_nuc(:,:) = Zero
  ! do i = 1, 3
  !   h_rms_nuc = h_rms_nuc + real(h_HFC(i,:,:))**2 + aimag(h_HFC(i,:,:))**2
  ! end do

  ! OPTIMIZED VERSION :: FULLY VECTORIZED
  h_rms_nuc(:,:) =  sum(abs(h_HFC(1:3,:,:))**2, dim=1)

  ! Scale with factor || I_sq ||^2
  NSpin_I = NucSpin(iAtom)
  I_sq = NSpin_I * (NSpin_I + One) * (Two * NSpin_I + One) / Three
  fac  = I_sq * NucGFac(iAtom)**2
  h_rms_nuc(:,:) = fac * h_rms_nuc(:,:)

  ! Update h_RMS for this nuclei
  h_hfc_rms(:,:) = h_hfc_rms(:,:) + h_rms_nuc(:,:)

  if (iAtom == NAtoms) then
    ! Note: con_to_MHz = gElectron*beta_e*beta_n*auToHz*1e-6_wp
    fac = -gElectron * beta_e * beta_n
    h_hfc_rms(:,:) = fac * sqrt(h_hfc_rms(:,:))
  endif
end subroutine
!======================================================================


!======================================================================
subroutine calc_h_PSO(iAtom,PROP)
  real(kind=wp), intent(in)       :: PROP(NSTATE,NSTATE,NPROP)
  real(kind=wp), allocatable      :: Im_h_PSO(:,:,:)
  integer(kind=iwp),intent(in)    :: iAtom
  integer(kind=iwp)               :: u

  call mma_allocate(Im_h_PSO,3,NSS,NSS, Label="Im_h_PSO")

  Im_h_PSO(:,:,:) = Zero

  do u = 1, 3
    call SMMAT(PROP,Im_h_PSO(u,:,:),NSS,PSO_idx(iAtom,u),u)
  enddo

  h_PSO(:,:,:) = cmplx(0.0_wp, Im_h_PSO(:,:,:), kind=wp)

  call mma_deallocate(Im_h_PSO)

  call to_cmpl_SO_states(h_PSO)

end subroutine
!======================================================================


!======================================================================
subroutine save_h_rms()
  implicit none
  integer(kind=iwp) LU, JSTA, ISS
  integer(kind=iwp), External:: IsFreeUnit

#ifdef _HDF5_
  call mh5_put_dset(wfn_h_hfc_rms,h_hfc_rms)
#endif

  Lu = 88
  Lu = IsFreeUnit(Lu)
  OPEN(UNIT=Lu,FILE='h_RMS.txt',STATUS='REPLACE')
  write(Lu,*) "NSS= ",NSS
  write(Lu,*) "#NROW NCOL REAL"
  do JSTA=1,NSS
    do ISS=1,NSS
    write(Lu,'(I6,1X,I6,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ', h_hfc_rms(ISS,JSTA)
    end do
  end do
  CLOSE(Lu)
end subroutine save_h_rms
!======================================================================


!======================================================================
subroutine calc_h_Zeeman(PROP)
  real(kind=wp), intent(in)   :: PROP(NSTATE,NSTATE,NPROP)
  real(kind=wp), allocatable  :: Re_h_Zeeman(:,:,:) , Im_h_Zeeman(:,:,:),Angmom(:,:,:)
  integer(kind=iwp)           :: IAMX, IAMY, IAMZ, IPROP

  call mma_allocate(Re_h_Zeeman, 3, NSS, NSS, Label='Re_h_Zeeman')
  call mma_allocate(Im_h_Zeeman, 3, NSS, NSS, Label='Im_h_Zeeman')
  call mma_allocate(Angmom,      3, NSS, NSS, Label='Angmom')
  call mma_allocate(h_Zeeman,3,NSS,NSS,Label='h_Zeeman')

  h_Zeeman(:,:,:) = cmplx(Zero,Zero, kind=wp)
  Re_h_Zeeman(:,:,:) = Zero
  Im_h_Zeeman(:,:,:) = Zero

  DO IPROP=1,NPROP
    IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
      IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
      IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
      IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
    END IF
  END DO

  Angmom(:,:,:) = Zero

  CALL SMMAT(PROP,Angmom(1,:,:),NSS,IAMX,0)
  CALL SMMAT(PROP,Angmom(2,:,:),NSS,IAMY,0)
  CALL SMMAT(PROP,Angmom(3,:,:),NSS,IAMZ,0)

  CALL SMMAT(PROP,Re_h_Zeeman(1,:,:),NSS,0,1)
  CALL SMMAT(PROP,Im_h_Zeeman(2,:,:),NSS,0,2)
  CALL SMMAT(PROP,Re_h_Zeeman(3,:,:),NSS,0,3)

  Re_h_Zeeman(1,:,:) = -gElectron * Re_h_Zeeman(1,:,:)
  Im_h_Zeeman(2,:,:) = -gElectron * Im_h_Zeeman(2,:,:)
  Re_h_Zeeman(3,:,:) = -gElectron * Re_h_Zeeman(3,:,:)
  Im_h_Zeeman(:,:,:) = Im_h_Zeeman(:,:,:) + Angmom(:,:,:)

  h_Zeeman(:,:,:) = cmplx(Re_h_Zeeman(:,:,:),Im_h_Zeeman(:,:,:), kind=wp)
  call to_cmpl_SO_states(h_Zeeman)

  h_Zeeman(:,:,:) = h_Zeeman(:,:,:) / Two

  call mma_deallocate(Re_h_Zeeman)
  call mma_deallocate(Im_h_Zeeman)
  call mma_deallocate(Angmom)
  ! h_Zeeman needs to deallocated at the end [after pNMR calc. finishes]

endsubroutine calc_h_Zeeman
!======================================================================


!======================================================================
subroutine calc_pNMR_Tensor(iAtom,h_HFC,contrib)
! INPUT VARIABLES
  integer(kind=iwp), intent(in)   :: iAtom
  complex(kind=wp), intent(in)    ::  h_HFC(3,NSS,NSS)
  integer(kind=iwp), intent(in)   ::  contrib
! LOCAL VARIABLES
  real(kind=wp)                   ::  fac
  integer(kind=iwp)               ::  lmb_a, lmb_ap, lmb

  integer(kind=iwp)                ::  u, w, ISS, iT

! VARIABLES DESCRIPTION
! REF: DOI: 10.1021/acs.jctc.6b00462   [eq. 1]
! lmb                 : "lambda" is the group of all degenerate states.
! lmb_a & lmb_ap      : states between lmb_a - lmb_ap are degenerate
! Z_HFC_interaction   : <i |h_Zeeman| j> <j |h_HFC| i>
! Z_HFC_over_dE       : Z_HFC_interaction / delta_E(ij)
! LR_tens             : Linear Response tensor [first term,  eq 1]
! C_tens              : Curie tensor           [second term, eq 1]


  do u = 1, 3
    do w = 1, 3
      Z_HFC_int_oper(u,w,:,:) = real(h_Zeeman(u,:,:) * conjg(h_HFC(w,:,:)), kind=wp)
      Z_HFC_over_dE(u,w,:,:)  = Z_HFC_int_oper(u,w,:,:) * dE_inv(:,:)
    end do
  end do

  LR_tens(:,:,:) = Zero
  C_tens(:,:,:)  = Zero
!--------------------Calculate numerators for pNMR and Curie term----------------------
  do iT = 1, NTP
    ! Curie term contribution [degenerate states]
    do ISS = 1, NSS

      ! Lambda (lmb) is a group of degenereate states DOI: 10.1021/acs.jctc.6b00462   [eq. 1]
      ! All degenerate states have the SAME lambda
      lmb     = degen_group(ISS)
      lmb_a   = degen_start_idx(lmb)
      lmb_ap  = degen_end_idx(lmb)

      ! Vectorized Curie tensor is trivial, loop between lmb_a and lmb_ap
      C_tens(iT,:,:) = C_tens(iT,:,:) +  pBoltz(iT,ISS) * sum(Z_HFC_int_oper(:,:,ISS,lmb_a:lmb_ap), dim=3)

      ! Non-degenerate states will contribute to Linear Response term
      ! Two if-logic has been used to avoid empty array (lmb_a == 1 means 1:lmb_a-1 == 1:0)
      if(lmb_a > 1)    LR_tens(iT,:,:) = LR_tens(iT,:,:) + pBoltz(iT,ISS)* sum(Z_HFC_over_dE(:,:,ISS,1:lmb_a-1), dim=3)
      if(lmb_ap < NSS) LR_tens(iT,:,:) = LR_tens(iT,:,:) + pBoltz(iT,ISS)* sum(Z_HFC_over_dE(:,:,ISS,lmb_ap+1:NSS), dim=3)

    end do
    fac = One / (kBoltzman_in_cm * Temp_in_K(iT))
    C_tens(iT,:,:) = fac * C_tens(iT,:,:)
  enddo !Temperature

  C_tens(:,:,:) =  to_ppm * C_tens(:,:,:)
  fac = to_ppm * Two
  LR_tens(:,:,:) = fac * LR_tens(:,:,:)
!-----------------------------------------------
  Curie_ChemShift(ipNMR_Calc,contrib,:)  = -(C_tens(:,1,1)+C_tens(:,2,2)+C_tens(:,3,3))/Three
  LinRes_ChemShift(ipNMR_Calc,contrib,:) = -(LR_tens(:,1,1)+LR_tens(:,2,2)+LR_tens(:,3,3))/Three


  write(6,*) ''
  write(6,'(3X,A66)') repeat('-',66)
  write(6,'(3X,A10,A6,A19,A7)') '>>> ATOM: ', LAtomLbl(iAtom), 'LINEAR RESPONSE :: ', contrib_lab(contrib)
  write(6,'(3X,A66)') repeat('-',66)
  write(6,'(3X,A8,16X,A7,15X,A11)') 'Temp (K)', 'Tensor', 'Shift (ppm)'
  do iT=1,NTP
    call print_pNMR_tens(Temp_in_K(iT), LR_tens(iT,:,:))
  end do
  write(6,'(3X,A66)') repeat('-',66)
  write(6,'(3X,A10,A6,A9,A7)') '>>> ATOM: ', LAtomLbl(iAtom), 'CURIE :: ', contrib_lab(contrib)
  write(6,'(3X,A66)') repeat('-',66)
  write(6,'(3X,A8,16X,A7,15X,A11)') 'Temp (K)', 'Tensor', 'Chemical Shift (ppm)'
  do iT=1,NTP
    call print_pNMR_tens(Temp_in_K(iT), C_tens(iT,:,:))
  end do
  write(6,*) ''
  write(6,*) ''
  write(6,*) ''

end subroutine
!======================================================================


!======================================================================
  subroutine print_pNMR_tens(temp, tensor)
    real(kind=wp), intent(in) :: temp, tensor(3,3)
    real(kind=wp)             :: shift
    integer(kind=iwp)         :: u
    shift = -(tensor(1,1)+tensor(2,2)+tensor(3,3))/Three
    write(6,'(2X,F6.1,3(2x,ES12.3),2X,F15.3)') temp, tensor(1,1:3), shift
    do u = 2, 3
      write(6,'(8X,3(2x,ES12.3))') tensor(u,1:3)
    end do
    write(6,*) ""
  end subroutine print_pNMR_tens
!======================================================================


!======================================================================
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
    endif

    if (allocated(pNMR_req))  then
      call mma_deallocate(Curie_ChemShift)
      call mma_deallocate(LinRes_ChemShift)
      call mma_deallocate(pNMR_req)
      call mma_deallocate(pBoltz)
      call mma_deallocate(Temp_in_K)

      call mma_deallocate(Z_HFC_int_oper)
      call mma_deallocate(Z_HFC_over_dE)
      call mma_deallocate(h_Zeeman)

      call mma_deallocate(dE_inv)
      call mma_deallocate(LR_tens)
      call mma_deallocate(C_tens)
      call mma_deallocate(degen_group)
    endif


    if (allocated(NucSpin))   call mma_deallocate(NucSpin)
    if (allocated(NucGFac))   call mma_deallocate(NucGFac)
    if (allocated(LCSTATES))  call mma_deallocate(LCSTATES)

    if (allocated(ASD_idx))   call mma_deallocate(ASD_idx)
    if (allocated(PSO_idx))   call mma_deallocate(PSO_idx)

    call mma_deallocate(degen_start_idx, safe='*')
    call mma_deallocate(degen_end_idx, safe='*')

    if (allocated(Atens_Req)) then
      call mma_deallocate(Atens_Req)
      call mma_deallocate(prin_vals)
      call mma_deallocate(signs_resolved)
    endif
  end subroutine
!======================================================================


end module

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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!               2026, Arpan Choudhury                                  *
!***********************************************************************

module mspdft

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
private

character(len=8) :: mspdftmethod = ' Unknown'
real(kind=wp), allocatable :: heff(:,:)
real(kind=wp), allocatable :: w_rf(:,:)
real(kind=wp), allocatable :: dao_pair(:,:), dao_core(:)
integer(kind=iwp), allocatable :: pair_roots(:,:)
real(kind=wp) :: e_rf_core
real(kind=wp) :: e_rf_pot_scf, e_rf_self_scf

! RUNFILE records that hold the slow component of a non-equilibrium reaction field:
character(len=16), parameter :: slow_ras(4) = ['RCTFLD          ','W_or_el         ','W_or_Inf        ','E_0_NN          '], &
                                slow_lpdft(4) = ['LPDFT RCTFLD    ','LPDFT W_or_el   ','LPDFT W_or_Inf  ','LPDFT E_0_NN    ']

public :: heff, w_rf, e_rf_core, dao_pair, pair_roots, mspdft_finalize, mspdft_init, mspdftmethod, dao_core, &
          lpdft_rf_step, lpdft_rf_scf, rf_runold_push, rf_runold_pop
contains

!> @brief Initializes mspdft calculation
!>
!>   Loads effective Hamiltonian from disk and initializes gradient variables as needed
!>
!> @author Matthew R. Hennefarth
subroutine mspdft_init()

  use rasscf_global, only: lroots
  use mcpdft_input, only: mcpdft_options
  use mspdftgrad, only: mspdftgrad_init
  use Constants, only: Zero

  character(len=18) :: matrix_info

  call mma_allocate(heff,lroots,lroots,label='heff')
  call mma_allocate(w_rf,lroots,lroots,label='w_rf')
  w_rf(:,:) = Zero
  e_rf_core = Zero

  if (mcpdft_options%lpdft) then
    heff(:,:) = Zero
    mspdftmethod = 'LPDFT'
  else
    call readmat2('ROT_HAM',matrix_info,heff,size(heff,1),size(heff,2),7,len(matrix_info),'T')
    call determine_method(matrix_info)
  end if

  if (mcpdft_options%grad) call mspdftgrad_init()

end subroutine

!> @brief determines method by which effective Hamiltonian is constructed
!>
!>   Sets mspdftmethod and and writes to stdout
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] matrix_info string
subroutine determine_method(matrix_info)

  use PrintLevel, only: VERBOSE
  use mcpdft_output, only: iprloc

  character(len=*), intent(in) :: matrix_info
  integer(kind=iwp) :: print_level
  character(len=len(mspdftmethod)) :: buffer

  print_level = iPrLoc(1)
  buffer = trim(adjustl(matrix_info))

  select case (buffer)
    case ('XMS-PDFT','CMS-PDFT','FMS-PDFT','VMS-PDFT')
      mspdftmethod = buffer
      if (print_level > VERBOSE) write(u6,'(6X,A,A)') 'The MS-PDFT method is ',trim(mspdftmethod)
    case default
      if (print_level > VERBOSE) write(u6,'(6X,A)') 'The MS-PDFT calculation is based on a user-supplied rotation matrix'
  end select

end subroutine determine_method

!> @brief Performs the final parts of a MS-PDFT calculation
!>
!>   1. diagonalize heff
!>   2. print out the final results
!>   3. save to jobiph if requested
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in]   nroots number of roots, or dimension of heff
subroutine mspdft_finalize(nroots)

  use mcpdft_output, only: iprloc
  use mspdft_util, only: print_effective_ham, print_final_energies, print_mspdft_vectors
  use mcpdft_input, only: mcpdft_options
  use write_pdft_job, only: writejob
  use PrintLevel, only: TERSE, USUAL
  use mspdftgrad, only: mspdftgrad_free
  use general_data, only: ntot1
  use Constants, only: Zero, Two

  integer(kind=iwp), intent(in) :: nroots
  integer(kind=iwp) :: info, dim_scratch, iprlev, iroot
  real(kind=wp) :: e_mspdft(nroots), si_pdft(nroots,nroots), wgronk(1)
  real(kind=wp) :: wdiag(nroots), e_rf_pot, e_rf_self
  logical(kind=iwp) :: runold
  character(len=120) :: Line
  real(kind=wp), allocatable :: scratch(:)

  iprlev = iprloc(1)

  if (iprlev >= USUAL) then
    write(u6,*)
    write(Line,'(6X,2A)') trim(MSPDFTMethod),' FINAL RESULTS'
    call CollapseOutput(1,Line)
    write(u6,'(6X,A)') repeat('-',len_trim(Line)-3)
    write(u6,*)

    if (mcpdft_options%otfnal%is_hybrid()) then
      write(u6,'(6X,3A)') 'Hybrid ',trim(MSPDFTMethod),' Effective Hamiltonian'
    else
      write(u6,'(6X,2A)') trim(mspdftmethod),' Effective Hamiltonian'
    end if
    call print_effective_ham(heff,nroots,10)
  end if

  ! Now we diagonalize the final MS-PDFT H-matrix
  ! Eigenvectors will be stored in heff.
  ! Eigenvalues will be stored in e_mspdft

  ! Since the dsyev_ call will override the heff with the orthonormal
  ! eigenvectors, I am going to set the new variable now.
  si_pdft = heff

  ! This first call is to get how big of a scratch we need
  call dsyev_('V','U',nroots,si_pdft,nroots,e_mspdft,wgronk,-1,info)

  dim_scratch = int(wgronk(1))
  call mma_allocate(scratch,dim_scratch,label='XScratch')
  ! Now we actually do the diagonalization
  call dsyev_('V','U',nroots,si_pdft,nroots,e_mspdft,scratch,dim_scratch,info)
  call mma_deallocate(scratch)

  if (mcpdft_options%lpdft .and. mcpdft_options%rf_on) then
    if (mcpdft_options%rf_scf) then
      e_rf_pot = e_rf_pot_scf
      e_rf_self = e_rf_self_scf
    else
      call rf_runold_push(runold)
      call get_dscalar('RF Pot Energy',e_rf_pot)
      call get_dscalar('RF Self Energy',e_rf_self)
      call rf_runold_pop(runold)
    end if
    call rf_rot_diag(nroots,w_rf,si_pdft,wdiag)
    do iroot=1,nroots
      e_mspdft(iroot) = e_mspdft(iroot)-wdiag(iroot)-e_rf_core+e_rf_pot+e_rf_self
    end do
  end if

  if (iprlev >= TERSE) call print_final_energies(e_mspdft,nroots)

  ! Update information on the runfile for possible gradient calculations.
  call put_dArray('Last energies',e_mspdft,nroots)
  if (mcpdft_options%grad) call Put_dScalar('Last energy',e_mspdft(mcpdft_options%rlxroot))

  ! Add info the checkfile for testing!
  call Add_Info('MSPDFTE',e_mspdft,nroots,8)

  if (iprlev >= USUAL) then
    if (mcpdft_options%otfnal%is_hybrid()) then
      write(u6,'(6X,3A)') 'Hybrid ',trim(MSPDFTMethod),' Eigenvectors:'
    else
      write(u6,'(6X,2A)') trim(MSPDFTMethod),' Eigenvectors:'
    end if
    call print_mspdft_vectors(si_pdft,nroots)
  end if

  ! Added by Chen to write energies and states of MS-PDFT into JOBIPH
  if (mcpdft_options%wjob) call writejob(e_mspdft,nroots,si_pdft=si_pdft)

  if (iprlev >= USUAL) then
    call CollapseOutput(0,Line)
    write(u6,*)
  end if

  if (mcpdft_options%grad) then
    call put_lscalar('CalcNAC_Opt     ',.false.)
    call put_lscalar('MECI_via_SLAPAF',.false.)
    if (mcpdft_options%nac) then
      call mspdftgrad_misc(si_pdft,mcpdft_options%nac_states)
    else
      call put_iscalar('Relax CASSCF root',mcpdft_options%rlxroot)
      call mspdftgrad_misc(si_pdft,[mcpdft_options%rlxroot,mcpdft_options%rlxroot])
    end if
  end if

  ! Deallocate here!
  call mma_deallocate(heff)
  call mma_deallocate(w_rf)
  if (allocated(dao_pair)) call mma_deallocate(dao_pair)
  if (allocated(pair_roots)) call mma_deallocate(pair_roots)
  if (allocated(dao_core)) call mma_deallocate(dao_core)
  if (mcpdft_options%grad) call mspdftgrad_free()

end subroutine mspdft_finalize

!> @brief One pass of the L-PDFT reaction-field loop (only called when rf_scf).
!>
!> @param[in]  nroots       dimension of heff
!> @param[out] e_states     RF-corrected L-PDFT energies
!> @param[out] d_rf         folded AO density (ntot1) of the RF root
!> @param[out] k            the RF root
subroutine lpdft_rf_step(nroots,e_states,d_rf,k)

  use general_data, only: ntot1
  use Constants, only: Two
  use mcpdft_input, only: mcpdft_options

  integer(kind=iwp), intent(in) :: nroots
  real(kind=wp), intent(out) :: e_states(nroots), d_rf(ntot1)
  integer(kind=iwp), intent(out) :: k
  integer(kind=iwp) :: info, lwork, iroot, iPair, ir, jr
  real(kind=wp) :: si(nroots,nroots), wdiag(nroots), wq(1), e_rf_pot, e_rf_self, s1
  real(kind=wp), allocatable :: work(:)

  si = heff
  call dsyev_('V','U',nroots,si,nroots,e_states,wq,-1,info)
  lwork = int(wq(1))
  call mma_allocate(work,lwork,label='RFWork')
  call dsyev_('V','U',nroots,si,nroots,e_states,work,lwork,info)
  call mma_deallocate(work)
  if (info /= 0) then
    call WarningMessage(2,'lpdft_rf_step: diagonalization of heff failed')
    call Abend()
  end if

  call get_dscalar('RF Pot Energy',e_rf_pot)
  call get_dscalar('RF Self Energy',e_rf_self)
  call rf_rot_diag(nroots,w_rf,si,wdiag)
  do iroot=1,nroots
    e_states(iroot) = e_states(iroot)-wdiag(iroot)-e_rf_core+e_rf_pot+e_rf_self
  end do

  k = mcpdft_options%rfroot
  if ((k <= 0) .or. (k > nroots)) then
    call WarningMessage(2,'LPDFT+RF: only a single RF root (1 <= RFRoot <= nroots) is supported')
    call Quit_OnUserError()
  end if

  d_rf(:) = dao_core(:)
  do iPair=1,size(dao_pair,2)
    ir = pair_roots(1,iPair)
    jr = pair_roots(2,iPair)
    s1 = si(ir,k)*si(jr,k)
    if (ir /= jr) s1 = Two*s1
    call daxpy_(ntot1,s1,dao_pair(:,iPair),1,d_rf,1)
  end do

end subroutine lpdft_rf_step

!> @brief Self-consistent reaction field for L-PDFT (only when rf_scf).

subroutine lpdft_rf_scf(nroots)

  use general_data, only: ntot1, ispin
  use rasscf_global, only: potnuc, exfac, weight
  use Constants, only: Zero
  use mcpdft_input, only: mcpdft_options
  use mcpdft_output, only: iPrLoc
  use PrintLevel, only: TERSE, USUAL


  integer(kind=iwp), intent(in) :: nroots
  integer(kind=iwp), parameter :: maxit = 50
  real(kind=wp), parameter :: thr = 1.0e-8_wp
  integer(kind=iwp) :: it, iroot, iPair, ir, jr, k, icharge
  real(kind=wp) :: heff0(nroots,nroots), e(nroots)
  real(kind=wp) :: repnuc, exf, de, e_sa, e_sa_old
  real(kind=wp), allocatable :: d_rf(:,:), v(:), h1d(:), g2d(:)
  logical(kind=iwp) :: conv, lFirst, lDff, lNonEq, lRFon, lDoDFT
  real(kind=wp), allocatable :: q_ras(:), q_lp(:)
  real(kind=wp) :: s_ras(3), s_lp(3)
  logical(kind=iwp) :: hq_ras, hs_ras(3), hq_lp, hs_lp(3)
  real(kind=wp), allocatable :: v_ras(:), c_ras(:)
  real(kind=wp) :: erfself_ras, erfpot_ras
  logical(kind=iwp) :: hv_ras, hc_ras, hself_ras, hpot_ras
  integer(kind=iwp) :: nv_ras, nc_ras
  character(len=4) :: dftfock
  integer(kind=iwp), external :: get_charge
  integer(kind=iwp) :: iprlev
  real(kind=wp), external :: dDot_

  call mma_allocate(d_rf,ntot1,2,label='D_RF_SCF')
  call mma_allocate(v,ntot1,label='V_RF_SCF')
  call mma_allocate(h1d,ntot1,label='H1_RF_SCF')
  call mma_allocate(g2d,ntot1,label='G2_RF_SCF')
  d_rf(:,:) = Zero

  heff0(:,:) = heff(:,:)-w_rf(:,:)
  do iroot=1,nroots
    heff0(iroot,iroot) = heff0(iroot,iroot)-e_rf_core
  end do

  lFirst = .true.; lDff = .false.; lRFon = .true.; lDoDFT = .false.
  lNonEq = mcpdft_options%noneq

  ! Keep the RASSCF reaction field untouched by L-PDFT, so that later modules
  ! (e.g. CASPT2 or RASSI with RFPERT) see the same records as without the
  ! L-PDFT run: save the field, its energies and the PCM charges here and
  ! restore them after the loop
  call Qpg_dArray('Reaction field',hv_ras,nv_ras)
  hv_ras = hv_ras .and. (nv_ras > 0)
  if (hv_ras) then
    call mma_allocate(v_ras,nv_ras,label='V_RAS')
    call Get_dArray('Reaction field',v_ras,nv_ras)
  end if
  call Qpg_dArray('PCM Charges',hc_ras,nc_ras)
  hc_ras = hc_ras .and. (nc_ras > 0)
  if (hc_ras) then
    call mma_allocate(c_ras,nc_ras,label='C_RAS')
    call Get_dArray('PCM Charges',c_ras,nc_ras)
  end if
  erfself_ras = Zero
  erfpot_ras = Zero
  call Qpg_dScalar('RF Self Energy',hself_ras)
  if (hself_ras) call Get_dScalar('RF Self Energy',erfself_ras)
  call Qpg_dScalar('RF Pot Energy',hpot_ras)
  if (hpot_ras) call Get_dScalar('RF Pot Energy',erfpot_ras)

  ! Same for the slow-component records
  call get_rf_slow(slow_ras,q_ras,hq_ras,s_ras,hs_ras)
  if (lNonEq) then
    ! Non-equilibrium L-PDFT uses only the slow component of a previous
    ! equilibrium L-PDFT calculation
    call get_rf_slow(slow_lpdft,q_lp,hq_lp,s_lp,hs_lp)
    if (.not. hq_lp) then
      call WarningMessage(2,'NONEQ: no slow charges from a previous equilibrium L-PDFT calculation')
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' Run an equilibrium L-PDFT with a  '
      write(u6,*) ' reaction field first              '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    ! DrvPCM/RctFld read the slow component from the standard records
    call put_rf_slow(slow_ras,q_lp,hq_lp,s_lp,hs_lp)
  end if
  dftfock = 'SCF '
  conv = .false.
  e_sa_old = Zero

  iprlev = iPrLoc(1)
  if (iprlev >= USUAL) then
    write(u6,*)
    write(u6,'(41X,A)') 'LPDFT RF iterations: Energy and convergence statistics'
    write(u6,'(41X,A)') repeat('-',55)
    write(u6,*)
    write(u6,'(6X,A4,A22,A16)') 'Iter','Energy','Energy change'
  end if
  do it=1,maxit
    heff(:,:) = heff0(:,:)+w_rf(:,:)
    do iroot=1,nroots
      heff(iroot,iroot) = heff(iroot,iroot)+e_rf_core
    end do

    call lpdft_rf_step(nroots,e,d_rf(:,1),k)

    e_sa = sum(weight(1:nroots)*e(1:nroots))
    if (it == 1) then
      de = Zero
    else
      de = e_sa-e_sa_old
    end if
    if (iprlev >= USUAL) write(u6,'(6X,I4,F22.8,ES16.2)') it,e_sa,de
    if ((it > 1) .and. (abs(de) < thr)) then
      conv = .true.
      exit
    end if
    e_sa_old = e_sa

    ! new field from the RF-root L-PDFT density
    h1d(:) = Zero
    g2d(:) = Zero
    repnuc = potnuc
    exf = exfac
    icharge = get_charge()
    call DrvXV(h1d,g2d,d_rf,repnuc,ntot1,lFirst,lDff,lNonEq,lRFon,'SCF',exf,icharge,ispin,dftfock,lDoDFT)

    call get_dArray('Reaction field',v,ntot1)
    e_rf_core = dDot_(ntot1,v,1,dao_core,1)
    do iPair=1,size(dao_pair,2)
      ir = pair_roots(1,iPair)
      jr = pair_roots(2,iPair)
      w_rf(ir,jr) = dDot_(ntot1,v,1,dao_pair(:,iPair),1)
      w_rf(jr,ir) = w_rf(ir,jr)
    end do
  end do

  if (conv) then
    if (iprlev >= TERSE) write(u6,'(6X,A,I3,A)') 'Convergence after',it,' iterations'
  else
    if (iprlev >= TERSE) write(u6,'(6X,A,I3,A)') 'No convergence after',maxit,' iterations'
    call WarningMessage(1,'L-PDFT reaction field did not converge')
  end if

  if (.not. lNonEq) then
    ! Save the slow component of the L-PDFT RF root for a later NONEQ L-PDFT run
    call get_rf_slow(slow_ras,q_lp,hq_lp,s_lp,hs_lp)
    call put_rf_slow(slow_lpdft,q_lp,hq_lp,s_lp,hs_lp)
  end if
  call Get_dScalar('RF Pot Energy',e_rf_pot_scf)
  call Get_dScalar('RF Self Energy',e_rf_self_scf)

  ! Restore the RASSCF records
  call put_rf_slow(slow_ras,q_ras,hq_ras,s_ras,hs_ras)
  if (hv_ras) call Put_dArray('Reaction field',v_ras,nv_ras)
  if (hc_ras) call Put_dArray('PCM Charges',c_ras,nc_ras)
  if (hself_ras) call Put_dScalar('RF Self Energy',erfself_ras)
  if (hpot_ras) call Put_dScalar('RF Pot Energy',erfpot_ras)
  if (allocated(q_lp)) call mma_deallocate(q_lp)
  if (allocated(q_ras)) call mma_deallocate(q_ras)
  if (allocated(v_ras)) call mma_deallocate(v_ras)
  if (allocated(c_ras)) call mma_deallocate(c_ras)

  call mma_deallocate(g2d)
  call mma_deallocate(h1d)
  call mma_deallocate(v)
  call mma_deallocate(d_rf)

end subroutine lpdft_rf_scf

!> @brief Read the slow-component records from RUNFILE.
subroutine get_rf_slow(lab,q,has_q,s,has_s)

  use Constants, only: Zero

  character(len=16), intent(in) :: lab(4)
  real(kind=wp), allocatable, intent(inout) :: q(:)
  logical(kind=iwp), intent(out) :: has_q, has_s(3)
  real(kind=wp), intent(out) :: s(3)
  integer(kind=iwp) :: i, n

  if (allocated(q)) call mma_deallocate(q)
  call Qpg_dArray(lab(1),has_q,n)
  has_q = has_q .and. (n > 0)
  if (has_q) then
    call mma_allocate(q,n,label='RF_SLOW')
    call Get_dArray(lab(1),q,n)
  else
    call mma_allocate(q,0,label='RF_SLOW')
  end if
  s(:) = Zero
  do i=1,3
    call Qpg_dScalar(lab(i+1),has_s(i))
    if (has_s(i)) call Get_dScalar(lab(i+1),s(i))
  end do

end subroutine get_rf_slow

!> @brief Write the slow-component records that exist under the labels in lab.
subroutine put_rf_slow(lab,q,has_q,s,has_s)

  character(len=16), intent(in) :: lab(4)
  real(kind=wp), intent(in) :: q(:), s(3)
  logical(kind=iwp), intent(in) :: has_q, has_s(3)
  integer(kind=iwp) :: i

  if (has_q) call Put_dArray(lab(1),q,size(q))
  do i=1,3
    if (has_s(i)) call Put_dScalar(lab(i+1),s(i))
  end do

end subroutine put_rf_slow

!> @brief With RFPERT, read the fixed reaction field from RUNOLD if that file
!>   exists, as RASSCF, CASPT2 and MC-PDFT do; otherwise (and for the
!>   self-consistent L-PDFT field, which is written to RUNFILE) nothing changes.
subroutine rf_runold_push(switched)

  use mcpdft_input, only: mcpdft_options

  logical(kind=iwp), intent(out) :: switched

  switched = .false.
  if (mcpdft_options%rfpert) call f_Inquire('RUNOLD',switched)
  if (switched) call NameRun('RUNOLD')

end subroutine rf_runold_push

!> @brief Switch back to RUNFILE after rf_runold_push.
subroutine rf_runold_pop(switched)

  logical(kind=iwp), intent(in) :: switched

  if (switched) call NameRun('#Pop')

end subroutine rf_runold_pop

!> @brief Diagonal of S^T W S, the interaction of each rotated state with the reaction field.
subroutine rf_rot_diag(n,w,s,d)

  use Constants, only: Zero, One

  integer(kind=iwp), intent(in) :: n
  real(kind=wp), intent(in) :: w(n,n), s(n,n)
  real(kind=wp), intent(out) :: d(n)
  integer(kind=iwp) :: i
  real(kind=wp) :: ws(n,n)
  real(kind=wp), external :: dDot_

  call dgemm_('N','N',n,n,n,One,w,n,s,n,Zero,ws,n)
  do i=1,n
    d(i) = dDot_(n,s(:,i),1,ws(:,i),1)
  end do

end subroutine rf_rot_diag

end module mspdft

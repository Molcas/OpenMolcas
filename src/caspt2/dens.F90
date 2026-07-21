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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DENS(IVEC,NDMAT,NSTATE,DMAT,UEFF,U0)

use Index_Functions, only: iTri, nTri_Elem
use CHOVEC_IO, only: nvloc_chobatch
use PrintLevel, only: DEBUG, VERBOSE
use EQSOLV, only: IVECC, IVECC2, IVECR, IVECW, IVECX
use ChoCASPT2, only: iALGO, MaxVec_PT2
use sguga_states, only: SGS
use caspt2_global, only: CLag, CLagFull, CMOPT2, DMIX, do_csf, do_grad, DPT2_AO_tot, DPT2_tot, DPT2C_AO_tot, DPT2C_tot, &
                         DPT2Canti_tot, DREF, FIFA, FIFA_all, FIMO, FIMO_all, IDCIEX, IDTCEX, if_invar, if_invaria, if_SSDM, &
                         imag_shift, iPrGlb, iRoot1, iRoot2, jStLag, NDREF, nOLag, OLag, OMGDER, real_shift, sigma_p_epsilon, &
                         SLag, TORB, Weight
use caspt2_module, only: DENORM, HZERO, IfChol, IFDENS, IFDW, IFMSCOUP, IFSADREF, iRlxRoot, JSTATE, MAXIT, NAES, NASH, NASHT, &
                         NBAS, NBAST, NBSQT, NCONF, NFROT, NISH, NORB, NOSQT, NRAS1T, NRAS2T, NRAS3T, NROOTS, NSYM, ORBIN, ZETA
use BDerNEV, only: BDerNEV_initial, BDerNEV_final1, BDerNEV_final2
use SC_NEVPT2, only: SC_prop
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
use caspt2_global, only: nCLag
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC, NDMAT, NSTATE
real(kind=wp), intent(inout) :: DMAT(NDMAT)
real(kind=wp), intent(in) :: UEFF(nState,nState), U0(nState,nState)
integer(kind=iwp) :: I, iBasI, iBasSq, iBasTr, ibk, IDM, IDMOFF, IDRF, IDSOFF, IDSUM, IFF, II, IP, IQ, iSQ, iState, iStLag, ISYM, &
                     IT, ITABS, iTR, ITTOT, IU, IUABS, IUTOT, J, jBasI, liBasSq, liBasTr, lT2AO, NA, nBasI, nch, NDPT, nDPTAO, NI, &
                     NLEV, NO, nOcc, nOrbI
real(kind=wp) :: CPE, CPTF0, CPTF10, CPUT, Scal, TIOE, TIOTF0, TIOTF10, val, WALLT, wgt, X
integer(kind=iwp), allocatable :: ISAV(:)
real(kind=wp), allocatable :: A_PT2(:), CI1(:), CLagT(:,:), DEPSA(:,:), DEPSA_diag(:), DI(:), DIA(:), DPT(:), DPT2(:), DPT2_AO(:), &
                              DPT2C_AO(:), DSUM(:), EigT(:,:), FPT2(:), FPT2_AO(:), FPT2C(:), FPT2C_AO(:), OMGT(:,:), RDMEIG(:,:), &
                              RDMSA(:,:), T2AO(:), Trf(:), VECROT(:), WRK1(:), WRK2(:)
real(kind=wp), allocatable, target :: DPT2Canti_(:), DPT2C(:)
real(kind=wp), pointer :: DPT2Canti(:)
integer(kind=iwp), parameter :: kstate=1

if (do_grad) then
  !! Set indices for densities and partial derivatives
  call mma_allocate(VECROT,nState,Label='VECROT')
  call GradPrep(nState,UEFF,VECROT)

  ! Compute total density matrix as symmetry-blocked array of
  ! triangular matrices in DMAT. Size of a triangular submatrix is
  !  nTri_Elem(NORB(ISYM)).
  NDPT = sum(NORB(1:nSym)**2)
  nDPTAO = sum(nBas(1:nSym)**2)
  ! shouldn't be necessary, is already done outside
  DMAT(1:NDMAT) = Zero
  ! First, put in the reference density matrix.
  IDMOFF = 0
  do ISYM=1,NSYM
    NI = NISH(ISYM)
    NA = NASH(ISYM)
    NO = NORB(ISYM)
    do II=1,NI
      IDM = IDMOFF+nTri_Elem(II)
      DMAT(IDM) = Two
    end do
    do IT=1,NA
      ITABS = NAES(ISYM)+IT
      ITTOT = NI+IT
      do IU=1,IT
        IUABS = NAES(ISYM)+IU
        IUTOT = NI+IU
        IDRF = iTri(ITABS,IUABS)
        IDM = IDMOFF+iTri(ITTOT,IUTOT)
        DMAT(IDM) = DREF(IDRF)
      end do
    end do
    IDMOFF = IDMOFF+nTri_Elem(NO)
  end do
  !write(u6,*) ' DENS. Initial DMAT:'
  !write(u6,'(1x,8f16.8)') (dmat(i),i=1,ndmat)
  ! Add the 1st and 2nd order density matrices:
  call mma_allocate(DPT,NDPT,Label='DPT')
  call mma_allocate(DSUM,NDPT,Label='DSUM')
  DPT(:) = Zero
  DSUM(:) = Zero

  !! Modify the solution (T; amplitude), if the real- or
  !! imaginary- shift is utilized. We need both the unmodified (T)
  !! and modified (T+\lambda) amplitudes. \lambda can be obtained
  !! by solving the CASPT2 equation, but it can alternatively
  !! obtained by a direct summation only if CASPT2-D.
  !! iVecX remains unchanged (iVecX = T)
  !! iVecR will be 2\lambda

  !! For MS-CASPT2, calling this subroutine is required.
  !! The lambda-equation is solved without iteration only when
  !! MS-CASPT2-D (shift?). Otherwise, solved iteratively.
  !! After this subroutine, iVecR has multi-state weighted (?)
  !! contributions.
  call TIMING(CPTF0,CPE,TIOTF0,TIOE)
  call CASPT2_Res(VECROT,nState)
  call TIMING(CPTF10,CPE,TIOTF10,TIOE)
  if (IPRGLB >= VERBOSE) then
    CPUT = CPTF10-CPTF0
    WALLT = TIOTF10-TIOTF0
    write(u6,'(a,2f10.2)') ' Lambda  : CPU/WALL TIME=',cput,wallt
  end if

  call TIMING(CPTF0,CPE,TIOTF0,TIOE)
  !! Diagonal part
  if (SC_prop) then
    call TRDNS2D(iVecW,iVecC,DPT,NDPT,VECROT(JSTATE))
  else
    call TRDNS2D(iVecX,iVecR,DPT,NDPT,VECROT(JSTATE))
  end if
  !! Remove the off-diagonal elements in inactive/secondary
  if (.not. if_invaria) call caspt2_grad_invaria1(NDPT,DPT)
  DSUM(:) = DSUM(:)+DPT(:)
  !write(u6,*) ' DPT after TRDNS2D.'
  !write(u6,'(1x,8f16.8)') (dpt(i),i=1,ndpt)
  !! Off-diagonal part, if full-CASPT2
  if (MAXIT /= 0) then
    !! off-diagonal are ignored for CASPT2-D
    DPT(:) = Zero
    call TRDNS2O(iVecX,iVecR,DPT,size(DPT),NDPT,VECROT(JSTATE))
    DSUM(:) = DSUM(:)+DPT(:)
  end if
  !write(u6,*) ' DPT after TRDNS2O.'
  !write(u6,'(1x,8f16.8)') (dpt(i),i=1,ndpt)
  call TIMING(CPTF10,CPE,TIOTF10,TIOE)
  if (IPRGLB >= VERBOSE) then
    CPUT = CPTF10-CPTF0
    WALLT = TIOTF10-TIOTF0
    write(u6,'(a,2f10.2)') ' TRDNS2DO: CPU/WALL TIME=',cput,wallt
  end if

  !! D^PT2 in MO
  call mma_allocate(DPT2,NBSQT,Label='DPT2')
  !! D^PT2(C) in MO
  call mma_allocate(DPT2C,NBSQT,Label='DPT2C')
  !! DPTAO1 (D^PT in AO, but not DPTA-01) couples with
  !! the CASSCF density (assume state-averaged) through ERIs.
  !! This density corresponds to the eigenvalue derivative.
  !! This is sometimes referred to as DPT2(AO) else where.
  call mma_allocate(DPT2_AO,NBSQT,Label='DPT2_AO')
  !! DPTAO2 couples with the inactive density.
  !! This density comes from derivative of the generalized
  !! Fock matrix (see for instance Eq. (24) in the 1990 paper).
  !! This is sometimes referred to as DPT2C(AO) else where.
  call mma_allocate(DPT2C_AO,NBSQT,Label='DPT2C_AO')
  !! DPTAO,DPTCAO,FPTAO,FPTCAO are in a block-squared form
  call mma_allocate(FPT2,NBSQT,Label='FPT2')
  call mma_allocate(FPT2C,NBSQT,Label='FPT2C')
  call mma_allocate(FPT2_AO,NBSQT,Label='FPT2_AO')
  call mma_allocate(FPT2C_AO,NBSQT,Label='FPT2C_AO')
  !! Transformation matrix
  call mma_allocate(Trf,NBSQT,Label='TRFMAT')
  nch = 0
  if (IfChol) nch = nvloc_chobatch(1)
  call mma_allocate(WRK1,max(nBasT**2,nch),Label='WRK1')
  call mma_allocate(WRK2,max(nBasT**2,nch),Label='WRK2')
  !! state-averaged density
  call mma_allocate(RDMSA,nAshT,nAshT,Label='RDMSA')
  !! Derivative of state-averaged density
  call mma_allocate(RDMEIG,nAshT,nAshT,Label='RDMEIG')
  NLEV = SGS(kstate)%NLEV
  if (nAshT /= SGS(kstate)%NLEV) then
    write(u6,*) 'Analytical gradients for nAshT /= SGS(kstate)%NLEV (GASPT2?) does not work'
    call abend()
  end if
  !write(u6,*) 'olag before'
  !call sqprt(olag,nbast)

  DPT2(:) = Zero
  DPT2C(:) = Zero
  DPT2_AO(:) = Zero
  DPT2C_AO(:) = Zero
  FPT2(:) = Zero
  FPT2C(:) = Zero
  FPT2_AO(:) = Zero
  FPT2C_AO(:) = Zero
  if (.not. IfChol) then
    FIMO_all(:) = Zero
    FIFA_all(:) = Zero
  end if
  RDMSA(:,:) = Zero
  RDMEIG(:,:) = Zero

  CLag(:,:) = Zero
  OLag(:) = Zero

  if ((nFroT /= 0) .or. (.not. if_invaria)) then
    call mma_allocate(DIA,NBSQT,Label='DIA')
    call mma_allocate(DI,NBSQT,Label='DI')
  else
    call mma_allocate(DIA,1,Label='DIA')
    call mma_allocate(DI,1,Label='DI')
  end if

  if (do_csf) then
    call mma_allocate(DPT2Canti_,NBSQT,Label='DPT2Canti')
    DPT2Canti_(:) = Zero
    DPT2Canti => DPT2Canti_
  else
    DPT2Canti => DPT2C
  end if

  !! DPT -> DPT2
  !! Note that DPT2 has the index of frozen orbitals.
  !! Note also that unrelaxed (w/o Z-vector) dipole moments with
  !! frozen orbitals must be wrong.
  !dpt(:) = Zero
  if ((nFroT == 0) .and. if_invaria) then
    DPT2(1:nOsqT) = DSUM(1:nOsqT)
  else
    call OLagFro0(NOSQT,NBSQT,DSUM,DPT2)
  end if

  !! Construct the transformation matrix
  !! It seems that we have to transform quasi-canonical
  !! to CASSCF orbitals. The forward transformation has been
  !! done in ORBCTL.
  !!   C(PT2) = C(CAS)*X    ->    C(CAS) = C(PT2)*X^T
  !!   -> L(CAS) = X*L(PT2)*X^T
  !! inactive and virtual orbitals are not affected.
  Trf(:) = Zero
  call CnstTrf(NBSQT,TOrb,Trf)
  !call sqprt(trf,nbast)

  !! Construct the density matrix used in the Fock operator
  if (IFSADREF) then
    WRK1(1:nDRef) = Zero
    do iState=1,nState
      WRK1(1:nDRef) = WRK1(1:nDRef)+Weight(iState)*DMix(1:nDRef,iState)
    end do
  else
    WRK1(1:nDRef) = DMix(1:nDRef,jState)
  end if
  call SQUARE(WRK1,RDMSA,1,nAshT,nAshT)
  !write(u6,*) 'state-averaged density matrix'
  !call sqprt(rdmsa,nasht)

  ! ----- Construct configuration Lagrangian -----

  !! For CI coefficient derivatives (CLag)
  !! Calculate the configuration Lagrangian
  !! This is done in the quasi-canonical basis
  call mma_allocate(DEPSA,nAshT,nAshT,Label='DEPSA')
  DEPSA(:,:) = Zero
  !! Derivative of off-diagonal H0 of <Psi1|H0|Psi1>
  if (MAXIT /= 0) then
    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    call SIGDER(iVecX,iVecR,VECROT(jState))
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    if (IPRGLB >= VERBOSE) then
      CPUT = CPTF10-CPTF0
      WALLT = TIOTF10-TIOTF0
      write(u6,'(a,2f10.2)') ' SIGDER  : CPU/WALL TIME=',cput,wallt
    end if
  end if
  IFF = 1
  if (HZero == 'DYALL') then
    IFF = 0
    call BDerNEV_initial()
  end if
  call CLagX(IFF,nConf,nRoots,nState,nAshT,CLag,DEPSA,VECROT)
  !call test3_dens(clag)
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) call GADGOP(DEPSA,nAshT**2,'+')
# endif
  if (HZERO == 'DYALL') call BDerNEV_final1(NBSQT,DPT2C)
  !write(u6,*) 'original depsa'
  !call sqprt(depsa,nasht)
  !write(u6,*) 'original depsa (sym)'
  do i=1,nasht
    do j=1,i-1
      val = (DEPSA(i,j)+DEPSA(j,i))*Half
      DEPSA(i,j) = val
      DEPSA(j,i) = val
    end do
  end do
  !call sqprt(depsa,nasht)

  if (NRAS1T+NRAS3T /= 0) then
    !! The density of the independent pairs (off-diagonal blocks)
    !! should be determined by solving Z-vector, so these blocks
    !! should be removed...?
    ! write(u6,*) 'removing DEPSA of off-diagonal blocks'
    ! write(u6,*) 'before'
    ! call sqprt(depsa,nasht)
    DEPSA(1:nRAS1T,nRAS1T+1:) = Zero
    DEPSA(nRAS1T+1:,1:nRAS1T) = Zero
    DEPSA(nRAS1T+1:nRAS1T+nRAS2T,nRAS1T+nRAS2T+1:) = Zero
    DEPSA(nRAS1T+nRAS2T+1:,nRAS1T+1:nRAS1T+nRAS2T) = Zero
    !write(u6,*) 'after'
    !call sqprt(depsa,nasht)
    if (IPRGLB >= DEBUG) write(u6,*) 'depsa (sym) after removing off-diagonal blocks'
  else
    if (IPRGLB >= DEBUG) write(u6,*) 'depsa (sym)'
  end if
  if (IPRGLB >= VERBOSE) call sqprt(depsa,nasht)

  !! Configuration Lagrangian for MS-CASPT2
  !! This is the partial derivative of the transition reduced
  !! density matrices
  if (IFMSCOUP) then
    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    call DerHEff(nConf,nRoots,nState,CLag,VECROT)
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    if (IPRGLB >= VERBOSE) then
      CPUT = CPTF10-CPTF0
      WALLT = TIOTF10-TIOTF0
      write(u6,'(a,2f10.2)') ' DerHEff : CPU/WALL TIME=',cput,wallt
      write(u6,*)
    end if
  end if

  !! I need to add the derivative of the effective Hamiltonian
  !! for MS-CASPT2, but this is done after orbital Lagrangian.
  !! I just have to have IVECC = T + lambda.

  !! If CASPT2 energy is not invariant to rotations in active
  !! orbitals, off-diagonal elements of the density obtained
  !! as DEPSA is incorrect, so remove them. The true density
  !! is computed after everything.
  if (.not. if_invar) then
    !! But, save the diagonal elements
    call mma_allocate(DEPSA_diag,nAshT,Label='DEPSA_diag')
    do i=1,nAshT
      DEPSA_diag(i) = DEPSA(i,i)
    end do
    !! Clear
    DEPSA(:,:) = Zero
  end if
  !write(u6,*) 'depsad'
  !call sqprt(depsa,nasht)

  !! Transform the quasi-variational amplitude (T+\lambda/2?)
  !! in SR (iVecX) to C (iVecC2)
  !! Note that the contribution is multiplied by two
  !! somewhere else (maybe in olagns?)
  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) then
    !! Have to weight the T-amplitude for MS-CASPT2
    if (IFMSCOUP) then
      if (.not. SC_prop) then
        !! add lambda
        call PLCVEC(VECROT(jState),Half,IVECX,IVECR)
        call PTRTOC(1,IVECR,IVECC2)
        !! T-amplitude
        do iStLag=1,nState
          if (iStLag == jState) cycle
          Scal = VECROT(iStLag)
          if (abs(Scal) <= 1.0e-12_wp) cycle
          call MS_Res(2,jStLag,iStLag,Scal*Half)
        end do
      end if
      if (do_csf) then
        !! Prepare for something <\Phi_K^{(1)}|Ers|L>
        ibk = IVECC2
        IVECC2 = 7
        call RHS_ZERO(IVECC2)
        do iStLag=1,nState
          if (iStLag == jState) cycle
          Scal = UEFF(iStLag,iRoot1)*UEFF(jStLag,iRoot2)-UEFF(jStLag,iRoot1)*UEFF(iStLag,iRoot2)
          Scal = Scal*Half
          if (abs(Scal) <= 1.0e-12_wp) cycle
          call MS_Res(2,jStLag,iStLag,Scal)
        end do
        IVECC2 = ibk
      end if
    else
      !! Add lambda to the T-amplitude
      call PLCVEC(Half,One,IVECR,IVECX)
      call PTRTOC(1,IVECX,IVECC2)
    end if
  end if

  !ipTrfL = 1+nAshT*nBasT+nAshT
  !Call DGemm_('n','N',nAshT,nAshT,nAshT,One,Trf(ipTrfL),nBasT,DEPSA,nAshT,Zero,dpt2c_ao,nAshT)
  !Call DGemm_('N','t',nAshT,nAshT,nAshT,One,dpt2c_ao,nAshT,Trf(ipTrfL),nBasT,Zero,DEPSA,nAshT)

  !! Just add DEPSA to DPT2
  call AddDEPSA(NBSQT,nAshT,DPT2,DEPSA)
  !! Just transform the density in MO to AO
  call DPT2_Trf(NBSQT,nAshT,DPT,DPT2_AO,CMOPT2,DEPSA,DSUM)
  !call mma_deallocate(DEPSA)
  !! Save the AO density
  !! ... write

  ! ----- Construct orbital Lagrangian -----

  if ((nFroT /= 0) .or. (.not. if_invaria)) then
    !! If frozen orbitals exist, we need to obtain
    !! electron-repulsion integrals with frozen orbitals to
    !! construct the orbital Lagrangian.
    if (.not. IfChol) call TRAFRO(1)

    !! Get density matrix (DIA) and inactive density
    !! matrix (DI) to compute FIFA and FIMO.
    call OLagFroD(NBSQT,nAshT,DIA,DI,RDMSA,Trf)
  end if

  !! Construct orbital Lagrangian that comes from the derivative
  !! of ERIs. Also, do the Fock transformation of the DPT2 and
  !! DPT2C densities.
  if (IfChol) then
    call mma_allocate(A_PT2,MaxVec_PT2**2,Label='A_PT2')
    A_PT2(:) = Zero
  else
    call mma_allocate(A_PT2,1,Label='A_PT2')
  end if
  do iSym=1,nSym
    nOcc = nIsh(iSym)+nAsh(iSym)
    lT2AO = 1
    if ((.not. IfChol) .or. (iALGO /= 1)) then
      lT2AO = nOcc*nOcc*nBasT*nBasT
      call mma_allocate(T2AO,lT2AO,Label='T2AO')
      T2AO(:) = Zero
    else
      call mma_allocate(T2AO,lT2AO,Label='T2AO')
    end if

    !! Orbital Lagrangian that comes from the derivative of ERIs.
    !! OLagNS computes only the particle orbitals.
    !write(u6,*) 'ialgo = ',ialgo
    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    if (IfChol .and. (iALGO == 1)) then
      call OLagNS_RI(iSym,NBSQT,MaxVec_PT2,DPT2C,DPT2Canti,A_PT2)
    else
      call OLagNS2(iSym,NBSQT,lT2AO,DPT2C,T2AO)
    end if
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    if (IPRGLB >= VERBOSE) then
      CPUT = CPTF10-CPTF0
      WALLT = TIOTF10-TIOTF0
      write(u6,'(a,2f10.2)') ' OLagNS  : CPU/WALL TIME=',cput,wallt
    end if
    !write(u6,*) 'DPT2C'
    !call sqprt(dpt2c,nbast)

    !! MO -> AO transformations for DPT2 and DPT2C
    if (((.not. IfChol) .or. (iALGO /= 1)) .or. ((nFroT == 0) .and. if_invaria)) then
      call OLagTrf(1,iSym,NBSQT,CMOPT2,DPT2,DPT2_AO,WRK1)
      call OLagTrf(1,iSym,NBSQT,CMOPT2,DPT2C,DPT2C_AO,WRK1)
      !write(u6,*) 'dpt2'
      !call sqprt(dpt2,nbast)
      !write(u6,*) 'dpt2ao'
      !call sqprt(dpt2_ao,nbast)
    end if

    !! Do some transformations relevant to avoiding (VV|VO)
    !! integrals. Orbital Lagrangian for the hole orbitals are
    !! computed. At the same time, F = G(D) transformations are
    !! also performed for D = DPT2 and DPT2C
    !! The way implemented (what?) is just a shit. I cannot find
    !! FIFA and FIMO for frozen orbitals, so I have to construct
    !! them. Here is the transformation of G(D^inact) and G(D).
    !! FIFA_all and FIMO_all computed in this subroutine
    !! is not yet correct. They are just two-electron after this
    !! subroutine.
    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    call OLagVVVO(iSym,NBSQT,lT2AO,MaxVec_PT2,DPT2_AO,DPT2C_AO,FPT2_AO,FPT2C_AO,T2AO,DIA,DI,FIFA_all,FIMO_all,A_PT2)
    !write(u6,*) 'olag after vvvo'
    !call sqprt(olag,nbast)
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    if (IPRGLB >= VERBOSE) then
      CPUT = CPTF10-CPTF0
      WALLT = TIOTF10-TIOTF0
      write(u6,'(a,2f10.2)') ' OLagVVVO: CPU/WALL TIME=',cput,wallt
    end if
    !write(u6,*) 'OLag'
    !do i=1,144
    !  write(u6,'(i3,f20.10)') i,olag(i)
    !end do
    !write(u6,*) 'fpt2ao'
    !call sqprt(fpt2_ao,12)
    !call abend()

    !! AO -> MO transformations for FPT2AO and FPT2CAO
    if (((.not. IfChol) .or. (iALGO /= 1)) .or. ((nFroT == 0) .and. if_invaria)) then
      call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2,FPT2_AO,WRK1)
      call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2C,FPT2C_AO,WRK1)
    end if

    call mma_deallocate(T2AO)
  end do
  call mma_deallocate(A_PT2)
  !! Add DPTC to DSUM for the correct unrelaxed density
  !! Also, symmetrize DSUM
  call AddDPTC(NBSQT,NDPT,DPT2C,DSUM)
  if (HZERO == 'DYALL') call BDerNEV_final2()

  !write(u6,*) 'fptao after olagns'
  !call sqprt(fpt2_ao,nbast)
  !write(u6,*) 'fptcao after olagns'
  !call sqprt(fpt2c_ao,nbast)
  !write(u6,*) 'olag after olagns'
  !call sqprt(olag,nbast)

  !! If frozen orbitals exist, frozen-inactive part of the
  !! unrelaxed PT2 density matrix is computed using the orbital
  !! Lagrangian. Additionally, Fock transformation is also
  !! required.
  if ((nFroT /= 0) .or. (.not. if_invaria)) then
    !! Compute DPT2 density for frozen-inactive
    !write(u6,*) 'dpt2 before frozen'
    !call sqprt(dpt2,nbast)
    if (.not. ifchol) then
      !! Construct FIFA and FIMO
      call OLagFro3(NBSQT,FIFA_all,FIMO_all,WRK1,WRK2)
      !! if possible, canonicalize frozen orbitals, and update
      !! FIMO and Trf
    end if
    !! Save DPT in order to subtract later
    WRK1(1:nDPTAO) = DPT2(1:nDPTAO)
    !! Add explicit FIMO and FIFA contributions. Implicit
    !! contributions are all symmetric in frozen + inactive
    !! orbitals, so they do not contribute to frozen density
    call DGEMM_('N','T',nBasT,nBasT,nBasT,One,FIMO_all,nBasT,DPT2C,nBasT,One,OLAG,nBasT)
    call DGEMM_('T','N',nBasT,nBasT,nBasT,One,FIMO_all,nBasT,DPT2C,nBasT,One,OLAG,nBasT)
    call DGEMM_('N','T',nBasT,nBasT,nBasT,One,FIFA_all,nBasT,DPT2,nBasT,One,OLAG,nBasT)
    call DGEMM_('T','N',nBasT,nBasT,nBasT,One,FIFA_all,nBasT,DPT2,nBasT,One,OLAG,nBasT)

    !! non-invariant in inactive and secondary
    if (.not. if_invaria) then
      !! Construct the density from orbital Lagrangian
      call caspt2_grad_invaria2(NBSQT,nOLag,DPT2,OLag)
      !! FIFA contributions from the non-invariant density
      DPT2(1:nDPTAO) = DPT2(1:nDPTAO)-WRK1(1:nDPTAO)
      !! Add the non-invariant contribution to unrelaxed density
      call AddDPTC(NBSQT,NDPT,DPT2,DSUM)
      call DGEMM_('N','T',nBasT,nBasT,nBasT,One,FIFA_all,nBasT,DPT2,nBasT,One,OLAG,nBasT)
      call DGEMM_('T','N',nBasT,nBasT,nBasT,One,FIFA_all,nBasT,DPT2,nBasT,One,OLAG,nBasT)
      !! Restore the second-order correlated density
      DPT2(1:nDPTAO) = DPT2(1:nDPTAO)+WRK1(1:nDPTAO)
      WRK1(1:nDPTAO) = DPT2(1:nDPTAO)
    end if
    !! Now, compute pseudo-density using orbital Lagrangian
    !! DSUM does not contain frozen orbitals,
    !! so the properties using this density may be inaccurate
    if (nFroT /= 0) call OLagFro1(NBSQT,nOLag,DPT2,OLag)

    !! Subtract the orbital Lagrangian added above.
    !! It is computed again in EigDer
    call DGEMM_('N','T',nBasT,nBasT,nBasT,-One,FIMO_all,nBasT,DPT2C,nBasT,One,OLAG,nBasT)
    call DGEMM_('T','N',nBasT,nBasT,nBasT,-One,FIMO_all,nBasT,DPT2C,nBasT,One,OLAG,nBasT)
    call DGEMM_('N','T',nBasT,nBasT,nBasT,-One,FIFA_all,nBasT,WRK1,nBasT,One,OLAG,nBasT)
    call DGEMM_('T','N',nBasT,nBasT,nBasT,-One,FIFA_all,nBasT,WRK1,nBasT,One,OLAG,nBasT)
    !write(u6,*) 'dpt after frozen'
    !call sqprt(dpt2,nbast)

    !! Fock transformation for frozen-inactive density
    if (IfChol) then
      iSym = 1
      !! MO -> AO transformations for DPT2 and DPT2C
      call OLagTrf(1,iSym,NBSQT,CMOPT2,DPT2,DPT2_AO,WRK1)
      call OLagTrf(1,iSym,NBSQT,CMOPT2,DPT2C,DPT2C_AO,WRK1)
      !! For DF-CASPT2, Fock transformation of DPT2, DPT2C, DIA,
      !! DA is done here, but not OLagVVVO
      !! It seems that it is not possible to do this
      !! transformation in OLagVVVO, because the frozen-part of
      !! the DPT2 is obtained after OLagVVVO.
      FPT2_AO(:) = Zero
      FPT2C_AO(:) = Zero
      call OLagFro4(NBSQT,1,1,1,1,1,DPT2_AO,DPT2C_AO,FPT2_AO,FPT2C_AO,WRK1)
      !! AO -> MO transformations for FPT2AO and FPT2CAO
      call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2,FPT2_AO,WRK1)
      call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2C,FPT2C_AO,WRK1)
    else
      !write(u6,*) 'dpt'
      !call sqprt(dpt2,nbast)
      call OLagFro2(NBSQT,DPT2,FPT2,WRK1,WRK2)
      !write(u6,*) 'fpt'
      !call sqprt(dpt2,nbast)
    end if
    !write(u6,*) 'fifa_all'
    !call sqprt(fifa_all,12)
    !write(u6,*) 'fimo_all'
    !call sqprt(fimo_all,12)
    !! Construct FIFA and FIMO
    !call OLagFro3(NBSQT,FIFA_all,FIMO_all,WRK1,WRK2)
  else ! there are no frozen orbitals
    iSQ = 0
    iTR = 0
    do iSym=1,nSym
      nOrbI = nOrb(iSym)
      call SQUARE(FIFA(1+iTR),FIFA_all(1+iSQ),1,nOrbI,nOrbI)
      call SQUARE(FIMO(1+iTR),FIMO_all(1+iSQ),1,nOrbI,nOrbI)
      iSQ = iSQ+nOrbI*nOrbI
      iTR = iTR+nTri_Elem(nOrbI)
    end do
  end if
  call mma_deallocate(DIA)
  call mma_deallocate(DI)
  !write(u6,*) 'fifa_all in dens'
  !call sqprt(fifa_all,nbast)
  !write(u6,*) 'fimo_all in dens'
  !call sqprt(fimo_all,nbast)
  !write(u6,*) 'FIFA in natural'
  !call DGemm_('N','N',nBasT,nBasT,nBasT,One,Trf,nBasT,fifa_all,nBasT,Zero,WRK1,nBasT)
  !call DGemm_('N','T',nBasT,nBasT,nBasT,One,WRK1,nBasT,Trf,nBasT,Zero,WRK2,nBasT)
  !call sqprt(wrk2,12)

  !! Do some post-process for the contributions that comes from
  !! the above two densities.
  !write(u6,*) 'olag before eigder'
  !call sqprt(olag,nbast)
  !write(u6,*) 'fpt2'
  !call sqprt(fpt2,nbast)
  call EigDer(NBSQT,nAshT,DPT2,DPT2C,FPT2_AO,FPT2C_AO,RDMEIG,CMOPT2,Trf,FPT2,FPT2C,FIFA_all,FIMO_all,RDMSA)
  !call test2_dens(olag,depsa)
  !write(u6,*) 'olag after eigder'
  !call sqprt(olag,nbast)
  !write(u6,*) 'Wlag after eigder'
  !call sqprt(wlag,nbast)
  !write(u6,*) 'rdmeig'
  !call sqprt(rdmeig,nasht)
  !call abend()

  !! Calculate the configuration Lagrangian again.
  !! The contribution comes from the derivative of eigenvalues.
  !! It seems that TRACI_RPT2 uses CI coefficients of RASSCF,
  !! so canonical -> natural transformation is required.
  !ipTrfL = 1+nAshT*nBasT+nAshT
  !call DGemm_('N','N',nAshT,nAshT,nAshT,One,Trf(ipTrfL),nBasT,RDMEIG,nAshT,Zero,WRK1,nAshT)
  !call DGemm_('N','T',nAshT,nAshT,nAshT,One,WRK1,nAshT,Trf(ipTrfL),nBasT,Zero,RDMEIG,nAshT)
  if (.not. if_invar) then !test
    call mma_allocate(CLagT,nConf,nState,Label='CLagT')
    call mma_allocate(EigT,nAshT,nAshT,Label='EigT')
    CLagT(:,:) = CLag(:,:)
    EigT(:,:) = RDMEIG(:,:)
    if (IFDW .and. (zeta >= Zero)) then
      call mma_allocate(OMGT,nState,nState,Label='OMGT')
      OMGT(:,:) = OMGDER(:,:)
    end if
  end if
  !! Use canonical CSFs rather than natural CSFs in CLagEig
  call mma_allocate(ISAV,size(IDCIEX),Label='ISAV')
  ISAV(:) = IDCIEX(:)
  IDCIEX(:) = IDTCEX(:)
  !! Now, compute the configuration Lagrangian
  call CLagEig(if_SSDM,.false.,nConf,nRoots,nState,NLEV,CLag,RDMEIG)
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) call GADGOP(CLag,nCLag,'+')
# endif

  !! Now, here is the best place to compute the true off-diagonal
  !! active density for non-invariant CASPT2
  if (.not. if_invar) then
    SLag(:,:) = Zero
    !! Add the density that comes from CI Lagrangian
    call DEPSAOffC(nConf,nState,nAshT,nBasT,CLag,DEPSA,FIFA_all,FIMO_all,WRK1,WRK2,U0)
    !! Add the density that comes from orbital Lagrangian
    call DEPSAOffO(nOLag,nAshT,NBSQT,OLag,DEPSA,FIFA_all)
    !! Restore the diagonal elements
    do i=1,nAshT
      DEPSA(i,i) = DEPSA_diag(i)
    end do
    call mma_deallocate(DEPSA_diag)
    if (IPRGLB >= VERBOSE) then
      write(u6,*) 'DEPSA computed again'
      call sqprt(depsa,nasht)
    end if
    if (NRAS1T+NRAS3T /= 0) then
      !! Remove the off-diagonal blocks for RASPT2
      DEPSA(1:nRAS1T,nRAS1T+1:) = Zero
      DEPSA(nRAS1T+1:,1:nRAS1T) = Zero
      DEPSA(nRAS1T+1:nRAS1T+nRAS2T,nRAS1T+nRAS2T+1:) = Zero
      DEPSA(nRAS1T+nRAS2T+1:,nRAS1T+1:nRAS1T+nRAS2T) = Zero
    end if
    !depsa(:,:) = Zero

    !! We have to do many things again...
    !! Just add DEPSA to DPT2
    call AddDEPSA(NBSQT,nAshT,DPT2,DEPSA)
    !! Just transform the density in MO to AO
    call DPT2_Trf(NBSQT,nAshT,DPT,DPT2_AO,CMOPT2,DEPSA,DSUM)
    !! For IPEA shift with state-dependent density
    if (if_SSDM .and. ((jState == iRlxRoot) .or. IFMSCOUP)) then
      iSym = 1
      call OLagTrf(1,iSym,NBSQT,CMOPT2,DPT2,DPT2_AO,WRK1)
    end if
    !! Some transformations similar to EigDer
    call EigDer2(NBSQT,nAshT,RDMEIG,Trf,FIFA_all,RDMSA,DEPSA,WRK1,WRK2)

    CLag(:,:) = CLagT(:,:) !test
    ! test
    RDMEIG(:,:) = RDMEIG(:,:)+EigT(:,:)
    SLag(:,:) = Zero
    call mma_deallocate(CLagT)
    call mma_deallocate(EigT)
    if (IFDW .and. (zeta >= Zero)) then
      OMGDER(:,:) = OMGT(:,:)
      call mma_deallocate(OMGT)
    end if
    !! RDMEIG contributions
    !! Use canonical CSFs rather than natural CSFs
    !! Now, compute the configuration Lagrangian
    call CLagEig(if_SSDM,.false.,nConf,nRoots,nState,NLEV,CLag,RDMEIG)
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) call GADGOP(CLag,nCLag,'+')
#   endif
    !! Now, compute the state Lagrangian and do some projections
    !call CLagFinal(nConf,nState,CLag,SLag)
  end if

  !! Restore integrals without frozen orbitals, although not sure
  !! this operation is required.
  if (((nFroT /= 0) .or. (.not. if_invaria)) .and. (.not. IfChol)) call TRAFRO(2)

  IDCIEX(:) = ISAV(:)
  call mma_deallocate(ISAV)
  !! Canonical -> natural transformation
  if (ORBIN == 'TRANSFOR') then
    do iState=1,nState
      call CLagX_TrfCI(nConf,CLag(1,iState))
    end do
  end if
  ! accumulate configuration Lagrangian only for MS,XMS,XDW,RMS,
  ! but not for SS-CASPT2
  if ((jState == iRlxRoot) .or. IFMSCOUP) CLagFull(1:nConf,1:nState) = CLagFull(1:nConf,1:nState)+CLag(1:nConf,1:nState)
  !call CLagFinal(nConf,nState,CLag,SLag)

  !! Transformations of DPT2 in quasi-canonical to natural orbital
  !! basis and store the transformed density so that the MCLR
  !! module can use them.
  ! accumulate only if MS,XMS,XDW or RMS calculation
  !call RecPrt('DPT2 before', '', DPT2_tot, nBast, nBast)
  if ((jState == iRlxRoot) .or. IFMSCOUP) then
    call DPT2_TrfStore(One,NBSQT,DPT2,DPT2_tot,Trf,WRK1)
    call DPT2_TrfStore(Two,NBSQT,DPT2C,DPT2C_tot,Trf,WRK1)
    if (do_csf) call DPT2_TrfStore(One,NBSQT,DPT2Canti,DPT2Canti_tot,Trf,WRK1)
  end if
  !call RecPrt('DPT2 after','',DPT2_tot,nBast,nBast)
  !! Save MO densities for post MCLR
  !call DGemm_('N','N',nBasT,nBasT,nBasT,One,Trf,nBasT,DPT,nBasT,Zero,WRK1,nBasT)
  !call DGemm_('N','T',nBasT,nBasT,nBasT,One,WRK1,nBasT,Trf,nBasT,Zero,WRK2,nBasT)
  !iSQ = 0
  !do iSym=1,nSym
  !  nOrbI = nBas(iSym)-nDel(iSym)
  !  nSQ = nOrbI*nOrbI
  !  DPT2_tot(iSQ+1:iSQ+nSQ) = DPT2_tot(iSQ+1:iSQ+nSQ)+WRK2(iSQ+1:iSQ+nSQ)
  !  iSQ = iSQ+nSQ
  !end do

  !! Do the same for DPT2C Save MO densities for post MCLR
  !call DGemm_('N','N',nBasT,nBasT,nBasT,One,Trf,nBasT,DPT2C,nBasT,Zero,WRK1,nBasT)
  !call DGemm_('N','T',nBasT,nBasT,nBasT,One,WRK1,nBasT,Trf,nBasT,Zero,WRK2,nBasT)
  !iSQ = 0
  !do iSym=1,nSym
  !  nOrbI = nBas(iSym)-nDel(iSym)
  !  nSQ = nOrbI*nOrbI
  !  DPT2C_tot(iSQ+1:iSQ+nSQ) = DPT2C_tot(iSQ+1:iSQ+nSQ)+Two*WRK2(iSQ+1:iSQ+nSQ)
  !  iSQ = iSQ+nSQ
  !end do
  !call abend()
  !call sqprt(RDMEIG,nAshT)

  !! square -> triangle so that the MCLR module can use the AO
  !! densities. Do this for DPT2AO and DPT2CAO (defined in
  !! caspt2_grad).
  ! accumulate only if MS,XMS,XDW or RMS calculation
  if ((jState == iRlxRoot) .or. IFMSCOUP) then
    iBasTr = 1
    iBasSq = 1
    do iSym=1,nSym
      nBasI = nBas(iSym)
      liBasTr = iBasTr
      liBasSq = iBasSq
      do iBasI=1,nBasI
        do jBasI=1,iBasI
          liBasSq = iBasSq+iBasI-1+nBasI*(jBasI-1)
          if (iBasI == jBasI) then
            DPT2_AO_tot(liBasTr) = DPT2_AO(liBasSq)
            DPT2C_AO_tot(liBasTr) = DPT2C_AO(liBasSq)
          else
            DPT2_AO_tot(liBasTr) = DPT2_AO(liBasSq)*Two
            DPT2C_AO_tot(liBasTr) = DPT2C_AO(liBasSq)*Two
          end if
          liBasTr = liBasTr+1
        end do
      end do
      iBasTr = iBasTr+nTri_Elem(nBasI)
      iBasSq = iBasSq+nBasI**2
    end do
  end if

  !! If the density matrix used in the Fock operator is different
  !! from the averaged density in the SCF calculation, we need an
  !! additional term for electron-repulsion integral.
  !! Here prepares such densities.
  !! The first one is just DPT2AO, while the second one is the
  !! difference between the SS and SA density matrix. because the
  !! SA density-contribution will be added and should be
  !! subtracted
  ! This should be done only for iRlxRoot
  if (if_SSDM .and. ((jState == iRlxRoot) .or. IFMSCOUP)) then
    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    !if (.not. if_invar) then
    !  write(u6,*) 'SS density matrix with IPEA not implemented'
    !  call abend()
    !end if

    !! Construct the SCF density
    !! We first need to construct the density averaged over all
    !! roots involved in SCF.
    WRK1(1:nDRef) = Zero
    call mma_allocate(CI1,nConf,Label='CI1')
    do iState=1,nState
      call LoadCI_XMS('N',1,nConf,nState,CI1,iState,U0)
      call POLY1(CI1,nConf)
      call GETDREF(WRK2,nDRef)
      wgt = Weight(iState)
      WRK1(1:nDRef) = WRK1(1:nDRef)+Wgt*WRK2(1:nDRef)
    end do
    call mma_deallocate(CI1)
    !! WRK2 is the SCF density (for nstate=nroots)
    call SQUARE(WRK1,WRK2,1,nAshT,nAshT)
    RDMSA(:,:) = RDMSA(:,:)-reshape(WRK2(1:nAshT**2),[nAshT,nAshT])

    !! Construct the SS minus SA density matrix in WRK1
    call OLagFroD(NBSQT,nAshT,WRK1,WRK2,RDMSA,Trf)
    !! Subtract the inactive part
    WRK1(1:nBasT**2) = WRK1(1:nBasT**2)-WRK2(1:nBasT**2)
    !! Here we should use DPT2_AO??
    !! Save
    if (IfChol) then
      call CnstAB_SSDM(NBSQT,DPT2_AO,WRK1)
    else
      !! Well, it is not working any more. I need to use
      !! Position='APPEND', but it is not possible if I need to
      !! use molcas_open or molcas_open_ext2
      write(u6,*) 'It is not possible to perform this calculation'
      write(u6,*) '(non-state averaged density without'
      write(u6,*) 'density-fitting or Cholesky decomposition)'
      write(u6,*) 'Please use DF or CD'
      call abend()
    end if
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    if (IPRGLB >= VERBOSE) then
      CPUT = CPTF10-CPTF0
      WALLT = TIOTF10-TIOTF0
      write(u6,'(a,2f10.2)') ' SSDM    : CPU/WALL TIME=',cput,wallt
    end if
  end if
  !write(u6,*) 'pt2ao'
  !call sqprt(DPT2_AO,12)
  call mma_deallocate(DEPSA)

  call mma_deallocate(DPT2)
  call mma_deallocate(DPT2C)
  call mma_deallocate(DPT2_AO)
  call mma_deallocate(DPT2C_AO)
  call mma_deallocate(FPT2)
  call mma_deallocate(FPT2C)
  call mma_deallocate(FPT2_AO)
  call mma_deallocate(FPT2C_AO)
  if (do_csf) call mma_deallocate(DPT2Canti_)
  nullify(DPT2Canti)

  !! Finalize OLag (anti-symmetrize) and construct WLag
  call OLagFinal(nOLag,NBSQT,OLag,Trf)

  call mma_deallocate(TRF)
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)
  call mma_deallocate(RDMSA)
  call mma_deallocate(RDMEIG)
  call mma_deallocate(VECROT)
  DENORM = One
  !! end of with gradient
else
  !! without gradient
  ! Compute total density matrix as symmetry-blocked array of
  ! triangular matrices in DMAT. Size of a triangular submatrix is
  ! nTri_Elem(NORB(ISYM)).
  NDPT = sum(NORB(1:nSym)**2)
  DMAT(1:NDMAT) = Zero
  ! First, put in the reference density matrix.
  IDMOFF = 0
  do ISYM=1,NSYM
    NI = NISH(ISYM)
    NA = NASH(ISYM)
    NO = NORB(ISYM)
    do II=1,NI
      IDM = IDMOFF+nTri_Elem(II)
      DMAT(IDM) = Two
    end do
    do IT=1,NA
      ITABS = NAES(ISYM)+IT
      ITTOT = NI+IT
      do IU=1,IT
        IUABS = NAES(ISYM)+IU
        IUTOT = NI+IU
        IDRF = iTri(ITABS,IUABS)
        IDM = IDMOFF+iTri(ITTOT,IUTOT)
        DMAT(IDM) = DREF(IDRF)
      end do
    end do
    IDMOFF = IDMOFF+nTri_Elem(NO)
  end do
  !write(u6,*) ' DENS. Initial DMAT:'
  !write(u6,'(1x,8f16.8)') (dmat(i),i=1,ndmat)
  ! Add the 1st and 2nd order density matrices:
  call mma_allocate(DPT,NDPT,Label='DPT')
  call mma_allocate(DSUM,NDPT,Label='DSUM')
  DSUM(1:NDPT) = Zero

  ! The 1st order contribution to the density matrix
  DPT(1:NDPT) = Zero
  call TRDNS1(IVEC,DPT,NDPT)
  DSUM(1:NDPT) = DSUM(1:NDPT)+DPT(1:NDPT)
  !write(u6,*) ' DPT after TRDNS1.'
  !write(u6,'(1x,8f16.8)') (dpt(i),i=1,ndpt)
  DPT(1:NDPT) = Zero
  call TRDNS2D(IVEC,IVEC,DPT,NDPT,One)
  if (IFDENS) then
    ! The exact density matrix evaluation:
    call TRDTMP(DPT,NDPT)
  else
    ! The approximate density matrix evaluation:
    call TRDNS2A(IVEC,IVEC,DPT,NDPT)
  end if
  DSUM(1:NDPT) = DSUM(1:NDPT)+DPT(1:NDPT)
  !write(u6,*) ' DPT after TRDNS2D.'
  !write(u6,'(1x,8f16.8)') (dpt(i),i=1,ndpt)
  DPT(1:NDPT) = Zero
  call TRDNS2O(IVEC,IVEC,DPT,size(DPT),NDPT,One)
  DSUM(1:NDPT) = DSUM(1:NDPT)+DPT(1:NDPT)
  !write(u6,*) ' DPT after TRDNS2O.'
  !write(u6,'(1x,8f16.8)') (dpt(i),i=1,ndpt)
end if

call mma_deallocate(DPT)
IDMOFF = 0
IDSOFF = 0
do ISYM=1,NSYM
  NO = NORB(ISYM)
  do IP=1,NO
    do IQ=1,IP
      IDM = IDMOFF+iTri(IP,IQ)
      IDSUM = IDSOFF+IP+NO*(IQ-1)
      DMAT(IDM) = DMAT(IDM)+DSUM(IDSUM)
    end do
  end do
  IDMOFF = IDMOFF+nTri_Elem(NO)
  IDSOFF = IDSOFF+NO**2
end do
call mma_deallocate(DSUM)
! Scale with 1/DENORM to normalize
X = One/DENORM
if (do_grad) X = One
DMAT(1:NDMAT) = X*DMAT(1:NDMAT)

!SVC: For true parallel calculations, replicate the DMAT array
! so that the slaves have the same density matrix as the master.
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (.not. KING()) DMAT(1:NDMAT) = Zero
  call GADGOP(DMAT,NDMAT,'+')
end if
#endif

end subroutine DENS

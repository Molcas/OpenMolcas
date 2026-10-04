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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

! Strongly contracted NEVPT2 (SC-NEVPT2)
! See Angeli, C.; Cimiraglia, R.; Malrieu, J.-P. J. Chem. Phys. 2002, 117, 9138. for equations

module SC_NEVPT2

use Index_Functions, only: iTri, nTri_Elem
use caspt2_global, only: imag_shift, iPrGlb, jStLag, LUSBT, real_shift, sigma_p_epsilon
use caspt2_module, only: CASES, EREF, IfChol, JSTATE, MXCASE, NASUP, NISUP, NSTATE, NSYM
use EQSOLV, only: IDBMAT, IDSMAT, IVECC, IVECC2, IVECW
use fake_GA, only: GA_Arrays
# ifdef _MOLCAS_MPP_
use GA_Wrapper, only: DBL_MB, GA_NodeId
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
private

!! Do_FIC       :: Compute partially contracted or fully internally contracted NEVPT2 (PC- or FIC-NEVPT2) energy
!! Do_SC        :: Compute strongly contracted NEVPT2 (SC-NEVPT2) energy
!! SC_prop      :: Do property (energy and gradients) calculations using SC-NEVPT2 wavefunction
!! SC_amplitude :: Compute the T-amplitude for SC-NEVPT2 wavefunction
!! SC_thres     :: Threshold for avoiding vanishing denominators: the main purpose is to avoid numerical blow up, not for efficiency
!! SC_merge     :: Choose perturber functions (.true. -> SCPF = MERGE, .false. -> SCPF = SPLIT)

!----- SC-NEVPT2
! ECORR_SC                         :: Correlation energy
! OVLAPS_SC                        :: Wavefunction (correction) overlap
! IDBMAT_NEVPT2(1:nSym,1:MXCASE,1) :: The original B matrix
! IDBMAT_NEVPT2(1:nSym,1:MXCASE,2) :: denominator of the amplitude

!----- QD-SC-NEVPT2
! DVALUE_SC :: temporary storage of DVALUE
! ENERGY_SC :: energy
! HEFF_SC   :: Effective Hamiltonian
! UEFF_SC   :: Eigenvector

integer(kind=iwp) :: IDBMAT_NEVPT2(8,13,2)
real(kind=wp) :: DVALUE_SC, ECORR_SC(0:8,0:MXCASE), OVLAPS_SC(0:8,0:MXCASE), SC_thres = 1.0e-9_wp
logical(kind=iwp) :: Do_FIC, Do_SC, SC_amplitude = .false., SC_merge = .false., SC_prop = .false.
real(kind=wp), allocatable :: ENERGY_SC(:), HEFF_SC(:,:), UEFF_SC(:,:)

public :: Do_FIC, Do_SC, DVALUE_SC, ECORR_SC, ENERGY_SC, HEFF_SC, IDBMAT_NEVPT2, OVLAPS_SC, SC_amplitude, SC_merge, &
          SC_NEVPT2_amplitude, SC_NEVPT2_CLagD, SC_NEVPT2_final, SC_NEVPT2_initial, SC_NEVPT2_Print, SC_NEVPT2_res, SC_prop, &
          SC_thres, UEFF_SC

contains

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_initial()

  call mma_allocate(ENERGY_SC,nState,Label='ENERGY_SC')
  call mma_allocate(HEFF_SC,nState,nState,Label='HEFF_SC')
  call mma_allocate(UEFF_SC,nState,nState,Label='UEFF_SC')

  ENERGY_SC(:) = Zero
  HEFF_SC(:,:) = Zero
  UEFF_SC(:,:) = Zero

end subroutine SC_NEVPT2_initial

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_final()

  call mma_deallocate(ENERGY_SC)
  call mma_deallocate(HEFF_SC)
  call mma_deallocate(UEFF_SC)

end subroutine SC_NEVPT2_final

!-----------------------------------------------------------------------

function SC_NEVPT2_Coupling(iCase)

  integer(kind=iwp), intent(in) :: iCase
  integer(kind=iwp) :: SC_NEVPT2_Coupling

  SC_NEVPT2_Coupling = 0
  if (.not. SC_merge) return

  select case (iCase)
    case (2,6,8,10)
      SC_NEVPT2_Coupling = +1
    case (1,4,5,12,13)
      SC_NEVPT2_Coupling = 0
    case (3,7,9,11)
      SC_NEVPT2_Coupling = -1
    case default
      write(u6,*) 'SC_NEVPT2_Coupling: illegal case number', iCase
      call abend()
  end select

end function SC_NEVPT2_Coupling

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_Map(iCase,iSym,NIS,Map)

  use Symmetry_Info, only: Mul
  use SUPERINDEX, only: KAGEB, KIGEJ, MAGTB, MIGTJ
  use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NISH, NSSH

  integer(kind=iwp), intent(in) :: iCase, iSym, NIS
  integer(kind=iwp), intent(out) :: Map(NIS)

  integer(kind=iwp) :: ia, iab, iaq, ibq, ii, iij, iiq, ijq, ioff, isp, iSym1, iSym2, jab, jij

  ! Construct the map array for converting inactive/virtual superindex of minus to plus

  select case (iCase)
    case (3)
      ! VJTIM -> VJTIP: i>j --> i>=j
      do isp=1,NIGTJ(iSym)
        iiq = MIGTJ(1,isp+NIGTJES(iSym))
        ijq = MIGTJ(2,isp+NIGTJES(iSym))
        Map(isp) = KIGEJ(iiq,ijq)-NIGEJES(iSym)
      end do
    case (7)
      ! VJAIM --> VJAIP: a*(i>j) --> a*(i>=j)
      isp = 0
      ioff = 0
      do iSym1=1,nSym
        iSym2 = Mul(iSym1,iSym)
        do iij=1,NIGTJ(iSym2)
          iiq = MIGTJ(1,iij+NIGTJES(iSym2))
          ijq = MIGTJ(2,iij+NIGTJES(iSym2))
          jij = KIGEJ(iiq,ijq)-NIGEJES(iSym2)
          do ia=1,NSSH(iSym1)
            isp = isp+1
            Map(isp) = ioff+ia+NSSH(iSym1)*(jij-1)
          end do
        end do
        ioff = ioff+NSSH(iSym1)*NIGEJ(iSym2)
      end do
    case (9)
      ! BVATM --> BVATP: a>b --> a>=b
      do isp=1,NAGTB(iSym)
        iaq = MAGTB(1,isp+NAGTBES(iSym))
        ibq = MAGTB(2,isp+NAGTBES(iSym))
        Map(isp) = KAGEB(iaq,ibq)-NAGEBES(iSym)
      end do
    case (11)
      ! BJATM -> BJATP: i*(a>b) --> i*(a>=b)
      isp = 0
      ioff = 0
      do iSym1=1,nSym
        iSym2 = Mul(iSym1,iSym)
        do iab=1,NAGTB(iSym2)
          iaq = MAGTB(1,iab+NAGTBES(iSym2))
          ibq = MAGTB(2,iab+NAGTBES(iSym2))
          jab = KAGEB(iaq,ibq)-NAGEBES(iSym2)
          do ii=1,NISH(iSym1)
            isp = isp+1
            Map(isp) = ioff+ii+NISH(iSym1)*(jab-1)
          end do
        end do
        ioff = ioff+NISH(iSym1)*NAGEB(iSym2)
      end do
    case default
      write(u6,*) 'SC_NEVPT2_Map: illegal case number', iCase
      call abend()
  end select

end subroutine SC_NEVPT2_Map

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_Print()

  use PrintLevel, only: TERSE, USUAL

  integer(kind=iwp) :: IC, iCase, IS, iSym, LAXITY
  real(kind=wp) :: DENORM_SC, E2CORR_SC, E2NONV_SC, E2TOT_SC, EAIVX_SC, EATVX_SC, EBJAI_SC, EBJAT_SC, EBVAT_SC, EVJAI_SC, &
                   EVJTI_SC, EVJTU_SC, REFWGT_SC, RNORM_SC
  character(len=8) :: LABEL
  integer(kind=iwp), external :: Cho_X_GetTol

  if (.not. Do_SC) return

  call SC_NEVPT2_Energy()

  do iCase=1,MXCASE
    ECORR_SC(0,iCase) = sum(ECORR_SC(1:nSym,iCase))
    OVLAPS_SC(0,iCase) = sum(OVLAPS_SC(1:nSym,iCase))
  end do

  do iSym=1,nSym
    ECORR_SC(iSym,0) = sum(ECORR_SC(iSym,0:MXCASE))
    OVLAPS_SC(iSym,0) = sum(OVLAPS_SC(iSym,0:MXCASE))
  end do

  ECORR_SC(0,0) = sum(ECORR_SC(1:nSym,1:MXCASE))
  OVLAPS_SC(0,0) = sum(OVLAPS_SC(1:nSym,1:MXCASE))

  EVJTU_SC = ECORR_SC(0,1)
  EVJTI_SC = ECORR_SC(0,2)+ECORR_SC(0,3)
  EATVX_SC = ECORR_SC(0,4)
  EAIVX_SC = ECORR_SC(0,5)
  EVJAI_SC = ECORR_SC(0,6)+ECORR_SC(0,7)
  EBVAT_SC = ECORR_SC(0,8)+ECORR_SC(0,9)
  EBJAT_SC = ECORR_SC(0,10)+ECORR_SC(0,11)
  EBJAI_SC = ECORR_SC(0,12)+ECORR_SC(0,13)

  E2NONV_SC = ECORR_SC(0,0)
  E2CORR_SC = E2NONV_SC !! Two*E2NONV_SC + OVLAPS_SC(0,0)
  E2TOT_SC = EREF+E2CORR_SC
  DENORM_SC = One+OVLAPS_SC(0,0)
  REFWGT_SC = One/DENORM_SC
  RNORM_SC = Zero

  if (IPRGLB > USUAL) then
    write(u6,*)
    write(u6,*) ' Correlation energy /case, /Symm, and sums (SC-NEVPT2):'
    do IC=1,13
      select case (SC_NEVPT2_Coupling(IC))
        case (+1)
          LABEL = CASES(IC)(1:len_trim(CASES(IC))-1)
        case (0)
          LABEL = CASES(IC)
        case (-1)
          cycle
      end select
      write(u6,'(1X,A8,9F12.8)') LABEL,(ECORR_SC(IS,IC),IS=1,NSYM),ECORR_SC(0,IC)
    end do
    write(u6,'(1X,A8,9F12.8)') 'Summed: ',(ECORR_SC(IS,0),IS=1,NSYM),ECORR_SC(0,0)
  end if

  if (IPRGLB >= TERSE) then
    write(u6,*)
    write(u6,*) ' FINAL SC-NEVPT2 RESULT:'
    write(u6,*)
    write(u6,'(6x,a,f18.10)') 'Reference energy:     ',EREF
    write(u6,'(6x,a,f18.10)') 'E2 (Non-variational): ',E2NONV_SC
    write(u6,'(6x,a,f18.10)') 'E2 (Variational):     ',E2CORR_SC
    write(u6,'(6x,a,f18.10)') 'Total energy:         ',E2TOT_SC
    write(u6,'(6x,a,f18.10)') 'Residual norm:        ',RNORM_SC
    write(u6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT_SC

    write(u6,*)
    write(u6,'(6x,a)') 'Contributions to the SC-NEVPT2 correlation energy'
    write(u6,'(6x,a,F18.10)') 'Active & Virtual only:    ',EATVX_SC+EBVAT_SC
    write(u6,'(6x,a,F18.10)') 'One Inactive Excited:     ',EVJTU_SC+EAIVX_SC+EBJAT_SC
    write(u6,'(6x,a,F18.10)') 'Two Inactive Excited:     ',EVJTI_SC+EVJAI_SC+EBJAI_SC
    write(u6,*)
  end if

  ENERGY_SC(JSTATE) = E2TOT_SC

  !! for verification
  LAXITY = 8
  if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
  call Add_Info('E_CASPT2',[E2TOT_SC],1,LAXITY)

end subroutine SC_NEVPT2_Print

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_Energy()

  integer(kind=iwp) :: ias, iCase, iCoupling, idisk, IHI1, ILO1, isp, iSym, jCase, JHI1, JLO1, lg_V, NAS, NAS2, NIS, NIS2
  real(kind=wp) :: Dummy(1), etmp, otmp, vale
  logical(kind=iwp) :: do_H

  integer(kind=iwp), allocatable :: Map(:)
  real(kind=wp), allocatable :: E_value(:), LBD(:), LID(:), ValHN(:,:), ValHN2(:,:)
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: LDV, MV, myRank
# endif

  ECORR_SC(:,0:11) = Zero
  OVLAPS_SC(:,0:11) = Zero

  if (.not. Do_SC) return

  ! Compute the SC-NEVPT2 energy
  ! Additionally, compute the strongly contracted zeroth-order active energy

  ! Hmm, I wish I could compute the T-amplitude in the SR basis...
  ! There is no linear dependence in SC-NEVPT2, so NIN = NAS.
  ! However, it is not easy to modify the variable later in the program.

  do_H = .false.
  if (.not. Do_FIC) do_H = .true.
  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero)) do_H = .true.

  do iSym=1,nSym
    do iCase=1,MXCASE
      !cycle
      !if ((icase /= 12) .and. (icase /= 13)) cycle ! H
      !if ((icase /= 10) .and. (icase /= 11)) cycle ! G
      !if ((icase /= 6) .and. (icase /= 7)) cycle ! E
      !if ((icase /= 8) .and. (icase /= 9)) cycle ! F
      !if ((icase /= 2) .and. (icase /= 3)) cycle ! B
      !if (icase /= 5) cycle ! D
      !if (icase /= 4) cycle ! C
      !if (icase /= 1) cycle ! A

      NAS = NASUP(iSym,iCase)
      NIS = NISUP(iSym,iCase)
      if (NAS*NIS == 0) cycle

      if ((iCase == 12) .or. (iCase == 13)) then
        ! Compute the H subspace only when real/imaginary shift has been used
        ! The H subspace energy is not vulnerable to inactive/virtual rotations
        if (do_H) then
          call RHS_ALLO(NAS,NIS,lg_V)
          call RHS_READ(NAS,NIS,lg_V,iCase,iSym,IVECW)

          call mma_allocate(LBD,NAS,Label='LBD')
          call mma_allocate(LID,NIS,Label='LID')
          idisk = IDBMAT(iSym,iCase)
          call dDaFile(LUSBT,2,LBD,NAS,idisk)
          call dDaFile(LUSBT,2,LID,NIS,idisk)

          ILO1 = 1
          IHI1 = NAS
          JLO1 = 1
          JHI1 = NIS
          etmp = Zero
          otmp = Zero
#         ifdef _MOLCAS_MPP_
          if (Is_Real_Par()) then
            myRank = GA_NodeID()
            call GA_DISTRIBUTION(lg_V,myRank,ILO1,IHI1,JLO1,JHI1)
            if ((ILO1 > 0) .and. (JLO1 > 0)) then
              call GA_Access(lg_V,ILO1,IHI1,JLO1,JHI1,MV,LDV)
              if (IHI1-ILO1+1 /= NAS) then
                write(u6,'(1x,A)') 'Something is wrong in SC_NEVPT2_Energy'
                call abend()
              end if
              call SC_NEVPT2_Energy_H_MPP(DBL_MB(MV),LDV,JLO1,JHI1,NAS,NIS,LBD,LID,etmp,otmp)
              call GA_RELEASE(lg_V,ILO1,IHI1,JLO1,JHI1)
            end if
          else
#         endif
            do isp=JLO1,JHI1
              do ias=ILO1,IHI1
                etmp = etmp-GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))/(lid(isp)+lbd(ias))
                otmp = otmp+GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))/(lid(isp)+lbd(ias))**2
              end do
            end do
#         ifdef _MOLCAS_MPP_
          end if
#         endif
          ECORR_SC(iSym,iCase) = etmp
          OVLAPS_SC(iSym,iCase) = otmp
          call mma_deallocate(LBD)
          call mma_deallocate(LID)
          call RHS_FREE(lg_V)
        end if
        cycle
      end if

      iCoupling = SC_NEVPT2_Coupling(iCase)

      if (iCoupling == -1) cycle ! minus coupling has been computed already (SC_Merge)
      if (SC_merge) then
        jCase = iCase+iCoupling
      else
        jCase = 0
      end if

      ! 1. compute VALH (= ValHN(:,1)) and VALN (= ValHN(:,2))
      ! If SC_merge is true, VALN and VALH, and thus E_value, are calculated using both plus and minus couplings

      call mma_allocate(ValHN,NIS,2,Label='ValHN')
      call SC_NEVPT2_ValHN(iCase,iSym,NAS,NIS,ValHN,2)

      NAS2 = 0
      NIS2 = 0
      if (iCoupling == +1) then
        ! minus coupling for SC_merge
        NAS2 = NASUP(iSym,jCase)
        NIS2 = NISUP(iSym,jCase)
        if (NAS2*NIS2 > 0) then
          call mma_allocate(ValHN2,NIS2,2,Label='ValHN2')
          call SC_NEVPT2_ValHN(jCase,iSym,NAS2,NIS2,ValHN2,2)

          call mma_allocate(Map,NIS2,Label='Map')
          call SC_NEVPT2_Map(jCase,iSym,NIS2,Map)
          do isp=1,NIS2
            ValHN(Map(isp),:) = ValHN(Map(isp),:)+ValHN2(isp,:)
          end do
          call mma_deallocate(ValHN2)
        end if
      end if

      ! 2. compute E_value using VALN and VALH
      ! For SC_merge, the minus coupling is included in VALN

      call mma_allocate(LID,NIS,Label='LID')
      idisk = IDBMAT(iSym,iCase)
      call dDaFile(LUSBT,0,Dummy,NAS,idisk)
      call dDaFile(LUSBT,2,LID,NIS,idisk)

      call mma_allocate(E_value,NIS,Label='E_value')

      etmp = Zero
      otmp = Zero
      do isp=1,NIS
        !! compute epsilon (the equation after Eq. (28))
        if (abs(ValHN(isp,2)) <= SC_thres) then
          vale = Zero
        else
          vale = ValHN(isp,1)/ValHN(isp,2)
        end if
        if (abs(lid(isp)+vale) > SC_thres) then
          etmp = etmp-ValHN(isp,2)/(lid(isp)+vale)
          otmp = otmp+ValHN(isp,2)/(lid(isp)+vale)**2
        end if
        E_value(isp) = vale
      end do
      ECORR_SC(iSym,iCase) = etmp
      OVLAPS_SC(iSym,iCase) = otmp

      call mma_deallocate(ValHN)

      ! 3. Save the denominator of the T-amplitude for future use

      ! Add the epsilon term to LID
      LID(1:NIS) = LID(1:NIS)+E_value(1:NIS)
      idisk = IDBMAT_NEVPT2(iSym,iCase,2)
      call dDaFile(LUSBT,1,LID,NIS,idisk)

      call mma_deallocate(LID)

      if ((iCoupling == +1) .and. (NAS2*NIS2 > 0)) then
        ! minus coupling for SC_merge
        call mma_allocate(LID,NIS2,Label='LID')
        idisk = IDBMAT(iSym,jCase)
        call dDaFile(LUSBT,0,Dummy,NAS2,idisk)
        call dDaFile(LUSBT,2,LID,NIS2,idisk)

        do isp=1,NIS2
          LID(isp) = LID(isp)+E_value(Map(isp))
        end do
        idisk = IDBMAT_NEVPT2(iSym,jCase,2)
        call dDaFile(LUSBT,1,LID,NIS2,idisk)

        call mma_deallocate(LID)
        call mma_deallocate(Map)
      end if
      call mma_deallocate(E_value)
    end do
  end do

# ifdef _MOLCAS_MPP_
  !! Gather the contributions for H if needed
  if (is_real_par()) then
    if (do_H) then
      call GADGOP(ECORR_SC(:,12:13),9*2,'+')
      call GADGOP(OVLAPS_SC(:,12:13),9*2,'+')
    end if
  end if
# endif

  return

end subroutine SC_NEVPT2_Energy

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_ReadSB(iCase,iSym,NAS,SMAT,BMAT)

  integer(kind=iwp), intent(in) :: iCase, iSym, NAS
  real(kind=wp), intent(out) :: SMAT(nTri_Elem(NAS)), BMAT(nTri_Elem(NAS))

  integer(kind=iwp) :: idisk
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: ii, jj, lg_S
  real(kind=wp), allocatable :: WRK(:,:)
# endif

  ! Read the original B matrix
  idisk = IDBMAT_NEVPT2(iSym,iCase,1)
  call dDaFile(LUSBT,2,BMAT,nTri_Elem(NAS),idisk)

  ! Read the S matrix
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par() .and. ((iCase == 1) .or. (iCase == 4))) then
    call mma_allocate(WRK,NAS,NAS,Label='WRK')
    call PSBMAT_GETMEM('S',lg_S,NAS)
    call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
    call RHS_GET(NAS,NAS,lg_S,WRK)
    do ii=1,NAS
      do jj=1,ii
        SMAT(iTri(jj,ii)) = WRK(jj,ii)
      end do
    end do
    call PSBMAT_FREEMEM(lg_S)
    call mma_deallocate(WRK)
  else
# endif
    idisk = IDSMAT(iSym,iCase)
    call dDaFile(LUSBT,2,SMAT,nTri_Elem(NAS),idisk)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine SC_NEVPT2_ReadSB

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_ValHN(iCase,iSym,NAS,NIS,ValHN,nVal)

  integer(kind=iwp), intent(in) :: iCase, iSym, NAS, NIS, nVal
  real(kind=wp), intent(out) :: ValHN(NIS,nVal)

  integer(kind=iwp) :: ias, IHI1, ILO1, isp, jas, JHI1, JLO1, lg_C, lg_V
  real(kind=wp), external :: ddot_
  real(kind=wp), allocatable :: BMAT(:), SMAT(:)
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: LDC, LDV, MC, MV, myRank
# endif

  ! Compute H and N
  ! ValHN(:,1): <Psi|V[H,V]|Psi> (Eq. (A2), (A6), ...)
  ! ValHN(:,2): norm (Eq. (11)--(18))
  ! ValHN(:,3): overlap of the RHS with IVECC (only for the gradient)

  ! Read the S matrix and the original B matrix
  call mma_allocate(SMAT,nTri_Elem(NAS),Label='SMAT')
  call mma_allocate(BMAT,nTri_Elem(NAS),Label='BMAT')
  call SC_NEVPT2_ReadSB(iCase,iSym,NAS,SMAT,BMAT)

  call RHS_ALLO(NAS,NIS,lg_V)
  call RHS_READ(NAS,NIS,lg_V,iCase,iSym,IVECW)
  if (nVal == 3) then
    call RHS_ALLO(NAS,NIS,lg_C)
    call RHS_READ(NAS,NIS,lg_C,iCase,iSym,IVECC)
  end if

  ILO1 = 1
  IHI1 = NAS
  JLO1 = 1
  JHI1 = NIS

  ValHN(:,:) = Zero
# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    myRank = GA_NodeID()
    call GA_DISTRIBUTION(lg_V,myRank,ILO1,IHI1,JLO1,JHI1)
    if ((ILO1 > 0) .and. (JLO1 > 0)) then
      call GA_Access(lg_V,ILO1,IHI1,JLO1,JHI1,MV,LDV)
      if (IHI1-ILO1+1 /= NAS) then
        write(u6,'(1x,A)') 'Something is wrong in SC_NEVPT2_ValHN'
        call abend()
      end if
      !! The local patch is passed by argument association on purpose:
      !! DBL_MB is declared with size 2 (mafdecls.fh), so indexing it directly in loops is out-of-bounds and formally undefined;
      !! aggressive loop optimizers may derive wrong trip counts from it (e.g., GCC13)
      call SC_NEVPT2_ValHN_MPP(DBL_MB(MV),LDV,JLO1,JHI1,NAS,NIS,BMAT,SMAT,ValHN,nVal)
      if (nVal == 3) then
        call GA_Access(lg_C,ILO1,IHI1,JLO1,JHI1,MC,LDC)
        call SC_NEVPT2_ValHN_MPP2(DBL_MB(MV),DBL_MB(MC),JLO1,JHI1,LDV,LDC,NAS,NIS,ValHN(1,3))
        call GA_Release(lg_C,ILO1,IHI1,JLO1,JHI1)
      end if
      call GA_RELEASE(lg_V,ILO1,IHI1,JLO1,JHI1)
    end if
    call GADGOP(ValHN,NIS*nVal,'+')
  else
# endif
    do isp=JLO1,JHI1
      do ias=ILO1,IHI1
        do jas=ias+1,IHI1
          ValHN(isp,1) = ValHN(isp,1)+Two*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(jas+NAS*(isp-1))*BMAT(iTri(jas,ias))
          ValHN(isp,2) = ValHN(isp,2)+Two*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(jas+NAS*(isp-1))*SMAT(iTri(jas,ias))
        end do
        ValHN(isp,1) = ValHN(isp,1)+GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*BMAT(nTri_Elem(ias))
        ValHN(isp,2) = ValHN(isp,2)+GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*SMAT(nTri_Elem(ias))
      end do
      if (nVal == 3) ValHN(isp,3) = DDot_(NAS,GA_Arrays(lg_V)%A(1+NAS*(isp-1)),1,GA_Arrays(lg_C)%A(1+NAS*(isp-1)),1)
    end do
# ifdef _MOLCAS_MPP_
  end if
# endif

  call mma_deallocate(BMAT)
  call mma_deallocate(SMAT)
  call RHS_FREE(lg_V)
  if (nVal == 3) call RHS_FREE(lg_C)

end subroutine SC_NEVPT2_ValHN

!-----------------------------------------------------------------------

# ifdef _MOLCAS_MPP_
subroutine SC_NEVPT2_Energy_H_MPP(V,LDV,JLO,JHI,NAS,NIS,LBD,LID,etmp,otmp)

  integer(kind=iwp), intent(in) :: LDV, JLO, JHI, NAS, NIS
  real(kind=wp), intent(in) :: V(LDV,JHI-JLO+1), LBD(NAS), LID(NIS)
  real(kind=wp), intent(inout) :: etmp, otmp

  integer(kind=iwp) :: ias, isp

  do isp=JLO,JHI
    do ias=1,NAS
      etmp = etmp-V(ias,isp-JLO+1)*V(ias,isp-JLO+1)/(LID(isp)+LBD(ias))
      otmp = otmp+V(ias,isp-JLO+1)*V(ias,isp-JLO+1)/(LID(isp)+LBD(ias))**2
    end do
  end do

end subroutine SC_NEVPT2_Energy_H_MPP

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_ValHN_MPP(V,LDV,JLO,JHI,NAS,NIS,BMAT,SMAT,ValHN,nVal)

  integer(kind=iwp), intent(in) :: LDV, JLO, JHI, NAS, NIS, nVal
  real(kind=wp), intent(in) :: V(LDV,JHI-JLO+1), BMAT(nTri_Elem(NAS)), SMAT(nTri_Elem(NAS))
  real(kind=wp), intent(inout) :: ValHN(NIS,nVal)

  integer(kind=iwp) :: ias, isp, jas

  do isp=JLO,JHI
    do ias=1,NAS
      do jas=1,ias-1
        ValHN(isp,1) = ValHN(isp,1)+Two*V(ias,isp-JLO+1)*V(jas,isp-JLO+1)*BMAT(iTri(ias,jas))
        ValHN(isp,2) = ValHN(isp,2)+Two*V(ias,isp-JLO+1)*V(jas,isp-JLO+1)*SMAT(iTri(ias,jas))
      end do
      ValHN(isp,1) = ValHN(isp,1)+V(ias,isp-JLO+1)*V(ias,isp-JLO+1)*BMAT(nTri_Elem(ias))
      ValHN(isp,2) = ValHN(isp,2)+V(ias,isp-JLO+1)*V(ias,isp-JLO+1)*SMAT(nTri_Elem(ias))
    end do
  end do

end subroutine SC_NEVPT2_ValHN_MPP

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_ValHN_MPP2(V1,V2,JLO,JHI,LDV,LDC,NAS,NIS,ValHN)

  integer(kind=iwp), intent(in) :: JLO, JHI, LDV, LDC, NAS, NIS
  real(kind=wp), intent(in) :: V1(LDV,JHI-JLO+1), V2(LDC,JHI-JLO+1)
  real(kind=wp), intent(inout) :: ValHN(NIS)

  integer(kind=iwp) :: isp
  real(kind=wp), external :: ddot_

  do isp=JLO,JHI
    ValHN(isp) = ddot_(NAS,V1(1,isp-JLO+1),1,V2(1,isp-JLO+1),1)
  end do

end subroutine SC_NEVPT2_ValHN_MPP2
#endif

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_Amplitude(NAS,NIS,iCase,iSym,lg_V)

  integer(kind=iwp), intent(in) :: NAS, NIS, iCase, iSym, lg_V
  integer(kind=iwp) :: ias, idisk, IHI, ILO, isp, JHI, JLO
  real(kind=wp) :: scal
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: LDV, MV, myRank
# endif
  real(kind=wp), allocatable :: LBD(:), LID(:)

  ! Compute the T-amplitude from the right hand side (lg_V = <Phi|H|Psi>)
  ! Different from CASPT2 and PC-NEVPT2, the amplitude and RHS must be in the contravariant MO basis

  if (NAS*NIS == 0) return

  call mma_allocate(LBD,NAS,Label='LBD')
  call mma_allocate(LID,NIS,Label='LID')

  idisk = IDBMAT(iSym,iCase)
  call dDaFile(LUSBT,2,LBD,NAS,idisk)
  call dDaFile(LUSBT,2,LID,NIS,idisk)

  if (iCase <= 11) then
    LBD(1:NAS) = Zero
    idisk = IDBMAT_NEVPT2(iSym,iCase,2)
    call dDaFile(LUSBT,2,LID,NIS,idisk)
  end if

  ILO = 1
  IHI = NAS
  JLO = 1
  JHI = NIS

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    myRank = GA_NodeID()
    call GA_DISTRIBUTION(lg_V,myRank,ILO,IHI,JLO,JHI)
    if ((ILO > 0) .and. (JLO > 0)) then
      call GA_Access(lg_V,ILO,IHI,JLO,JHI,MV,LDV)
      !! The local patch is passed by argument association on purpose
      if (iCase <= 11) then
        call SC_NEVPT2_Amplitude_MPP(DBL_MB(MV),LDV,ILO,IHI,JLO,JHI,LBD,LID,.false.)
      else if ((iCase == 12) .or. (iCase == 13)) then
        call SC_NEVPT2_Amplitude_MPP(DBL_MB(MV),LDV,ILO,IHI,JLO,JHI,LBD,LID,.true.)
      end if
      call GA_Release_Update(lg_V,ILO,IHI,JLO,JHI)
    end if
  else
# endif
    if (iCase <= 11) then
      do isp=JLO,JHI
        if (abs(LID(isp)) <= SC_thres) then
          scal = Zero
        else
          scal = -One/LID(isp)
        end if
        GA_Arrays(lg_V)%A(NAS*(isp-1)+ILO:NAS*(isp-1)+IHI) = scal*GA_Arrays(lg_V)%A(NAS*(isp-1)+ILO:NAS*(isp-1)+IHI)
      end do
    else if ((iCase == 12) .or. (iCase == 13)) then
      do isp=JLO,JHI
        do ias=ILO,IHI
          if (abs(LBD(ias)+LID(isp)) <= SC_thres) then
            !! it is the MP2-like term, this should not happen...
            scal = Zero
          else
            scal = -One/(LBD(ias)+LID(isp))
          end if
          GA_Arrays(lg_V)%A(ias+NAS*(isp-1)) = scal*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))
        end do
      end do
    end if
# ifdef _MOLCAS_MPP_
  end if
# endif

  call mma_deallocate(LBD)
  call mma_deallocate(LID)

end subroutine SC_NEVPT2_Amplitude

!-----------------------------------------------------------------------

# ifdef _MOLCAS_MPP_
subroutine SC_NEVPT2_Amplitude_MPP(V,LDV,ILO,IHI,JLO,JHI,LBD,LID,use_bd)

  integer(kind=iwp), intent(in) :: LDV, ILO, IHI, JLO, JHI
  real(kind=wp), intent(inout) :: V(LDV,JHI-JLO+1)
  real(kind=wp), intent(in) :: LBD(IHI), LID(JHI)
  logical(kind=iwp), intent(in) :: use_bd
  integer(kind=iwp) :: ias, isp
  real(kind=wp) :: scal

  ! Scale the local patch by the (negative inverse) denominators:
  ! T = -RHS/denominator. With use_bd, the denominator is LBD+LID
  ! (cases 12 and 13); otherwise it is LID alone (cases 1 to 11).

  do isp=JLO,JHI
    if (use_bd) then
      do ias=ILO,IHI
        if (abs(LBD(ias)+LID(isp)) <= SC_thres) then
          scal = Zero
        else
          scal = -One/(LBD(ias)+LID(isp))
        end if
        V(ias-ILO+1,isp-JLO+1) = V(ias-ILO+1,isp-JLO+1)*scal
      end do
    else
      if (abs(LID(isp)) <= SC_thres) then
        scal = Zero
      else
        scal = -One/LID(isp)
      end if
      V(1:IHI-ILO+1,isp-JLO+1) = V(1:IHI-ILO+1,isp-JLO+1)*scal
    end if
  end do

end subroutine SC_NEVPT2_Amplitude_MPP
#endif

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,VECROT)

  use BDerNEV, only: BDN_G3, BDNA, BDNB, BDNC, BDND, BDNE, BDNF, BDNG

  integer(kind=iwp), intent(in) :: NASHT, NG3, NSTATE
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT), G3(NG3), VECROT(NSTATE)
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT), DG3(NG3)

  integer(kind=iwp) :: ias, iCase, iCoupling, idisk, isp, iSym, jas, jCase, JHI1, JLO1, lg_V1, lg_V2, lg_V3, NAS, NAS2, NIS, NIS2, &
                       NVEC

  integer(kind=iwp), allocatable :: Map(:)
  real(kind=wp), allocatable :: BDER(:,:), BMAT(:), derHNS(:,:), LID(:), SDER(:,:), SMAT(:), ValHN(:,:), ValHN2(:,:)

# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHI1, IHI2, IHI3, ILO1, ILO2, ILO3, JHI2, JHI3, JLO2, JLO3, LDV, LDV1, LDV2, MV1, MV2, MV3, myRank
# endif

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) call RHS_ZERO(iVecC2)
# endif

  do iCase=1,MXCASE
    !cycle
    !if ((icase /= 12) .and. (icase /= 13)) cycle ! H
    !if ((icase /= 10) .and. (icase /= 11)) cycle ! G
    !if ((icase /= 6) .and. (icase /= 7)) cycle ! E
    !if ((icase /= 8) .and. (icase /= 9)) cycle ! F
    !if ((icase /= 2) .and. (icase /= 3)) cycle ! B
    !if (icase /= 5) cycle ! D
    !if (icase /= 4) cycle ! C
    !if (icase /= 1) cycle ! A
    do iSym=1,nSym
      NIS = NISUP(iSym,iCase)
      NAS = NASUP(iSym,iCase)
      NVEC = NAS*NIS
      if (NVEC == 0) cycle

      if ((iCase == 12) .or. (iCase == 13)) then
        call RHS_ALLO(NAS,NIS,lg_V1)
        call RHS_READ(NAS,NIS,lg_V1,iCase,iSym,iVecC)
        call SC_NEVPT2_amplitude(NAS,NIS,iCase,iSym,lg_V1)
        call RHS_SAVE(NAS,NIS,lg_V1,iCase,iSym,iVecC2)
        call RHS_FREE(lg_V1)
        cycle
      end if

      ! Read the S matrix and the original B matrix
      call mma_allocate(SMAT,nTri_Elem(NAS),Label='SMAT')
      call mma_allocate(BMAT,nTri_Elem(NAS),Label='BMAT')
      call SC_NEVPT2_ReadSB(iCase,iSym,NAS,SMAT,BMAT)

      iCoupling = SC_NEVPT2_Coupling(iCase)

      if (SC_merge) then
        jCase = iCase + iCoupling
      else
        jCase = 0
      end if

      !! Construct the T-amplitude in lg_V1
      call RHS_ALLO(NAS,NIS,lg_V1)
      call RHS_READ(NAS,NIS,lg_V1,iCase,iSym,iVecW)
      call SC_NEVPT2_Amplitude(NAS,NIS,iCase,iSym,lg_V1)
      !! Construct T+lambda in lg_V2 with (S+S')/2
      call RHS_ALLO(NAS,NIS,lg_V2)
      call RHS_READ(NAS,NIS,lg_V2,iCase,iSym,iVecC)
      call SC_NEVPT2_amplitude(NAS,NIS,iCase,iSym,lg_V2)

      call mma_allocate(LID,max(NAS,NIS),Label='LID')
      idisk = IDBMAT_NEVPT2(iSym,iCase,2)
      call dDaFile(LUSBT,2,LID,NIS,idisk)

      call mma_allocate(derHNS,NIS,3,Label='derHNS')
      derHNS(:,:) = Zero

      !! Derivative factors of H and N from the norm, the Dyall expectation value, and the overlap with T+lambda.
      !! The quantities of the amplitudes are those of the RHS divided by the squared denominator, which cancels in the factors.
      !! The two cases of a merged pair share the factors made from the sums over the pair.
      call mma_allocate(ValHN,NIS,3,Label='ValHN')
      call SC_NEVPT2_ValHN(iCase,iSym,NAS,NIS,ValHN,3)

      NAS2 = 0
      NIS2 = 0
      if (iCoupling /= 0) then
        ! the other coupling for SC_merge
        NAS2 = NASUP(iSym,jCase)
        NIS2 = NISUP(iSym,jCase)
        if (NAS2*NIS2 /= 0) then
          call mma_allocate(ValHN2,NIS2,3,Label='ValHN2')
          call SC_NEVPT2_ValHN(jCase,iSym,NAS2,NIS2,ValHN2,3)

          if (iCoupling == +1) then
            ! add the minus coupling to the plus coupling
            call mma_allocate(Map,NIS2,Label='Map')
            call SC_NEVPT2_Map(jCase,iSym,NIS2,Map)
            do isp=1,NIS2
              ValHN(Map(isp),:) = ValHN(Map(isp),:)+ValHN2(isp,:)
            end do
          else if (iCoupling == -1) then
            ! add the plus coupling to the minus coupling
            call mma_allocate(Map,NIS,Label='Map')
            call SC_NEVPT2_Map(iCase,iSym,NIS,Map)
            do isp=1,NIS
              ValHN(isp,:) = ValHN(isp,:)+ValHN2(Map(isp),:)
            end do
          end if
          call mma_deallocate(ValHN2)
          call mma_deallocate(Map)
        end if
      end if
      do isp=1,NIS
        if (abs(ValHN(isp,2)) <= SC_thres) cycle
        ! derivative of H (through vale)
        derHNS(isp,1) = ValHN(isp,3)/ValHN(isp,2)
        ! derivative of N (through vale)
        derHNS(isp,2) = -ValHN(isp,3)*ValHN(isp,1)/ValHN(isp,2)**2
        ! derivative of N: the explicit derivative (numerator)
        ! Compensate the division by lid
        derHNS(isp,3) = -lid(isp)
      end do
      call mma_deallocate(ValHN)

      JLO1 = 1
      JHI1 = NIS

      call mma_allocate(BDER,NAS,NAS,Label='BDER')
      call mma_allocate(SDER,NAS,NAS,Label='SDER')
      BDER(:,:) = Zero
      SDER(:,:) = Zero

      !! Construct BDER and SDER
#     ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        myRank = GA_NodeID()
        call GA_DISTRIBUTION(lg_V1,myRank,ILO1,IHI1,JLO1,JHI1)
        if ((ILO1 > 0) .and. (JLO1 > 0)) then
          call GA_Access(lg_V1,ILO1,IHI1,JLO1,JHI1,MV1,LDV1)
          if (IHI1-ILO1+1 /= NAS) then
            write(u6,'(1x,A)') 'Something is wrong in SC_NEVPT2_CLagD'
            call abend()
          end if
          call SC_NEVPT2_CLagD1_MPP(DBL_MB(MV1),LDV1,JLO1,JHI1,NAS,NIS,derHNS,VECROT(jStLag),BDER,SDER)
          call GA_RELEASE(lg_V1,ILO1,IHI1,JLO1,JHI1)
        end if
      else
#     endif
        do isp=JLO1,JHI1
          do ias=1,NAS
            do jas=1,NAS
              BDER(ias,jas) = BDER(ias,jas)+derHNS(isp,1)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
              SDER(ias,jas) = SDER(ias,jas)+derHNS(isp,2)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
              SDER(ias,jas) = SDER(ias,jas)+ &
                              VECROT(jStLag)*derHNS(isp,3)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
            end do
          end do
        end do
#     ifdef _MOLCAS_MPP_
      end if
#     endif

#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call GADGOP(BDER,NAS**2,'+')
        call GADGOP(SDER,NAS**2,'+')
      end if
#     endif

      !! BDER and SDER --> DG1, DG2, DG3
      if (iCase == 1) call BDNA(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
      if ((iCase == 2) .or. (iCase == 3)) call BDNB(iSym,iCase,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
      if (iCase == 4) call BDNC(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
      if (iCase == 5) call BDND(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
      if ((iCase == 6) .or. (iCase == 7)) call BDNE(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)
      if ((iCase == 8) .or. (iCase == 9)) call BDNF(iSym,iCase,NAS,NG3,BDER,SDER,G2,G3,DG2,DG3)
      if ((iCase == 10) .or. (iCase == 11)) call BDNG(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)

      call mma_deallocate(BDER)
      call mma_deallocate(SDER)

      !! Construct integral derivatives
      call RHS_ALLO(NAS,NIS,lg_V3)
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call RHS_READ(NAS,NIS,lg_V3,iCase,iSym,IVECC2) ! dummy read
        call GA_Sync()
        myRank = GA_NodeID()
        call GA_DISTRIBUTION(lg_V1,myRank,ILO1,IHI1,JLO1,JHI1)
        call GA_DISTRIBUTION(lg_V2,myRank,ILO2,IHI2,JLO2,JHI2)
        call GA_DISTRIBUTION(lg_V3,myRank,ILO3,IHI3,JLO3,JHI3)
        if ((ILO1 > 0) .and. (JLO1 > 0)) then
          call GA_Access(lg_V1,ILO1,IHI1,JLO1,JHI1,MV1,LDV1)
          call GA_Access(lg_V2,ILO2,IHI2,JLO2,JHI2,MV2,LDV2)
          call GA_Access(lg_V3,ILO3,IHI3,JLO3,JHI3,MV3,LDV)
          !! The local patches are passed by argument association on purpose
          call SC_NEVPT2_CLagD2_MPP(DBL_MB(MV1),LDV1,DBL_MB(MV2),LDV2,DBL_MB(MV3),LDV,JLO1,JHI1,NAS,NIS,BMAT,SMAT,LID,derHNS)
          call GA_Release(lg_V1,ILO1,IHI1,JLO1,JHI1)
          call GA_Release(lg_V2,ILO2,IHI2,JLO2,JHI2)
          call GA_Release_Update(lg_V3,ILO3,IHI3,JLO3,JHI3)
        end if
      else
#     endif
        GA_Arrays(lg_V3)%A(:) = Zero
        do isp=1,NIS
          if (abs(LID(isp)) <= SC_thres) cycle
          derHNS(isp,1) = -derHNS(isp,1)/LID(isp)
          derHNS(isp,2) = -derHNS(isp,2)/LID(isp)
          derHNS(isp,3) = One
          do ias=1,NAS
            do jas=1,NAS
              ! <lambda|H0-E0|Psi1> terms
              GA_Arrays(lg_V3)%A(ias+NAS*(isp-1)) = GA_Arrays(lg_V3)%A(ias+NAS*(isp-1))+ &
                                                    derHNS(isp,1)*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))*BMAT(iTri(ias,jas))
              ! <lambda|H0-E0|Psi1> terms
              GA_Arrays(lg_V3)%A(ias+NAS*(isp-1)) = GA_Arrays(lg_V3)%A(ias+NAS*(isp-1))+ &
                                                    derHNS(isp,2)*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))*SMAT(iTri(ias,jas))
            end do
            ! <lambda|H|Psi0> term(s)
            GA_Arrays(lg_V3)%A(ias+NAS*(isp-1)) = GA_Arrays(lg_V3)%A(ias+NAS*(isp-1))+ &
                                                  derHNS(isp,3)*GA_Arrays(lg_V2)%A(ias+NAS*(isp-1))
          end do
        end do
#     ifdef _MOLCAS_MPP_
      end if
#     endif
      call RHS_SAVE(NAS,NIS,lg_V3,iCase,iSym,iVecC2)

      call RHS_FREE(lg_V1)
      call RHS_FREE(lg_V2)
      call RHS_FREE(lg_V3)

      call mma_deallocate(LID)
      call mma_deallocate(BMAT)
      call mma_deallocate(SMAT)
      call mma_deallocate(derHNS)
    end do
  end do

  !! Compute correct DG3 contributions
  !! E4 contributions are evaluated elsewhere
  call BDN_G3(DG1,DG2,DG3)

end subroutine SC_NEVPT2_CLagD

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_res(VECROT)

  real(kind=wp), intent(in) :: VECROT(*)
  integer(kind=iwp) :: iCase, iStLag, iSym, lg_V1, lg_V2, NAS, NIS, NVEC
  real(kind=wp) :: scal

  !! Construct the residual vector
  !! IVECC will contain (S+S')/2*T (S' = S transpose) but icase = 12 and 13 is omitted

  !! derivative of the Jth state wrt the T-amplitude of Jth state
  !! (residue) = VECROT(jStLag)*(\partial <Psi1|H|Psi0>)/(\partial T)
  !! S matrix is symmetric
  call RHS_ZERO(IVECC)
  do iSym=1,nSym
    do iCase=1,MXCASE
      NAS = NASUP(iSym,iCase)
      NIS = NISUP(iSym,iCase)
      NVEC = NAS*NIS
      if (NVEC == 0) cycle

      if ((iCase == 12) .or. (iCase == 13)) then
        call RHS_ALLO(NAS,NIS,lg_V1)
        !! Read RHS (residual of the target state)
        call RHS_READ(NAS,NIS,lg_V1,iCase,iSym,IVECW)
        !! Compute appropriate contributions
        call RHS_SCAL(NAS,NIS,lg_V1,VECROT(jStLag))
        !! Save the residual vector
        call RHS_SAVE(NAS,NIS,lg_V1,iCase,iSym,IVECC)
        call RHS_FREE(lg_V1)
      else
        call RHS_ALLO(NAS,NIS,lg_V1)
        call RHS_ALLO(NAS,NIS,lg_V2)
        !! Read RHS (residual of the target state)
        call RHS_READ(NAS,NIS,lg_V1,iCase,iSym,IVECW)
        call RHS_SCAL(NAS,NIS,lg_V2,Zero)
        !! Compute appropriate contributions
        call RHS_STRANS(NAS,NIS,VECROT(jStLag),lg_V1,lg_V2,ICASE,ISYM)
        !! Save the residual vector
        call RHS_SAVE(NAS,NIS,lg_V2,iCase,iSym,IVECC)
        call RHS_FREE(lg_V1)
        call RHS_FREE(lg_V2)
      end if
    end do
  end do

  !! derivative of the non-Jth states
  do iStLag=1,nState
    scal = VECROT(iStLag)
    if (iStLag == jStLag) cycle
    if (abs(scal) <= 1.0e-12_wp) cycle
    call MS_Res(1,iStLag,jStLag,scal)
  end do

end subroutine SC_NEVPT2_res

!-----------------------------------------------------------------------

# ifdef _MOLCAS_MPP_
subroutine SC_NEVPT2_CLagD1_MPP(V1,LDV1,JLO,JHI,NAS,NIS,derHNS,scal,BDER,SDER)

  integer(kind=iwp), intent(in) :: LDV1, JLO, JHI, NAS, NIS
  real(kind=wp), intent(in) :: V1(LDV1,JHI-JLO+1), derHNS(NIS,3), scal
  real(kind=wp), intent(inout) :: BDER(NAS,NAS), SDER(NAS,NAS)
  integer(kind=iwp) :: ias, icol, isp, jas

  ! Accumulate the local-column contributions to BDER and SDER

  do isp=JLO,JHI
    icol = isp-JLO+1
    do jas=1,NAS
      do ias=1,NAS
        BDER(ias,jas) = BDER(ias,jas)+derHNS(isp,1)*V1(ias,icol)*V1(jas,icol)
        SDER(ias,jas) = SDER(ias,jas)+derHNS(isp,2)*V1(ias,icol)*V1(jas,icol)
        SDER(ias,jas) = SDER(ias,jas)+scal*derHNS(isp,3)*V1(ias,icol)*V1(jas,icol)
      end do
    end do
  end do

end subroutine SC_NEVPT2_CLagD1_MPP

!-----------------------------------------------------------------------

subroutine SC_NEVPT2_CLagD2_MPP(V1,LDV1,V2,LDV2,V3,LDV3,JLO,JHI,NAS,NIS,BMAT,SMAT,LID,derHNS)

  integer(kind=iwp), intent(in) :: LDV1, LDV2, LDV3, JLO, JHI, NAS, NIS
  real(kind=wp), intent(in) :: V1(LDV1,JHI-JLO+1), V2(LDV2,JHI-JLO+1), BMAT(nTri_Elem(NAS)), SMAT(nTri_Elem(NAS)), LID(NIS)
  real(kind=wp), intent(inout) :: V3(LDV3,JHI-JLO+1), derHNS(NIS,3)
  integer(kind=iwp) :: ias, icol, isp, jas

  ! Construct the integral derivative vector V3 from V1, V2, and derHNS

  do isp=JLO,JHI
    icol = isp-JLO+1
    V3(1:NAS,icol) = Zero
    if (abs(LID(isp)) <= SC_thres) cycle
    derHNS(isp,1) = -derHNS(isp,1)/LID(isp)
    derHNS(isp,2) = -derHNS(isp,2)/LID(isp)
    derHNS(isp,3) = One
    do ias=1,NAS
      do jas=1,ias-1
        ! <lambda|H0-E0|Psi1> terms
        V3(jas,icol) = V3(jas,icol)+derHNS(isp,1)*V1(ias,icol)*BMAT(iTri(ias,jas))
        V3(jas,icol) = V3(jas,icol)+derHNS(isp,2)*V1(ias,icol)*SMAT(iTri(ias,jas))
        V3(ias,icol) = V3(ias,icol)+derHNS(isp,1)*V1(jas,icol)*BMAT(iTri(ias,jas))
        V3(ias,icol) = V3(ias,icol)+derHNS(isp,2)*V1(jas,icol)*SMAT(iTri(ias,jas))
      end do
      V3(ias,icol) = V3(ias,icol)+derHNS(isp,1)*V1(ias,icol)*BMAT(nTri_Elem(ias))
      V3(ias,icol) = V3(ias,icol)+derHNS(isp,2)*V1(ias,icol)*SMAT(nTri_Elem(ias))
      ! <lambda|H|Psi0> term(s)
      V3(ias,icol) = V3(ias,icol)+V2(ias,icol)
    end do
  end do

end subroutine SC_NEVPT2_CLagD2_MPP
#endif

end module SC_NEVPT2

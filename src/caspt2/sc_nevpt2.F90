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
logical(kind=iwp) :: Do_FIC, Do_SC, SC_amplitude = .false., SC_prop = .false.
real(kind=wp), allocatable :: ENERGY_SC(:), HEFF_SC(:,:), UEFF_SC(:,:)

public :: Do_FIC, Do_SC, DVALUE_SC, ECORR_SC, ENERGY_SC, HEFF_SC, IDBMAT_NEVPT2, OVLAPS_SC, SC_amplitude, SC_NEVPT2_amplitude, &
          SC_NEVPT2_CLagD, SC_NEVPT2_final, SC_NEVPT2_initial, SC_NEVPT2_Print, SC_NEVPT2_res, SC_prop, SC_thres, UEFF_SC

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

subroutine SC_NEVPT2_Print()

  use PrintLevel, only: TERSE, USUAL

  integer(kind=iwp) :: IC, iCase, IS, iSym, LAXITY
  real(kind=wp) :: DENORM_SC, E2CORR_SC, E2NONV_SC, E2TOT_SC, EAIVX_SC, EATVX_SC, EBJAI_SC, EBJAT_SC, EBVAT_SC, EVJAI_SC, &
                   EVJTI_SC, EVJTU_SC, REFWGT_SC, RNORM_SC
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
      write(u6,'(1X,A8,9F12.8)') CASES(IC),(ECORR_SC(IS,IC),IS=1,NSYM),ECORR_SC(0,IC)
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

  integer(kind=iwp) :: ias, iCase, idisk, IHI1, ILO1, isp, iSym, jas, JHI1, JLO1, lg_V, NAS, NIS
  real(kind=wp) :: etmp, otmp, vale, valh, valn
  logical(kind=iwp) :: do_H
  real(kind=wp), allocatable :: BMAT(:), E_value(:), LBD(:), LID(:), SMAT(:)
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: ii, jj, LDV, lg_S, MV, myRank
  real(kind=wp), allocatable :: WRK(:,:)
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
        !! Compute the H subspace only when real/imaginary shift has been used
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
              do isp=JLO1,JHI1
                valn = Zero
                do ias=1,NAS
                  etmp = etmp-DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+ias-1+NAS*(isp-JLO1))/(lid(isp)+lbd(ias))
                  otmp = otmp+DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+ias-1+NAS*(isp-JLO1))/(lid(isp)+lbd(ias))**2
                end do
              end do
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

      ! Read the original B matrix
      call mma_allocate(BMAT,nTri_Elem(NAS),Label='BMAT')
      idisk = IDBMAT_NEVPT2(iSym,iCase,1)
      call dDaFile(LUSBT,2,BMAT,nTri_Elem(NAS),idisk)

      ! Read the S matrix
      call mma_allocate(SMAT,nTri_Elem(NAS),Label='SMAT')
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par() .and. ((iCase == 1) .or. (iCase == 4))) then
        call mma_allocate(WRK,NAS,NAS,Label='WRK')
        call PSBMAT_GETMEM('S',lg_S,NAS)
        call PSBMAT_READ('S',iCase,iSYM,lg_S,NAS)
        call RHS_GET(NAS,NAS,lg_S,WRK)
        do ii=1,NAS
          do jj=ii,NAS
            SMAT(iTri(jj,ii)) = WRK(jj,ii)
          end do
        end do
        call PSBMAT_FREEMEM(lg_S)
        call mma_deallocate(WRK)
      else
#     endif
        idisk = IDSMAT(iSym,iCase)
        call dDaFile(LUSBT,2,SMAT,nTri_Elem(NAS),idisk)
#     ifdef _MOLCAS_MPP_
      end if
#     endif

      call RHS_ALLO(NAS,NIS,lg_V)
      call RHS_READ(NAS,NIS,lg_V,iCase,iSym,IVECW)

      call mma_allocate(LBD,NAS,Label='LBD')
      call mma_allocate(LID,NIS,Label='LID')
      idisk = IDBMAT(iSym,iCase)
      call dDaFile(LUSBT,2,LBD,NAS,idisk)
      call dDaFile(LUSBT,2,LID,NIS,idisk)

      call mma_allocate(E_value,NIS,Label='E_value')
      E_value(1:NIS) = Zero

      ILO1 = 1
      IHI1 = NAS
      JLO1 = 1
      JHI1 = NIS

      etmp = Zero
      otmp = Zero
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        myRank = GA_NodeID()
        call GA_DISTRIBUTION(lg_V,myRank,ILO1,IHI1,JLO1,JHI1)
        if ((ILO1 > 0) .and. (JLO1 > 0)) then
          call GA_Access(lg_V,ILO1,IHI1,JLO1,JHI1,MV,LDV)
          if (IHI1-ILO1+1 /= NAS) then
            write(u6,'(1x,A)') 'Something is wrong in SC_NEVPT2_Energy'
            call abend()
          end if
          do isp=JLO1,JHI1
            valh = Zero
            valn = Zero
            do ias=1,NAS
              do jas=1,ias-1 ! ias+1, NAS
                valh = valh+Two*DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+jas-1+NAS*(isp-JLO1))*BMAT(iTri(ias,jas))
                valn = valn+Two*DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+jas-1+NAS*(isp-JLO1))*SMAT(iTri(ias,jas))
              end do
              valh = valh+DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+ias-1+NAS*(isp-JLO1))*BMAT(nTri_Elem(ias))
              valn = valn+DBL_MB(MV+ias-1+NAS*(isp-JLO1))*DBL_MB(MV+ias-1+NAS*(isp-JLO1))*SMAT(nTri_Elem(ias))
            end do
            if (abs(valn) <= SC_thres) then
              vale = Zero
            else
              vale = valh/valn
            end if
            if (abs(lid(isp)+vale) > SC_thres) then
              etmp = etmp-valn/(lid(isp)+vale)
              otmp = otmp+valn/(lid(isp)+vale)**2
            end if
            E_value(isp) = vale
          end do
          call GA_RELEASE(lg_V,ILO1,IHI1,JLO1,JHI1)
        end if
        call GADGOP(E_value,NIS,'+')
      else
#     endif
        do isp=JLO1,JHI1
          valh = Zero !! <Psi|V[H,V]|Psi> (Eq. (A2), (A6), ...)
          valn = Zero !! N: norm (Eq. (11)--(18))
          do ias=ILO1,IHI1
            do jas=ias+1,IHI1
              valh = valh+Two*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(jas+NAS*(isp-1))*BMAT(iTri(jas,ias))
              valn = valn+Two*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(jas+NAS*(isp-1))*SMAT(iTri(jas,ias))
            end do
            valh = valh+GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*BMAT(nTri_Elem(ias))
            valn = valn+GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V)%A(ias+NAS*(isp-1))*SMAT(nTri_Elem(ias))
          end do
          !! compute epsilon (the equation after Eq. (28))
          if (abs(valn) <= SC_thres) then
            vale = Zero
          else
            vale = valh/valn
          end if
          if (abs(lid(isp)+vale) > SC_thres) then
            etmp = etmp-valn/(lid(isp)+vale)
            otmp = otmp+valn/(lid(isp)+vale)**2
          end if
          E_value(isp) = vale
        end do
#     ifdef _MOLCAS_MPP_
      end if
#     endif
      ECORR_SC(iSym,iCase) = etmp
      OVLAPS_SC(iSym,iCase) = otmp

      call mma_deallocate(BMAT)
      call mma_deallocate(SMAT)

      !! Save the denominator of the T-amplitude for future use
      !! LID gets the epsilon term
      LID(1:NIS) = LID(1:NIS)+E_value(1:NIS)
      idisk = IDBMAT_NEVPT2(iSym,iCase,2)
      call dDaFile(LUSBT,1,LID,NIS,idisk)

      call mma_deallocate(LBD)
      call mma_deallocate(LID)
      call mma_deallocate(E_value)

      call RHS_FREE(lg_V)
    end do
  end do

# ifdef _MOLCAS_MPP_
  !! Gather the contributions A~G (not H)
  if (is_real_par()) then
    if (do_H) then
      call GADGOP(ECORR_SC(:,0:13),9*14,'+')
      call GADGOP(OVLAPS_SC(:,0:13),9*14,'+')
    else
      call GADGOP(ECORR_SC(:,0:11),9*12,'+')
      call GADGOP(OVLAPS_SC(:,0:11),9*12,'+')
    end if
  end if
# endif

  return

end subroutine SC_NEVPT2_Energy

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
      if (iCase <= 11) then
        do isp=JLO,JHI
          if (abs(LID(isp)) <= SC_thres) then
            scal = Zero
          else
            scal = -One/LID(isp)
          end if
          do ias=ILO,IHI
            DBL_MB(MV+ias-ILO+LDV*(isp-JLO)) = DBL_MB(MV+ias-ILO+LDV*(isp-JLO))*scal
          end do
        end do
      else if ((iCase == 12) .or. (iCase == 13)) then
        do isp=JLO,JHI
          do ias=ILO,IHI
            if (abs(LBD(ias)+LID(isp)) <= SC_thres) then
              scal = Zero
            else
              scal = -One/(LBD(ias)+LID(isp))
            end if
            DBL_MB(MV+ias-ILO+LDV*(isp-JLO)) = DBL_MB(MV+ias-ILO+LDV*(isp-JLO))*scal
          end do
        end do
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

subroutine SC_NEVPT2_CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,VECROT)

  use BDerNEV, only: BDN_G3, BDNA, BDNB, BDNC, BDND, BDNE, BDNF, BDNG

  integer(kind=iwp), intent(in) :: NASHT, NG3, NSTATE
  real(kind=wp), intent(in) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT), G3(NG3), VECROT(NSTATE)
  real(kind=wp), intent(inout) :: DG1(NASHT,NASHT), DG2(NASHT,NASHT,NASHT,NASHT), DG3(NG3)
  integer(kind=iwp) :: ias, iCase, idisk, isp, iSym, jas, JHI1, JLO1, lg_V1, lg_V2, lg_V3, NAS, NIS, NVEC
  real(kind=wp) :: vale, valh, valn
  real(kind=wp), allocatable :: BDER(:,:), BMAT(:), derHNS(:,:), LID(:), SDER(:,:), SMAT(:)
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: idx_bs, idx_v1, idx_v2, idx_v3, IHI1, IHI2, IHI3, ii, ILO1, ILO2, ILO3, JHI2, JHI3, jj, JLO2, JLO3, LDV, &
                       LDV1, LDV2, lg_S, MV1, MV2, MV3, myRank
  real(kind=wp) :: tmph, tmpn

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

      ! Read the original B matrix
      call mma_allocate(BMAT,nTri_Elem(NAS),Label='BMAT')
      idisk = IDBMAT_NEVPT2(iSym,iCase,1)
      call DDAFILE(LUSBT,2,BMAT,nTri_Elem(NAS),idisk)

      ! Read the S matrix
      call mma_allocate(SMAT,nTri_Elem(NAS),Label='SMAT')
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par() .and. ((iCase == 1) .or. (iCase == 4))) then
        call mma_allocate(SDER,NAS,NAS,Label='SDER')
        call PSBMAT_GETMEM('S',lg_S,NAS)
        call PSBMAT_READ('S',iCase,iSYM,lg_S,NAS)
        call RHS_GET(NAS,NAS,lg_S,SDER)
        do ii=1,NAS
          do jj=1,ii
            SMAT(iTri(ii,jj)) = SDER(jj,ii)
          end do
        end do
        call PSBMAT_FREEMEM(lg_S)
        call mma_deallocate(SDER)
      else
#     endif
        idisk = IDSMAT(iSym,iCase)
        call DDAFILE(LUSBT,2,SMAT,nTri_Elem(NAS),idisk)
#     ifdef _MOLCAS_MPP_
      end if
#     endif

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

#     ifdef _MOLCAS_MPP_
      ILO1 = 1
      IHI1 = NAS
#     endif
      JLO1 = 1
      JHI1 = NIS

      !! Compute derivative of H and N
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        myRank = GA_NodeID()
        call GA_DISTRIBUTION(lg_V1,myRank,ILO1,IHI1,JLO1,JHI1)
        call GA_DISTRIBUTION(lg_V2,myRank,ILO2,IHI2,JLO2,JHI2)
        if ((ILO1 > 0) .and. (JLO1 > 0) .and. (ILO2 > 0) .and. (JLO2 > 0)) then
          call GA_Access(lg_V1,ILO1,IHI1,JLO1,JHI1,MV1,LDV1)
          call GA_Access(lg_V2,ILO2,IHI2,JLO2,JHI2,MV2,LDV2)
          if ((IHI1-ILO1+1 /= NAS) .or. (ILO1 /= ILO2) .or. (IHI1 /= IHI2) .or. (JLO1 /= JLO2) .or. (JHI1 /= JHI2) .or. &
              (LDV1 /= LDV2)) then
            write(u6,'(1x,A)') 'Something is wrong in SC_NEVPT2_CLagD'
            call abend()
          end if
          do isp=JLO1,JHI1
            vale = Zero
            valh = Zero
            valn = Zero
            idx_v1 = MV1-1+NAS*(isp-JLO1)
            idx_v2 = MV2-1+NAS*(isp-JLO2)
            vale = vale+sum(DBL_MB(idx_v1+1:idx_v1+NAS)*DBL_MB(idx_v2+1:idx_v2+NAS))
            do ias=1,NAS
              idx_bs = nTri_Elem(ias-1)
              tmph = Zero
              tmpn = Zero
              if (ias >= 2) then
                tmph = sum(DBL_MB(idx_v1+1:idx_v1+ias-1)*BMAT(idx_bs+1:idx_bs+ias-1))
                tmpn = sum(DBL_MB(idx_v1+1:idx_v1+ias-1)*SMAT(idx_bs+1:idx_bs+ias-1))
              end if
              valh = valh+DBL_MB(idx_v1+ias)*(Two*tmph+DBL_MB(idx_v1+ias)*BMAT(idx_bs+ias))
              valn = valn+DBL_MB(idx_v1+ias)*(Two*tmpn+DBL_MB(idx_v1+ias)*SMAT(idx_bs+ias))
            end do
            !do ias=1,NAS
            !  vale = vale+DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*DBL_MB(MV2+ias-1+NAS*(isp-JLO2))
            !  do jas=1,ias-1
            !    valh = valh+Two*DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*DBL_MB(MV1+jas-1+NAS*(isp-JLO1))*BMAT(iTri(ias,jas))
            !    valn = valn+Two*DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*DBL_MB(MV1+jas-1+NAS*(isp-JLO1))*SMAT(iTri(ias,jas))
            !  end do
            !  valh = valh+DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*BMAT(nTri_Elem(ias))
            !  valn = valn+DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*DBL_MB(MV1+ias-1+NAS*(isp-JLO1))*SMAT(nTri_Elem(ias))
            !end do
            if (abs(valn) <= SC_thres) then
              derHNS(isp,1) = Zero
              derHNS(isp,2) = Zero
              derHNS(isp,3) = Zero
            else
              derHNS(isp,1) = vale/valn
              derHNS(isp,2) = -vale*valh/(valn*valn)
              derHNS(isp,3) = -One*lid(isp)
            end if
          end do
          call GA_RELEASE(lg_V1,ILO1,IHI1,JLO1,JHI1)
          call GA_RELEASE(lg_V2,ILO2,IHI2,JLO2,JHI2)
        end if
        call GADGOP(derHNS,NIS*3,'+')
      else
#     endif
        do isp=1,NIS
          vale = Zero ! derivative contributions of e
          valh = Zero
          valn = Zero
          do ias=1,NAS
            vale = vale+GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V2)%A(ias+NAS*(isp-1))
            do jas=1,ias-1
              valh = valh+Two*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))*BMAT(iTri(ias,jas))
              valn = valn+Two*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))*SMAT(iTri(ias,jas))
            end do
            valh = valh+GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*BMAT(nTri_Elem(ias))
            valn = valn+GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*SMAT(nTri_Elem(ias))
          end do
          if (abs(valn) <= SC_thres) cycle
          ! derivative of H (through vale)
          derHNS(isp,1) = vale/valn
          ! derivative of N (through vale)
          derHNS(isp,2) = -vale*valh/(valn*valn)
          ! derivative of N: the explicit derivative (numerator)
          ! Compensate the division by lid
          derHNS(isp,3) = -One*lid(isp)
        end do
#     ifdef _MOLCAS_MPP_
      end if
#     endif

      call mma_allocate(BDER,NAS,NAS,Label='BDER')
      call mma_allocate(SDER,NAS,NAS,Label='SDER')
      BDER(:,:) = Zero
      SDER(:,:) = Zero

      !! Construct BDER and SDER
      do isp=JLO1,JHI1
        do ias=1,NAS
          do jas=1,NAS
#           ifdef _MOLCAS_MPP_
            if (is_real_par()) then
              BDER(ias,jas) = BDER(ias,jas)+derHNS(isp,1)*DBL_MB(MV1+ias-ILO1+NAS*(isp-JLO1))*DBL_MB(MV1+jas-ILO1+NAS*(isp-JLO1))
              SDER(ias,jas) = SDER(ias,jas)+derHNS(isp,2)*DBL_MB(MV1+ias-ILO1+NAS*(isp-JLO1))*DBL_MB(MV1+jas-ILO1+NAS*(isp-JLO1))
              SDER(ias,jas) = SDER(ias,jas)+ &
                              VECROT(jStLag)*derHNS(isp,3)*DBL_MB(MV1+ias-ILO1+NAS*(isp-JLO1))*DBL_MB(MV1+jas-ILO1+NAS*(isp-JLO2))
            else
#           endif
              BDER(ias,jas) = BDER(ias,jas)+derHNS(isp,1)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
              SDER(ias,jas) = SDER(ias,jas)+derHNS(isp,2)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
              SDER(ias,jas) = SDER(ias,jas)+ &
                              VECROT(jStLag)*derHNS(isp,3)*GA_Arrays(lg_V1)%A(ias+NAS*(isp-1))*GA_Arrays(lg_V1)%A(jas+NAS*(isp-1))
#           ifdef _MOLCAS_MPP_
            end if
#           endif
          end do
        end do
      end do

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
          call GA_Access(lg_V1,ILO1,IHI1,JLO1,JHI1,MV1,LDV)
          call GA_Access(lg_V2,ILO2,IHI2,JLO2,JHI2,MV2,LDV)
          call GA_Access(lg_V3,ILO3,IHI3,JLO3,JHI3,MV3,LDV)
          do isp=JLO1,JHI1
            DBL_MB(MV3+NAS*(isp-JLO1):MV3+NAS+NAS*(isp-JLO1)) = Zero
            if (abs(LID(isp)) <= SC_thres) cycle
            derHNS(isp,1) = -derHNS(isp,1)/LID(isp)
            derHNS(isp,2) = -derHNS(isp,2)/LID(isp)
            derHNS(isp,3) = One
            idx_v1 = MV1-1+NAS*(isp-JLO1)
            idx_v2 = MV2-1+NAS*(isp-JLO1)
            idx_v3 = MV3-1+NAS*(isp-JLO1)
            do ias=1,NAS
              do jas=1,ias-1
                ! <lambda|H0-E0|Psi1> terms
                DBL_MB(idx_v3+jas) = DBL_MB(idx_v3+jas)+derHNS(isp,1)*DBL_MB(idx_v1+ias)*BMAT(iTri(ias,jas))
                DBL_MB(idx_v3+jas) = DBL_MB(idx_v3+jas)+derHNS(isp,2)*DBL_MB(idx_v1+ias)*SMAT(iTri(ias,jas))
                DBL_MB(idx_v3+ias) = DBL_MB(idx_v3+ias)+derHNS(isp,1)*DBL_MB(idx_v1+jas)*BMAT(iTri(ias,jas))
                DBL_MB(idx_v3+ias) = DBL_MB(idx_v3+ias)+derHNS(isp,2)*DBL_MB(idx_v1+jas)*SMAT(iTri(ias,jas))
              end do
              DBL_MB(idx_v3+ias) = DBL_MB(idx_v3+ias)+derHNS(isp,1)*DBL_MB(idx_v1+ias)*BMAT(nTri_Elem(ias))
              DBL_MB(idx_v3+ias) = DBL_MB(idx_v3+ias)+derHNS(isp,2)*DBL_MB(idx_v1+ias)*SMAT(nTri_Elem(ias))
              ! <lambda|H|Psi0> term(s)
              DBL_MB(idx_v3+ias) = DBL_MB(idx_v3+ias)+DBL_MB(idx_v2+ias)
            end do
            !do ias=1,NAS
            !  do jas=1,NAS
            !    ! <lambda|H0-E0|Psi1> terms
            !    DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1)) = DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1))+ &
            !                                          derHNS(isp,1)*DBL_MB(MV1+jas-ILO1+LDV*(isp-JLO1))*BMAT(iTri(ias,jas))
            !    DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1)) = DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1))+ &
            !                                          derHNS(isp,2)*DBL_MB(MV1+jas-ILO1+LDV*(isp-JLO1))*SMAT(iTri(ias,jas))
            !  end do
            !  ! <lambda|H|Psi0> term(s)
            !  DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1)) = DBL_MB(MV3+ias-ILO1+LDV*(isp-JLO1))+ &
            !                                        derHNS(isp,3)*DBL_MB(MV2+ias-ILO1+LDV*(isp-JLO1))
            !end do
          end do
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

end module SC_NEVPT2

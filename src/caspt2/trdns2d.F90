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
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine TRDNS2D(IVEC,JVEC,DPT2,NDPT2,SCAL)
! Add to the diagonal blocks of transition density matrix,
!    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>,
! i.e. inactive/inactive, active/active, and virt/virt
! submatrices.IVEC, JVEC stands for the 1st-order perturbed
! CASPT2 wave functions in vectors nr IVEC, JVEC on LUSOLV.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use EQSOLV, only: iDBMat
use fake_GA, only: GA_Arrays
use caspt2_global, only: do_grad, imag_shift, LISTS, LUSBT, sigma_p_epsilon
use caspt2_module, only: MXCASE, nASup, nInDep, nISup, nSym
use SC_NEVPT2, only: SC_NEVPT2_Amplitude, SC_Prop
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, NDPT2
real(kind=wp), intent(inout) :: DPT2(NDPT2)
real(kind=wp), intent(in) :: SCAL
integer(kind=iwp) :: ICASE, ISYM, jD, lg_v1, lg_v2, NAS, NIN, NIS, nVec
integer(kind=iwp), allocatable :: NINDEP_tmp(:,:)
real(kind=wp), allocatable :: BD(:), ID(:)
#ifdef _MOLCAS_MPP_
real(kind=wp), allocatable :: VEC1(:), VEC2(:)
#endif

!! SC-NEVPT2 does not use linearly-dependence basis,
!! so replace NINDEP with NASUP
if (SC_Prop) then
  !! Ensure that the original NASUP is not used
  if (sigma_p_epsilon /= Zero) then
    call warningMessage(2,'sigma regularization cannot be used for SC-NEVPT2 properties.')
    call quit_onUserError()
  end if
  call mma_allocate(NINDEP_tmp,8,MXCASE,Label='NINDEP_tmp')
  NINDEP_tmp(:,:) = NINDEP(:,:)
  NINDEP(:,:) = NASUP(:,:)
end if

! Inact/Inact and Virt/Virt blocks:
do ICASE=1,13
  !if ((icase /= 12) .and. (icase /= 13)) cycle ! H
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NIS = NISUP(ISYM,ICASE)
    NVEC = NIN*NIS
    if (NVEC == 0) cycle
    !! lg_V1: T+lambda
    !! lg_V2: T
    !! IVEC = iVecX
    !! JVEC = iVecR
    call RHS_ALLO(NIN,NIS,lg_V1)
    call RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
    !! lg_V1 is IVECW if SC-NEVPT2
    if (SC_prop) call SC_NEVPT2_Amplitude(NIN,NIS,iCase,iSym,lg_V1)
    if ((IVEC == JVEC) .and. (.not. SC_prop)) then
      lg_V2 = lg_V1
    else
      call RHS_ALLO(NIN,NIS,lg_V2)
      call RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
      if (SC_prop) then
        if (do_grad) then
          !! lg_V2 = IVECW
          call SC_NEVPT2_Amplitude(NIN,NIS,iCase,iSym,lg_V2)
        else
          call RHS_SCAL(NIN,NIS,lg_V2,Zero)
          call RHS_STRANS(NIN,NIS,+One,lg_V1,lg_V2,ICASE,ISYM)
        end if
      else
        if (do_grad) then
          call RHS_SCAL(NIN,NIS,lg_V1,SCAL)
          if (sigma_p_epsilon /= Zero) then
            !! derivative of the numerator
            !! multiply the lambda part (lg_V2) only
            nAS = nASUP(iSym,iCase)
            call mma_allocate(BD,nAS,Label='BD')
            call mma_allocate(ID,nIS,Label='ID')
            jD = iDBMat(iSym,iCase)
            call dDaFile(LUSBT,2,BD,nAS,jD)
            call dDaFile(LUSBT,2,ID,nIS,jD)
            call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,BD,ID)
            call mma_deallocate(BD)
            call mma_deallocate(ID)
          end if
          call RHS_DAXPY(NIN,NIS,One,lg_V2,lg_V1)
          call RHS_READ_SR(lg_V2,ICASE,ISYM,IVEC)
        end if
      end if
    end if

    !SVC: DIADNS can currently not handle pieces of RHS, so pass the
    ! full array in case we are running in parallel
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      if (KING()) then
        ! copy global array to local buffer
        call mma_allocate(VEC1,NVEC,Label='VEC1')
        call GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
        if (IVEC == JVEC) then
          call DIADNS(ISYM,ICASE,VEC1,nVEC,VEC1,nVEC,DPT2,nDPT2,LISTS,size(LISTS))
        else
          call mma_allocate(VEC2,NVEC,Label='VEC2')
          call GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
          call DIADNS(ISYM,ICASE,VEC1,nVEC,VEC2,nVec,DPT2,nDPT2,LISTS,size(LISTS))
          call mma_deallocate(VEC2)
        end if
        call mma_deallocate(VEC1)
      end if
      call GASYNC()
    else
#   endif
      call DIADNS(ISYM,ICASE,GA_Arrays(lg_V1)%A,NVEC,GA_Arrays(lg_V2)%A,NVEC,DPT2,nDPT2,LISTS,size(LISTS))
#   ifdef _MOLCAS_MPP_
    end if
#   endif
    if (do_grad .and. ((.not. SC_prop) .and. ((imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero)))) then
      !! for sigma-p CASPT2, derivative of the denominator
      nAS = nASUP(iSym,iCase)
      call mma_allocate(BD,nAS,Label='BD')
      call mma_allocate(ID,nIS,Label='ID')
      jD = iDBMat(iSym,iCase)
      call dDaFile(LUSBT,2,BD,nAS,jD)
      call dDaFile(LUSBT,2,ID,nIS,jD)

      call RHS_READ_SR(lg_V1,ICASE,ISYM,IVEC)
      call RHS_READ_SR(lg_V2,ICASE,ISYM,JVEC)
      call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,BD,ID)

      DPT2(:) = -DPT2(:)
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        if (KING()) then
          ! copy global array to local buffer
          call mma_allocate(VEC1,NVEC,Label='VEC1')
          call GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
          if (IVEC == JVEC) then
            call DIADNS(ISYM,ICASE,VEC1,nVEC,VEC1,nVec,DPT2,nDPT2,LISTS,size(LISTS))
          else
            call mma_allocate(VEC2,NVEC,Label='VEC2')
            call GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
            call DIADNS(ISYM,ICASE,VEC1,nVec,VEC2,nVec,DPT2,nDPT2,LISTS,size(LISTS))
            call mma_deallocate(VEC2)
          end if
          call mma_deallocate(VEC1)

        end if
        call GASYNC()
      else
#     endif
        call DIADNS(ISYM,ICASE,GA_Arrays(lg_V1)%A,NVEC,GA_Arrays(lg_V2)%A,NVEC,DPT2,nDPT2,LISTS,size(LISTS))
#     ifdef _MOLCAS_MPP_
      end if
#     endif
      DPT2(:) = -DPT2(:)
      call mma_deallocate(BD)
      call mma_deallocate(ID)
    end if

    call RHS_FREE(lg_V1)
    if ((IVEC /= JVEC) .or. SC_prop) call RHS_FREE(lg_V2)
  end do
end do
#ifdef _MOLCAS_MPP_
if (Is_Real_Par() .and. do_grad) call gadgop(DPT2,NDPT2,'+')
#endif

if (SC_Prop) then
  NINDEP(:,:) = NINDEP_tmp(:,:)
  call mma_deallocate(NINDEP_tmp)
end if

end subroutine TRDNS2D

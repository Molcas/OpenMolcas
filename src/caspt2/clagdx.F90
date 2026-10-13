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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CLagDX(iSym,iCase,nIN,nIS,nAS,nState,VECROT,BDERmat,SDERmat)

use EQSOLV, only: IDBMAT, IDTMAT, IRHS, IVECR, IVECW, IVECX
use fake_GA, only: GA_Arrays
use caspt2_global, only: do_lindep, idSCMat, idSDMat, imag_shift, iVecG, iVecL, LUSBT, LUSTD, real_shift, sigma_p_epsilon
use caspt2_module, only: IFMSCOUP, JSTATE, MAXIT
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, iCase, nIN, nIS, nAS, nState
real(kind=wp), intent(in) :: VECROT(nState)
real(kind=wp), intent(inout) :: BDERmat(NAS*NAS), SDERmat(NAS*NAS)
integer(kind=iwp) :: i, idB, idSD, idT, j, jICB, lg_V1, lg_V2
real(kind=wp) :: Dummy(1), EigJ, SCAL, Tmp
logical(kind=iwp) :: do_lambda
real(kind=wp), allocatable :: DIS(:), EIG(:), TRANS(:), WRK1(:), WRK2(:), WRK3(:), XC(:)
#ifdef _MOLCAS_MPP_
real(kind=wp), allocatable :: BUF1(:), BUF2(:)
#endif

! At most two vectors are read in this subroutine:
!   T      (iVecX): solution in the IC basis
!   lambda (iVecR): solution of the lambda equation in the IC basis
!   RHS    (iRHS) : RHS in the IC basis
!   RHS    (iVecL): RHS in the IC basis of the coupled states, with the multistate coupling instead of iRHS
!   RHS    (iVecW): RHS in the contravariant MO basis
!   RHS    (iVecG): RHS of the other states with the multistate coupling and linear dependence in the contravariant MO basis
! In parallel, the master copies the two vectors to local buffers and does everything.
! Work arrays:
!   WRK1 (nIN,nIN): density, and then the derivative of the overlap, in the IC basis
!   WRK2 (nIN,nAS): (T+lambda/2)*RHS with the RHS in the MO basis. It is also the scratch array before that.
!                   WRK2 may also be used as a working array of WRK2(nAS,nAS)
!   WRK3 (nAS,nAS): explicit derivative of the overlap in the MO basis

do_lambda = (real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP

call mma_allocate(EIG,nIN,Label='EIG')
idB = idBMAT(iSym,iCase)
call DDAFILE(LUSBT,2,EIG,nIN,idB)
call mma_allocate(WRK1,nIN**2,Label='WRK1')
call mma_allocate(WRK2,nAS**2,Label='WRK2')
call mma_allocate(WRK3,nAS**2,Label='WRK3')
call mma_allocate(TRANS,nAS*nIN,Label='TRANS')
if (do_lindep .and. (nAS /= nIN) .and. IFMSCOUP) call mma_allocate(XC,nIN*nAS,Label='XC')
WRK1(:) = Zero

SCAL = One
if (IFMSCOUP) SCAL = VECROT(jState)

! Form the density in internally contracted basis
! The G subspace is employed in the following comments as an example.
! i  : inactive
! a,b: secondary
! t,u: active
! o,p: internally contracted configuration (basis)
! WRK1(o,p) = \sum_{iab} T_{o,i}^{ab}*T_{p,i}^{ab}
! WRK1 is the effective density in the IC basis, and will be the B derivative contribution.

! WRK1 = T*T + (T*lambda+lambda*T)/2, and the terms specific to sigma^p and the imaginary shift
call RHS_ALLO(nIN,nIS,lg_V1)
call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
if (do_lambda) then
  call RHS_ALLO(nIN,nIS,lg_V2)
  call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
  if (sigma_p_epsilon /= Zero) call Contract(1,lg_V1,nIN,lg_V2,nIN)
  if (imag_shift /= Zero) then
    ! WRK1 = WRK1 - (T'*lambda'+lambda'*T')/2, with T and lambda scaled in place
    call mma_allocate(DIS,nIS,Label='DIS')
    idB = idBMAT(iSym,iCase)
    call DDAFILE(LUSBT,0,Dummy,nAS,idB)
    call DDAFILE(LUSBT,2,DIS,nIS,idB)
    call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,EIG,DIS)
    call mma_deallocate(DIS)
    call Contract(2,lg_V1,nIN,lg_V2,nIN)
    call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
    call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
  end if
  ! lg_V2 = (T+lambda)/2, so that T*T + (T*lambda+lambda*T)/2 is lg_V2*T symmetrized
  call RHS_DAXPY(nIN,nIS,SCAL,lg_V1,lg_V2)
  call RHS_SCAL(nIN,nIS,lg_V2,Half)
  call Contract(3,lg_V2,nIN,lg_V1,nIN)
else
  call Contract(3,lg_V1,nIN,lg_V1,nIN)
end if

#ifdef _MOLCAS_MPP_
if ((.not. Is_Real_Par()) .or. King()) then
#endif
  idT = idTMAT(iSym,iCase)
  call DDAFILE(LUSBT,2,TRANS,nAS*nIN,idT)

  ! Transform the internally contracted density to active MO basis
  ! BDERmat(t,u) = ST(t,o)*WRK1(o,p)*ST(u,p)
  ! BDERmat is the derivative contribution of the B matrix in the MO basis
  call DGEMM_('N','N',nAS,nIN,nIN,One,TRANS,nAS,WRK1,nIN,Zero,WRK2,nAS)
  call DGEMM_('N','T',nAS,nAS,nIN,One,WRK2,nAS,TRANS,nAS,Zero,BDERmat,nAS)

  ! Implicit derivative of the IC vector. This derivative comes from the derivative of the eigenvalue only.
  ! Other contributions of the derivative of the IC vector is considered later.
  ! -(e_o + e_p)*dS/da
  do jICB=1,nIN
    EigJ = EIG(jICB)
    WRK1(nIN*(jICB-1)+1:nIN*jICB) = -WRK1(nIN*(jICB-1)+1:nIN*jICB)*(EIG(:)+EigJ)*Half
  end do
#ifdef _MOLCAS_MPP_
end if
#endif

! Derivative of the overlap. The shift-specific terms do not contribute to the RHS terms.
! 1) Implicit: WRK1(o,p) = WRK1(o,p) - RHS(o,i,a,b)*(T+lambda/2)_{p,i}^{ab}
! 2) Explicit: WRK2 = (T+lambda/2)*RHS, with the RHS in the MO basis
! E = 2<1|H|0> + <1|H0-E0|1>
if (IFMSCOUP) then
  ! The RHS that couples with T (iVecL) and that with lambda (iRHS) are different.
  ! WRK1 = WRK1 - RHS*T
  call RHS_READ_SR(lg_V2,iCase,iSym,iVecL)
  call Contract(4,lg_V2,nIN,lg_V1,nIN)
  call RHS_FREE(lg_V2)
  ! WRK2 = T*RHS
  call RHS_ALLO(nAS,nIS,lg_V2)
  call RHS_READ(nAS,nIS,lg_V2,iCase,iSym,iVecW)
  call Contract(5,lg_V1,nIN,lg_V2,nAS)
  if (do_lindep .and. (nAS /= nIN) .and. IFMSCOUP) then
    ! XC = T*RHS, with the RHS of the other states
    call RHS_READ(nAS,nIS,lg_V2,iCase,iSym,iVecG)
    call Contract(8,lg_V1,nIN,lg_V2,nAS)
    call RHS_READ(nAS,nIS,lg_V2,iCase,iSym,iVecW)
  end if
  ! WRK2 = WRK2 + lambda*RHS/2
  call RHS_READ_SR(lg_V1,iCase,iSym,iVecR)
  call Contract(6,lg_V1,nIN,lg_V2,nAS)
  call RHS_FREE(lg_V2)
  ! WRK1 = WRK1 - RHS*lambda/2
  call RHS_ALLO(nIN,nIS,lg_V2)
  call RHS_READ_SR(lg_V2,iCase,iSym,iRHS)
  call Contract(7,lg_V2,nIN,lg_V1,nIN)
  call RHS_FREE(lg_V1)
else
  ! lg_V2 = T+lambda/2
  if (do_lambda) then
    call RHS_DAXPY(nIN,nIS,Half*SCAL,lg_V1,lg_V2)
    call RHS_FREE(lg_V1)
  else
    lg_V2 = lg_V1
  end if
  !! WRK2 = (T+lambda/2)*RHS
  call RHS_ALLO(nAS,nIS,lg_V1)
  call RHS_READ(nAS,nIS,lg_V1,iCase,iSym,iVecW)
  call Contract(5,lg_V2,nIN,lg_V1,nAS)
  call RHS_FREE(lg_V1)
  !! WRK1 = WRK1 - RHS*(T+lambda/2)
  call RHS_ALLO(nIN,nIS,lg_V1)
  call RHS_READ_SR(lg_V1,iCase,iSym,iRHS)
  call Contract(4,lg_V1,nIN,lg_V2,nIN)
  call RHS_FREE(lg_V1)
end if
call RHS_FREE(lg_V2)

#ifdef _MOLCAS_MPP_
if ((.not. Is_Real_Par()) .or. King()) then
#endif
  ! Explicit overlap derivative in the MO basis
  WRK3(:) = Zero
  ! Add the contributions from the off-diagonal coupling (i.e., CASPT2-N).
  if (MAXIT /= 0) then
    idSD = idSDMat(iSym,iCase)
    call DDAFILE(LuSTD,2,WRK3,nAS*nAS,idSD)
    ! T*(T+lambda) + (T+lambda)*T is saved, so 1/2
    WRK3(:) = Half*WRK3(:)
  end if
  call DGEMM_('N','N',nAS,nAS,nIN,Two,TRANS,nAS,WRK2,nIN,One,WRK3,nAS) ! WRK3 = GDer

  ! Convert the implicit overlap derivative to the MO basis
  call DGEMM_('N','N',nAS,nIN,nIN,One,TRANS,nAS,WRK1,nIN,Zero,WRK2,nAS)
  call DGEMM_('N','T',nAS,nAS,nIN,One,WRK2,nAS,TRANS,nAS,Zero,SDERmat,nAS)

  ! Add the contributions due to the dependence on the linearly independent space
  if (do_lindep .and. (nAS /= nIN)) then
    WRK2(:) = Zero
    ! Coupling of H0 with the lower cases for CASPT2-N
    if (MAXIT /= 0) then
      idSD = idSCMat(iSym,iCase)
      call DDAFILE(LuSTD,2,WRK2,nAS*nAS,idSD)
      WRK2(:) = Half*WRK2(:)
    end if
    ! RHS of the other states
    if (IFMSCOUP) call DGEMM_('N','N',nAS,nAS,nIN,Two,TRANS,nAS,XC,nIN,One,WRK2,nAS) ! WRK2 = GCDer
    call LinDepLag(BDERmat,SDERmat,WRK3,WRK2,nAS,nIN,iSym,iCase)
  end if

  SDERmat(1:NAS**2) = SDERmat(1:NAS**2)+WRK3(1:NAS**2)
#ifdef _MOLCAS_MPP_
end if
#endif

call mma_deallocate(WRK1)
call mma_deallocate(WRK2)
if (do_lindep .and. (nAS /= nIN) .and. IFMSCOUP) call mma_deallocate(XC)
call mma_deallocate(WRK3)
call mma_deallocate(TRANS)
call mma_deallocate(EIG)

contains

subroutine Contract(iStep,lg_A,nRowA,lg_B,nRowB)

integer(kind=iwp), intent(in) :: iStep, lg_A, nRowA, lg_B, nRowB

! Contract the two vectors on lg_A and lg_B for the step iStep

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (King()) then
    ! copy global array to local buffer
    call mma_allocate(BUF1,nRowA*nIS,Label='BUF1')
    call mma_allocate(BUF2,nRowB*nIS,Label='BUF2')
    call GA_GET(lg_A,1,nRowA,1,nIS,BUF1,nRowA)
    call GA_GET(lg_B,1,nRowB,1,nIS,BUF2,nRowB)
    call Contract_(iStep,nRowA,nRowB,BUF1,BUF2)
    call mma_deallocate(BUF1)
    call mma_deallocate(BUF2)
  end if
  call GASync()
else
#endif
  call Contract_(iStep,nRowA,nRowB,GA_Arrays(lg_A)%A,GA_Arrays(lg_B)%A)
#ifdef _MOLCAS_MPP_
end if
#endif

end subroutine Contract

subroutine Contract_(iStep,nRowA,nRowB,VEC1,VEC2)

integer(kind=iwp), intent(in) :: iStep, nRowA, nRowB
real(kind=wp), intent(in) :: VEC1(nRowA*nIS), VEC2(nRowB*nIS)

select case (iStep)
  case (1)
    ! VEC1 = T, VEC2 = lambda
    call CLagDX_Sigma(iSym,iCase,nIN,nIS,nAS,1,nIS,VEC1,VEC2,EIG,WRK1)
  case (2,3)
    ! 2: VEC1 = scaled T, VEC2 = scaled lambda, WRK1 = WRK1 - (VEC1*VEC2+VEC2*VEC1)/2
    ! 3: VEC1 = (T+lambda)/2 (or T), VEC2 = T,  WRK1 = WRK1 + VEC1*VEC2+VEC2*VEC1
    Tmp = One
    if (iStep == 2) Tmp = -Half
    if ((iStep == 3) .and. (.not. do_lambda)) Tmp = Half*SCAL
    call DGEMM_('N','T',nIN,nIN,nIS,Tmp,VEC1,nIN,VEC2,nIN,Zero,WRK2,nIN)
    do j=1,nIN
      do i=1,nIN
        WRK1(i+nIN*(j-1)) = WRK1(i+nIN*(j-1))+WRK2(i+nIN*(j-1))+WRK2(j+nIN*(i-1))
      end do
    end do
  case (4)
    ! VEC1 = RHS, VEC2 = T (or T+lambda/2)
    call DGEMM_('N','T',nIN,nIN,nIS,-One,VEC1,nIN,VEC2,nIN,One,WRK1,nIN)
  case (5)
    ! VEC1 = T (or T+lambda/2), VEC2 = RHS in the MO basis
    call DGEMM_('N','T',nIN,nAS,nIS,SCAL,VEC1,nIN,VEC2,nAS,Zero,WRK2,nIN)
  case (6)
    ! VEC1 = lambda, VEC2 = RHS in the MO basis
    call DGEMM_('N','T',nIN,nAS,nIS,Half,VEC1,nIN,VEC2,nAS,One,WRK2,nIN)
  case (7)
    ! VEC1 = RHS, VEC2 = lambda
    call DGEMM_('N','T',nIN,nIN,nIS,-Half,VEC1,nIN,VEC2,nIN,One,WRK1,nIN)
  case (8)
    ! VEC1 = T, VEC2 = RHS of the other states in the MO basis
    call DGEMM_('N','T',nIN,nAS,nIS,One,VEC1,nIN,VEC2,nAS,Zero,XC,nIN)
end select

end subroutine Contract_

end subroutine CLagDX

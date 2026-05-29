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

subroutine CLagDX(Mode,iSym,iCase,VEC1,VEC2,VEC3,VEC4,nIN,nIS,nAS,nState,VECROT,VEC5,lg_V2,BDERmat,SDERmat)

use EQSOLV, only: IDBMAT, IDTMAT, IVECR
use caspt2_global, only: do_lindep, idSDMat, imag_shift, LUSBT, LUSTD, real_shift, sigma_p_epsilon
use caspt2_module, only: IFMSCOUP, JSTATE, MAXIT
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mode, iSym, iCase, nIN, nIS, nAS, nState, lg_V2
real(kind=wp), intent(in) :: VEC1(NIN*NIS), VEC3(NIN*NIS), VEC4(NAS*NIS), VECROT(nState), VEC5(NIN*NIS)
real(kind=wp), intent(inout) :: VEC2(NIN*NIS), BDERmat(NAS*NAS), SDERmat(NAS*NAS)
integer(kind=iwp) :: idB, idSD, idT, iICB, jICB
real(kind=wp) :: EigI, EigJ, SCAL, tmp
logical(kind=iwp) :: invar_act
real(kind=wp), allocatable :: EIG(:), TRANS(:), WRK1(:), WRK2(:), WRK3(:)

!! sigma^P may not introduce non-invariance, so the name may be
!! simply confusing. I just do not know how to apply the
!! non-canonical approach for the derivative of the exponential
!! that appears in the denominator of the sigma^P regularization
invar_act = .true.
if (sigma_p_epsilon /= Zero) invar_act = .false.

call mma_allocate(WRK1,nAS**2,Label='WRK1')
call mma_allocate(WRK2,max(nAS**2,nAS*nIS),Label='WRK2')
call mma_allocate(WRK3,nAS**2,Label='WRK3')
call mma_allocate(TRANS,nAS*nIN,Label='TRANS')
call mma_allocate(EIG,nIN,Label='EIG')

idT = idTMAT(iSym,iCase)
call DDAFILE(LUSBT,2,TRANS,nAS*nIN,idT)
idB = idBMAT(iSym,iCase)
call DDAFILE(LUSBT,2,EIG,nIN,idB)

SCAL = One
if (IFMSCOUP) SCAL = VECROT(jState)

!! VEC1: solution in IC basis
!! VEC2: lambda   in IC basis
!! VEC3: RHS      in IC basis
!! VEC4: RHS      in MO basis

!! Form the density in internally contracted basis
!! The G subspace is employed in the following comments
!! as an example.
!! i  : inactive
!! a,b: secondary
!! t,u: active
!! o,p: internally contracted configuration (basis)
!! WRK1(o,p) = \sum_{iab} T_{o,i}^{ab}*T_{p,i}^{ab}
!! WRK1 is the effective density in the IC basis,
!! and will be the B derivative contribution.

if (Mode == 0) then
  !! WRK1 = T*T
  call DGEMM_('N','T',nIN,nIN,nIS,SCAL,VEC1,nIN,VEC1,nIN,Zero,WRK1,nIN)
else
  WRK1(1:nIN*nIN) = Zero
end if

if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) then
  !! WRK1 = T*T + (T*lambda+lambda*T)/2
  !! For sigma-p CASPT2, this if branch computes the pseudo-
  !! density that comes from the numerator of the shift.
  call DGEMM_('N','T',nIN,nIN,nIS,Half,VEC2,nIN,VEC1,nIN,One,WRK1,nIN)
  call DGEMM_('N','T',nIN,nIN,nIS,Half,VEC1,nIN,VEC2,nIN,One,WRK1,nIN)
end if
if ((sigma_p_epsilon /= Zero) .and. (mode == 0)) then
  !! the remaining is the derivative of 2<1|H|0>, so the unscaled
  !! lambda is loaded
# ifdef _MOLCAS_MPP_
  !! --- To Do ---
  !! Parallel sigma^P gradient does not work at the moment. I know
  !! that the next part has to be fixed, but it is not possible to
  !! load the unscaled lambda in the KING environment, so this
  !! subroutine should be completely rewritten or an additional
  !! (relatively large) array has to be allocated. The number of
  !! allocated VEC arrays is already large (5?), so it is better
  !! to rewrite this subroutine...
  if (Is_Real_Par()) then
    if (KING()) call GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
  else
# endif
    call RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
# ifdef _MOLCAS_MPP_
  end if
# endif
end if

if (invar_act) then
  !! Transform the internally contracted density to
  !! active MO basis
  !! WRK3(t,u) = ST(t,o)*WRK1(o,p)*ST(u,p)
  !! WRK3 is the derivative contribution of the B matrix
  !! in the MO basis
  call DGEMM_('N','N',nAS,nIN,nIN,One,TRANS,nAS,WRK1,nIN,Zero,WRK2,nAS)
  call DGEMM_('N','T',nAS,nAS,nIN,One,WRK2,nAS,TRANS,nAS,Zero,WRK3,nAS)
  !write(u6,*) 'B derivative in MO'
  !call sqprt(WRK3,nas)

  !! Implicit derivative of the IC vector. This derivative
  !! comes from the derivative of the eigenvalue only. Other
  !! contributions of the derivative of the IC vector is
  !! considered later.
  !! -(e_o + e_p)*dS/da
  do jICB=1,nIN
    EigJ = EIG(jICB)
    WRK1(nIN*(jICB-1)+1:nIN*jICB) = -WRK1(nIN*(jICB-1)+1:nIN*jICB)*(EIG(:)+EigJ)*Half
  end do
else
  WRK3(1:NIN**2) = WRK1(1:NIN**2)
  WRK1(1:NIN**2) = Zero
end if

!! Derivative of the overlap in the IC basis.
!! WRK1(o,p) = WRK1(o,p) - T_{o,i}^{ab}*RHS(p,i,a,b)
!! This contribution should not be done for the imaginary
!! shift-specific term
!  1) Implicit overlap derivative of the 2<1|H|0> part
if (Mode == 0) then
  !! WRK1 = -RHS*T
  call DGEMM_('N','T',nIN,nIN,nIS,-One,VEC5,nIN,VEC1,nIN,One,WRK1,nIN)
  !! WRK1 = -RHS*(T+lambda/2)
  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) &
    call DGEMM_('N','T',nIN,nIN,nIS,-Half,VEC3,nIN,VEC2,nIN,One,WRK1,nIN)
end if

if (.not. invar_act) then
  !! it seems that only the diagonal elements are correct?
  !! construct the off-diagonal
  do iICB=1,NIN
    EigI = EIG(iICB)
    do jICB=1,iICB-1 !NIN
      EigJ = EIG(jICB)
      tmp = (WRK1(iICB+NIN*(jICB-1))-WRK1(jICB+NIN*(iICB-1)))/(EigI-EigJ)
      WRK3(iICB+NIN*(jICB-1)) = tmp
      WRK3(jICB+NIN*(iICB-1)) = tmp
    end do
  end do
  !! -(e_o + e_p)*dS/da
  do jICB=1,nIN
    EigJ = EIG(jICB)
    WRK1(nIN*(jICB-1)+1:nIN*jICB) = WRK1(nIN*(jICB-1)+1:nIN*jICB)-WRK3(nIN*(jICB-1)+1:nIN*jICB)*(EIG(:)+EigJ)*Half
  end do
  !! IC -> MO (B matrix)
  call DGEMM_('N','N',nAS,nIN,nIN,One,TRANS,nAS,WRK3,nIN,Zero,WRK2,nAS)
  call DGEMM_('N','T',nAS,nAS,nIN,One,WRK2,nAS,TRANS,nAS,Zero,WRK3,nAS)
end if
!! Convert the IC basis to the MO basis
call DGEMM_('N','N',nAS,nIN,nIN,One,TRANS,nAS,WRK1,nIN,Zero,WRK2,nAS)
call DGEMM_('N','T',nAS,nAS,nIN,One,WRK2,nAS,TRANS,nAS,Zero,WRK1,nAS)

!! Add some trivial contributions due to the dependence
!! on the linearly independent space
if (do_lindep .and. (nAS /= nIN)) call LinDepLag(WRK3,WRK1,nAS,nIN,iSym,iCase)

!  2) Explicit overlap derivative of the 2<1|H|0> part
!     Again, not for imaginary shift-specific terms
if (Mode == 0) then
  !! E = 2<1|H|0> + <1|H0-E0|1>
  call DGEMM_('N','N',nAS,nIS,nIN,SCAL,TRANS,nAS,VEC1,nIN,Zero,WRK2,nAS)
  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) &
    call DGEMM_('N','N',nAS,nIS,nIN,Half,TRANS,nAS,VEC2,nIN,One,WRK2,nAS)
  call DGEMM_('N','T',nAS,nAS,nIS,Two,WRK2,nAS,VEC4,nAS,One,WRK1,nAS)
end if

!! Add the contributions from the off-diagonal coupling
!! (i.e., CASPT2-N). Of course, this is not for imaginary shift-
!! specific terms.
if ((MAXIT /= 0) .and. (Mode == 0)) then
  idSD = idSDMat(iSym,iCase)
  call DDAFILE(LuSTD,2,WRK2,nAS*nAS,idSD)
  !! T*(T+lambda) + (T+lambda)*T is saved, so 1/2
  WRK1(1:NAS**2) = WRK1(1:NAS**2)+Half*WRK2(1:NAS**2)
end if

if (mode == 0) then
  BDERmat(1:NAS**2) = BDERmat(1:NAS**2)+WRK3(1:NAS**2)
  SDERmat(1:NAS**2) = SDERmat(1:NAS**2)+WRK1(1:NAS**2)
else
  BDERmat(1:NAS**2) = BDERmat(1:NAS**2)-WRK3(1:NAS**2)
  SDERmat(1:NAS**2) = SDERmat(1:NAS**2)-WRK1(1:NAS**2)
end if
!! Now, convert the above contributions to derivatives of RDM,
!! weighted Fock, etc.
!! WRK3 is the derivative of B in the MO basis
!! WRK1 is the derivative of S in the MO basis
!! See and be consistent with mkbmat and mksmat
!! Notice that F2 and G2 in mkbmat and mksmat are halved
!! (see getdpref).
!! they are moved to CLagDXA, ... CLagDXG

call mma_deallocate(WRK1)
call mma_deallocate(WRK2)
call mma_deallocate(WRK3)
call mma_deallocate(TRANS)
call mma_deallocate(EIG)

end subroutine CLagDX

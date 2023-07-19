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
! Copyright (C) 2005, Christian Ander                                  *
!***********************************************************************

subroutine EU(dq,dg,gi,H,nH)
!                                                                      *
!     Implemented by Christian Ander, 2005, christian@eximius.se       *
!                                                                      *
!***********************************************************************
!                                                                      *
!     Hessian update method; EU ;from Bofill "Remarks on the Updated   *
!     Hessian Matrix Methods" 2003.                                    *
!                                                                      *
!     Vectors are column vectors.                                      *
!                                                                      *
!     Stuff to do:                                                     *
!                  Double definition for gi, check again if gi is more *
!                  correct.                                            *
!                                                                      *
!                  the expression for dE is making the surface weird   *
!***********************************************************************
!                                                                      *
!     gi    :  gradient from previous it.(nH)
!     dg    :  gradient difference       (nH)
!     de    :  energy constant           (real*8)
!     dq    :  Perturbation Geometry     (nH)
!     Eval  :  EigenValues of Hessian    (nH*(nH+1)/2)
!     Evec  :  EigenVector(s) of Hessian (nH,nH)
!     M     :  M-Matrix                  (nH,nH)
!     E     :  Error matrix              (nH,nH)
!     WorkM :  Temporary working matrix  (nH,nH)
!     WorkV :  Temporary working vector  (nH)
!     WorkR :  Temporary working varii   (real*8)
!     H     :  Hessian                   (nH,nH)
!     nH    :  Hessian size              (integer)
!     p,f   :  Multi-used vectors        (nh)
!     mi    :  Used varii                (real*8)

use Index_Functions, only: nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(in) :: dq(nH), dg(nH), gi(nH)
real(kind=wp), intent(inout) :: H(nH,nH)
integer(kind=iwp) :: i, ij, j
real(kind=wp) :: lim, WorkR
real(kind=wp), allocatable :: E(:,:), Eval(:), EVec(:,:), f(:), M(:,:), p(:), u(:), v(:), WorkM(:,:), WorkV(:)
real(kind=wp), external :: ddot_
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ii
real(kind=wp) :: de, mi

! Make a comment in logfile
write(u6,*) 'hello from EU_'
call RecPrt('H matrix',' ',H,nH,nH)
#endif

call mma_allocate(M,nH,nH,Label='M')
call mma_allocate(WorkM,nH,nH,Label='WorkM')
call mma_allocate(E,nH,nH,Label='E')
call mma_allocate(EVec,nH,nH,Label='EVec')
call mma_allocate(u,nH,Label='u')
call mma_allocate(v,nH,Label='v')
call mma_allocate(EVal,nTri_Elem(nH),Label='Eval')
call mma_allocate(p,nH,Label='p')
call mma_allocate(f,nH,Label='f')
call mma_allocate(WorkV,nH,Label='WorkV')

! Calculate Eval and Evec

! Make a triangular form of H that Jacob/NIDiag accepts

ij = 0
do i=1,nH
  do j=1,i
    ij = ij+1
    EVal(ij) = (H(i,j)+H(j,i))*Half
  end do
end do

! Make a unit Matrix

call unitmat(Evec,nH)

! Get the eigenvalues and eigenvectors

call NIDiag_new(Eval,Evec,nH,nH)
#ifdef _DEBUGPRINT_
call RecPrt('Evec matrix',' ',Evec,nH,nH)
#endif

! Calculate mi, diagonal elements of M matrix (eq. 18)

M(:,:) = Zero
#ifdef _DEBUGPRINT_
call RecPrt('gi-vector',' ',gi,1,nH)
call RecPrt('dq-vector',' ',dq,1,nH)
#endif

! p = <t|dq>

call DGEMM_('N','N',1,nH,nH,One,dq,1,Evec,nH,Zero,p,1)

! f = <t|gi>

call DGEMM_('N','N',1,nH,nH,One,gi,1,Evec,nH,Zero,f,1)

#ifdef _DEBUG
call RecPrt('p-vector',' ',p,1,nH)
call RecPrt('f-vector',' ',f,1,nH)
#endif

! The mi calculation loop

do i=1,nH

  ! If p(i) = 0 we ignore this mode
  lim = 5.0e-7_wp
  if (abs(p(i)) > lim) then
#   ifdef _DEBUGPRINT_

    ! The triangular indexation, ii.
    ii = nTri_Elem(i)

    ! Negative sign for the TS-reaction coordinate.
    ! WorkR = (-)f*p

    if (Eval(ii) > Zero) then
      WorkR = f(i)*p(i)
      de = abs(One/(Two*Eval(ii)*p(i)**2)+WorkR)
    else
      WorkR = -f(i)*p(i)
      de = One
    end if

    mi = exp(-(One/Two*abs(Eval(ii))*p(i)**2+WorkR+f(i)**2/(Two*abs(Eval(ii))))/de)*sqrt(abs(Eval(ii))/(Two*Pi*de))
#   endif

    ! Experimental! (worked well)

    !mi = One-exp(-mi)

    ! Store mi value in M
    M(i,i) = One
#   ifdef _DEBUGPRINT_
    write(u6,*) 'mi ',mi,'  Energy ',de,'  Eigenvalue ',Eval(ii)
#   endif

  else
    ! p(i) = 0, so mi is set to zero.
    M(i,i) = One
#   ifdef _DEBUGPRINT_
    ii = nTri_Elem(i)
    write(u6,*) 'mi = p 1.0  Eigenvalue ',Eval(ii)
#   endif
  end if
end do

! Building of M-matrix (equation 16)

! M = WorkM * T^T = T * M * T^T

call DGEMM_('N','N',nH,nH,nH,One,Evec,nH,M,nH,Zero,WorkM,nH)
call DGEMM_('N','T',nH,nH,nH,One,WorkM,nH,Evec,nH,Zero,M,nH)
#ifdef _DEBUGPRINT_
call RecPrt('M-matrix',' ',M,nH,nH)
! correct

! Building of error matrix E (equation 5)

! u = ( dq^T * M * dq )^-1 *  M * dq = ( dq^T * WorkV )^-1 * WorkV
!   = WorkR * WorkV

! WorkV = M * dq

WorkR = DDot_(nH,dq,1,dq,1)
write(u6,*) WorkR,' = <dq|dq>, should be one?'
#endif
call DGEMM_('N','N',nH,1,nH,One,M,nH,dq,nH,Zero,WorkV,nH)

! WorkR = ( dq^T * WorkV )^-1

WorkR = DDot_(nH,dq,1,WorkV,1)
if (WorkR /= Zero) WorkR = One/WorkR

! u = u + WorkR * WorkV

u(:) = WorkR*WorkV(:)

! v = dg - H * dq (equation 3, quasi-Newton condition)

! according to the paper, we want the Next iterationstep, but
! who doesn't? So we take the current and the previously.

v(:) = dg(:)
call DGEMM_('N','N',nH,1,nH,-One,H,nH,dq,nH,One,v,nH)

! E = v*u^T + u*v^T - ( v^T*dq ) u*u^T (equation 5)

! E = v*u^T

call DGEMM_('N','N',nH,nH,1,One,v,nH,u,1,Zero,E,nH)

! WorkR = -v^T*dq

WorkR = -DDot_(nH,v,1,dq,1)

! E = E + WorkR * u*u^T = v*u^T - ( v^T*dq ) * u*u^T

call DGEMM_('N','N',nH,nH,1,WorkR,u,nH,u,1,One,E,nH)

! E = E + u*v^T = v*u^T - ( v^T*dq ) u*u^T + u*v^T

call DGEMM_('N','N',nH,nH,1,One,u,nH,v,1,One,E,nH)
#ifdef _DEBUGPRINT_
call RecPrt('Error matrix',' ',E,nH,nH)
#endif

! The new Hessian (H = H + E, equation 1)

H(:,:) = H(:,:)+E(:,:)
#ifdef _DEBUGPRINT_
call RecPrt('new Hessian',' ',H,nH,nH)
#endif

! Checking the Quasi-Newton condition.

! WorkV = H * dq

call DGEMM_('N','N',nH,1,nH,One,H,nH,dq,nH,Zero,WorkV,nH)

! WorkV = WorkV - dg = 0.0

#ifdef _DEBUGPRINT_
WorkV(:) = WorkV(:)-dg(:)
call RecPrt('Quasi-Newton',' ',WorkV,1,nH)

write(u6,*) 'goodbye from EU_'
#endif

call mma_deallocate(WorkV)
call mma_deallocate(f)
call mma_deallocate(p)
call mma_deallocate(EVal)
call mma_deallocate(v)
call mma_deallocate(u)
call mma_deallocate(EVec)
call mma_deallocate(E)
call mma_deallocate(WorkM)
call mma_deallocate(M)

return

end subroutine EU

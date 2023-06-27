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

implicit real*8(a-h,o-z)
real*8 H(nH,nH), dq(nH), dg(nH), gi(nH)
#include "stdalloc.fh"
real*8, allocatable :: M(:), WorkM(:), E(:), EVec(:), u(:), v(:), Eval(:), p(:), f(:), WorkV(:)

call mma_allocate(M,nH**2,Label='M    ')
call mma_allocate(WorkM,nH**2,Label='WorkM')
call mma_allocate(E,nH**2,Label='E    ')
call mma_allocate(EVec,nH**2,Label='EVec ')
call mma_allocate(u,nH,Label='u    ')
call mma_allocate(v,nH,Label='v    ')
call mma_allocate(EVal,nH*(nH+1)/2,Label='Eval ')
call mma_allocate(p,nH,Label='p    ')
call mma_allocate(f,nH,Label='f    ')
call mma_allocate(WorkV,nH,Label='WorkV')

call EU_(dq,dg,gi,H,nH,M,WorkM,E,EVec,u,v,Eval,p,f,WorkV)

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

subroutine EU_(dq,dg,gi,H,nH,M,WorkM,E,EVec,u,v,Eval,p,f,WorkV)

implicit none
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
!                  the expression for dE is makeing the surface wierd  *
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

#include "real.fh"
integer nH, i, j, ij
real*8 M(nH,nH), WorkM(nH,nH), E(nH,nH), Evec(nH,nH), H(nH,nH)
real*8 dq(nH), u(nH), v(nH), dg(nH), gi(nH), Eval(nH*(nH+1)/2)
real*8 p(nH), f(nH), WorkV(nH)
real*8 WorkR, ddot_, lim
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer ii
real*8 mi, de

! Make a comment in logfile
write(6,*) 'hello from EU_'
call RecPrt('H matrix',' ',H,nH,nH)
#endif

! Calculate Eval and Evec

! Make a triangular form of H that Jacob/NIDiag accepts

do i=1,nH
  do j=1,i
    ij = i*(i-1)/2+j
    EVal(ij) = (H(i,j)+H(j,i))/2.0d0
  end do
end do

! Make a unit Matrix

call FZero(Evec,nH**2)
call dcopy_(nH,[1.0d0],0,Evec,nH+1)

! Get the eigenvalues and eigenvectors

call NIDiag_new(Eval,Evec,nH,nH)
#ifdef _DEBUGPRINT_
call RecPrt('Evec matrix',' ',Evec,nH,nH)
#endif

! Calculate mi, diagonal elements of M matrix (eq. 18)

call FZero(M,nH**2)
#ifdef _DEBUGPRINT_
call RecPrt('gi-vector',' ',gi,1,nH)
call RecPrt('dq-vector',' ',dq,1,nH)
#endif

! p = <t|dq>

call DGEMM_('N','N',1,nH,nH,1.0d0,dq,1,Evec,nH,0.0d0,p,1)

! f = <t|gi>

call DGEMM_('N','N',1,nH,nH,1.0d0,gi,1,Evec,nH,0.0d0,f,1)

#ifdef _DEBUG
call RecPrt('p-vector',' ',p,1,nH)
call RecPrt('f-vector',' ',f,1,nH)
#endif

! The mi calculation loop

do i=1,nH

  ! If p(i) = 0 we ignore this mode
  lim = 5.0d-7
  if (abs(p(i)) > lim) then
#   ifdef _DEBUGPRINT_

    ! The triangular indexation, ii.
    ii = i*(i+1)/2

    ! Negative sign for the TS-reaction coordinate.
    ! WorkR = (-)f*p

    if (Eval(ii) > 0.0d0) then
      WorkR = f(i)*p(i)
      de = abs(1.0d0/(2.0d0*Eval(ii)*p(i)**2)+WorkR)
    else
      WorkR = -f(i)*p(i)
      de = 1.0d0
    end if

    mi = exp(-(1.0d0/2.0d0*abs(Eval(ii))*p(i)**2+WorkR+f(i)**2/(2.0d0*abs(Eval(ii))))/de)*sqrt(abs(Eval(ii))/(2.0d0*Pi*de))
#   endif

    ! Experimental! (worked well)

    !mi = (1.0d0-exp(-mi/1.0d0))*1.0d0

    ! Store mi value in M
    M(i,i) = 1.0d0
#   ifdef _DEBUGPRINT_
    write(6,*) 'mi ',mi,'  Energy ',de,'  Eigenvalue ',Eval(ii)
#   endif

  else
    ! p(i) = 0, so mi is set to zero.
    M(i,i) = 1.0d0
#   ifdef _DEBUGPRINT_
    ii = i*(i+1)/2
    write(6,*) 'mi = p 1.0  Eigenvalue ',Eval(ii)
#   endif
  end if
end do

! Building of M-matrix (equation 16)

! M = WorkM * T^T = T * M * T^T

call DGEMM_('N','N',nH,nH,nH,1.0d0,Evec,nH,M,nH,0.0d0,WorkM,nH)
call DGEMM_('N','T',nH,nH,nH,1.0d0,WorkM,nH,Evec,nH,0.0d0,M,nH)
#ifdef _DEBUGPRINT_
call RecPrt('M-matrix',' ',M,nH,nH)
! correct

! Building of error matrix E (equation 5)

! u = ( dq^T * M * dq )^-1 *  M * dq = ( dq^T * WorkV )^-1 * WorkV
!   = WorkR * WorkV

! WorkV = M * dq

WorkR = DDot_(nH,dq,1,dq,1)
write(6,*) WorkR,' = <dq|dq>, should be one?'
#endif
call DGEMM_('N','N',nH,1,nH,1.0d0,M,nH,dq,nH,0.0d0,WorkV,nH)

! WorkR = ( dq^T * WorkV )^-1

WorkR = DDot_(nH,dq,1,WorkV,1)
if (WorkR /= 0.0d0) then
  WorkR = 1.0d0/WorkR
end if

! u = u + WorkR * WorkV

call FZero(u,nH)
call DaxPy_(nH,WorkR,WorkV,1,u,1)

! v = dg - H * dq (equation 3, quasi-Newton condition)

! according to the paper, we want the Next iterationstep, but
! who doesn't? So we take the current and the previously.

call dcopy_(nH,dg,1,v,1)
call DGEMM_('N','N',nH,1,nH,-1.0d0,H,nH,dq,nH,1.0d0,v,nH)

! E = v*u^T + u*v^T - ( v^T*dq ) u*u^T (equation 5)

! E = v*u^T

call DGEMM_('N','N',nH,nH,1,1.0d0,v,nH,u,1,0.0d0,E,nH)

! WorkR = -v^T*dq

WorkR = -DDot_(nH,v,1,dq,1)

! E = E + WorkR * u*u^T = v*u^T - ( v^T*dq ) * u*u^T

call DGEMM_('N','N',nH,nH,1,WorkR,u,nH,u,1,1.0d0,E,nH)

! E = E + u*v^T = v*u^T - ( v^T*dq ) u*u^T + u*v^T

call DGEMM_('N','N',nH,nH,1,1.0d0,u,nH,v,1,1.0d0,E,nH)
#ifdef _DEBUGPRINT_
call RecPrt('Error matrix',' ',E,nH,nH)
#endif

! The new Hessian (H = H + E, equation 1)

call DaxPy_(nH*nH,1.0d0,E,1,H,1)
#ifdef _DEBUGPRINT_
call RecPrt('new Hessian',' ',H,nH,nH)
#endif

! Checking the Quasi-Newton condition.

! WorkV = H * dq

call DGEMM_('N','N',nH,1,nH,1.0d0,H,nH,dq,nH,0.0d0,WorkV,nH)

! WorkV = WorkV - dg = 0.00

#ifdef _DEBUGPRINT_
call DaxPy_(nH,-1.0d0,dg,1,WorkV,1)
call RecPrt('Quasi-Newton',' ',WorkV,1,nH)

write(6,*) 'goodbye from EU_'
#endif

return

end subroutine EU_

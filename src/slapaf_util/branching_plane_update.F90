!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!  Branching_Plane_Update
!
!> @brief
!>   Apply the Branching plane update method of Maeda et al.
!> @author Roland Lindh
!> @modified_by Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> The branching plane of conical intersection is defined by two vectors,
!> the gradient difference vector (GDV) and the coupling derivative vector
!> (CDV). Equivalently, it can be defined with two orthonormal vectors
!> \f$ x \f$ and \f$ y \f$, such that \f$ x \f$ has the same direction as
!> GDV and \f$ \{x,y\} \f$ spans the same space as GDV and CDV.
!>
!> This routine obtains an approximate \f$ y \f$ vector, by assuming that
!> the branching plane defined by \f$ \{x,y\} \f$ changes the least from
!> iteration to iteration \cite Mae2010-JCTC-6-1538. Thus:
!>
!> \f[ y_k = \beta y_{k-1} - \alpha x_{k-1} \\
!>     \alpha = \frac{y_{k-1}\cdot x_k}{w} \\
!>     \beta  = \frac{x_{k-1}\cdot x_k}{w} \\
!>     w      = \sqrt{(y_{k-1}\cdot x_k)^2+(x_{k-1}\cdot x_k)^2} \f]
!>
!> As an initial guess for the CDV, the average gradient vector (AGV) at
!> the first iteration is used.
!>
!> @param[in]     AGV   Average gradient vector(s)
!> @param[in]     DGV   Difference gradient vector(s)
!> @param[in,out] CDV   Approximate coupling derivative vector, \f$ y \f$
!> @param[in]     n     Size of the vectors
!> @param[in]     nIter Iteration number
!***********************************************************************

subroutine Branching_Plane_Update(AGV,DGV,CDV,n,nIter)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n, nIter
real(kind=wp), intent(in) :: AGV(n,nIter), DGV(n,nIter)
real(kind=wp), intent(inout) :: CDV(n)
#include "print.fh"
integer(kind=iwp) :: iPrint, iRout, iter
real(kind=wp) :: alpha, beta, r, xx, yx, yx_xx
real(kind=wp), allocatable :: x0(:), x1(:)
real(kind=wp), external :: DDot_

iRout = 31
iPrint = nPrint(iRout)

if (iPrint >= 6) then
  write(u6,*) 'Branching plane'
  write(u6,*) 'n,nIter=',n,nIter
  call RecPrt('AGV',' ',AGV,n,nIter)
  call RecPrt('DGV',' ',DGV,n,nIter)
  call RecPrt('CDV (init)',' ',CDV,n,1)
end if
call mma_allocate(x0,n,Label='x0')
call mma_allocate(x1,n,Label='x1')

! Get the directional vector for the first difference gradient vector (DGV).

x0(:) = DGV(:,1)/sqrt(DDot_(n,DGV(:,1),1,DGV(:,1),1))
x1(:) = x0(:)

! Initial guess for the coupling derivative vector (CDV) is the
! average gradient vector (AGV).
! ... or rather its component perpendicular to DGV
! so, remove the projection of CDV along DGV and renormalize

CDV(:) = AGV(:,1)-DDot_(n,AGV(:,1),1,x0,1)*x0(:)
CDV(:) = CDV(:)/sqrt(DDot_(n,CDV,1,CDV,1))
if (iPrint >= 6) call RecPrt('CDV(0)',' ',CDV,n,1)

! Apply the MOM update method to correct the guessed CDV.

! ipx0: xk-1, ~DGV of previous iteration
! ipx1: xk,   ~DGV of current iteration

do iter=2,nIter
  x1(:) = DGV(:,iter)
  r = sqrt(DDot_(n,x1,1,x1,1))
  x1(:) = DGV(:,iter)/r

  xx = DDot_(n,x0,1,x1,1)
  yx = DDot_(n,CDV,1,x1,1)
  yx_xx = sqrt(yx**2+xx**2)

  ! different signs from the paper, should not matter,
  ! but will keep the y vector sign if x does not change
  alpha = -yx/yx_xx
  beta = xx/yx_xx
  CDV(:) = beta*CDV(:)+alpha*x0(:)

  ! remove the projection of CDV along DGV and renormalize

  if (iPrint >= 6) then
    write(u6,*)
    write(u6,*) 'iter=',iter
    write(u6,*) 'r(DGV)=',r
    write(u6,*) 'xx=',xx
    write(u6,*) 'yx=',yx
    write(u6,*) 'alpha,beta=',alpha,beta
  end if

  CDV(:) = CDV(:)-DDot_(n,CDV,1,x1,1)*x1(:)
  r = sqrt(DDot_(n,CDV,1,CDV,1))
  CDV(:) = CDV(:)/r

  if (iPrint >= 6) write(u6,*) 'r(CDV)=',r

  if (iter /= nIter) x1(:) = x0(:)

end do

call mma_deallocate(x1)
call mma_deallocate(x0)
if (iPrint >= 6) call RecPrt('CDV',' ',CDV,n,1)

return

end subroutine Branching_Plane_Update

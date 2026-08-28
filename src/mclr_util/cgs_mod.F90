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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************
!
! (Preconditioned) Conjugate Gradient Squared (CGS) method
!
! Ref: Sonneveld, P. SIAM J. Sci. Stat. Comput. 1989, 10, 36-52.
!      Itoh, S.; Sugihara, M. Trans. Jpn. Soc. Ind. Appl. Math. 2013, 23, 253-286.
! The preconditioning comes from Algorithm 4 in the second ref (Japanese literature).
!
! CGS does not assume the electronic Hessian is symmetric positive definite.
! The electronic Hessian is non-symmetric if PCM or dynamically weighted things are used, so CGS is preferred.
! Note, however, that CGS requires twice more computation per iteration and memory.
! Unfortunately, CGS requires more iteration than the halve with CG.
! The actual calculation is done in cgs_x.F90.

module cgs_mod

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none

! Activate the CGS method
! This option can perhaps be useful for non-PCM runs as well in some difficult cases
logical(kind=iwp) :: CGS = .false.

! Algorithm number written in Itoh, S. et al (above)
! BiCG algorithms are not implemented mainly because the transpose of the Hessian is not available
! PCGS = 4: the standard preconditioning
! PCGS = 5: an additional preconditioning before the loop
! Other algorithms are not implemented. Algorithm 5 may improve the convergence but I do not know.
integer(kind=iwp), parameter :: PCGS = 4

type CGS_type
  !! CI rotations
  integer(kind=iwp) :: ipPvec, ipQvec, ipUvec, ipR0
  !! orbital rotations
  real(kind=wp), allocatable :: Pvec(:), Qvec(:), Uvec(:), R0(:)
end type CGS_type

type(CGS_type) :: CGSvec

contains

!-----------------------------------------------------------------------

subroutine CGS_init(nconf1,nRoots,nDens2)

  use ipPage, only: ipget

  integer(kind=iwp), intent(in) :: nconf1, nRoots, nDens2

  if ((PCGS /= 4) .and. (PCGS /= 5)) then
    write(u6,*) 'The selected CGS algorithm is not implemented'
    call abend()
  end if

  call mma_allocate(CGSvec%Pvec,nDens2+6,Label='Pvec')
  call mma_allocate(CGSvec%Qvec,nDens2+6,Label='Qvec')
  call mma_allocate(CGSvec%Uvec,nDens2+6,Label='Uvec')
  call mma_allocate(CGSvec%R0,nDens2+6,Label='R0')

  CGSvec%ipPvec = ipGet(nconf1*nroots)
  CGSvec%ipQvec = ipGet(nconf1*nroots)
  CGSvec%ipUvec = ipGet(nconf1*nroots)
  CGSvec%ipR0 = ipGet(nconf1*nroots)

end subroutine CGS_init

!-----------------------------------------------------------------------

subroutine CGS_final()

  implicit none

  call mma_deallocate(CGSvec%Pvec)
  call mma_deallocate(CGSvec%Qvec)
  call mma_deallocate(CGSvec%Uvec)
  call mma_deallocate(CGSvec%R0)

end subroutine CGS_final

end module cgs_mod

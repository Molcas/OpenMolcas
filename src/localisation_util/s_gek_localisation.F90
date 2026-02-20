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
! Copyright (C) 2026, Lila Zapp                                        *
!***********************************************************************

subroutine S_GEK_localisation(nOrb2Loc,lastmiter,coordinates,gradients)
use Definitions, only: wp, iwp, u6
use Constants, only: One, Zero
use Localisation_globals, only: Debug

implicit none

integer(kind=iwp),intent(in) :: nOrb2Loc,lastmiter
real(kind=wp),intent(in) :: coordinates(nOrb2Loc,nOrb2Loc,nDiis), gradients(nOrb2Loc,nOrb2Loc,nDiis)

real(kind=wp), allocatable :: e_f(:,:), e_r(:,:)
integer(kind=iwp) :: k, fsdim, subspdim
logical,parameter :: fullspace = .true., sgek_debug=.true.

! full space dimensionality
fsdim = nOrb2Loc*(nOrb2Loc-1)/2

if (sgek_debug) then
    write(u6,*) 'Enter S-GEK Optimizer'
end if

if (fullspace) then
    subspdim = fsdim
!else
!    subspdim = !you get that from GS orthogonalizsation = number of linearly independent vecs
end if

! basis vectors in subspace representation
call mma_allocate(e_f,fsdim,subsdim,Label='e_f')

! basis vectors in subspace representation
call mma_allocate(e_r,subsdim,subsdim,Label='e_r')
e_r(:,:) = Zero
do k=1,subsdim
    e_r(k,k) = One
end do



end subroutine S_GEK_localisation

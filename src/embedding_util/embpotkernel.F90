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
! Copyright (C) Thomas Dresselhaus                                     *
!***********************************************************************

subroutine EmbPotKernel( &
#                       define _CALLING_
#                       include "int_interface.fh"
                       )
!***********************************************************************
!                                                                      *
! Object: kernel routine to calculate integrals over an embedding      *
!         potential on a grid.                                         *
!                                                                      *
! Called from: OneEl                                                   *
!                                                                      *
!     Author: Thomas Dresselhaus                                       *
!                                                                      *
!***********************************************************************

use Embedding_Global, only: embDebug, nEmbGridPoints, embGridCoord, embPotVal, embWeight
use Index_Functions, only: nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
integer(kind=iwp) :: i, j, ia, ib, m, ix, iy, iz, nShellA, nShellB
! dRA, dRB: distance of the current grid point to A/B
real(kind=wp) :: dRA(3), dRB(3), prefactor
real(kind=wp), allocatable :: radA(:), radB(:), sphA(:), sphB(:)
real(kind=wp), external :: gaussRad

#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(nIC)
unused_var(nHer)
unused_var(Array)
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(lOper)
unused_var(iCho)
unused_var(iStabM)
unused_var(PtChrg)
unused_var(nGrid)
unused_var(iAddPot)

!***** Initialization **************************************************

nShellA = (la+1)*(la+2)/2
nShellB = (lb+1)*(lb+2)/2

! Intermediate arrays, radial and angular part of gaussians
call mma_allocate(radA,nAlpha,label='radA')
call mma_allocate(radB,nBeta,label='radB')
call mma_allocate(sphA,nShellA,label='sphA')
call mma_allocate(sphB,nShellB,label='sphB')

!***** Calculation *****************************************************

! Init result var
rFinal(:,:,:,:) = Zero

! Now loop over grid points first
do m=1,nEmbGridPoints
  ! Calculate distances from gaussian centers to current grid point
  do i=1,3
    dRA(i) = embGridCoord(i,m)-A(i)
    dRB(i) = embGridCoord(i,m)-RB(i)
  end do

  prefactor = embPotVal(m)*embWeight(m)

  ! Precalculate radial factors
  do ia=1,nAlpha
    radA(ia) = gaussRad(Alpha(ia),dRA)
  end do

  do ib=1,nBeta
    radB(ib) = gaussRad(Beta(ib),dRB)
  end do

  ! Precalculate spherial factors
  ! shell A
  i = 0
  do ix=la,0,-1
    do iy=la-ix,0,-1
      iz = la-ix-iy
      i = i+1
      sphA(i) = dRA(1)**ix*dRA(2)**iy*dRA(3)**iz
    end do
  end do
  ! shell B
  i = 0
  do ix=lb,0,-1
    do iy=lb-ix,0,-1
      iz = lb-ix-iy
      i = i+1
      sphB(i) = dRB(1)**ix*dRB(2)**iy*dRB(3)**iz
    end do
  end do

  ! Put them together
  do ia=1,nAlpha
    do ib=1,nBeta
      do i=1,nShellA
        do j=1,nShellB
          rFinal(ia,ib,i,j) = rFinal(ia,ib,i,j)+(prefactor*radA(ia)*radB(ib)*sphA(i)*sphB(j))
        end do
      end do
    end do
  end do

  if (embDebug) then
    if (mod(m,587) == 0) then
      write(u6,*) '-------------------------------------------------'
      write(u6,*) 'm=',m
      write(u6,*) 'A=',A(:)
      write(u6,*) 'r=',embGridCoord(:,m)
      write(u6,*) 'dRA=',dRA(1),dRA(2),dRA(3)
      write(u6,*) 'Alpha(1)=',Alpha(1),'; gaussRad(...)=',gaussRad(Alpha(1),dRA)
      write(u6,*) '-------------------------------------------------'
      write(u6,*) 'm=',m,'prefactor=',prefactor
      write(u6,*) 'radA(1)=',radA(1),'; sphA(1)=',sphA(1),', radB(2)=',radB(2),'; sphB(2)=',sphB(2)
      write(u6,*) '-------------------------------------------------'
    end if
  end if

end do ! m, grid point

if (embDebug) then
  write(u6,*) '-------------------------------------------------'
  do ia=1,nAlpha
    do ib=1,nBeta
      do i=1,nShellA
        do j=1,nShellB
          write(u6,*) rFinal(ia,ib,i,j)
        end do
      end do
    end do
  end do
  !write(u6,*) 'rFinal(1,1,1,1)=',rFinal(1,1,1,1),'; rFinal(1,2,1,2)=',rFinal(1,2,1,2)
  write(u6,*) '-------------------------------------------------'
  write(u6,*) '-------------------------------------------------'
end if

!***** Done. Tidy up. **************************************************

call mma_deallocate(radA)
call mma_deallocate(radB)
call mma_deallocate(sphA)
call mma_deallocate(sphB)

return

end subroutine EmbPotKernel

!***********************************************************************
! Returns the radial part of the value of a GTO with given exponent,   *
! centered at the origin.                                              *
!***********************************************************************
function gaussRad(alpha,r)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: gaussRad
real(kind=wp), intent(in) :: alpha, r(3)
integer(kind=iwp) :: i
real(kind=wp) :: rSquare

rSquare = Zero
do i=1,3
  rSquare = rSquare+r(i)**2
end do

! Radial part
gaussRad = exp(-alpha*rSquare)

return

end function gaussRad

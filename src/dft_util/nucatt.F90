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
! Copyright (C) 2000, Roland Lindh                                     *
!***********************************************************************

subroutine NucAtt(mGrid,iSpin)
!***********************************************************************
!      Author:Roland Lindh, Department of Chemical Physics, University *
!             of Lund, SWEDEN. November 2000                           *
!***********************************************************************

use nq_Grid, only: F_xc, Grid, Rho, vRho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mGrid, iSpin
integer(kind=iwp) :: i, iGrid, iOff, mCenter, n, nCenter, nSym
real(kind=wp) :: Attr, d_alpha, d_beta, DTot, Fact, x, y, z
integer(kind=iwp), allocatable :: nStab(:)
real(kind=wp), allocatable :: Eff(:), RA(:,:), ZA(:)

call Get_nAtoms_All(mCenter)
call mma_allocate(RA,3,mCenter,Label='RA')
call Get_Coord_All(RA,mCenter)

call Get_iScalar('Unique atoms',nCenter)
call mma_allocate(nStab,nCenter,Label='nStab')
call Get_iArray('nStab',nStab,nCenter)
call mma_allocate(Eff,nCenter,Label='Eff')
call Get_dArray('Effective Nuclear Charge',Eff,nCenter)

call Get_iScalar('nSym',nSym)

call mma_allocate(ZA,mCenter,Label='ZA')
iOff = 0
do i=1,nCenter
  n = nSym/nStab(i)
  ZA(iOff+1:iOff+n) = Eff(i)
  iOff = iOff+n
end do

call mma_deallocate(Eff)
call mma_deallocate(nStab)

!                                                                      *
!***********************************************************************
!                                                                      *

vRho(:,:) = Zero
if (iSpin == 1) then
  ! iSpin=1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    DTot = Two*d_alpha

    ! Accumulate contributions to the nuclear attraction energy

    Attr = Zero
    do i=1,mCenter
      x = Grid(1,iGrid)-RA(1,i)
      y = Grid(2,iGrid)-RA(2,i)
      z = Grid(3,iGrid)-RA(3,i)
      Fact = ZA(i)/sqrt(x**2+y**2+z**2)
      Attr = Attr+Fact
    end do
    F_xc(iGrid) = F_xc(iGrid)-Attr*DTot

    vRho(1,iGrid) = -Attr

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  ! iSpin=/=1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iGrid=1,mGrid

    d_alpha = Rho(1,iGrid)
    d_beta = Rho(2,iGrid)
    DTot = d_alpha+d_beta

    ! Accumulate contributions to the nuclear attraction energy

    Attr = Zero
    do i=1,mCenter
      x = Grid(1,iGrid)-RA(1,i)
      y = Grid(2,iGrid)-RA(2,i)
      z = Grid(3,iGrid)-RA(3,i)
      Fact = ZA(i)/sqrt(x**2+y**2+z**2)
      Attr = Attr+Fact
    end do
    F_xc(iGrid) = F_xc(iGrid)-Attr*DTot

    vRho(1,iGrid) = -Attr
    vRho(2,iGrid) = -Attr

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

call mma_deallocate(ZA)
call mma_deallocate(RA)

return

end subroutine NucAtt

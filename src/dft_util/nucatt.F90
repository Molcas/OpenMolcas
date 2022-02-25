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

use nq_Grid, only: Grid
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: mGrid, iSpin
integer(kind=iwp) :: i, iOff, mCenter, n, nCenter, nSym
integer(kind=iwp), allocatable :: nStab(:)
real(kind=wp), allocatable :: Eff(:), RA(:,:), ZA(:)

call Get_nAtoms_All(mCenter)
call mma_allocate(RA,3,mCenter,Label='RA')
call Get_Coord_All(RA,mCenter)

call mma_allocate(ZA,mCenter,Label='ZA')
call Get_iScalar('Unique atoms',nCenter)

call mma_allocate(nStab,nCenter,Label='nStab')
call Get_iArray('nStab',nStab,nCenter)
call mma_allocate(Eff,nCenter,Label='Eff')
call Get_dArray('Effective Nuclear Charge',Eff,nCenter)
call Get_iScalar('nSym',nSym)

iOff = 1
do i=1,nCenter
  n = nSym/nStab(i)
  call dcopy_(n,[Eff(i)],0,ZA(iOff),1)
  iOff = iOff+n
end do

call mma_deallocate(Eff)
call mma_deallocate(nStab)

call Do_NucAtt_(mGrid,iSpin,Grid,RA,ZA,mCenter)

call mma_deallocate(ZA)
call mma_deallocate(RA)

return

end subroutine NucAtt

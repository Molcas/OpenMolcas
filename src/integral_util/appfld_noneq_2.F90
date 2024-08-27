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

subroutine AppFld_NonEq_2(Cavxyz,radius,Eps,lmax,EpsInf,NonEq)

use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer lmax
real*8 Cavxyz((lMax+1)*(lMax+2)*(lMax+3)/6)
real*8, allocatable :: CavSph(:,:)
real*8 radius, EPS, EPSInf
logical NonEq

call mma_allocate(CavSph,lmax+1,lmax+1,Label='CavSph')
call AppFld_2(Cavxyz,CavSph,radius,Eps,lmax,EpsInf,NonEq)
call mma_deallocate(CavSph)

return

end subroutine AppFld_NonEq_2

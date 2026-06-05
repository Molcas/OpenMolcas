
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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

Subroutine sg_d1mat(SGS,CIS,EXS,CI,lCI,ISYCI,D1MAT,lD1MAT)

use stdalloc, only: mma_allocate, mma_deallocate
use sguga, only: CIStruct, EXStruct, SGStruct
use constants, only: Zero, One
use definitions, only: iwp, wp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: lCI, ISYCI, lD1MAT
real(kind=wp), intent(in) :: CI(lCI)
real(kind=wp), intent(out) :: D1MAT(lD1MAT)

real(kind=wp), allocatable:: SGM(:)
real(kind=wp), parameter :: CPQ=One
integer(kind=iwp) iOrb, jOrb, iSym, jSym

Call mma_allocate(SGM,lCI,Label='SGM')

D1MAT(:)=Zero
Do iOrb =1, SGS%nLev
   iSym = SGS%ISM(iOrb)
   Do jOrb =1, iOrb
      jSym = SGS%ISM(jOrb)
      If (jSym/=iSym) Cycle

      SGM(:)=Zero
      Call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,ISYCI,CI,SGM)

      D1MAT(iOrb*(iOrb-1)/2+jOrb)=Dot_Product(CI,SGM)

   End Do
End Do

Call mma_deallocate(SGM)

End Subroutine sg_d1mat


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

subroutine sg_one_tdm(SGS,CIS,EXS,Bra,Ket,lCI,ISYCI,TD1MAT,lTD1MAT)

use Index_Functions, only: iTri
use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: lCI, ISYCI, lTD1MAT
real(kind=wp), intent(in) :: Bra(lCI), Ket(lCI)
real(kind=wp), intent(out) :: TD1MAT(lTD1MAT)
integer(kind=iwp) iOrb, jOrb, iSym, jSym
real(kind=wp), allocatable :: SGM(:)
real(kind=wp), parameter :: CPQ = One

call mma_allocate(SGM,lCI,Label='SGM')

TD1MAT(:) = Zero
do iOrb=1,SGS%nLev
  iSym = SGS%ISM(iOrb)
  do jOrb=1,SGS%nLev
    jSym = SGS%ISM(jOrb)
    if (jSym /= iSym) cycle

    SGM(:) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,ISYCI,Ket,SGM)

    TD1MAT(iOrb+(jOrb-1)*SGS%nLev) = dot_product(Bra,SGM)

  end do
end do

call mma_deallocate(SGM)

end subroutine sg_one_tdm

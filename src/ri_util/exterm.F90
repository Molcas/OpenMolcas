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

module exterm

use Data_Structures, only: Alloc1DiArray_Type, DSBA_Type
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iMP2prpt, iOff_Ymnij(8,5), LuAVector(2), LuBVector(2), nAuxVe, nYmnij(8,5)
type(Alloc1DiArray_Type) :: Ymnij(5)
type(DSBA_Type), target :: CMOi(5), DMLT(5)
real(kind=wp), allocatable :: A(:), AMP2(:,:), BMP2(:,:)
real(kind=wp), allocatable, target :: BklK(:), CijK(:), CilK(:), VJ(:), Yij(:,:,:)

public :: A, AMP2, BklK, BMP2, CijK, CilK, CMOi, DMLT, iMP2prpt, iOff_Ymnij, LuAVector, LuBVector, nAuxVe, nYmnij, VJ, Yij, Ymnij

end module exterm

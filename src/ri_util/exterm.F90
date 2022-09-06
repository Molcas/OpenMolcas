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

use Data_Structures, only: DSBA_type

private

public :: CijK, VJ, CilK, BklK
public :: Ymnij, ipYmnij, nYmnij, iOff_Ymnij
public :: Yij
public :: A, AMP2, BMP2
public :: iMP2prpt, nAuxVe
public :: LuAVector, LuBVector
public :: CMOi, DMLT

real*8, allocatable, target :: CijK(:), VJ(:), CilK(:), BklK(:)
integer, allocatable :: Ymnij(:)
integer ipYmnij(5), nYmnij(8,5), iOff_Ymnij(8,5)
real*8, allocatable, target :: Yij(:,:,:)
real*8, allocatable :: A(:)
! Cholesky Mp2-gradients
real*8, allocatable :: AMP2(:,:)
real*8, allocatable :: BMP2(:,:)
integer :: iMP2prpt, nAuxVe
integer :: LuAVector(2), LuBVector(2)
type(DSBA_Type), target :: CMOi(5), DMLT(5)

end module exterm

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
Module exterm
Private
Public :: CijK, VJ, CilK, BklK
Public :: Ymnij, ipYmnij, nYmnij, iOff_Ymnij
Public :: Yij

Public :: AMP2, BMP2
Public :: iMP2prpt, nAuxVe

Real*8, Allocatable, Target:: CijK(:), VJ(:), CilK(:), BklK(:)
Integer, Allocatable:: Ymnij(:)
Integer ipYmnij(5), nYmnij(8,5), iOff_Ymnij(8,5)
Real*8, Allocatable:: Yij(:,:,:)

!  Cholesky Mp2-gradients
Real*8, Allocatable:: AMP2(:,:)
Real*8, Allocatable:: BMP2(:,:)
Integer :: iMP2prpt, nAuxVe
End Module exterm

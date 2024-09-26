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
! Copyright (C) 2016, Roland Lindh                                     *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     This module contains the most-important globally-allocated arrays*
!     in the SCF program.                                              *
!                                                                      *
!***********************************************************************

module SCF_Arrays
!   Dens    : density matrix - vector containing some (NumDT) last     *
!             (optimized) density matrix differences - (nDT,nD,NumDT)  *
!   TwoHam  : two-el. part of the Fock matrix - vector containing      *
!             corresponding 2-el. contributions - (nDT,nD,NumDT)       *
!   Vxc     : Vxc     part of the Fock matrix - vector containing      *
!             corresponding 2-el. contributions - (nDT,nD,NumDT)       *
!   CMO     : molecular orbitals of length nCMO                        *
!   OccNo   : occupation numbers of length lthO                        *

use Definitions, only: wp

implicit none
private

real(kind=wp), allocatable :: CMO(:,:), CMO_ref(:,:), Darwin(:), EDFT(:), EOrb(:,:), FockAO(:,:), HDiag(:), KntE(:), MssVlc(:), &
                              OccNo(:,:), OneHam(:), Ovrlp(:), TrDD(:,:,:), TrDh(:,:,:), TrDP(:,:,:), TrM(:,:)
real(kind=wp), allocatable, target :: Dens(:,:,:), FockMO(:,:), TwoHam(:,:,:), Vxc(:,:,:)

public :: CMO, CMO_ref, Darwin, Dens, EDFT, EOrb, FockAO, FockMO, HDiag, KntE, MssVlc, OccNo, OneHam, Ovrlp, TrDD, TrDh, TrDP, &
          TrM, TwoHam, Vxc

end module SCF_Arrays

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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************

module Arrays

implicit none
private
public :: Hss, FAMO, FIMO, F0SQMO, FAMO_spinp, FAMO_spinm, SFock, G2mp, G2pp, G2mm, Fm, Fp, G1p, G1m, CMO_Inv, CMO, INT1, pINT1, &
          INT2, pInt2, G2t, G2sq, G1t, KAIN1, KINT2, KINT2A, TI1, TI2

real*8, allocatable :: Hss(:)
real*8, allocatable :: FAMO(:), FIMO(:), FAMO_spinp(:), FAMO_spinm(:), SFock(:)
real*8, allocatable :: Fm(:), Fp(:)
real*8, allocatable :: F0SQMO(:)
! Various one- and two-particle densities
real*8, allocatable :: G1t(:)
real*8, allocatable :: G1p(:), G1m(:)
real*8, allocatable :: G2t(:)
real*8, allocatable :: G2sq(:)
real*8, allocatable :: G2mp(:), G2pp(:), G2mm(:)
! MO coefficients
real*8, allocatable, target :: CMO(:)
real*8, allocatable :: CMO_Inv(:)

!INT1  : 1-electron integrals
!INT2  : 2-electron integrals
!PINT1 : Offsets to symmetry blocks
!PINT2 : Offsets to symmetry blocks

real*8, target, allocatable :: INT1(:)
integer, allocatable :: pINT1(:)
real*8, allocatable :: INT2(:)
integer, allocatable :: pINT2(:)

real*8, pointer :: KAIN1(:), KINT2(:), KINT2A(:)
real*8, target, allocatable :: TI1(:), TI2(:)

end module Arrays

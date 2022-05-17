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
Module Arrays

Implicit None
Private
Public :: Hss, FAMO, FIMO, F0SQMO, FAMO_spinp, FAMO_spinm, SFock, G2mp, G2pp, G2mm, Fm, Fp, &
          G1p, G1m, CMO_Inv, CMO, INT1, pINT1, INT2, pInt2, &
          G2t, G2sq, G1t,  &
          KAIN1, KINT2, KINT2A, TI1, TI2

Real*8, Allocatable:: Hss(:)
Real*8, Allocatable:: FAMO(:), FIMO(:), FAMO_spinp(:), FAMO_spinm(:), SFock(:)
Real*8, Allocatable:: Fm(:), Fp(:)
Real*8, Allocatable:: F0SQMO(:)
!     Various one- and two-particle densities
Real*8, Allocatable:: G1t(:)
Real*8, Allocatable:: G1p(:), G1m(:)
Real*8, Allocatable:: G2t(:)
Real*8, Allocatable:: G2sq(:)
Real*8, Allocatable:: G2mp(:), G2pp(:), G2mm(:)
!     MO coefficients
Real*8, Allocatable:: CMO(:)
Real*8, Allocatable:: CMO_Inv(:)

!         INT1        :  1-electron integrals
!         INT2        :  2-electron integrals
!         PINT1       :  Offsets to symmetry blocks
!         PINT2       :  Offsets to symmetry blocks

Real*8,  Target, Allocatable::  INT1(:)
Integer, Allocatable:: pINT1(:)
Real*8,  Allocatable::  INT2(:)
Integer, Allocatable:: pINT2(:)

Real*8, Pointer:: KAIN1(:), KINT2(:), KINT2A(:)
Real*8, Target, Allocatable:: TI1(:), TI2(:)
End Module Arrays

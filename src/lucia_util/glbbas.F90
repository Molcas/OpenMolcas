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

module GLBBAS

use lucia_data, only: MXPOBS

private

real*8, allocatable :: INT1(:), INT1O(:)
integer, allocatable :: PINT1(:), PINT2(:)
real*8, allocatable :: VEC3(:)

! DETERMINE BASE ADDRESSES
! DFTP        : OPEN SHELL DETERMINANTS OF PROTO TYPE
! CFTP        : BRANCHING DIAGRAMS FOR PROTO TYPES
! DTOC        : CSF-DET TRANSFORMATION FOR PROTO TYPES
! CONF_OCC(I) : SPACE FOR STORING NCNSM CONFIGURATION EXPANSIONS
integer, allocatable :: DFTP(:)
integer, allocatable :: CFTP(:)
real*8, allocatable :: DTOC(:)
type iArray
  integer, allocatable :: I(:)
end type iArray
type(iArray), target :: SDREO_I(8)
integer, pointer :: SDREO(:)
type(iArray) :: CONF_OCC(8), CONF_REO(8), PGINT1A(MXPOBS), PGINT1(MXPOBS)
type(iArray), allocatable :: Z_PTDT(:), REO_PTDT(:)
integer, allocatable :: LSM1(:), LSM2(:)
real*8, allocatable :: RHO1(:), SRHO1(:)
integer, allocatable :: KINH1_NOCCSYM(:), KINH1(:)
real*8, allocatable, target :: CI_VEC(:)
real*8, allocatable :: SIGMA_VEC(:)

public INT1, PINT1, PINT2, LSM1, LSM2, RHO1, VEC3, KINH1, PGINT1, INT1O, PGINT1A, SRHO1, KINH1_NOCCSYM, CONF_OCC, CONF_REO, DFTP, &
       CFTP, DTOC, SDREO_I, Z_PTDT, REO_PTDT, CI_VEC, SIGMA_VEC, SDREO

end module GLBBAS

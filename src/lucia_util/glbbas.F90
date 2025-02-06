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

! DETERMINE BASE ADDRESSES
! DFTP        : OPEN SHELL DETERMINANTS OF PROTO TYPE
! CFTP        : BRANCHING DIAGRAMS FOR PROTO TYPES
! DTOC        : CSF-DET TRANSFORMATION FOR PROTO TYPES
! CONF_OCC(I) : SPACE FOR STORING NCNSM CONFIGURATION EXPANSIONS

use Data_Structures, only: Alloc1DiArray_Type
use lucia_data, only: MXPOBS
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), allocatable :: CFTP(:), DFTP(:), KINH1(:), KINH1_NOCCSYM(:), LSM1(:), LSM2(:), PINT1(:), PINT2(:)
integer(kind=iwp), pointer :: SDREO(:)
real(kind=wp), allocatable :: DTOC(:), INT1(:), INT1O(:), RHO1(:), SIGMA_VEC(:), SRHO1(:), VEC3(:)
real(kind=wp), allocatable, target :: CI_VEC(:)
type(Alloc1DiArray_Type) :: CONF_OCC(8), CONF_REO(8), PGINT1(MXPOBS), PGINT1A(MXPOBS)
type(Alloc1DiArray_Type), target :: SDREO_I(8)
type(Alloc1DiArray_Type), allocatable :: REO_PTDT(:), Z_PTDT(:)

public :: CFTP, CI_VEC, CONF_OCC, CONF_REO, DFTP, DTOC, INT1, INT1O, KINH1, KINH1_NOCCSYM, LSM1, LSM2, PGINT1, PGINT1A, PINT1, &
          PINT2, REO_PTDT, RHO1, SDREO, SDREO_I, SIGMA_VEC, SRHO1, VEC3, Z_PTDT

end module GLBBAS

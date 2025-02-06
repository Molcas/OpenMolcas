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

module STRBAS

use Data_Structures, only: Alloc1DiArray_Type
use lucia_data, only: MXPSTT
use Definitions, only: iwp

implicit none
private

type(Alloc1DiArray_Type) :: ISTSO(MXPSTT), NSTSO(MXPSTT), OCSTR(MXPSTT), STREO(MXPSTT), STSTM(MXPSTT,2), Zmat(MXPSTT)
integer(kind=iwp), allocatable :: IOCLS(:), ISTSGP(:), NSTSGP(:), SPGPAN(:), SPGPCR(:)

public :: IOCLS, ISTSGP, ISTSO, NSTSGP, NSTSO, OCSTR, SPGPAN, SPGPCR, STREO, STSTM, Zmat

end module STRBAS

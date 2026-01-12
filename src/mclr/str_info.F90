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
!               2025, Ignacio Fdez. galvan                             *
!***********************************************************************

module Str_Info

use Definitions, only: wp, iwp

implicit none
private

! OCSTR  : Offsets for occupation of strings
! STREO  : reordering array
! STSM   : Symmetry of each string
! STCL   : Class of each string
! NSTSO  : Number of strings per symmetry and occupation
! ISTSO  : Offset of strings per symmetry and occupation
! EL1    : Number of electrons in RAS1 per sub type
! EL3    : Number of electrons in RAS3 per sub type
! Z      : Lexical addressing of arrays
! EL123  : -"-    But array
! STSTMI : Explicit offsets and lengths
! STSTMN :           "
! STSTM  : ?

type String_Info
  integer(kind=iwp), pointer :: OCSTR(:) => null(), STREO(:) => null(), STSM(:) => null(), STCL(:) => null(), NSTSO(:) => null(), &
                                ISTSO(:) => null(), EL1(:) => null(), EL3(:) => null(), Z(:) => null(), EL123(:) => null(), &
                                STSTMI(:) => null(), STSTMN(:) => null(), STSTM(:,:) => null()
end type String_Info

type String_Hidden
  integer(kind=iwp), allocatable :: OCSTR(:), STREO(:), STSM(:), STCL(:), NSTSO(:), ISTSO(:), EL1(:), EL3(:), Z(:), EL123(:), &
                                    STSTMI(:), STSTMN(:), STSTM(:,:)
end type String_Hidden

! ICONF : NCNSM  CONFIGURATION EXPANSIONS
! ICTS  : address of determinant I in STRING ordering for
!         determinant I in CSF ordering
!         reference symmetry IREFSM.
type Storage
  integer(kind=iwp), allocatable :: ICONF(:), ICTS(:)
end type Storage

integer(kind=iwp), parameter :: MXCNSM = 8, NSTTYP_MAX = 6+1  ! "+1" is the dummy layer

! INITIALIZED IN STRTYP
!  NSTTYP  : Number of string types
!  MNRS1   : Min ras1
!  MXRS1   : Max ras1
!  MNRS3   : Min ras3
!  MXRS3   : Max ras3
!  NELEC   : Number of electrons
!  IAZTP   : Pointer to alpha types
!  IBZTP   : Pointer to beta types
!  iuniqmp : Unique types (not necessary here just 0order space)
!  ISTAC   : Stringtype maping; a(or a+) i -> istac(j,1(2))
!
!  NOCTYP  : Number of occupation classes for given type
!  NSTFTP  : Number of strings of this type
!  DFTP    : OPEN SHELL DETERMINANTS OF PROTO TYPE
!  CFTP    : BRANCHING DIAGRAMS FOR PROTO TYPES
!  DTOC    : CSF-DET TRANSFORMATION FOR PROTO TYPES

integer(kind=iwp) :: IATPM1, IATPM2, IAZTP, IBTPM1, IBTPM2, IBZTP, ISTAC(NSTTYP_MAX,2), ITYP_DUMMY = 0, iuniqmp(NSTTYP_MAX) = 0, &
                     iuniqtp(NSTTYP_MAX) = 0, MNRS1(NSTTYP_MAX) = 0, MNRS3(NSTTYP_MAX) = 0, MXRS1(NSTTYP_MAX) = 0, &
                     MXRS3(NSTTYP_MAX) = 0, NELEC(NSTTYP_MAX) = 0, NOCTYP(NSTTYP_MAX) = 0, NSTFTP(NSTTYP_MAX) = 0, NSTTYP
type(Storage) :: CNSM(MXCNSM)
integer(kind=iwp), allocatable :: CFTP(:), DFTP(:)
real(kind=wp), allocatable :: DTOC(:)
type(String_Info) :: Str(NSTTYP_MAX)
type(String_Hidden), target :: Str_Hidden(NSTTYP_MAX)

public :: CFTP, CNSM, DFTP, DTOC, IATPM1, IATPM2, IAZTP, IBTPM1, IBTPM2, IBZTP, ISTAC, ITYP_DUMMY, iuniqmp, iuniqtp, MNRS1, MNRS3, &
          MXRS1, MXRS3, NELEC, NOCTYP, NSTFTP, NSTTYP, Str, Str_Hidden

end module Str_Info

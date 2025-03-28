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

module Str_Info
!  OCSTR  : Offsets for occupation of strings
!  STREO  : reordering array
!  STSM   : Symmetry of each string
!  STCL   : Class of each string
!  NSTSO  : Number of strings per symmetry and occupation
!  ISTSO  : Offset of strings per symmetry and occupation
!  EL1    : Number of electrons in RAS1 per sub type
!  EL3    : Number of electrons in RAS3 per sub type
!  ACTP   : is sub-type active
!  Z      : Lexical addressing of arrays
!  EL123  : -"-    But array
!  STSTMI : Explicit offsets and lengths
!  STSTMN :           "
!  STSTM  : ?
!  NDMAP  : Down mappings of strings containing the same number of electrons
!  NUMAP  :  Up mappings of strings containing the same number of electrons

! Not used
!  COBSM : Symmetry of conjugated orbitals
!  NIFSJ :
!  IFSJ  :
!  IFSJO :
!  STSTX : Symmetry of excitation connecting strings of given symmetry

use MCLR_Data, only: MXCNSM

implicit none
private
public :: String_Info, Str, NSTTYP_MAX, ITYP_DUMMY, NSTTYP, NELEC, MNRS1, MXRS1, MNRS3, MXRS3, IZORR, ISTTP, iuniqmp, iuniqtp, &
          IAZTP, IBZTP, IARTP, IBRTP, NZSTTP, NRSTTP, IATPM1, IATPM2, IBTPM1, IBTPM2, ISTAC, NOCTYP, NSTFTP, INUMAP, INDMAP, &
          MXNSTR, DFTP, CFTP, DTOC, Storage, CNSM

type String_Info
  sequence
  integer, pointer :: OCSTR(:) => null()
  integer, allocatable :: OCSTR_hidden(:)
  integer, pointer :: STREO(:) => null()
  integer, allocatable :: STREO_hidden(:)
  integer, pointer :: STSM(:) => null()
  integer, allocatable :: STSM_hidden(:)
  integer, pointer :: STCL(:) => null()
  integer, allocatable :: STCL_hidden(:)
  integer, pointer :: NSTSO(:) => null()
  integer, allocatable :: NSTSO_hidden(:)
  integer, pointer :: ISTSO(:) => null()
  integer, allocatable :: ISTSO_hidden(:)
  integer, pointer :: EL1(:) => null()
  integer, allocatable :: EL1_hidden(:)
  integer, pointer :: EL3(:) => null()
  integer, allocatable :: EL3_hidden(:)
  integer, pointer :: ACTP(:) => null()
  integer, allocatable :: ACTP_hidden(:)
  integer, pointer :: Z(:) => null()
  integer, allocatable :: Z_hidden(:)
  integer, pointer :: EL123(:) => null()
  integer, allocatable :: EL123_hidden(:)
  integer, allocatable :: STSTMI(:)
  integer, allocatable :: STSTMN(:)
  integer, pointer :: STSTM(:,:) => null()
  integer, allocatable :: STSTM_hidden(:,:)
  integer, allocatable :: NUMAP(:)
  integer, allocatable :: NDMAP(:)
end type String_Info

type(String_Info), allocatable, target :: Str(:)
!integer, allocatable :: COBSM(:)
!integer, allocatable :: NIFSJ(:)
!integer, allocatable :: IFSJ(:)
!integer, allocatable :: IFSJO(:)
!integer, allocatable :: STSTX(:)

! INITITIALIZED IN STRTYP
!  NSTTYP : Number of string types
!  MNRS1  : Min ras1
!  MXRS1  : Max ras1
!  MNRS3  : Min ras3
!  MXRS3  : Max ras3
!  NELEC  : Number of electrons
!  IZORR  : Zero orde y/n
!  IAZTP  : Pointer to alpha types
!  IBZTP  : Pointer to beta types
!
! not in use (just zero order space)
!  IARTP   : Give type nr to certain exc
!  IBRTP   : Give type nr to certain exc
!
!  NZSTTP  : Not in use
!  NRSTTP  : Not in use
!
!  ISTTP   : Space (0=zero order)
!  iuniqmp : Unique types (not necessary here just 0order space)
integer, parameter :: NSTTYP_MAX = 6+1   ! "+1" is the dummy layer
integer :: ITYP_DUMMY = 0
integer :: NSTTYP, NELEC(NSTTYP_MAX) = 0, MNRS1(NSTTYP_MAX) = 0, MXRS1(NSTTYP_MAX) = 0, MNRS3(NSTTYP_MAX) = 0, &
           MXRS3(NSTTYP_MAX) = 0, IZORR(NSTTYP_MAX) = 0, ISTTP(NSTTYP_MAX) = 0, iuniqmp(NSTTYP_MAX) = 0, iuniqtp(NSTTYP_MAX) = 0, &
           IAZTP, IBZTP, IARTP(3,10), IBRTP(3,10), NZSTTP, NRSTTP, IATPM1, IATPM2, IBTPM1, IBTPM2

!  ISTAC  : Stringtype maping; a(or a+) i -> istac(j,1(2))
!  NOCTYP : Number of occupation classes for given type
!  NSTFTP : Number of strings of this type
!  INUMAP : Mapping of string type to next more general type
!  INDMAP : Mapping of string type to next more restricted type
!  MXNSTR : Largest number of strings of given sym and type

integer :: ISTAC(NSTTYP_MAX,2), NOCTYP(NSTTYP_MAX) = 0, NSTFTP(NSTTYP_MAX) = 0, INUMAP(NSTTYP_MAX) = 0, INDMAP(NSTTYP_MAX) = 0, &
           MXNSTR

!  DFTP          : OPEN SHELL DETERMINANTS OF PROTO TYPE
!  CFTP          : BRANCHING DIAGRAMS FOR PROTO TYPES
!  DTOC          : CSF-DET TRANSFORMATION FOR PROTO TYPES
!  CNSM(:)%ICONF : NCNSM  CONFIGURATION EXPANSIONS
!  CNSM(I)%ICTS  : address of determinant I in STRING ordering for
!                  determinant I in CSF ordering
!                  reference symmetry IREFSM.
integer, allocatable :: DFTP(:)
integer, allocatable :: CFTP(:)
real*8, allocatable :: DTOC(:)
type Storage
  integer, allocatable :: ICONF(:)
  integer, allocatable :: ICTS(:)
end type Storage
type(Storage) :: CNSM(MXCNSM)

end module Str_Info

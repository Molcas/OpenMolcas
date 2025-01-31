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
! Copyright (C) 1994, Jeppe Olsen                                      *
!***********************************************************************

subroutine MEMSTR_GAS()
! Construct pointers for saving information about strings and
! their mappings
!
! GAS version
!
!========
! Input :
!========
! Number and groups of strings defined by /GASSTR/
! Symmetry information stored in         /CSM/
! String information stored in           /STINF/
!=========
! Output
!=========
! Pointers stored in Module STRBAS
!
! Jeppe Olsen, Winter of 1994

use stdalloc, only: mma_allocate
use strbas, only: OCSTR, STREO, NSTSGP, ISTSGP, NSTSO, ISTSO, ZMAT, STSTM, IOCLS, SPGPAN, SPGPCR
use lucia_data, only: NMXOCCLS, NGAS
use lucia_data, only: NGRP, NSTTP, NTSPGP, IGSFGP, NELFGP, NSPGPFTP, NSTFGP
use lucia_data, only: NACOB, NOBPT
use lucia_data, only: ISTAC
use csm_data, only: NSMST

implicit none
integer IGRP, NSTRIN, LSTRIN, ITP, IEL, IGAS, IORB, ISTRIN, LENGTH

! Start of string information

!  Offsets for occupation and reorder array of strings

do IGRP=1,NGRP
  NSTRIN = NSTFGP(IGRP)
  LSTRIN = NSTRIN*NELFGP(IGRP)
  call mma_allocate(OCSTR(IGRP)%I,LSTRIN,Label='OCSTR()')
  call mma_allocate(STREO(IGRP)%I,NSTRIN,Label='STREO()')
end do

! Number of strings per symmetry and offset for strings of given sym
! for groups

call mma_allocate(NSTSGP(1)%I,NSMST*NGRP,Label='NSTSGP(1)')
call mma_allocate(ISTSGP(1)%I,NSMST*NGRP,Label='ISTSGP(1)')

! Number of strings per symmetry and offset for strings of given sym
! for types

do ITP=1,NSTTP
  call mma_allocate(NSTSO(ITP)%I,NSPGPFTP(ITP)*NSMST,Label='NSTSO(ITP)')
  call mma_allocate(ISTSO(ITP)%I,NSPGPFTP(ITP)*NSMST,Label='ISTSO(ITP)')
end do

! Lexical addressing of arrays : use array indices for complete active space

! Not in use so
do IGRP=1,NGRP
  call mma_allocate(Zmat(IGRP)%I,NACOB*NELFGP(IGRP),Label='ZMat()')
end do

! Mappings between different groups

do IGRP=1,NGRP
  IEL = NELFGP(IGRP)
  IGAS = IGSFGP(IGRP)
  IORB = NOBPT(IGAS)
  ISTRIN = NSTFGP(IGRP)
  ! If creation is involve : Use full orbital notation
  ! If only annihilation is involved, compact form will be used
  LENGTH = 1
  if (ISTAC(IGRP,2) /= 0) then
    LENGTH = IORB*ISTRIN
  else if (ISTAC(IGRP,1) /= 0) then
    ! Only annihilation map so
    LENGTH = IEL*ISTRIN
  end if
  call mma_allocate(STSTM(IGRP,1)%I,LENGTH,LABEL='STSTM(IGRP,1)')
  call mma_allocate(STSTM(IGRP,2)%I,LENGTH,LABEL='STSTM(IGRP,2)')
end do

! Occupation classes

call mma_allocate(IOCLS,NMXOCCLS*NGAS,Label='IOCLS')
! Annihilation/Creation map of supergroup types
call mma_allocate(SPGPAN,NTSPGP*NGAS,Label='SPGPAN')
call mma_allocate(SPGPCR,NTSPGP*NGAS,Label='SPGPCR')

end subroutine MEMSTR_GAS

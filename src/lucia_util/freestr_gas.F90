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

subroutine FREESTR_GAS()

use stdalloc, only: mma_deallocate
use strbas, only: OCSTR, STREO, NSTSGP, ISTSGP, NSTSO, ISTSO, ZMAT, STSTM, IOCLS, SPGPAN, SPGPCR
! allocations during strinf_gas
use distsym, only: ISMDFGP, ISMSCR, NACTSYM
! Deallocate the memory that was set up in MEMSTR_GAS
use lucia_data, only: NGRP, NSTTP

implicit none
integer IGRP, ITP

! Offsets for occupation and reorder array of strings

do IGRP=1,NGRP
  call mma_deallocate(OCSTR(IGRP)%I)
  call mma_deallocate(STREO(IGRP)%I)
end do

! Number of strings per symmetry and offset for strings of given sym
! for groups

call mma_deallocate(NSTSGP(1)%I)
call mma_deallocate(ISTSGP(1)%I)

! Number of strings per symmetry and offset for strings of given sym
! for types

do ITP=1,NSTTP
  call mma_deallocate(NSTSO(ITP)%I)
  call mma_deallocate(ISTSO(ITP)%I)
end do

! Lexical addressing of arrays : use array indices for complete active space

! Not in use so
do IGRP=1,NGRP
  call mma_deallocate(Zmat(IGRP)%I)
end do

! Mappings between different groups

do IGRP=1,NGRP
  ! IF creation is involve : Use full orbital notation
  ! If only annihilation is involved, compact form will be used
  call mma_deallocate(STSTM(IGRP,1)%I)
  call mma_deallocate(STSTM(IGRP,2)%I)
end do

! Occupation classes

call mma_deallocate(IOCLS)
! Annihilation/Creation map of supergroup types
call mma_deallocate(SPGPAN)
call mma_deallocate(SPGPCR)

! Allocated during strinf_gas call
call mma_deallocate(ISMDFGP)
call mma_deallocate(NACTSYM)
call mma_deallocate(ISMSCR)

end subroutine FREESTR_GAS

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

use lucia_data, only: IOCLS, ISMDFGP, ISMSCR, ISTSGP, ISTSO, NACTSYM, NGRP, NSTSGP, NSTSO, NSTTP, OCCSTR, SPGPAN, SPGPCR, STREO, &
                      STSTM, ZMAT
use stdalloc, only: mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IGRP, ITP

! Offsets for occupation and reorder array of strings

do IGRP=1,NGRP
  call mma_deallocate(OCCSTR(IGRP)%A)
  call mma_deallocate(STREO(IGRP)%A)
end do

! Number of strings per symmetry and offset for strings of given sym
! for groups

call mma_deallocate(NSTSGP)
call mma_deallocate(ISTSGP)

! Number of strings per symmetry and offset for strings of given sym
! for types

do ITP=1,NSTTP
  call mma_deallocate(NSTSO(ITP)%A)
  call mma_deallocate(ISTSO(ITP)%A)
end do

! Lexical addressing of arrays : use array indices for complete active space

! Not in use so
do IGRP=1,NGRP
  call mma_deallocate(Zmat(IGRP)%A)
end do

! Mappings between different groups

do IGRP=1,NGRP
  ! IF creation is involve : Use full orbital notation
  ! If only annihilation is involved, compact form will be used
  call mma_deallocate(STSTM(IGRP,1)%A)
  call mma_deallocate(STSTM(IGRP,2)%A)
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

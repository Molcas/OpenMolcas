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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************

subroutine NEWDIA()
! Post-diagonalization modification of diagonal energy
! denominator terms for active and for non-active superindex.

use EQSOLV, only: IDBMAT
use caspt2_global, only: LUSBT
use caspt2_module, only: nASup, nInDep, nISup, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: I, ICASE, ISYM, JD, NAS, NIN, NIS
real(kind=wp), allocatable :: BD(:), C1(:), C2(:), ID(:)

do ICASE=1,13
  do ISYM=1,NSYM
    NIN = NINDEP(ISYM,ICASE)
    if (NIN == 0) cycle
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NIS == 0) cycle
    ! Remember: NIN values in BDIAG, but must read NAS for correct
    ! positioning.
    call mma_allocate(BD,NAS,LABEL='BD')
    call mma_allocate(ID,NIS,LABEL='ID')
    call mma_allocate(C1,NAS,LABEL='C1')
    call mma_allocate(C2,NIS,LABEL='C2')
    JD = IDBMAT(ISYM,ICASE)
    ! Active, and non-active, energy denominators:
    call DDAFILE(LUSBT,2,BD,NAS,JD)
    call DDAFILE(LUSBT,2,ID,NIS,JD)
    ! Active, and non-active, corrections:
    ! (Replace this strange example with something sensible)
    C1(:) = Zero
    C2(:) = Zero
    ! Modifications are added to the usual diagonal energies:
    do I=1,NAS
      BD(I) = BD(I)+C1(I)
    end do
    do I=1,NIS
      ID(I) = ID(I)+C2(I)
    end do
    JD = IDBMAT(ISYM,ICASE)
    call DDAFILE(LUSBT,1,BD,NAS,JD)
    call DDAFILE(LUSBT,1,ID,NIS,JD)
    ! Added modifications are saved on LUSBT.
    call DDAFILE(LUSBT,1,C1,NAS,JD)
    call DDAFILE(LUSBT,1,C2,NIS,JD)

    call mma_deallocate(BD)
    call mma_deallocate(ID)
    call mma_deallocate(C1)
    call mma_deallocate(C2)
    JD = IDBMAT(ISYM,ICASE)
  end do
end do

end subroutine NEWDIA

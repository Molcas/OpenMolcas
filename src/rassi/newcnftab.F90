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

subroutine NEWCNFTAB(NEL,NORB,MINOP,MAXOP,LSYM,NGAS,NGASORB,NGASLIM,IFORM,ICASE)

use Index_Functions, only: nTri_Elem, nTri_Elem1
use rassi_global_arrays, only: CnfTab1, CnfTab2
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NEL, NORB, MINOP, MAXOP, LSYM, NGAS, NGASORB(nIrrep,NGAS), NGASLIM(2,NGAS), IFORM, ICASE
integer(kind=iwp) :: IFPOSS, IGAS, IPOS, ISUM, ISYM, KCNFEND, KCNFSTA, KINFO, L, LENCNF, MXO, NCLS, NCNF, NNCNF1, NNCNF2, NOCC, &
                     NOPN, NTAB
integer(kind=iwp), allocatable :: NCNF1(:), NCNF2(:)
integer(kind=iwp), pointer :: CnfTab(:)

! Note how input parameter LSYM is used: If non-zero, only those configurations
! with symmetry label LSYM are selected. But if LSYM=0, they are all selected.
! We must figure out sizes before allocating the new configuration table.
! Set up a table NCNF1(nIrrep,NPOS) with NPOS=nTri_Elem1(NEL)
NNCNF1 = nIrrep*nTri_Elem1(NEL)
call mma_allocate(NCNF1,NNCNF1,Label='NCNF1')
! We need also a table NCNF2, temporarily. Need to know mx nr of active orbitals
! in any GAS subspace:
MXO = 0
do IGAS=1,NGAS
  ISUM = sum(NGASORB(1:nIrrep,IGAS))
  MXO = max(MXO,ISUM)
end do
NNCNF2 = nIrrep*nTri_Elem1(MXO)
call mma_allocate(NCNF2,NNCNF2,Label='NCNF2')
call NRCNF1(NEL,NORB,NGAS,NGASLIM,NGASORB,NCNF1,MXO,NCNF2)
call mma_deallocate(NCNF2)

! NCNF1(ISYM,IPOS) contains the number of possible configurations for symmetry
! label ISYM, nr of closed shells NCLS, and nr of open shells NOPN. The latter
! are combined as pair index IPOS=1+NOPN+nTri_Elem(NOCC) with NOCC=NCLS+NOPN

! Header (See below for contents):
NTAB = 10
! NGASORB array:
NTAB = NTAB+(nIrrep+1)*(NGAS+1)
! NGASLIM array:
NTAB = NTAB+2*NGAS
! Save offset to INFO table for later use:
KINFO = NTAB+1
! INFO array:
NTAB = NTAB+3*nIrrep*(MAXOP-MINOP+1)
! Save offset to configuration arrays for later use:
KCNFSTA = NTAB+1
! Configuration arrays:
do NOPN=MINOP,min(2*NORB-NEL,NEL,MAXOP)
  NCLS = (NEL-NOPN)/2
  if (NCLS < 0) cycle
  if (2*NCLS+NOPN /= NEL) cycle
  NOCC = NCLS+NOPN
  if (NOCC > NORB) cycle
  do ISYM=1,nIrrep
    NCNF = 0
    if ((LSYM >= 1) .and. (LSYM <= nIrrep)) then
      IPOS = 1+NOPN+nTri_Elem(NOCC)
      NCNF = NCNF1(ISYM+nIrrep*(IPOS-1))
    end if
    LENCNF = NOCC
    if (IFORM == 2) LENCNF = NORB
    if (IFORM == 3) LENCNF = (NOCC+3)/4
    if (IFORM == 4) LENCNF = (NORB+14)/15
    NTAB = NTAB+NCNF*LENCNF
  end do
end do

! Sizes and offsets are known. Now, we can allocate the table:
select case (ICASE)
  case (1)
    call mma_allocate(CnfTab1,NTAB,Label='CnfTab1')
    CnfTab => CnfTab1(:)
  case (2)
    call mma_allocate(CnfTab2,NTAB,Label='CnfTab2')
    CnfTab => CnfTab2(:)
  case default
    call ABEND()
    CnfTab => CnfTab1(:) ! dummy
end select

! Enter header:
CnfTab(1) = NTAB
CnfTab(2) = 37
CnfTab(3) = NEL
CnfTab(4) = NORB
CnfTab(5) = MINOP
CnfTab(6) = MAXOP
CnfTab(7) = nIrrep
CnfTab(8) = LSYM
CnfTab(9) = NGAS
CnfTab(10) = IFORM
! Enter copy of NGASORB array:
do IGAS=1,NGAS
  ISUM = 0
  do ISYM=1,nIrrep
    L = 11+ISYM+(nIrrep+1)*IGAS
    CnfTab(L) = NGASORB(ISYM,IGAS)
    ISUM = ISUM+NGASORB(ISYM,IGAS)
  end do
  CnfTab(11+(nIrrep+1)*IGAS) = ISUM
end do
do ISYM=0,nIrrep
  ISUM = 0
  do IGAS=1,NGAS
    L = 11+ISYM+(nIrrep+1)*IGAS
    ISUM = ISUM+CnfTab(L)
  end do
  CnfTab(11+ISYM) = ISUM
end do
L = 10+(nIrrep+1)*(NGAS+1)
! Enter copy of NGASLIM array:
do IGAS=1,NGAS
  L = L+1
  CnfTab(L) = NGASLIM(1,IGAS)
  L = L+1
  CnfTab(L) = NGASLIM(2,IGAS)
end do
! Construct and enter INFO table.
! The INFO table has a relative pointer to configuration arrays:
KCNFEND = KCNFSTA-1
do NOPN=MINOP,MAXOP
  NCLS = (NEL-NOPN)/2
  IFPOSS = 1
  if (NCLS < 0) IFPOSS = 0
  if (2*NCLS+NOPN /= NEL) IFPOSS = 0
  if (NCLS+NOPN > NORB) IFPOSS = 0
  if (IFPOSS == 0) then
    do ISYM=1,nIrrep
      ! No such configuration is possible.
      !INFO(1,ISYM,NOPN) = NCNF
      CnfTab(KINFO+0+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = 0
      !INFO(2,ISYM,NOPN) = NTAB+1
      CnfTab(KINFO+1+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = -1
      !INFO(3,ISYM,NOPN) = LENCNF
      CnfTab(KINFO+2+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = 0
    end do
  else
    NOCC = NCLS+NOPN
    do ISYM=1,nIrrep
      NCNF = 0
      ! If LSYM=0, all symmetry labels will be accepted. Else, only
      ! those with ISYM=LSYM.
      if ((LSYM == 0) .or. (ISYM == LSYM)) then
        IPOS = 1+NOPN+nTri_Elem(NOCC)
        NCNF = NCNF1(ISYM+nIrrep*(IPOS-1))
      end if
      if (NCNF == 0) then
        ! INFO(1,ISYM,NOPN) = NCNF
        CnfTab(KINFO+0+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = 0
        ! INFO(2,ISYM,NOPN) = NTAB+1
        CnfTab(KINFO+1+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = -1
        ! INFO(3,ISYM,NOPN) = LENCNF
        CnfTab(KINFO+2+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = 0
      else
        LENCNF = NOCC
        if (IFORM == 2) LENCNF = NORB
        if (IFORM == 3) LENCNF = (NOCC+3)/4
        if (IFORM == 4) LENCNF = (NORB+14)/15
        ! The relative pointer into this configuration array:
        KCNFSTA = KCNFEND+1
        KCNFEND = KCNFEND+NCNF*LENCNF
        ! INFO(1,ISYM,NOPN) = NCNF
        CnfTab(KINFO+0+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = NCNF
        ! INFO(2,ISYM,NOPN) = NTAB+1
        CnfTab(KINFO+1+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = KCNFSTA
        ! INFO(3,ISYM,NOPN) = LENCNF
        CnfTab(KINFO+2+3*(ISYM-1+nIrrep*(NOPN-MINOP))) = LENCNF
        CnfTab(KCNFSTA:KCNFSTA+NCNF*LENCNF-1) = 0
      end if
    end do
  end if
end do
! The NCNF1 array is no longer needed.
call mma_deallocate(NCNF1)

! Finally, only now when we know where to store each (ISYM,NOPN) block of
! configurations, can we compute the actual configuration arrays:
call MKCONF(CnfTab)
nullify(CnfTab)

end subroutine NEWCNFTAB

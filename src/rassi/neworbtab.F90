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

subroutine NEWORBTAB(IPRTTAB)

use rassi_global_arrays, only: OrbTab
use Cntrl, only: MORSBITS
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IPRTTAB(*)
integer(kind=iwp) :: I, IEXTNUM, INPART, INSBP, INSYM(8), IPAC, IPART, IPFR, IPIN, IPSE, ISMLAB, ISOIND, ISORB, ISPART, ISUM, &
                     ISYM, KOINFO, KSPART, N, NAPART, NASPO, NASPRT, NOES(8), NORBT, NPART, NSORBT, NSPART, NSYM, NTAB

ISOIND = 0 ! dummy initialize
! Executable statements
NPART = IPRTTAB(3)
NSYM = IPRTTAB(4)
! NAPART: Nr of partitions of active orbitals.
NAPART = NPART-4
! Partitions for inactive, secondary,frozen, and deleted orbitals:
IPIN = NAPART+1
IPSE = NAPART+2
IPFR = NAPART+3
! Total nr of orbitals:
NORBT = IPRTTAB(5)
NSORBT = 2*NORBT
! Table words 1--10 contain some header info.
! Table words 11--18 contain start index of each CMO symmetry block
! Table words 19-- contain info for each separate spin orbital
! Presently 8 table entries for each spin orbital.
KOINFO = 19
! Table words 19+8*NSORBT-- contain nr of sp-orbs for each subpartition
KSPART = 19+8*NSORBT
! Nr of sub-partitions:
NSPART = 0
do IPART=1,NPART
  N = 2*IPRTTAB(5+(NSYM+1)*IPART)
  if (N > 0) NSPART = NSPART+(N+MORSBITS-1)/MORSBITS
end do
NTAB = KSPART+NSPART-1
! Allocate the orbital table.
call mma_allocate(ORBTAB,NTAB,Label='OrbTab')
OrbTab(1) = NTAB
OrbTab(2) = 1
OrbTab(3) = NSORBT
! Nr of active spin-orbitals:
ISUM = 0
do IPART=1,NAPART
  ISUM = ISUM+IPRTTAB(5+(NSYM+1)*IPART)
end do
NASPO = 2*ISUM
OrbTab(4) = NASPO
OrbTab(5) = NSYM
OrbTab(6) = NAPART+4
! Accumulated nr of orbital functions/symm:
ISUM = 0
do ISYM=1,NSYM
  NOES(ISYM) = ISUM
  ISUM = ISUM+IPRTTAB(5+ISYM)
end do
! Relative pointers to CMO symmetry blocks:
ISUM = 0
do ISYM=1,NSYM
  OrbTab(9+ISYM) = ISUM+1
  ISUM = ISUM+IPRTTAB(5+ISYM)**2
end do
! Spin orbital number:
ISORB = 0
! First, active orbitals by partition, and by symmetry
! Previous MO indices within each symmetry.
INSYM(1:NSYM) = IPRTTAB(6+(NSYM+1)*IPFR:5+NSYM+(NSYM+1)*IPFR)+IPRTTAB(6+(NSYM+1)*IPIN:5+NSYM+(NSYM+1)*IPIN)
! Increase subpartition index as needed.
ISPART = 0
do IPART=1,NAPART
  N = IPRTTAB(5+(NSYM+1)*IPART)
  if (N == 0) cycle
  ISPART = ISPART+1
  INSBP = 0
  INPART = 0
  do ISYM=1,NSYM
    ISMLAB = ISYM
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPART)
    do I=1,N
      ISOIND = 1+INSYM(ISYM)
      INSYM(ISYM) = ISOIND
      IEXTNUM = NOES(ISYM)+ISOIND
      ! Next spin orbital, with alpha spin:
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      ! Fill in table:
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = 1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
      ! Next spin orbital, same, but beta spin:
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      ! Fill in table:
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = -1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
    end do
  end do
end do
NASPRT = ISPART
! Inactive:
! Must set up start index within each symmetry.
do ISYM=1,NSYM
  INSYM(1:NSYM) = IPRTTAB(6+(NSYM+1)*IPFR:5+NSYM+(NSYM+1)*IPFR)
end do
IPART = NAPART+1
N = IPRTTAB(5+(NSYM+1)*IPART)
if (N /= 0) then
  ISPART = ISPART+1
  INPART = 0
  INSBP = 0
  do ISYM=1,NSYM
    ISMLAB = ISYM
    !N = NISH(ISYM)
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPART)
    do I=1,N
      ISOIND = ISOIND+1
      IEXTNUM = NOES(ISYM)+ISOIND
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = 1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = -1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
    end do
  end do
end if
! Secondary:
IPART = NAPART+2
N = IPRTTAB(5+(NSYM+1)*IPART)
if (N /= 0) then
  ISPART = ISPART+1
  INPART = 0
  INSBP = 0
  ! Must set up start index within each symmetry.
  do ISYM=1,NSYM
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPFR)
    N = N+IPRTTAB(5+ISYM+(NSYM+1)*IPIN)
    do IPAC=1,NAPART
      N = N+IPRTTAB(5+ISYM+(NSYM+1)*IPAC)
    end do
    INSYM(ISYM) = N
  end do

  do ISYM=1,NSYM
    ISMLAB = ISYM
    ISOIND = INSYM(ISYM)
    !N = NSSH(ISYM)
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPART)
    do I=1,N
      ISOIND = ISOIND+1
      IEXTNUM = NOES(ISYM)+ISOIND
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = 1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = -1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
    end do
  end do
end if
! Frozen:
IPART = NAPART+3
N = IPRTTAB(5+(NSYM+1)*IPART)
if (N /= 0) then
  ISPART = ISPART+1
  INPART = 0
  INSBP = 0
  do ISYM=1,NSYM
    ISMLAB = ISYM
    !N = NFRO(ISYM)
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPART)
    ISOIND = 0
    do I=1,N
      ISOIND = ISOIND+1
      IEXTNUM = NOES(ISYM)+ISOIND
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = 1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = -1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
    end do
  end do
end if
! Deleted:
IPART = NAPART+4
N = IPRTTAB(5+(NSYM+1)*IPART)
if (N /= 0) then
  ISPART = ISPART+1
  INPART = 0
  INSBP = 0
  ! Must set up start index within each symmetry.
  do ISYM=1,NSYM
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPFR)
    N = N+IPRTTAB(5+ISYM+(NSYM+1)*IPIN)
    do IPAC=1,NAPART
      N = N+IPRTTAB(5+ISYM+(NSYM+1)*IPAC)
    end do
    N = N+IPRTTAB(5+ISYM+(NSYM+1)*IPSE)
    INSYM(ISYM) = N
  end do
  do ISYM=1,NSYM
    ISMLAB = ISYM
    ISOIND = INSYM(ISYM)
    !N = NDEL(ISYM)
    N = IPRTTAB(5+ISYM+(NSYM+1)*IPART)
    do I=1,N
      ISOIND = ISOIND+1
      IEXTNUM = NOES(ISYM)+ISOIND
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = 1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
      ISORB = ISORB+1
      INPART = INPART+1
      INSBP = INSBP+1
      if (INSBP > MORSBITS) then
        ISPART = ISPART+1
        INSBP = INSBP-MORSBITS
      end if
      OrbTab(KOINFO+0+(ISORB-1)*8) = IEXTNUM
      OrbTab(KOINFO+1+(ISORB-1)*8) = ISMLAB
      OrbTab(KOINFO+2+(ISORB-1)*8) = ISOIND
      OrbTab(KOINFO+3+(ISORB-1)*8) = -1
      OrbTab(KOINFO+4+(ISORB-1)*8) = IPART
      OrbTab(KOINFO+5+(ISORB-1)*8) = INPART
      OrbTab(KOINFO+6+(ISORB-1)*8) = ISPART
      OrbTab(KOINFO+7+(ISORB-1)*8) = INSBP
    end do
  end do
end if
if (ISPART /= NSPART) then
  write(u6,*) 'NEWORBTAB Error: Nr of subpartitions'
  write(u6,*) 'generated does not match number of'
  write(u6,*) 'subpartitions computed!'
  write(u6,*) 'Generated ISPART=',ISPART
  write(u6,*) 'Computed  NSPART=',NSPART
end if
OrbTab(7) = NSPART
OrbTab(8) = NAPART
OrbTab(9) = NASPRT
OrbTab(10) = KSPART
OrbTab(KSPART:KSPART+NSPART-1) = 0
do ISORB=1,NSORBT
  ISPART = OrbTab(KOINFO+6+(ISORB-1)*8)
  N = 1+OrbTab(KSPART-1+ISPART)
  OrbTab(KSPART-1+ISPART) = N
end do

end subroutine NEWORBTAB

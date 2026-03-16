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
! Copyright (C) 2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine WRVEC_DYSON(filename,LUNIT,NSYM,NBAS,ORBNUM,CMO,AMPS,DYSEN,TITLE,NZ)
! Sort Dyson orbitals according to symmetry so that they can be
! written with WrVec

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One

implicit none
character(len=*) :: filename, TITLE
integer :: LUNIT, NSYM, NBAS(NSYM), ORBNUM, NZ
real*8 :: CMO(NZ,ORBNUM), AMPS(ORBNUM), DYSEN(ORBNUM)
integer :: DUMMY(7,8), I, J, NB(0:NSYM), NBAS_(NSYM), NBT, NORB(NSYM), NSYM_, OFF, ORBSYM(ORBNUM)
real*8 :: RSUM
logical :: DODESYM
real*8, allocatable :: DESYM(:,:), REORD(:)

! First count how many orbitals in each symmetry
NB(0) = 1
do I=1,NSYM
  NB(I) = NB(I-1)+NBAS(I)
end do
DODESYM = .false.
NORB(:) = 0
ORBSYM(:) = 0
outer: do I=1,ORBNUM
  do J=1,NSYM
    RSUM = sum(abs(CMO(NB(J-1):NB(J)-1,I)))
    if (RSUM > Zero) then
      ! If there is any orbital with mixture of symmetries,
      ! we have to desymmetrize the whole thing
      if (ORBSYM(I) /= 0) then
        DODESYM = .true.
        exit outer
      end if
      ORBSYM(I) = J
      NORB(J) = NORB(J)+1
    end if
  end do
end do outer
if (DODESYM) then
  NSYM_ = 1
  NBAS_(1) = sum(NBAS(1:NSYM))
  NB(1) = NBAS_(1)
  NORB(1) = ORBNUM
  ORBSYM(:) = 1
else
  NSYM_ = NSYM
  NBAS_(:) = NBAS(:)
end if
NBT = 0
do I=1,NSYM_
  NBT = NBT+NORB(I)*NBAS_(I)
end do
call mma_allocate(REORD,NBT,Label='REORD')
if (DODESYM) then
  ! Here do a plain desymmetrization
  NBT = NBAS_(1)
  call mma_allocate(DESYM,NBT,NBT,Label='DESYM')
  call get_dArray('SM',DESYM,NBT**2)
  call DGEMM_('N','N',NBT,NORB(1),NBT,One,DESYM,NBT,CMO,NBT,Zero,REORD,NBT)
  call mma_deallocate(DESYM)
else
  ! Here distribute the coefficients
  OFF = 0
  do I=1,NSYM_
    do J=1,ORBNUM
      if (ORBSYM(J) /= I) cycle
      REORD(OFF+1:OFF+NBAS_(I)) = CMO(NB(I-1):NB(I)-1,J)
      OFF = OFF+NBAS_(I)
    end do
  end do
end if

! And call WrVec with the reordered data
call WRVEC(filename,LUNIT,'COE',NSYM_,NBAS_,NORB,REORD,AMPS,DYSEN,DUMMY,TITLE)

call mma_deallocate(REORD)

return

end subroutine WRVEC_DYSON

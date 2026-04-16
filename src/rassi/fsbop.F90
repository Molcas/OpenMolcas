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

subroutine FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB,ICASE)

use rassi_global_arrays, only: FSBANN1, FSBANN2, FSBANN3, FSBANN4
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IMODE, ISORB, IORBTAB(*), ISSTAB(*), IFSBTAB(*), ICASE
integer(kind=iwp) :: IDET2, IERR, IFSB1, IFSB2, INSBPT, ISPART, ISST1, ISST2, ITYPE, KHSH2, KOINFO, KPOS1, KPOS2, KSSTAN, KSSTCR, &
                     KSSTOP, KSSTTB, KSTARR1, KSTARR2, LHSH2, LPOS, LSSTARR2, MORSBITS, NASPRT, NASPRT1, NASPRT2, NDET1, NDET2, &
                     NDETS2, NFSB1, NFSB2, NHEAD, NHSH2, NSBS1, NSBS2, NSSTARR2, NTAB2, NTAB2NEW, NULLV
integer(kind=iwp), pointer :: FSBANN(:)

! The orbital table:
NASPRT = IORBTAB(9)
KOINFO = 19
! Needed properties of spin-orbital ISORB
ISPART = IORBTAB(KOINFO+6+8*(ISORB-1))
INSBPT = IORBTAB(KOINFO+7+8*(ISORB-1))
! The substring table:
MORSBITS = ISSTAB(6)
KSSTTB = 15
KSSTAN = ISSTAB(9)
KSSTCR = ISSTAB(10)
if (IMODE == 1) then
  KSSTOP = KSSTCR
else
  KSSTOP = KSSTAN
end if

! The FS blocks of the PSI wave function:
ITYPE = IFSBTAB(2)
NFSB1 = IFSBTAB(3)
NASPRT1 = IFSBTAB(4)
KSTARR1 = 8
! Count how big the new table will be:
NASPRT2 = NASPRT1
KSTARR2 = KSTARR1
IFSB2 = 0
IDET2 = 0
do IFSB1=1,NFSB1
  KPOS1 = KSTARR1+(NASPRT1+2)*(IFSB1-1)
  NDET1 = IFSBTAB(KPOS1+NASPRT)
  ! The substring type to be annihilated from or created in:
  ISST1 = IFSBTAB(KPOS1-1+ISPART)
  ! The resulting substring type:
  ISST2 = ISSTAB(KSSTOP-1+INSBPT+MORSBITS*(ISST1-1))
  if (ISST2 == 0) cycle
  IFSB2 = IFSB2+1
  KPOS2 = KSTARR2+(NASPRT2+2)*(IFSB2-1)
  ! Replace substring type:
  ! Old vs. new nr of substrings:
  NSBS1 = ISSTAB(KSSTTB+5*(ISST1-1))
  NSBS2 = ISSTAB(KSSTTB+5*(ISST2-1))
  NDET2 = (NDET1*NSBS2)/NSBS1
  IDET2 = IDET2+NDET2
end do
NFSB2 = IFSB2
NDETS2 = IDET2
NHEAD = 7
!LSSTARR2 = 1+NHEAD
NSSTARR2 = (NASPRT2+2)*NFSB2
!LHSH2 = LSSTARR2+NSSTARR2
NHSH2 = 997+2*NFSB2
NTAB2 = NHEAD+NSSTARR2+2*NHSH2
! NTAB2 is now known. Make a new FSB table:
select case (ICASE)
  case (1)
    call mma_allocate(FSBANN1,NTAB2,Label='FSBANN1')
    FSBANN => FSBANN1(:)
  case (2)
    call mma_allocate(FSBANN2,NTAB2,Label='FSBANN2')
    FSBANN => FSBANN2(:)
  case (3)
    call mma_allocate(FSBANN3,NTAB2,Label='FSBANN3')
    FSBANN => FSBANN3(:)
  case (4)
    call mma_allocate(FSBANN4,NTAB2,Label='FSBANN4')
    FSBANN => FSBANN4(:)
  case DEFAULT
    write(u6,*) 'FSBOP: Illegal iCase value.'
    write(u6,*) 'ICASE=',ICASE
    call ABEND()
    FSBANN => FSBANN1(:) !dummy
end select

FSBANN(1) = NTAB2
FSBANN(2) = ITYPE
FSBANN(4) = NASPRT2
KSTARR2 = KSTARR1
IFSB2 = 0
IDET2 = 0
do IFSB1=1,NFSB1
  KPOS1 = KSTARR1+(NASPRT1+2)*(IFSB1-1)
  NDET1 = IFSBTAB(KPOS1+NASPRT)
  ! The substring type to be annihilated from or created in:
  ISST1 = IFSBTAB(KPOS1-1+ISPART)
  ! The resulting substring type:
  ISST2 = ISSTAB(KSSTOP-1+INSBPT+MORSBITS*(ISST1-1))
  if (ISST2 == 0) cycle
  IFSB2 = IFSB2+1
  KPOS2 = KSTARR2+(NASPRT2+2)*(IFSB2-1)
  FSBANN(KPOS2:KPOS2+NASPRT1-1) = IFSBTAB(KPOS1:KPOS1+NASPRT1-1)
  ! Replace substring type:
  FSBANN(KPOS2-1+ISPART) = ISST2
  ! Old vs. new nr of substrings:
  NSBS1 = ISSTAB(KSSTTB+5*(ISST1-1))
  NSBS2 = ISSTAB(KSSTTB+5*(ISST2-1))
  NDET2 = (NDET1*NSBS2)/NSBS1
  FSBANN(KPOS2+NASPRT) = NDET2
  FSBANN(KPOS2+NASPRT+1) = IDET2+1
  IDET2 = IDET2+NDET2
end do
! Store this block in the PSI2 hash structure.
NHEAD = 7
LSSTARR2 = 1+NHEAD
NSSTARR2 = (NASPRT2+2)*NFSB2
LHSH2 = LSSTARR2+NSSTARR2
NHSH2 = 997+2*NFSB2
NTAB2NEW = NHEAD+NSSTARR2+2*NHSH2
if (NTAB2 /= NTAB2NEW) then
  write(u6,*) ' FSBOP Error: NTAB2 /= NTAB2NEW!'
  write(u6,*) ' (This should be impossible!)'
  write(u6,*) ' Program RASSI is forced to stop, sorry!'
  call ABEND()
end if
! Make the hash map: NULLV is a null marker. Suggested value=-1.
NULLV = -1
call HSHINI(NHSH2,FSBANN(LHSH2:),NULLV)
KHSH2 = 1+NHEAD+NSSTARR2
! Header data:
FSBANN(1) = NTAB2
FSBANN(2) = ITYPE
FSBANN(3) = NFSB2
FSBANN(4) = NASPRT2
FSBANN(5) = NDETS2
FSBANN(6) = NHSH2
FSBANN(7) = KHSH2
! Store values in the map:
do IFSB2=1,NFSB2
  call HSHPUT(NASPRT2,NASPRT2+2,FSBANN(LSSTARR2:),NHSH2,FSBANN(LHSH2:),IFSB2)
end do

IERR = 0
do IFSB2=1,NFSB2
  LPOS = LSSTARR2+(NASPRT+2)*(IFSB2-1)
  do ISPART=1,NASPRT
    ISST2 = FSBANN(LPOS-1+ISPART)
    if (ISST2 < 1) IERR = IERR+1
  end do
end do
if (IERR > 0) then
  write(u6,*) ' Bad substrings in FSBOP!'
  write(u6,*) ' IERR=',IERR
  call PRFSBTAB(FSBANN(:))
  call ABEND()
end if

nullify(FSBANN)

end subroutine FSBOP

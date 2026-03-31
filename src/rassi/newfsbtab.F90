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

subroutine NEWFSBTAB(NACTEL,MSPIN2,LSYM,REST,SSTAB,ICASE)
! Purpose: Construct an FSB table and return its address in the
! FSBTAB1/FiSBTAB2 array.

use rassi_global_arrays, only: FSBARR, FSBTAB1, FSBTAB2
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: NACTEL, MSPIN2, LSYM, REST(*), SSTAB(*), ICASE
integer(kind=iwp) :: IERR, IFSB, ITYPE, JFSB, KHSHMAP, KORB, KREST, LSSTARR, NASPRT, NFSB, NFSB0, NHEAD, NHSHMAP, NLEN, NPART, &
                     NRDETS, NRDETS0, NSIZE, NSSTARR, NSYM
integer(kind=iwp), pointer :: FSBTAB(:)

! ITYPE=73 is the check code for this table.
ITYPE = 73
if (SSTAB(2) /= 19) then
  write(u6,*) ' NEWFSBTAB error: Not a Substring Table.'
  call ABEND()
end if
if (REST(2) /= 91) then
  write(u6,*) ' NEWFSBTAB error: Not a GAS Restriction Table.'
  call ABEND()
end if
NSYM = SSTAB(4)
NASPRT = SSTAB(5)
NPART = REST(3)
KORB = 5
KREST = KORB+(NSYM+1)*(NPART+1)
! Table consists of a 6-word header with data (see below), then
! an array dimensioned (NASPRT+2)*NFSB, and finally a hash map
! with suitable capacity e.g. 2*NFSB (50% usage). Each item in a
! hash map takes up 2 integers. Capacity must be at least NFSB+997.
call mkVERTAB(NACTEL,MSPIN2,LSYM,NPART,REST(KORB),REST(KREST),SSTAB,NFSB0,NRDETS0,NFSB,NRDETS)
NHEAD = 7
NSSTARR = (NASPRT+2)*NFSB
NHSHMAP = 997+2*NFSB
NSIZE = NHEAD+NSSTARR+2*NHSHMAP
select case (iCase)
  case (1)
    call mma_allocate(FSBTAB1,NSIZE,Label='FSBTAB1')
    FSBTAB => FSBTAB1(:)
  case (2)
    call mma_allocate(FSBTAB2,NSIZE,Label='FSBTAB2')
    FSBTAB => FSBTAB2(:)
  case default
    write(u6,*) 'NEWFSBTAB: Illegal ICASE value'
    write(u6,*) 'ICASE=',ICASE
    call Abend()
    FSBTAB => FSBTAB1(:) !dummy
end select
NLEN = (NASPRT+2)*NFSB
call ICOPY(NLEN,FSBARR(1:NLEN),1,FSBTAB(1+NHEAD:NLEN+NHEAD),1)
call mma_deallocate(FSBARR)

LSSTARR = 1+NHEAD

! Position of hash table ('Map')
KHSHMAP = 1+NHEAD+NSSTARR
FSBTAB(1) = NSIZE
FSBTAB(2) = ITYPE
FSBTAB(3) = NFSB
FSBTAB(4) = NASPRT
FSBTAB(5) = NRDETS
FSBTAB(6) = NHSHMAP
FSBTAB(7) = KHSHMAP
! Make the hash map: NULL is a null marker. Suggested value=-1.
!NULL = -1
!call HSHINI(NHSHMAP,FSBTAB(KHSHMAP),NULL)
! In conflict with the null() pointer
call HSHINI(NHSHMAP,FSBTAB(KHSHMAP:),-1)
! Store values in the map:
do IFSB=1,NFSB
  call HSHPUT(NASPRT,NASPRT+2,FSBTAB(LSSTARR:),NHSHMAP,FSBTAB(KHSHMAP:),IFSB)
end do
! Check that they can be obtained back:
IERR = 0
do IFSB=1,NFSB
  call HSHGET(FSBTAB(LSSTARR+(NASPRT+2)*(IFSB-1):),NASPRT,NASPRT+2,FSBTAB(LSSTARR:),NHSHMAP,FSBTAB(KHSHMAP:),JFSB)
  if (IFSB /= JFSB) IERR = IERR+1
end do
if (IERR > 0) then
  write(u6,*) 'NEWFSBTAB Hash index errors. IERR=',IERR
  call ABEND()
end if
nullify(FSBTAB)

end subroutine NEWFSBTAB

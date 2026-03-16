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

subroutine PrSSTab(SSTAB)

use cntrl, only: MORSBITS
use Definitions, only: u6

implicit none
integer SSTAB(*)
integer LPOS
integer NSSTP, ISSTP, KSSTP, NSBS, NPOP, ISYM, MS2, ISPART
integer KSSTANN, KSSTCRE, KSBSMRS, KMRSSBS, KSBSANN, KSBSCRE
integer I, ISBS, NSBSTOT, NMORS, IMRS, NASPRT
integer NRSBST, ISBSSTA, ISBSEND, IMRSSTA, IMRSEND

write(u6,*)
write(u6,*) '============================================='
write(u6,*) ' Substring table printout.'
write(u6,*) ' (A) Header and TOC:'
write(u6,'(a,i16)') '            Table size:',SSTAB(1)
write(u6,'(a,i16)') '       Table type code:',SSTAB(2)
write(u6,'(a,i16)') ' Orbital Table        :',SSTAB(3)
write(u6,'(a,i16)') ' Nr of symmetry labels:',SSTAB(4)
write(u6,'(a,i16)') ' Nr of active subpart :',SSTAB(5)
write(u6,'(a,i16)') ' Nr of bits/morsel    :',SSTAB(6)
write(u6,'(a,i16)') ' Nr of Substring Types:',SSTAB(7)
write(u6,'(a,i16)') ' Nr of Substrings     :',SSTAB(8)
write(u6,'(a,i16)') ' Substr Type Ann Table:',SSTAB(9)
write(u6,'(a,i16)') ' Substr Type Cre Table:',SSTAB(10)
write(u6,'(a,i16)') ' Substr/Morsel Table  :',SSTAB(11)
write(u6,'(a,i16)') ' Morsel/Substr Table  :',SSTAB(12)
write(u6,'(a,i16)') ' Substring Annih Table:',SSTAB(13)
write(u6,'(a,i16)') ' Substring Creat Table:',SSTAB(14)
write(u6,*) '---------------------------------------------'
write(u6,*) ' (B) Substring Types:'
write(u6,*) ' NSBS  = Nr of substrings of each type.'
write(u6,*) ' NPOP  = Electron population of each type.'
write(u6,*) ' ISYM  = Combined symmetry of each type.'
write(u6,*) ' MS2   = Twice spin proj of each type.'
write(u6,*) ' ISPART= Subpartition of each type.'
write(u6,*) '---------------------------------------------'
write(u6,*) '           NSBS   NPOP   ISYM    MS2  ISBPRT'
NSSTP = SSTAB(7)
NSBSTOT = SSTAB(8)
KSSTP = 15
do ISSTP=1,NSSTP
  NSBS = SSTAB(KSSTP+0+5*(ISSTP-1))
  NPOP = SSTAB(KSSTP+1+5*(ISSTP-1))
  ISYM = SSTAB(KSSTP+2+5*(ISSTP-1))
  MS2 = SSTAB(KSSTP+3+5*(ISSTP-1))
  ISPART = SSTAB(KSSTP+4+5*(ISSTP-1))
  write(u6,'(1x,6I7)') ISSTP,NSBS,NPOP,ISYM,MS2,ISPART
end do
write(u6,*) '---------------------------------------------'
write(u6,*) ' (C1) Morsels to Substrings:'
NASPRT = SSTAB(5)
KSBSMRS = SSTAB(11)
KMRSSBS = SSTAB(12)
NMORS = 2**MORSBITS
do IMRSSTA=0,NMORS-1,15
  IMRSEND = min(NMORS-1,IMRSSTA+14)
  write(u6,*)
  write(u6,'(1x,a,3x,15I5)') 'MRS:',(IMRS,IMRS=IMRSSTA,IMRSEND)
  do ISPART=1,NASPRT
    write(u6,'(1x,a,i3,15I5)') 'SBS:',ISPART,(SSTAB(KMRSSBS+2*(IMRS+NMORS*(ISPART-1))),IMRS=IMRSSTA,IMRSEND)
  end do
end do

NRSBST = SSTAB(8)
write(u6,*) ' (C2) Substrings to Morsels:'
do ISBSSTA=1,NRSBST,15
  ISBSEND = min(NRSBST,ISBSSTA+14)
  write(u6,*)
  write(u6,'(1x,a,15I5)') 'SBS:',(ISBS,ISBS=ISBSSTA,ISBSEND)
  write(u6,'(1x,a,15I5)') 'MRS:',(SSTAB(KSBSMRS+2*(ISBS-1)),ISBS=ISBSSTA,ISBSEND)
end do
write(u6,*) '---------------------------------------------'
write(u6,*) ' (D1) Substring Type Annihilation Table:'
KSSTANN = SSTAB(9)
KSSTCRE = SSTAB(10)
do ISSTP=1,NSSTP
  LPOS = KSSTANN+MORSBITS*(ISSTP-1)
  write(u6,'(1x,I8,5X,8I8)') ISSTP,(SSTAB(LPOS-1+I),I=1,MORSBITS)
end do
write(u6,*) ' (D2) Substring Type Creation Table:'
do ISSTP=1,NSSTP
  LPOS = KSSTCRE+MORSBITS*(ISSTP-1)
  write(u6,'(1x,I8,5X,8I8)') ISSTP,(SSTAB(LPOS-1+I),I=1,MORSBITS)
end do
write(u6,*) '---------------------------------------------'
write(u6,*) ' (E1) Substring Annihilation Table:'
KSBSANN = SSTAB(13)
KSBSCRE = SSTAB(14)
do ISBS=1,NSBSTOT
  LPOS = KSBSANN+MORSBITS*(ISBS-1)
  write(u6,'(1x,I8,5X,8I8)') ISBS,(SSTAB(LPOS-1+I),I=1,MORSBITS)
end do
write(u6,*) '---------------------------------------------'
write(u6,*) ' (E2) Substring Creation Table:'
do ISBS=1,NSBSTOT
  LPOS = KSBSCRE+MORSBITS*(ISBS-1)
  write(u6,'(1x,I8,5X,8I8)') ISBS,(SSTAB(LPOS-1+I),I=1,MORSBITS)
end do
write(u6,*) '============================================='

end subroutine PrSSTab

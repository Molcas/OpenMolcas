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

subroutine PrOrbTab(ORBTAB)

use Definitions, only: u6

implicit none
integer ORBTAB(*)
character(len=8) STRING8
character(len=8), external :: ORBNAM
integer IEXTNUM, INPART, INSBP, ISORB, IPART, ISMLAB, ISOIND
integer ISPART, ISPLAB, NSPORB, KOINFO, KSPART
integer NSPART

write(u6,*)
write(u6,*) '============================================='
write(u6,*) ' Orbital table printout.'
write(u6,'(a,i16)') '            Table size:',ORBTAB(1)
write(u6,'(a,i16)') '       Table type code:',ORBTAB(2)
write(u6,'(a,i16)') '   Nr of spin-orbitals:',ORBTAB(3)
write(u6,'(a,i16)') ' Nr of active sp-orbs :',ORBTAB(4)
write(u6,'(a,i16)') ' Nr of symmetry labels:',ORBTAB(5)
write(u6,'(a,i16)') ' Nr of partitions     :',ORBTAB(6)
write(u6,'(a,i16)') ' Nr of sub-partitions :',ORBTAB(7)
write(u6,'(a,i16)') ' Nr of active part    :',ORBTAB(8)
write(u6,'(a,i16)') ' Nr of active subpart :',ORBTAB(9)
write(u6,*) '---------------------------------------------'
write(u6,*) ' IEXTNUM  = Orbital index, in external order.'
write(u6,*) ' ISMLAB= Symmetry label'
write(u6,*) ' ISOIND= In-Symmetry orbital index, external order.'
write(u6,*) ' ISPLAB= Spin component label'
write(u6,*) ' IPART= Orbital partition label'
write(u6,*) ' INPART= In-Partition orbital index'
write(u6,*) ' ISPART= Orbital sub-partition label'
write(u6,*) ' INSBP = In-Subpartition orbital index'
write(u6,*) ' ORBNAM= Orbital name'
write(u6,*) '---------------------------------------------'
NSPORB = ORBTAB(3)
KOINFO = 19
NSPART = ORBTAB(7)
KSPART = ORBTAB(10)
write(u6,*) '         IEXTNUM ISMLAB ISOIND ISPLAB  IPART INPART ISPART  INSBP ORBNAM'
do ISORB=1,NSPORB
  IEXTNUM = ORBTAB(KOINFO+0+(ISORB-1)*8)
  ISMLAB = ORBTAB(KOINFO+1+(ISORB-1)*8)
  ISOIND = ORBTAB(KOINFO+2+(ISORB-1)*8)
  ISPLAB = ORBTAB(KOINFO+3+(ISORB-1)*8)
  IPART = ORBTAB(KOINFO+4+(ISORB-1)*8)
  INPART = ORBTAB(KOINFO+5+(ISORB-1)*8)
  ISPART = ORBTAB(KOINFO+6+(ISORB-1)*8)
  INSBP = ORBTAB(KOINFO+7+(ISORB-1)*8)
  STRING8 = ORBNAM(ISORB,ORBTAB)
  write(u6,'(1x,9I7,2X,A8)') ISORB,IEXTNUM,ISMLAB,ISOIND,ISPLAB,IPART,INPART,ISPART,INSBP,STRING8
end do
write(u6,*)
write(u6,*) ' Nr of spin orbitals in each subpartition:'
write(u6,'(1x,20i5)') (ORBTAB(KSPART-1+ISPART),ISPART=1,NSPART)
write(u6,*) '============================================='

end subroutine PrOrbTab

************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE PrOrbTab(ORBTAB)
      IMPLICIT NONE
      INTEGER ORBTAB(*)
      CHARACTER(LEN=8)  STRING8
      CHARACTER(LEN=8), EXTERNAL :: ORBNAM
      INTEGER IEXTNUM,INPART,INSBP,ISORB,IPART,ISMLAB,ISOIND
      INTEGER ISPART,ISPLAB,NSPORB,KOINFO,KSPART
      INTEGER NSPART
      WRITE(6,*)
      WRITE(6,*)'============================================='
      WRITE(6,*)' Orbital table printout.'
      WRITE(6,'(a,i16)')'            Table size:',ORBTAB(1)
      WRITE(6,'(a,i16)')'       Table type code:',ORBTAB(2)
      WRITE(6,'(a,i16)')'   Nr of spin-orbitals:',ORBTAB(3)
      WRITE(6,'(a,i16)')' Nr of active sp-orbs :',ORBTAB(4)
      WRITE(6,'(a,i16)')' Nr of symmetry labels:',ORBTAB(5)
      WRITE(6,'(a,i16)')' Nr of partitions     :',ORBTAB(6)
      WRITE(6,'(a,i16)')' Nr of sub-partitions :',ORBTAB(7)
      WRITE(6,'(a,i16)')' Nr of active part    :',ORBTAB(8)
      WRITE(6,'(a,i16)')' Nr of active subpart :',ORBTAB(9)
      WRITE(6,*)'---------------------------------------------'
      WRITE(6,*)' IEXTNUM  = Orbital index, in external order.'
      WRITE(6,*)' ISMLAB= Symmetry label'
      WRITE(6,*)' ISOIND= In-Symmetry orbital index, external order.'
      WRITE(6,*)' ISPLAB= Spin component label'
      WRITE(6,*)' IPART= Orbital partition label'
      WRITE(6,*)' INPART= In-Partition orbital index'
      WRITE(6,*)' ISPART= Orbital sub-partition label'
      WRITE(6,*)' INSBP = In-Subpartition orbital index'
      WRITE(6,*)' ORBNAM= Orbital name'
      WRITE(6,*)'---------------------------------------------'
      NSPORB= ORBTAB(3)
      KOINFO=19
      NSPART=ORBTAB(7)
      KSPART=ORBTAB(10)
      WRITE(6,*)'         IEXTNUM ISMLAB ISOIND ISPLAB  IPART INPART'
     &           //' ISPART  INSBP ORBNAM'
      DO ISORB=1,NSPORB
        IEXTNUM= ORBTAB(KOINFO+ 0+(ISORB-1)*8)
        ISMLAB = ORBTAB(KOINFO+ 1+(ISORB-1)*8)
        ISOIND = ORBTAB(KOINFO+ 2+(ISORB-1)*8)
        ISPLAB = ORBTAB(KOINFO+ 3+(ISORB-1)*8)
        IPART  = ORBTAB(KOINFO+ 4+(ISORB-1)*8)
        INPART = ORBTAB(KOINFO+ 5+(ISORB-1)*8)
        ISPART = ORBTAB(KOINFO+ 6+(ISORB-1)*8)
        INSBP  = ORBTAB(KOINFO+ 7+(ISORB-1)*8)
        STRING8= ORBNAM(ISORB,ORBTAB)
        WRITE(6,'(1x,9I7,2X,A8)') ISORB,IEXTNUM,ISMLAB,ISOIND,ISPLAB,
     &                           IPART,INPART,ISPART,INSBP,STRING8
      END DO
      WRITE(6,*)
      WRITE(6,*)' Nr of spin orbitals in each subpartition:'
      WRITE(6,'(1x,20i5)') (ORBTAB(KSPART-1+ISPART),ISPART=1,NSPART)
      WRITE(6,*)'============================================='

      END SUBROUTINE PrOrbTab

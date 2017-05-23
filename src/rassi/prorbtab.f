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
      SUBROUTINE PrOrbTab(LORB)
      IMPLICIT NONE
      INTEGER LORB
#include "WrkSpc.fh"
      CHARACTER*8 STRING8,ORBNAM
      EXTERNAL ORBNAM
      INTEGER IEXTNUM,INPART,INSBP,ISORB,IPART,ISMLAB,ISOIND
      INTEGER ISPART,ISPLAB,NSPORB,KOINFO,LOINFO,KSPART,LSPART
      INTEGER NSPART
      WRITE(6,*)
      WRITE(6,*)'============================================='
      WRITE(6,*)' Orbital table printout.'
      WRITE(6,'(a,i16)')'     Workspace pointer:',LORB
      WRITE(6,'(a,i16)')'            Table size:',IWORK(LORB)
      WRITE(6,'(a,i16)')'       Table type code:',IWORK(LORB+1)
      WRITE(6,'(a,i16)')'   Nr of spin-orbitals:',IWORK(LORB+2)
      WRITE(6,'(a,i16)')' Nr of active sp-orbs :',IWORK(LORB+3)
      WRITE(6,'(a,i16)')' Nr of symmetry labels:',IWORK(LORB+4)
      WRITE(6,'(a,i16)')' Nr of partitions     :',IWORK(LORB+5)
      WRITE(6,'(a,i16)')' Nr of sub-partitions :',IWORK(LORB+6)
      WRITE(6,'(a,i16)')' Nr of active part    :',IWORK(LORB+7)
      WRITE(6,'(a,i16)')' Nr of active subpart :',IWORK(LORB+8)
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
      NSPORB= IWORK(LORB+2)
      KOINFO=19
      LOINFO=LORB-1+KOINFO
      NSPART=IWORK(LORB+6)
      KSPART=IWORK(LORB+9)
      LSPART=LORB-1+KSPART
      WRITE(6,*)'         IEXTNUM ISMLAB ISOIND ISPLAB  IPART INPART'
     &           //' ISPART  INSBP ORBNAM'
      DO ISORB=1,NSPORB
        IEXTNUM= IWORK(LOINFO+ 0+(ISORB-1)*8)
        ISMLAB= IWORK(LOINFO+ 1+(ISORB-1)*8)
        ISOIND= IWORK(LOINFO+ 2+(ISORB-1)*8)
        ISPLAB= IWORK(LOINFO+ 3+(ISORB-1)*8)
        IPART= IWORK(LOINFO+ 4+(ISORB-1)*8)
        INPART= IWORK(LOINFO+ 5+(ISORB-1)*8)
        ISPART= IWORK(LOINFO+ 6+(ISORB-1)*8)
        INSBP= IWORK(LOINFO+ 7+(ISORB-1)*8)
        STRING8=ORBNAM(ISORB,LORB)
        WRITE(6,'(1x,9I7,2X,A8)') ISORB,IEXTNUM,ISMLAB,ISOIND,ISPLAB,
     &                           IPART,INPART,ISPART,INSBP,STRING8
      END DO
      WRITE(6,*)
      WRITE(6,*)' Nr of spin orbitals in each subpartition:'
      WRITE(6,'(1x,20i5)') (IWORK(LSPART-1+ISPART),ISPART=1,NSPART)
      WRITE(6,*)'============================================='
      RETURN
      END

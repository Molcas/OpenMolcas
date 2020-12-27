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
      INTEGER FUNCTION NEWORBTAB(IPRTTAB)
      IMPLICIT NONE
      INTEGER LORB
      INTEGER IPRTTAB(*)
#include "WrkSpc.fh"
#include "Morsel.fh"
      INTEGER I,N,IEXTNUM,INPART,INSBP,ISORB,IPART,ISMLAB,ISOIND
      INTEGER ISPART,ISUM,KSPART
C     INTEGER IPFR,IPIN,IPAC,IPSE,IPDE,IP
      INTEGER IPFR,IPIN,IPAC,IPSE
      INTEGER ISYM
      INTEGER NPART,NSPART,NORBT,NSORBT,NSYM
      INTEGER NAPART,NASPRT,NASPO,KOINFO,LOINFO,NTAB
      INTEGER NOES(8),INSYM(8)

      ISOIND = 0 ! dummy initialize
C Executable statements
      NPART= IPRTTAB(3)
      NSYM = IPRTTAB(4)
C NAPART: Nr of partitions of active orbitals.
      NAPART=NPART-4
C Partitions for inactive, secondary,frozen, and deleted orbitals:
      IPIN=NAPART+1
      IPSE=NAPART+2
      IPFR=NAPART+3
C Total nr of orbitals:
      NORBT=IPRTTAB(5)
      NSORBT=2*NORBT
C Table words 1--10 contain some header info.
C Table words 11--18 contain start index of each CMO symmetry block
C Table words 19-- contain info for each separate spin orbital
C Presently 8 table entries for each spin orbital.
      KOINFO=19
C Table words 19+8*NSORBT-- contain nr of sp-orbs for each subpartition
      KSPART=19+8*NSORBT
C Nr of sub-partitions:
      NSPART=0
      DO IPART=1,NPART
       N=2*IPRTTAB(5+(NSYM+1)*IPART)
       IF(N.GT.0) NSPART=NSPART+(N+MORSBITS-1)/MORSBITS
      END DO
      NTAB= KSPART+NSPART-1
C Allocate the orbital table.
      CALL GETMEM('ORBTAB','ALLOCATE','INTEGER',LORB,NTAB)
      IWORK(LORB)=NTAB
      IWORK(LORB+1)= 1
      IWORK(LORB+2)=NSORBT
C Nr of active spin-orbitals:
      ISUM=0
      DO IPART=1,NAPART
        ISUM=ISUM+IPRTTAB(5+(NSYM+1)*IPART)
      END DO
      NASPO=2*ISUM
      IWORK(LORB+3)=NASPO
      IWORK(LORB+4)=NSYM
      IWORK(LORB+5)=NAPART+4
C Accumulated nr of orbital functions/symm:
      ISUM=0
      DO ISYM=1,NSYM
        NOES(ISYM)=ISUM
        ISUM=ISUM+ IPRTTAB(5+ISYM)
      END DO
C Relative pointers to CMO symmetry blocks:
      ISUM=0
      DO ISYM=1,NSYM
        IWORK(LORB+8+ISYM)=ISUM+1
        ISUM=ISUM+ IPRTTAB(5+ISYM)**2
      END DO
C Spin orbital number:
      ISORB=0
C First, active orbitals by partition, and by symmetry
C Previous MO indices within each symmetry.
      DO ISYM=1,NSYM
        INSYM(ISYM)=IPRTTAB(5+ISYM+(NSYM+1)*IPFR)+
     &                 IPRTTAB(5+ISYM+(NSYM+1)*IPIN)
      END DO
C Increase subpartition index as needed.
      ISPART=0
      LOINFO=LORB-1+KOINFO
      DO IPART=1,NAPART
       N= IPRTTAB(5+(NSYM+1)*IPART)
       IF(N.EQ.0) GOTO 100
       ISPART=ISPART+1
       INSBP=0
       INPART=0
       DO ISYM=1,NSYM
        ISMLAB=ISYM
        N= IPRTTAB(5+ISYM+(NSYM+1)*IPART)
        DO I=1,N
         ISOIND=1+INSYM(ISYM)
         INSYM(ISYM)=ISOIND
         IEXTNUM=NOES(ISYM)+ISOIND
C Next spin orbital, with alpha spin:
         ISORB=ISORB+1
         INPART=INPART+1
         INSBP=INSBP+1
         IF(INSBP.GT.MORSBITS) THEN
           ISPART=ISPART+1
           INSBP=INSBP-MORSBITS
         END IF
C Fill in table:
          IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
          IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
          IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
          IWORK(LOINFO+ 3+(ISORB-1)*8)=1
          IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
          IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
          IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
          IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
C Next spin orbital, same, but beta spin:
         ISORB=ISORB+1
         INPART=INPART+1
         INSBP=INSBP+1
         IF(INSBP.GT.MORSBITS) THEN
           ISPART=ISPART+1
           INSBP=INSBP-MORSBITS
         END IF
C Fill in table:
          IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
          IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
          IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
          IWORK(LOINFO+ 3+(ISORB-1)*8)=-1
          IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
          IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
          IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
          IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
        END DO
       END DO
 100   CONTINUE
      END DO
      NASPRT=ISPART
C Inactive:
C Must set up start index within each symmetry.
      DO ISYM=1,NSYM
        INSYM(ISYM)=IPRTTAB(5+ISYM+(NSYM+1)*IPFR)
      END DO
      IPART=NAPART+1
      N= IPRTTAB(5+(NSYM+1)*IPART)
      IF(N.EQ.0) GOTO 200
      ISPART=ISPART+1
      INPART=0
      INSBP=0
      DO ISYM=1,NSYM
       ISMLAB=ISYM
C       N=NISH(ISYM)
       N= IPRTTAB(5+ISYM+(NSYM+1)*IPART)
       DO I=1,N
        ISOIND=ISOIND+1
        IEXTNUM=NOES(ISYM)+ISOIND
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=-1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
       END DO
      END DO
 200  CONTINUE
C Secondary:
      IPART=NAPART+2
      N= IPRTTAB(5+(NSYM+1)*IPART)
      IF(N.EQ.0) GOTO 300
      ISPART=ISPART+1
      INPART=0
      INSBP=0
C Must set up start index within each symmetry.
      DO ISYM=1,NSYM
        N=IPRTTAB(5+ISYM+(NSYM+1)*IPFR)
        N=N+IPRTTAB(5+ISYM+(NSYM+1)*IPIN)
        DO IPAC=1,NAPART
         N=N+IPRTTAB(5+ISYM+(NSYM+1)*IPAC)
        END DO
        INSYM(ISYM)=N
      END DO

      DO ISYM=1,NSYM
       ISMLAB=ISYM
       ISOIND=INSYM(ISYM)
C       N=NSSH(ISYM)
       N= IPRTTAB(5+ISYM+(NSYM+1)*IPART)
       DO I=1,N
        ISOIND=ISOIND+1
        IEXTNUM=NOES(ISYM)+ISOIND
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=-1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
       END DO
      END DO
 300  CONTINUE
C Frozen:
      IPART=NAPART+3
      N= IPRTTAB(5+(NSYM+1)*IPART)
      IF(N.EQ.0) GOTO 400
      ISPART=ISPART+1
      INPART=0
      INSBP=0
      DO ISYM=1,NSYM
       ISMLAB=ISYM
C        N=NFRO(ISYM)
       N= IPRTTAB(5+ISYM+(NSYM+1)*IPART)
       ISOIND=0
       DO I=1,N
        ISOIND=ISOIND+1
        IEXTNUM=NOES(ISYM)+ISOIND
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=-1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
       END DO
      END DO
 400  CONTINUE
C Deleted:
      IPART=NAPART+4
      N= IPRTTAB(5+(NSYM+1)*IPART)
      IF(N.EQ.0) GOTO 500
      ISPART=ISPART+1
      INPART=0
      INSBP=0
C Must set up start index within each symmetry.
      DO ISYM=1,NSYM
        N=IPRTTAB(5+ISYM+(NSYM+1)*IPFR)
        N=N+IPRTTAB(5+ISYM+(NSYM+1)*IPIN)
        DO IPAC=1,NAPART
         N=N+IPRTTAB(5+ISYM+(NSYM+1)*IPAC)
        END DO
        N=N+IPRTTAB(5+ISYM+(NSYM+1)*IPSE)
        INSYM(ISYM)=N
      END DO
      DO ISYM=1,NSYM
       ISMLAB=ISYM
       ISOIND= INSYM(ISYM)
C        N=NDEL(ISYM)
       N= IPRTTAB(5+ISYM+(NSYM+1)*IPART)
       DO I=1,N
        ISOIND=ISOIND+1
        IEXTNUM=NOES(ISYM)+ISOIND
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
        ISORB=ISORB+1
        INPART=INPART+1
        INSBP=INSBP+1
        IF(INSBP.GT.MORSBITS) THEN
          ISPART=ISPART+1
          INSBP=INSBP-MORSBITS
        END IF
         IWORK(LOINFO+ 0+(ISORB-1)*8)=IEXTNUM
         IWORK(LOINFO+ 1+(ISORB-1)*8)=ISMLAB
         IWORK(LOINFO+ 2+(ISORB-1)*8)=ISOIND
         IWORK(LOINFO+ 3+(ISORB-1)*8)=-1
         IWORK(LOINFO+ 4+(ISORB-1)*8)=IPART
         IWORK(LOINFO+ 5+(ISORB-1)*8)=INPART
         IWORK(LOINFO+ 6+(ISORB-1)*8)=ISPART
         IWORK(LOINFO+ 7+(ISORB-1)*8)=INSBP
       END DO
      END DO
 500  CONTINUE
      IF(ISPART.NE.NSPART) THEN
        WRITE(6,*)'NEWORBTAB Error: Nr of subpartitions'
        WRITE(6,*)'generated does not match number of'
        WRITE(6,*)'subpartitions computed!'
        WRITE(6,*)'Generated ISPART=',ISPART
        WRITE(6,*)'Computed  NSPART=',NSPART
      END IF
      IWORK(LORB+6)=NSPART
      IWORK(LORB+7)=NAPART
      IWORK(LORB+8)=NASPRT
      IWORK(LORB+9)=KSPART
      CALL ICOPY(NSPART,[0],0,IWORK(LORB-1+KSPART),1)
      DO ISORB=1,NSORBT
        ISPART=IWORK(LOINFO+ 6+(ISORB-1)*8)
        N=1+IWORK(LORB-1+KSPART-1+ISPART)
        IWORK(LORB-1+KSPART-1+ISPART)=N
      END DO
      NEWORBTAB=LORB
      RETURN
      END

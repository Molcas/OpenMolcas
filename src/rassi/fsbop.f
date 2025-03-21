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
      SUBROUTINE FSBOP(IMODE,ISORB,IORBTAB,ISSTAB,IFSBTAB,ICASE)
      use rassi_global_arrays, only: FSBANN1, FSBANN2, FSBANN3,
     &                               FSBANN4, FSBANN
      use stdalloc, only: mma_allocate
      IMPLICIT NONE
      INTEGER IMODE,ISORB
      INTEGER IORBTAB(*),NASPRT
      INTEGER ISSTAB(*)
      INTEGER IFSBTAB(*)
      INTEGER ICASE

      INTEGER KOINFO,ISPART,INSBPT,MORSBITS
      INTEGER KSSTAN,KSSTCR
      INTEGER KSSTOP,ITYPE,NFSB1,NASPRT1
      INTEGER KSTARR1,NTAB2,NTAB2NEW,NASPRT2
      INTEGER KSTARR2,IFSB2,IDET2,IFSB1,KPOS1
      INTEGER ISST1,ISST2,KPOS2,NDETS2,NHEAD
      INTEGER LSSTARR2,NSSTARR2,NFSB2,LHSH2,NHSH2
      INTEGER KSSTTB,NSBS1,NSBS2,NDET1,NDET2
      INTEGER NULL,KHSH2,IERR,LPOS

C The orbital table:
      NASPRT= IORBTAB(9)
      KOINFO=19
C Needed properties of spin-orbital ISORB
      ISPART=IORBTAB(KOINFO+6+8*(ISORB-1))
      INSBPT=IORBTAB(KOINFO+7+8*(ISORB-1))
C The substring table:
      MORSBITS=ISSTAB(6)
      KSSTTB=15
      KSSTAN=ISSTAB( 9)
      KSSTCR=ISSTAB(10)
      IF(IMODE.EQ.1) THEN
        KSSTOP=KSSTCR
      ELSE
        KSSTOP=KSSTAN
      END IF

C The FS blocks of the PSI wave function:
      ITYPE  =IFSBTAB(2)
      NFSB1  =IFSBTAB(3)
      NASPRT1=IFSBTAB(4)
      KSTARR1=8
C Count how big the new table will be:
      NASPRT2=NASPRT1
      KSTARR2=KSTARR1
      IFSB2=0
      IDET2=0
      DO IFSB1=1,NFSB1
        KPOS1=KSTARR1+(NASPRT1+2)*(IFSB1-1)
        NDET1=IFSBTAB(KPOS1+NASPRT)
C The substring type to be annihilated from or created in:
        ISST1=IFSBTAB(KPOS1-1+ISPART)
C The resulting substring type:
        ISST2=ISSTAB(KSSTOP-1+INSBPT+MORSBITS*(ISST1-1))
        if (ISST2 == 0) cycle
        IFSB2=IFSB2+1
        KPOS2=KSTARR2+(NASPRT2+2)*(IFSB2-1)
C Replace substring type:
C Old vs. new nr of substrings:
        NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))
        NSBS2=ISSTAB(KSSTTB+5*(ISST2-1))
        NDET2=(NDET1*NSBS2)/NSBS1
        IDET2=IDET2+NDET2
      END DO
      NFSB2=IFSB2
      NDETS2=IDET2
      NHEAD=7
*      LSSTARR2=1+NHEAD
      NSSTARR2=(NASPRT2+2)*NFSB2
*      LHSH2=LSSTARR2+NSSTARR2
      NHSH2=997+2*NFSB2
      NTAB2=NHEAD+NSSTARR2+2*NHSH2
C NTAB2 is now known. Make a new FSB table:
      SELECT CASE(ICASE)
      CASE(1)
         CALL mma_allocate(FSBANN1,NTAB2,Label='FSBANN1')
         FSBANN=>FSBANN1(:)
      CASE(2)
         CALL mma_allocate(FSBANN2,NTAB2,Label='FSBANN2')
         FSBANN=>FSBANN2(:)
      CASE(3)
         CALL mma_allocate(FSBANN3,NTAB2,Label='FSBANN3')
         FSBANN=>FSBANN3(:)
      CASE(4)
         CALL mma_allocate(FSBANN4,NTAB2,Label='FSBANN4')
         FSBANN=>FSBANN4(:)
      CASE DEFAULT
         Write (6,*) 'FSBOP: Illegal iCase value.'
         Write (6,*) 'ICASE=',ICASE
         Call ABEND()
      END SELECT

      FSBANN(1)=NTAB2
      FSBANN(2)=ITYPE
      FSBANN(4)=NASPRT2
      KSTARR2=KSTARR1
      IFSB2=0
      IDET2=0
      DO IFSB1=1,NFSB1
        KPOS1=KSTARR1+(NASPRT1+2)*(IFSB1-1)
        NDET1=IFSBTAB(KPOS1+NASPRT)
C The substring type to be annihilated from or created in:
        ISST1=IFSBTAB(KPOS1-1+ISPART)
C The resulting substring type:
        ISST2=ISSTAB(KSSTOP-1+INSBPT+MORSBITS*(ISST1-1))
        if (ISST2 == 0) cycle
        IFSB2=IFSB2+1
        KPOS2=KSTARR2+(NASPRT2+2)*(IFSB2-1)
        CALL ICOPY(NASPRT1,IFSBTAB(KPOS1:KPOS1+NASPRT1-1),1,
     &                      FSBANN(KPOS2:KPOS2+NASPRT1-1),1)
C Replace substring type:
        FSBANN(KPOS2-1+ISPART)=ISST2
C Old vs. new nr of substrings:
        NSBS1=ISSTAB(KSSTTB+5*(ISST1-1))
        NSBS2=ISSTAB(KSSTTB+5*(ISST2-1))
        NDET2=(NDET1*NSBS2)/NSBS1
        FSBANN(KPOS2+NASPRT  )=NDET2
        FSBANN(KPOS2+NASPRT+1)=IDET2+1
        IDET2=IDET2+NDET2
      END DO
C Store this block in the PSI2 hash structure.
      NHEAD=7
      LSSTARR2=1+NHEAD
      NSSTARR2=(NASPRT2+2)*NFSB2
      LHSH2=LSSTARR2+NSSTARR2
      NHSH2=997+2*NFSB2
      NTAB2NEW=NHEAD+NSSTARR2+2*NHSH2
      IF(NTAB2.NE.NTAB2NEW) THEN
        WRITE(6,*)' FSBOP Error: NTAB2.NE.NTAB2NEW!'
        WRITE(6,*)' (This should be impossible!)'
        WRITE(6,*)' Program RASSI is forced to stop, sorry!'
        CALL ABEND()
      END IF
C Make the hash map: NULL is a null marker. Suggested value=-1.
      NULL=-1
      CALL HSHINI(NHSH2,FSBANN(LHSH2:),NULL)
      KHSH2=1+NHEAD+NSSTARR2
C Header data:
      FSBANN(1)=NTAB2
      FSBANN(2)=ITYPE
      FSBANN(3)=NFSB2
      FSBANN(4)=NASPRT2
      FSBANN(5)=NDETS2
      FSBANN(6)=NHSH2
      FSBANN(7)=KHSH2
C Store values in the map:
      DO IFSB2=1,NFSB2
        CALL HSHPUT(NASPRT2,NASPRT2+2,FSBANN(LSSTARR2:),
     &                         NHSH2,FSBANN(LHSH2:),IFSB2)
      END DO

      IERR=0
      DO IFSB2=1,NFSB2
        LPOS=LSSTARR2+(NASPRT+2)*(IFSB2-1)
        DO ISPART=1,NASPRT
          ISST2=FSBANN(LPOS-1+ISPART)
          IF(ISST2.LT.1) IERR=IERR+1
        END DO
      END DO
      IF(IERR.GT.0) THEN
        WRITE(6,*)' Bad substrings in FSBOP!'
        WRITE(6,*)' IERR=',IERR
        CALL PRFSBTAB(FSBANN(:))
        CALL ABEND()
      END IF

      END SUBROUTINE FSBOP

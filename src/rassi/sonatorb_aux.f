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
      SUBROUTINE SONATORB_PLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB_PLOT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      CHARACTER*25 FNAME
      CHARACTER(LEN=*) FILEBASE
      CHARACTER*16 KNUM
      CHARACTER*16 FNUM,XNUM
      CHARACTER*8 CHARTYPE
      CHARACTER CDIR
      INTEGER ASS,BSS
      DIMENSION Dummy(1),iDummy(7,8)




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c LSZZ  - AO Overlap integral
c LVEC  - AO Overlap eigenvectors
c LEIG  - AO Overlap eigenvalues
c LVEC2 - Eigenvectors of density matrix
c LSCR  - Temporary for matrix multiplication
C NOTE: LSCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NBTRI)
      CALL GETMEM('VEC   ','ALLO','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','ALLO','REAL',LVEC2,NBMX2)
      CALL GETMEM('SCR   ','ALLO','REAL',LSCR,NBMX2)
      CALL GETMEM('EIG   ','ALLO','REAL',LEIG,NBST)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSZZ),1)
      CALL DCOPY_(NBSQ,[0.0D00],0,WORK(LVEC),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LVEC2),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
      CALL DCOPY_(NBST,[0.0D00],0,WORK(LEIG),1)

      CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
      CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
      CALL DCOPY_(NBSQ,[0.0D00],0,WORK(LVNAT),1)
      CALL DCOPY_(NBST,[0.0D00],0,WORK(LOCC),1)

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSZZ
      LV=LVEC
      LE=LEIG
      CALL FZERO(WORK(LVEC),NBSQ)
      DO 700 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 620 I=1,NB**2,(NB+1)
          WORK(LV-1+I)=1.0D00
620      CONTINUE
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 630 I=1,NB
          EIG=WORK(LS1)
          WORK(LE1)=EIG
          X=1.0D00/SQRT(MAX(EIG,1.0D-14))
          CALL DSCAL_(NB,X,WORK(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
630      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
700   CONTINUE

      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NBTRI)

      CALL GETMEM('TDMAT ','ALLO','REAL',LDMAT,NBMX2)

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

        CDIR='?'
        IF(IDIR.EQ.1) CDIR='X'
        IF(IDIR.EQ.2) CDIR='Y'
        IF(IDIR.EQ.3) CDIR='Z'

        INV=1
        II2=0
        IOCC=0
        LV=LVEC
        LE=LEIG
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) GOTO 1750

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LDMAT),1)
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)/2.0d0
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)/2.0d0
            ELSE
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMAT),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCR),NB,
     &                 0.0D0,WORK(LDMAT),NB)

          ID1=1
          ID2=1
          DO I=1,NB
            EIG=WORK(LE-1+I)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
          ISCR=LSCR
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              WORK(ISCR)=(WORK(LDMAT-1+IJ)+WORK(LDMAT-1+JI))/2.0d0

c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) WORK(ISCR)=WORK(ISCR)*2.0d0
              ISCR=ISCR+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LVEC2),1)
          CALL DCOPY_(NB,[1.0D0],0,WORK(LVEC2),NB+1)

          CALL JACOB(WORK(LSCR),WORK(LVEC2),NB,NB)
          CALL JACORD(WORK(LSCR),WORK(LVEC2),NB,NB)

C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=LSCR-1
          DO I=1,NB
            II=II+I
            WORK(LOCC-1+IOCC+NB+1-I)=WORK(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVALUES IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2),NB,
     &                 0.0D0,WORK(LSCR),NB)
          I1=LSCR
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,WORK(I1),1,WORK(LVNAT-1+I2),1)
            I1=I1+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
1750      CONTINUE
        END DO

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2)') ASS,".",BSS
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".",CDIR
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNAT), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_PLOT", WORK(LVNAT), 1, 4)

c    ONLYFOR NATURAL ORBITALS
      if(ITYPE.EQ.1)
     &       CALL ADD_INFO("SONATORB_NO_OCC", WORK(LOCC), SUM(NBASF), 4)

      END DO

      CALL GETMEM('TDMAT ','FREE','REAL',LDMAT,NBMX2)
      CALL GETMEM('VEC   ','FREE','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','FREE','REAL',LVEC2,NBMX2)
      CALL GETMEM('SCR   ','FREE','REAL',LSCR,NBMX2)
      CALL GETMEM('EIG   ','FREE','REAL',LEIG,NBST)
      CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
      CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)


      RETURN
      END




      SUBROUTINE SONATORB_CPLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SONATORB_CPLOT')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION DENS(6,NBTRI)
      CHARACTER*25 FNAME
      CHARACTER(LEN=*) FILEBASE
      CHARACTER*16 KNUM
      CHARACTER*16 FNUM,XNUM
      CHARACTER*8 CHARTYPE
      CHARACTER CDIR
      INTEGER ASS,BSS
      DIMENSION IDUM(1),Dummy(1),iDummy(7,8)




CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PLOTTING SECTION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Get the proper type of the property
      ITYPE=0
      IF(CHARTYPE.EQ.'HERMSING') ITYPE=1
      IF(CHARTYPE.EQ.'ANTISING') ITYPE=2
      IF(CHARTYPE.EQ.'HERMTRIP') ITYPE=3
      IF(CHARTYPE.EQ.'ANTITRIP') ITYPE=4
      IF(ITYPE.EQ.0) THEN
        WRITE(6,*)'RASSI/SONATORB internal error.'
        WRITE(6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c LSZZ  - AO Overlap integral
c LVEC  - AO Overlap eigenvectors
c LEIG  - AO Overlap eigenvalues
c LVEC2 - Eigenvectors of density matrix
c LSCR  - Temporary for matrix multiplication
C NOTE: LSCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL GETMEM('SZZ   ','ALLO','REAL',LSZZ,NBTRI)
      CALL GETMEM('VEC   ','ALLO','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','ALLO','REAL',LVEC2,NBMX2)
      CALL GETMEM('VEC2I  ','ALLO','REAL',LVEC2I,NBMX2)
      CALL GETMEM('SCR   ','ALLO','REAL',LSCR,NBMX2)
      CALL GETMEM('SCRI   ','ALLO','REAL',LSCRI,NBMX2)
      CALL GETMEM('EIG   ','ALLO','REAL',LEIG,NBST)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSZZ),1)
      CALL DCOPY_(NBSQ,[0.0D00],0,WORK(LVEC),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LVEC2),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LVEC2I),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
      CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCRI),1)
      CALL DCOPY_(NBST,[0.0D00],0,WORK(LEIG),1)

      CALL GETMEM('VNAT  ','ALLO','REAL',LVNAT,NBSQ)
      CALL GETMEM('VNATI  ','ALLO','REAL',LVNATI,NBSQ)
      CALL GETMEM('OCC   ','ALLO','REAL',LOCC,NBST)
      CALL DCOPY_(NBSQ,[0.0D00],0,WORK(LVNAT),1)
      CALL DCOPY_(NBSQ,[0.0D00],0,WORK(LVNATI),1)
      CALL DCOPY_(NBST,[0.0D00],0,WORK(LOCC),1)

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=6
      ICMP=1
      ISYLAB=1
      CALL RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSZZ),ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=LSZZ
      LV=LVEC
      LE=LEIG
      CALL FZERO(WORK(LVEC),NBSQ)
      DO 1700 ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO 1620 I=1,NB**2,(NB+1)
          WORK(LV-1+I)=1.0D00
1620      CONTINUE
        CALL JACOB(WORK(LS),WORK(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO 1630 I=1,NB
          EIG=WORK(LS1)
          WORK(LE1)=EIG
          X=1.0D00/SQRT(MAX(EIG,1.0D-14))
          CALL DSCAL_(NB,X,WORK(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
1630      CONTINUE
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
1700   CONTINUE

      CALL GETMEM('SZZ   ','FREE','REAL',LSZZ,NBTRI)

      CALL GETMEM('TDMAT ','ALLO','REAL',LDMAT,NBMX2)
      CALL GETMEM('TDMATI ','ALLO','REAL',LDMATI,NBMX2)

      IF(ITYPE.LE.2) THEN
        ISTART=3
        IEND=3
      ELSE
        ISTART=1
        IEND=3
      END IF

      DO IDIR=ISTART,IEND

          CDIR='?'
          IF(IDIR.EQ.1) CDIR='X'
          IF(IDIR.EQ.2) CDIR='Y'
          IF(IDIR.EQ.3) CDIR='Z'


cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
C read in ao matrix for angmom or mltpl
      CALL GETMEM('SANG  ','ALLO','REAL',LSANG,NBTRI)
      CALL DCOPY_(NBTRI,[0.0D00],0,WORK(LSANG),1)

      IRC=-1
      IOPT=6

      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        ICMP=1
        CALL iRDONE(IRC,   1,'MLTPL  0',ICMP,IDUM,       ISYLAB)
        IF (IRC.EQ.0) NSIZ=IDUM(1)
        CALL  RDONE(IRC,IOPT,'MLTPL  0',ICMP,WORK(LSANG),ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      MLTPL0 INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
        ICMP=3
        CALL iRDONE(IRC,   1,'ANGMOM  ',ICMP,IDUM,       ISYLAB)
        IF (IRC.EQ.0) NSIZ=IDUM(1)
        CALL  RDONE(IRC,IOPT,'ANGMOM  ',ICMP,WORK(LSANG),ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(6,*)
          WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(6,*)'      ANGMOM INTEGRALS ARE NOT AVAILABLE'
          WRITE(6,*)'      IRC:',IRC
          WRITE(6,*)
          CALL ABEND()
        END IF

      END IF

cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
cccccccccccccccccccccccc
        INV=1
        II2=0
        IOCC=0
        LV=LVEC
        LE=LEIG
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) GOTO 1800

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LDMAT),1)
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LDMATI),1)
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCRI),1)

          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)/2.0d0
              WORK(LDMAT-1+JI)=DENS(IDIR,II2)/2.0d0
              WORK(LDMATI-1+IJ)=-1.0d0*DENS(IDIR+3,II2)/2.0d0
              WORK(LDMATI-1+JI)=DENS(IDIR+3,II2)/2.0d0
            ELSE
              WORK(LDMAT-1+IJ)=DENS(IDIR,II2)
              WORK(LDMATI-1+JI)=DENS(IDIR+3,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMAT),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LDMATI),NB,WORK(LV),NB,
     &                 0.0D0,WORK(LSCRI),NB)



          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCR),NB,
     &                 0.0D0,WORK(LDMAT),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LSCRI),NB,
     &                 0.0D0,WORK(LDMATI),NB)

          ID1=1
          ID2=1
          DO I=1,NB
            EIG=WORK(LE-1+I)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMAT-1+ID2),1)
            CALL DSCAL_(NB,EIG,WORK(LDMATI-1+ID1),NB)
            CALL DSCAL_(NB,EIG,WORK(LDMATI-1+ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCR),1)
          CALL DCOPY_(NBMX2,[0.0D00],0,WORK(LSCRI),1)

          ISCR=LSCR
          ISCRI=LSCRI
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              WORK(ISCR)=(WORK(LDMAT-1+JI)+WORK(LDMAT-1+IJ))/2.0d0
              WORK(ISCRI)=(WORK(LDMATI-1+JI)-WORK(LDMATI-1+IJ))/2.0d0
c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) WORK(ISCR)=WORK(ISCR)*2.0d0
              IF(ITYPE.GE.3) WORK(ISCRI)=WORK(ISCRI)*2.0d0
              ISCR=ISCR+1
              ISCRI=ISCRI+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LVEC2),1)
          CALL DCOPY_(NBMX2,[0.0D0],0,WORK(LVEC2I),1)

          CALL CPLOT_DIAG(WORK(LSCR),WORK(LSCRI), NB,
     &                    WORK(LVEC2),WORK(LVEC2I))

C LAPACK ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=LSCR-1
          DO I=1,NB
            II=II+I
            WORK(LOCC-1+IOCC+NB+1-I)=WORK(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2),NB,
     &                 0.0D0,WORK(LSCR),NB)
          CALL DGEMM_('N','N',NB,NB,NB,1.0D0,
     &                 WORK(LV),NB,WORK(LVEC2I),NB,
     &                 0.0D0,WORK(LSCRI),NB)

          I1=LSCR
          I1I=LSCRI
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,WORK(I1),1,WORK(LVNAT-1+I2),1)
            CALL DCOPY_(NB,WORK(I1I),1,WORK(LVNATI-1+I2),1)
            I1=I1+NB
            I1I=I1I+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
1800      CONTINUE
        END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC TESTING
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IPGLOB.GE.DEBUG) THEN

      CALL GETMEM('SANGF ','ALLO','REAL',LSANGF,NBMX**2)
      CALL GETMEM('SANGTR  ','ALLO','REAL',LSANGTR,NBMX**2)
      CALL GETMEM('SANGTI  ','ALLO','REAL',LSANGTI,NBMX**2)
      CALL GETMEM('SANGTR2  ','ALLO','REAL',LSANGTR2,NBMX**2)
      CALL GETMEM('SANGTI2  ','ALLO','REAL',LSANGTI2,NBMX**2)
      CALL DCOPY_(NBMX**2,[0.0D00],0,WORK(LSANGF),1)
      CALL DCOPY_(NBMX**2,[0.0D00],0,WORK(LSANGTR),1)
      CALL DCOPY_(NBMX**2,[0.0D00],0,WORK(LSANGTI),1)
      CALL DCOPY_(NBMX**2,[0.0D00],0,WORK(LSANGTR2),1)
      CALL DCOPY_(NBMX**2,[0.0D00],0,WORK(LSANGTI2),1)

      INV=0
      INV2=0
      II=0
      SUM = 0.0d0
      SUMI = 0.0d0

      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        IF(NB.EQ.0) GOTO 1860

c       Expand integrals for this symmetry to full storage
        CALL DCOPY_(NBMX**2,[0.0d0],0,WORK(LSANGF),1)

        DO J=1,NB
        DO I=1,J
          IJ=NB*(J-1)+I-1
          JI=NB*(I-1)+J-1

          WORK(LSANGF+JI) = WORK(LSANG+II)

          IF(I.NE.J) THEN
            IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
              WORK(LSANGF+IJ) = 1.0d0 * WORK(LSANG+II)
            ELSE
              WORK(LSANGF+IJ) = WORK(LSANG+II)
            END IF
          END IF

          II=II+1

        END DO
        END DO

        IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNAT+INV),NB,0.0d0,WORK(LSANGTR),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &              WORK(LVNATI+INV),NB,0.0d0,WORK(LSANGTI),NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTR2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTR2),NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTI2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTI2),NB)

        ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNAT+INV),NB,0.0d0,WORK(LSANGTI),NB)
          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LSANGF),NB,
     &             WORK(LVNATI+INV),NB,0.0d0,WORK(LSANGTR),NB)

          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTR2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTR2),NB)

          CALL DGEMM_('T','N',NB,NB,NB,-1.0d0,WORK(LVNATI+INV),NB,
     &             WORK(LSANGTR),NB,0.0d0,WORK(LSANGTI2),NB)
          CALL DGEMM_('T','N',NB,NB,NB,1.0d0,WORK(LVNAT+INV),NB,
     &             WORK(LSANGTI),NB,1.0d0,WORK(LSANGTI2),NB)

        END IF

c Sum over the trace
        DO I = 1,NB
          IJ = I+(I-1)*NB-1
          SUM  = SUM  + WORK(LOCC-1+I+INV2) * WORK(LSANGTR2+IJ)
          SUMI = SUMI + WORK(LOCC-1+I+INV2) * WORK(LSANGTI2+IJ)
        END DO

1860    CONTINUE

        INV=INV+NB**2
        INV2=INV2+NB

      END DO

        WRITE(6,*) "Ben P TEST for JA:"
        WRITE(6,*) "REAL: ",SUM
        WRITE(6,*) "IMAG: ",SUMI

        CALL GETMEM('SANGF ','FREE','REAL',LSANGF,NBMX**2)
        CALL GETMEM('SANGTR  ','FREE','REAL',LSANGTR,NBMX**2)
        CALL GETMEM('SANGTI  ','FREE','REAL',LSANGTI,NBMX**2)
        CALL GETMEM('SANGTR2  ','FREE','REAL',LSANGTR2,NBMX**2)
        CALL GETMEM('SANGTI2  ','FREE','REAL',LSANGTI2,NBMX**2)
      END IF ! IPGLOB >= DEBUG

      CALL GETMEM('SANG  ','FREE','REAL',LSANG,NBTRI)

C WRITE OUT THIS SET OF NATURAL SPIN ORBITALS
C REAL PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","R"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","R"
       END IF
       WRITE(FNUM,'(I8)') BSS
       FNUM=ADJUSTL(FNUM)
       IF (ASS.NE.BSS) THEN
         WRITE(XNUM,'(I8,A)') ASS,'_'//TRIM(FNUM)
         FNUM=ADJUSTL(XNUM)
       END IF
       IF (ITYPE.GT.2) FNUM=CDIR//TRIM(FNUM)

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.R'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNAT), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

C IMAGINARY PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","I"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","I"
       END IF

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.I'
       IF(ITYPE.EQ.1)
     &        WRITE(6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        IFOCC=1
        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &     WORK(LVNATI), WORK(LOCC), Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_CPLOTR", WORK(LVNAT), 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTI", WORK(LVNATI), 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTO", WORK(LOCC), 1, 4)

      END DO

      CALL GETMEM('TDMAT ','FREE','REAL',LDMAT,NBMX2)
      CALL GETMEM('TDMATI ','FREE','REAL',LDMATI,NBMX2)
      CALL GETMEM('VEC   ','FREE','REAL',LVEC,NBSQ)
      CALL GETMEM('VEC2  ','FREE','REAL',LVEC2,NBMX2)
      CALL GETMEM('VEC2I  ','FREE','REAL',LVEC2I,NBMX2)
      CALL GETMEM('SCR   ','FREE','REAL',LSCR,NBMX2)
      CALL GETMEM('SCRI   ','FREE','REAL',LSCRI,NBMX2)
      CALL GETMEM('EIG   ','FREE','REAL',LEIG,NBST)
      CALL GETMEM('VNAT  ','FREE','REAL',LVNAT,NBSQ)
      CALL GETMEM('VNATI  ','FREE','REAL',LVNATI,NBSQ)
      CALL GETMEM('OCC   ','FREE','REAL',LOCC,NBST)

      RETURN
      END





      SUBROUTINE CPLOT_DIAG(MATR, MATI, DIM, EIGVECR, EIGVECI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER DIM
      REAL*8 MATR(DIM*(DIM+1)/2),MATI(DIM*(DIM+1)/2)
      REAL*8 EIGVECR(DIM,DIM),EIGVECI(DIM,DIM)
      REAL*8 CEIGVAL(DIM)
      COMPLEX*16 MATFULL((DIM*(DIM+1)/2))
      COMPLEX*16 CEIGVEC(DIM,DIM)
      COMPLEX*16 ZWORK(2*DIM-1)
      REAL*8 RWORK(3*DIM-2)
      INTEGER INFO

      DO J=1,(DIM*(DIM+1)/2)
          MATFULL(J) = DCMPLX(MATR(J),MATI(J))
c          MATFULL(J) = DCMPLX(MATR(J),0.0d0)
      END DO


      call zhpev_('V','U',DIM,MATFULL,CEIGVAL,
     &           CEIGVEC,DIM,ZWORK,RWORK,INFO)


      IF(INFO.NE.0) THEN
          WRITE(6,*) "Error in diagonalization"
          WRITE(6,*) "INFO: ",INFO
          CALL ABEND()
      END IF

      DO I=1,DIM
      DO J=1,DIM
          EIGVECR(I,J) = REAL(CEIGVEC(I,J))
          EIGVECI(I,J) = AIMAG(CEIGVEC(I,J))
      END DO
      END DO

      CALL DCOPY_(DIM*(DIM+1)/2,[0.0D00],0,MATR,1)
      CALL DCOPY_(DIM*(DIM+1)/2,[0.0D00],0,MATI,1)

      DO J=1,DIM
         MATR((J*(J-1)/2)+J) = CEIGVAL(J)
      END DO

      RETURN
      END

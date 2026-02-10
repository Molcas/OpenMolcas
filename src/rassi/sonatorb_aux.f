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
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, One, Two
      use OneDat, only: sNoNuc, sNoOri
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NBTRI,NBMX,NBASF,NBSQ,NBST

      IMPLICIT None
      real(kind=wp), intent(in):: DENS(6,NBTRI)
      CHARACTER(LEN=*), intent(in):: FILEBASE
      CHARACTER(LEN=8), intent(in):: CHARTYPE
      INTEGER(KIND=IWP), INTENT(IN):: ASS,BSS

      CHARACTER(LEN=25) FNAME
      CHARACTER(LEN=16) KNUM
      CHARACTER(LEN=16) FNUM,XNUM
      CHARACTER(LEN=8) LABEL
      CHARACTER CDIR
      Real(kind=wp) Dummy(1)
      Integer(kind=iwp) iDummy(7,8)
      Real(kind=wp), allocatable:: SZZ(:), VEC(:), VEC2(:), DMAT(:),
     &                             SCR(:)
      Real(kind=wp), allocatable:: VNAT(:), EIG(:), OCC(:)
      Integer(kind=iwp) ITYPE, NBMX2, IRC, IOPT, ICMP, ISYLAB, LS, LV,
     &                  LE, ISYM, NB, I, LS1, LV1, LE1, ISTART, IEND,
     &                  IDIR, INV, II2, IOCC, J, IJ, JI, ID1, ID2, ISCR,
     &                  II, I1, I2, LuXXVEC
      Integer(kind=iwp), External:: IsFreeUnit
      REAL(kind=wp) X

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
        WRITE(u6,*)'RASSI/SONATORB internal error.'
        WRITE(u6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c SZZ  - AO Overlap integral
c VEC  - AO Overlap eigenvectors
c EIG  - AO Overlap eigenvalues
c VEC2 - Eigenvectors of density matrix
c SCR  - Temporary for matrix multiplication
C NOTE: SCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL mma_allocate(SZZ,NBTRI,Label='SZZ')
      CALL mma_allocate(VEC,NBSQ,Label='VEC')
      CALL mma_allocate(VEC2,NBMX2,Label='VEC2')
      CALL mma_allocate(SCR,NBMX2,Label='SCR')
      CALL mma_allocate(EIG,NBST,Label='EIG')
      SZZ(:)=Zero
      VEC(:)=Zero
      VEC2(:)=Zero
      SCR(:)=Zero
      EIG(:)=Zero

      CALL mma_allocate(VNAT,NBSQ,Label='VNAT')
      VNAT(:)=Zero
      CALL mma_allocate(OCC,NBST,Label='OCC')
      OCC(:)=Zero

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(u6,*)
        WRITE(u6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(u6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(u6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      VEC(:)=Zero
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO I=1,NB**2,(NB+1)
          VEC(LV-1+I)=One
        END DO
        CALL JACOB(SZZ(LS),VEC,NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          EIG(LE1)=SZZ(LS1)
          X=One/SQRT(MAX(SZZ(LS1),1.0D-14))
          CALL DSCAL_(NB,X,VEC(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO

      CALL mma_deallocate(SZZ)

      CALL mma_allocate(DMAT,NBMX2,Label='DMAT')

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
        LV=1
        LE=1
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) CYCLE

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          DMAT(:)=Zero
          CALL DCOPY_(NBMX2,[Zero],0,SCR,1)
          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              DMAT(IJ)=DENS(IDIR,II2)/Two
              DMAT(JI)=DENS(IDIR,II2)/Two
            ELSE
              DMAT(IJ)=DENS(IDIR,II2)
              DMAT(JI)=DENS(IDIR,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 DMAT,NB,VEC(LV),NB,
     &                 Zero,SCR,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,
     &                 VEC(LV),NB,SCR,NB,
     &                 Zero,DMAT,NB)

          ID1=1
          ID2=1
          DO I=1,NB
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          SCR(:)=Zero
          ISCR=1
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              SCR(ISCR)=(DMAT(IJ)+DMAT(JI))/Two

c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) SCR(ISCR)=SCR(ISCR)*Two
              ISCR=ISCR+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          CALL DCOPY_(NBMX2,[Zero],0,VEC2,1)
          CALL DCOPY_(NB,[One],0,VEC2,NB+1)

          CALL JACOB(SCR,VEC2,NB,NB)
          CALL JACORD(SCR,VEC2,NB,NB)

C JACORD ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=0
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=SCR(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVALUES IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 VEC(LV),NB,VEC2,NB,
     &                 Zero,SCR,NB)
          I1=1
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
            I1=I1+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB

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
     &        WRITE(u6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(u6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(u6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(u6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNAT, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_PLOT", VNAT, 1, 4)

c    ONLYFOR NATURAL ORBITALS
      if(ITYPE.EQ.1)
     &       CALL ADD_INFO("SONATORB_NO_OCC", OCC, SUM(NBASF), 4)

      END DO

      CALL mma_deallocate(DMAT)
      CALL mma_deallocate(VEC)
      CALL mma_deallocate(VEC2)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(EIG)
      CALL mma_deallocate(VNAT)
      CALL mma_deallocate(OCC)

      END SUBROUTINE SONATORB_PLOT

      SUBROUTINE SONATORB_CPLOT (DENS, FILEBASE, CHARTYPE, ASS, BSS)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, One, Two
      use OneDat, only: sNoNuc, sNoOri, sOpSiz
      use rassi_aux, only: ipglob
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep
      use rassi_data, only: NBTRI,NBMX,NBASF,NBSQ,NBST
      IMPLICIT NONE
      Real(kind=wp) DENS(6,NBTRI)
      CHARACTER(LEN=*) FILEBASE
      CHARACTER(LEN=8) CHARTYPE
      INTEGER(kind=iwp) ASS,BSS

      CHARACTER(LEN=25) FNAME
      CHARACTER(LEN=16) KNUM
      CHARACTER(LEN=16) FNUM,XNUM
      CHARACTER(LEN=8) LABEL
      CHARACTER CDIR
      Real(kind=wp) Dummy(1)
      Integer(kind=iwp) IDUM(1),iDummy(7,8)
      Real(kind=wp), Allocatable:: SZZ(:), VEC(:), VEC2(:), VEC2I(:),
     &                             SCR(:)
      Real(kind=wp), Allocatable:: SCRI(:), EIG(:)
      Real(kind=wp), Allocatable:: VNAT(:), VNATI(:), OCC(:)
      Real(kind=wp), Allocatable:: DMAT(:), DMATI(:)
      Real(kind=wp), Allocatable:: SANG(:)
      Real(kind=wp), Allocatable:: SANGF(:), SANGTR(:), SANGTI(:)
      Real(kind=wp), Allocatable:: SANGTR2(:), SANGTI2(:)

      Integer(kind=iwp) ITYPE, NBMX2, IRC, IOPT, ICMP, ISYLAB, LS, LV,
     &                  LE, ISYM, NB, I, LS1, LV1, LE1, ISTART, IEND,
     &                  IDIR, INV, II2, IOCC, J, IJ, JI, ID1, ID2, ISCR,
     &                  II, I1, I2, LuXXVEC, JOPT, I1I, INV2, ISCRI
      Integer(kind=iwp), External:: IsFreeUnit
      REAL(kind=wp) X, SUM, SUMI


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
        WRITE(u6,*)'RASSI/SONATORB internal error.'
        WRITE(u6,*)'Erroneous property type:',CHARTYPE
        CALL ABEND()
      END IF

      NBMX2=NBMX**2

c SZZ  - AO Overlap integral
c VEC  - AO Overlap eigenvectors
c EIG  - AO Overlap eigenvalues
c VEC2 - Eigenvectors of density matrix
c SCR  - Temporary for matrix multiplication
C NOTE: SCR COULD PROBABLY BE SOMETHING LIKE NBMX*(NBMX+1)/2
C       ALTHOUGH IT PROBABLY DOESN'T SAVE MUCH
C       (JACOB TAKES A TRIANGULAR MATRIX LIKE ZHPEV DOES?)
      CALL mma_allocate(SZZ,NBTRI,Label='SZZ')
      SZZ(:)=Zero
      CALL mma_allocate(VEC,NBSQ,Label='VEC')
      VEC(:)=Zero
      CALL mma_allocate(VEC2,NBMX2,Label='VEC2')
      VEC2(:)=Zero
      CALL mma_allocate(VEC2I,NBMX2,Label='VEC2I')
      VEC2I(:)=Zero
      CALL mma_allocate(SCR,NBMX2,Label='SCR')
      SCR(:)=Zero
      CALL mma_allocate(SCRI,NBMX2,Label='SCRI')
      SCRI(:)=Zero
      CALL mma_allocate(EIG,NBST,Label='EIG')
      EIG(:)=Zero

      CALL mma_allocate(VNAT,NBSQ,Label='VNAT')
      VNAT(:)=Zero
      CALL mma_allocate(VNATI,NBSQ,Label='VNATI')
      VNATI(:)=Zero
      CALL mma_allocate(OCC,NBST,Label='OCC')
      OCC(:)=Zero

C READ ORBITAL OVERLAP MATRIX.
      IRC=-1

c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      CALL RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF ( IRC.NE.0 ) THEN
        WRITE(u6,*)
        WRITE(u6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(u6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(u6,*)
        CALL ABEND()
      ENDIF


C DIAGONALIZE EACH SYMMETRY BLOCK OF THE OVERLAP MATRIX.
      LS=1
      LV=1
      LE=1
      VEC(:)=Zero
      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        DO I=1,NB**2,(NB+1)
          VEC(LV-1+I)=One
        END DO
        CALL JACOB(SZZ(LS),VEC(LV),NB,NB)
C SCALE EACH VECTOR TO OBTAIN AN ORTHONORMAL BASIS.
        LS1=LS
        LV1=LV
        LE1=LE
        DO I=1,NB
          EIG(LE1)=SZZ(LS1)
          X=One/SQRT(MAX(SZZ(LS1),1.0D-14))
          CALL DSCAL_(NB,X,VEC(LV1),1)
          LS1=LS1+I+1
          LV1=LV1+NB
          LE1=LE1+1
        END DO
        LS=LS+(NB*(NB+1))/2
        LV=LV+NB**2
        LE=LE+NB
      END DO

      CALL mma_deallocate(SZZ)

      CALL mma_allocate(DMAT,NBMX2,Label='DMAT')
      CALL mma_allocate(DMATI,NBMX2,Label='DMATI')

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
      CALL mma_allocate(SANG,NBTRI,Label='SANG')
      SANG(:)=Zero

      IRC=-1
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      JOPT=ibset(0,sOpSiz)

      IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
        ICMP=1
        LABEL='MLTPL  0'
        CALL iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
        CALL  RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(u6,*)
          WRITE(u6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(u6,*)'      MLTPL0 INTEGRALS ARE NOT AVAILABLE'
          WRITE(u6,*)'      IRC:',IRC
          WRITE(u6,*)
          CALL ABEND()
        END IF

      ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
        ICMP=3
        LABEL='ANGMOM'
        CALL iRDONE(IRC,JOPT,LABEL,ICMP,IDUM,ISYLAB)
        CALL  RDONE(IRC,IOPT,LABEL,ICMP,SANG,ISYLAB)

        IF ( IRC.NE.0 ) THEN
          WRITE(u6,*)
          WRITE(u6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
          WRITE(u6,*)'      ANGMOM INTEGRALS ARE NOT AVAILABLE'
          WRITE(u6,*)'      IRC:',IRC
          WRITE(u6,*)
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
        LV=1
        LE=1
        DO ISYM=1,NSYM
          NB=NBASF(ISYM)
          IF(NB.EQ.0) cycle

C TRANSFORM TO ORTHONORMAL BASIS. THIS REQUIRES THE CONJUGATE
C BASIS, BUT SINCE WE USE CANONICAL ON BASIS THIS AMOUNTS TO A
C SCALING WITH THE EIGENVALUES OF THE OVERLAP MATRIX:

C expand the triangular matrix for this symmetry to a square matrix
          DMAT(:)=Zero
          DMATI(:)=Zero
          SCR(:)=Zero
          SCRI(:)=Zero

          DO J=1,NB
          DO I=1,J
            II2=II2+1
            IJ=NB*(J-1)+I
            JI=NB*(I-1)+J
            IF(I.NE.J) THEN
              DMAT(IJ)=DENS(IDIR,II2)/Two
              DMAT(JI)=DENS(IDIR,II2)/Two
              DMATI(IJ)=-DENS(IDIR+3,II2)/Two
              DMATI(JI)= DENS(IDIR+3,II2)/Two
            ELSE
              DMAT(IJ)=DENS(IDIR,II2)
              DMATI(JI)=DENS(IDIR+3,II2)
            END IF
          END DO
          END DO

          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 DMAT,NB,VEC(LV),NB,
     &                 Zero,SCR,NB)
          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 DMATI,NB,VEC(LV),NB,
     &                 Zero,SCRI,NB)



          CALL DGEMM_('T','N',NB,NB,NB,One,
     &                 VEC(LV),NB,SCR,NB,
     &                 Zero,DMAT,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,
     &                 VEC(LV),NB,SCRI,NB,
     &                 Zero,DMATI,NB)

          ID1=1
          ID2=1
          DO I=1,NB
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMAT(ID2),1)
            CALL DSCAL_(NB,EIG(LE-1+I),DMATI(ID1),NB)
            CALL DSCAL_(NB,EIG(LE-1+I),DMATI(ID2),1)
            ID1=ID1+1
            ID2=ID2+NB
          END DO


C SYMMETRIZE THIS BLOCK INTO SCRATCH AREA, TRIANGULAR STORAGE:
          SCR(:)=Zero
          SCRI(:)=Zero

          ISCR=1
          ISCRI=1
          DO I=1,NB
            DO J=1,I
              IJ=I+NB*(J-1)
              JI=J+NB*(I-1)
c simple averaging
              SCR(ISCR)=(DMAT(JI)+DMAT(IJ))/Two
              SCRI(ISCRI)=(DMATI(JI)-DMATI(IJ))/Two
c add a factor of two to convert spin -> sigma
              IF(ITYPE.GE.3) SCR(ISCR)=SCR(ISCR)*Two
              IF(ITYPE.GE.3) SCRI(ISCRI)=SCRI(ISCRI)*Two
              ISCR=ISCR+1
              ISCRI=ISCRI+1
            END DO
          END DO

C DIAGONALIZE THE DENSITY MATRIX BLOCK:
          VEC2(:)=Zero
          VEC2I(:)=Zero

          CALL CPLOT_DIAG(SCR,SCRI, NB,VEC2,VEC2I)

C LAPACK ORDERS BY INCREASING EIGENVALUE. REVERSE THIS ORDER.
          II=0
          DO I=1,NB
            II=II+I
            OCC(IOCC+NB+1-I)=SCR(II)
          END DO
          IOCC=IOCC+NB

C REEXPRESS THE EIGENVECTORS IN AO BASIS FUNCTIONS. REVERSE ORDER.
          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 VEC(LV),NB,VEC2,NB,
     &                 Zero,SCR,NB)
          CALL DGEMM_('N','N',NB,NB,NB,One,
     &                 VEC(LV),NB,VEC2I,NB,
     &                 Zero,SCRI,NB)

          I1=1
          I1I=1
          I2=INV+NB**2
          DO I=1,NB
            I2=I2-NB
            CALL DCOPY_(NB,SCR(I1),1,VNAT(I2),1)
            CALL DCOPY_(NB,SCRI(I1I),1,VNATI(I2),1)
            I1=I1+NB
            I1I=I1I+NB
          END DO
          INV=INV+NB**2
          LV=LV+NB**2
          LE=LE+NB
        END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCC TESTING
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF(IPGLOB.GE.4) THEN

      CALL mma_allocate(SANGF,NBMX**2,Label='SANGF')
      SANGF(:)=Zero
      CALL mma_allocate(SANGTR,NBMX**2,Label='SANGTR')
      CALL mma_allocate(SANGTI,NBMX**2,Label='SANGTI')
      SANGTR(:)=Zero
      SANGTI(:)=Zero
      CALL mma_allocate(SANGTR2,NBMX**2,Label='SANGTR2')
      CALL mma_allocate(SANGTI2,NBMX**2,Label='SANGTI2')
      SANGTR2(:)=Zero
      SANGTI2(:)=Zero

      INV=0
      INV2=0
      II=0
      SUM = Zero
      SUMI = Zero

      DO ISYM=1,NSYM
        NB=NBASF(ISYM)
        IF (NB/=0) THEN

c       Expand integrals for this symmetry to full storage
        SANGF(:)=Zero

        DO J=1,NB
        DO I=1,J
          IJ=NB*(J-1)+I-1
          JI=NB*(I-1)+J-1

          SANGF(1+JI) = SANG(1+II)

          IF(I.NE.J) THEN
            IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN
              SANGF(1+IJ) = -SANG(1+II)
            ELSE
              SANGF(1+IJ) =  SANG(1+II)
            END IF
          END IF

          II=II+1

        END DO
        END DO

        IF(ITYPE.EQ.1.OR.ITYPE.EQ.3) THEN
          CALL DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,
     &             VNAT(1+INV),NB,Zero,SANGTR,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,
     &              VNATI(1+INV),NB,Zero,SANGTI,NB)

          CALL DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,
     &             SANGTR,NB,Zero,SANGTR2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,VNATI(1+INV),NB,
     &             SANGTI,NB,One,SANGTR2,NB)

          CALL DGEMM_('T','N',NB,NB,NB,-One,VNATI(1+INV),NB,
     &             SANGTR,NB,Zero,SANGTI2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,
     &             SANGTI,NB,One,SANGTI,NB)

        ELSE IF(ITYPE.EQ.2.OR.ITYPE.EQ.4) THEN

          CALL DGEMM_('T','N',NB,NB,NB,One,SANGF,NB,
     &             VNAT(1+INV),NB,Zero,SANGTI,NB)
          CALL DGEMM_('T','N',NB,NB,NB,-One,SANGF,NB,
     &             VNATI(1+INV),NB,Zero,SANGTR,NB)

          CALL DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,
     &             SANGTR,NB,Zero,SANGTR2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,VNATI(1+INV),NB,
     &             SANGTI,NB,One,SANGTR2,NB)

          CALL DGEMM_('T','N',NB,NB,NB,-One,VNATI(1+INV),NB,
     &             SANGTR,NB,Zero,SANGTI2,NB)
          CALL DGEMM_('T','N',NB,NB,NB,One,VNAT(1+INV),NB,
     &             SANGTI,NB,One,SANGTI2,NB)

        END IF

c Sum over the trace
        DO I = 1,NB
          IJ = I+(I-1)*NB-1
          SUM  = SUM  + OCC(I+INV2) * SANGTR2(1+IJ)
          SUMI = SUMI + OCC(I+INV2) * SANGTI2(1+IJ)
        END DO

        END IF

        INV=INV+NB**2
        INV2=INV2+NB

      END DO

        WRITE(u6,*) "Ben P TEST for JA:"
        WRITE(u6,*) "REAL: ",SUM
        WRITE(u6,*) "IMAG: ",SUMI

        CALL mma_deallocate(SANGF)
        CALL mma_deallocate(SANGTR)
        CALL mma_deallocate(SANGTI)
        CALL mma_deallocate(SANGTR2)
        CALL mma_deallocate(SANGTI2)
      END IF ! IPGLOB >= 4

      CALL mma_deallocate(SANG)

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
     &        WRITE(u6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(u6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(u6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(u6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNAT, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

C IMAGINARY PART
       IF(ITYPE.LE.2) THEN
         WRITE(KNUM,'(I2.2,A,I2.2,A,A)') ASS,".",BSS,".","I"
       ELSE
         WRITE(KNUM,'(I2.2,A,I2.2,A,A,A,A)')ASS,".",BSS,".",CDIR,".","I"
       END IF

       FNAME=FILEBASE//'.'//TRIM(FNUM)//'.I'
       IF(ITYPE.EQ.1)
     &        WRITE(u6,'(A,A)')' NATURAL ORBITALS FOR ',KNUM
       IF(ITYPE.EQ.2)
     &        WRITE(u6,'(A,A)')' ANTISING NATURAL ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.3)
     &        WRITE(u6,'(A,A)')' NATURAL SPIN ORBITALS FOR  ',KNUM
       IF(ITYPE.EQ.4)
     &        WRITE(u6,'(A,A)')' ANTITRIP NATURAL ORBITALS FOR  ',KNUM

       WRITE(u6,'(A,A)') ' ORBITALS ARE WRITTEN ONTO FILE ',FNAME

        LuxxVec=50
        LuxxVec=isfreeunit(LuxxVec)

        CALL WRVEC(FNAME,LUXXVEC,'CO',NSYM,NBASF,NBASF,
     &             VNATI, OCC, Dummy, iDummy,
     &     '* DENSITY FOR PROPERTY TYPE ' // CHARTYPE // KNUM )

c       Test a few values
C        CALL ADD_INFO("SONATORB_CPLOTR", VNAT, 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTI", VNATI, 1, 4)
C        CALL ADD_INFO("SONATORB_CPLOTO", OCC, 1, 4)

      END DO

      CALL mma_deallocate(DMAT)
      CALL mma_deallocate(DMATI)
      CALL mma_deallocate(VEC)
      CALL mma_deallocate(VEC2)
      CALL mma_deallocate(VEC2I)
      CALL mma_deallocate(SCR)
      CALL mma_deallocate(SCRI)
      CALL mma_deallocate(VNAT)
      CALL mma_deallocate(VNATI)
      CALL mma_deallocate(OCC)

      END SUBROUTINE SONATORB_CPLOT




      SUBROUTINE CPLOT_DIAG(MATR, MATI, DIM, EIGVECR, EIGVECI)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero
      IMPLICIT NONE
      INTEGER(KIND=IWP), INTENT(IN):: DIM
      REAL(KIND=WP), INTENT(INOUT):: MATR(DIM*(DIM+1)/2),
     &                               MATI(DIM*(DIM+1)/2)
      REAL(KIND=WP), INTENT(OUT):: EIGVECR(DIM,DIM),EIGVECI(DIM,DIM)

      REAL(KIND=WP) CEIGVAL(DIM)
      COMPLEX(KIND=WP) MATFULL((DIM*(DIM+1)/2))
      COMPLEX(KIND=WP) CEIGVEC(DIM,DIM)
      COMPLEX(KIND=WP) ZWORK(2*DIM-1)
      REAL(KIND=WP) RWORK(3*DIM-2)
      INTEGER(KIND=IWP) INFO, I, J

      DO J=1,(DIM*(DIM+1)/2)
          MATFULL(J) = CMPLX(MATR(J),MATI(J),kind=WP)
      END DO


      call zhpev_('V','U',DIM,MATFULL,CEIGVAL,
     &           CEIGVEC,DIM,ZWORK,RWORK,INFO)


      IF(INFO.NE.0) THEN
          WRITE(u6,*) "Error in diagonalization"
          WRITE(u6,*) "INFO: ",INFO
          CALL ABEND()
      END IF

      DO I=1,DIM
      DO J=1,DIM
          EIGVECR(I,J) = REAL(CEIGVEC(I,J))
          EIGVECI(I,J) = AIMAG(CEIGVEC(I,J))
      END DO
      END DO

      CALL DCOPY_(DIM*(DIM+1)/2,[Zero],0,MATR,1)
      CALL DCOPY_(DIM*(DIM+1)/2,[Zero],0,MATI,1)

      DO J=1,DIM
         MATR((J*(J-1)/2)+J) = CEIGVAL(J)
      END DO

      END SUBROUTINE CPLOT_DIAG

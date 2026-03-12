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
      SUBROUTINE NEWSCTAB(MINOP,MAXOP,MLTPL,MS2,ICASE)
      use definitions, only: iwp, u6
      use stdalloc, only: mma_allocate
      use rassi_global_arrays, only: TRANS1, TRANS2, TRANS
      use rassi_global_arrays, only: SPNTAB1, SPNTAB2, SPNTAB
      IMPLICIT NONE
      INTEGER(KIND=IWP), INTENT(IN):: MINOP,MAXOP,MLTPL,MS2, ICASE

      INTEGER(KIND=IWP), PARAMETER :: ASPIN=1, BSPIN=0
      INTEGER(KIND=IWP), PARAMETER :: UPCPL=1,DWNCPL=0, NULLPTR=-1
      INTEGER(KIND=IWP) IBLK,IOPEN,KSPCPL,KSPDET,LTRANS,LTRANS0,NA,NBLK,&
     &                  NCP,ND,NSPCPL,NSPDET,NTAB,NTRANS,NB
      INTEGER(KIND=IWP), EXTERNAL:: NGENE,NOVERM

      IF(MLTPL-MS2.LT.1 .OR. MLTPL+MS2.LT.1) THEN
      WRITE(u6,*)'NewSCTab: Contradictory values of MLTPL vs. MS2.'
      WRITE(u6,*)                                                       &
     &        'The function was invoked with the following arguments:'
      WRITE(u6,'(1X,A,I9)')' MINOP:',MINOP
      WRITE(u6,'(1X,A,I9)')' MAXOP:',MAXOP
      WRITE(u6,'(1X,A,I9)')' MLTPL:',MLTPL
      WRITE(u6,'(1X,A,I9)')' MS2  :',MS2
      CALL ABEND()
      END IF

! Run through construction loop twice. First get size of table
! and total nr of transformation coefficients:
      NSPCPL=0
      NSPDET=0
      NTRANS=0
      IBLK=0
      DO IOPEN=MINOP,MAXOP
        IBLK=IBLK+1
        NCP=NGENE(IOPEN,MLTPL)
        IF(NCP.EQ.0) CYCLE
        NA=(IOPEN+MS2)/2
        ND=NOVERM(IOPEN,NA)
        NSPCPL=NSPCPL+IOPEN*NCP
        NSPDET=NSPDET+IOPEN*ND
        NTRANS=NTRANS+ND*NCP
      END DO
      NBLK=IBLK
      NTAB=8+6*NBLK+NSPCPL+NSPDET
! The spincoupling table will be NTAB integers long, and there
! will be NTRANS spin-coupling coefficients (real*8).
! The table consists of 8 header words, then an array (6,NBLK)
! with pointers and sizes to the spin coupling, spin determinant
! and spin-coupling coefficient arrays, and finally the
! spin coupling and spin determinant arrays themselves.
! The transformation coefficients are real*8 data and stored
! in a separate array.
      SELECT CASE(ICASE)
      CASE (1)
         CALL mma_allocate(SPNTAB1,NTAB,Label='SPNTAB1')
         SPNTAB=>SPNTAB1(:)
         CALL mma_allocate(TRANS1,NTRANS,Label='TRANS1')
         TRANS=>TRANS1(:)
      CASE (2)
         CALL mma_allocate(SPNTAB2,NTAB,Label='SPNTAB2')
         SPNTAB=>SPNTAB2(:)
         CALL mma_allocate(TRANS2,NTRANS,Label='TRANS2')
         TRANS=>TRANS2(:)
      CASE DEFAULT
      END SELECT
      KSPCPL=9+6*NBLK
      KSPDET=KSPCPL+NSPCPL
! Table size
      SPNTAB(1)=NTAB
! Table type identifier
      SPNTAB(2)=47
! Spin multiplicity
      SPNTAB(3)=MLTPL
! Spin projection
      SPNTAB(4)=MS2
! Min and max nr of open shells
      SPNTAB(5)=MINOP
      SPNTAB(6)=MAXOP
! Associated workspace array for Re*8 data (transf matrices)
      SPNTAB(7)=-1 ! not used
      SPNTAB(8)=NTRANS
! Individual information for each separate nr of open shells:
      NTAB=6
      NTRANS=0
      IBLK=0
      LTRANS0=1
      LTRANS=1
      DO IOPEN=MINOP,MAXOP
        IBLK=IBLK+1
        NCP=NGENE(IOPEN,MLTPL)
        NA=(IOPEN+MS2)/2
        NB= IOPEN-NA
        IF(NCP.GT.0) THEN
          ND=NOVERM(IOPEN,NA)
          SPNTAB(9+(IBLK-1)*6)=IOPEN
          SPNTAB(10+(IBLK-1)*6)=NCP
          SPNTAB(11+(IBLK-1)*6)=ND
! Compute spin couplings:
          CALL PROTOCSF(IOPEN,MLTPL,NCP,SPNTAB(KSPCPL:))
          SPNTAB(12+(IBLK-1)*6)=KSPCPL
! Compute spin determinants:
          CALL PROTOSD(NA,NB,ND,SPNTAB(KSPDET:))
          SPNTAB(13+(IBLK-1)*6)=KSPDET
! Compute spin coupling coefficients:
          CALL PROTOT(IOPEN,ND,SPNTAB(KSPDET:),NCP,                     &
     &                      SPNTAB(KSPCPL:),TRANS(LTRANS:))
          SPNTAB(14+(IBLK-1)*6)=LTRANS-LTRANS0+1
          KSPCPL=KSPCPL+IOPEN*NCP
          KSPDET=KSPDET+IOPEN*ND
          LTRANS=LTRANS+ND*NCP
        ELSE
          SPNTAB(9+(IBLK-1)*6)=IOPEN
          SPNTAB(10+(IBLK-1)*6)=0
          SPNTAB(11+(IBLK-1)*6)=0
          SPNTAB(12+(IBLK-1)*6)=NULLPTR
          SPNTAB(13+(IBLK-1)*6)=NULLPTR
          SPNTAB(14+(IBLK-1)*6)=NULLPTR
        END IF
      END DO

      nullify(TRANS,SPNTAB)

      END SUBROUTINE NEWSCTAB

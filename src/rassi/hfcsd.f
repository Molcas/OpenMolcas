************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
      SUBROUTINE HFCSD(LABEL,IC,BUFF,NBUFF,NSIZ,ISCHK)
      use definitions, only: iwp, wp, u6
      use constants, only: Zero, Two, Three, Four
      use stdalloc, only: mma_allocate, mma_deallocate
      use hfc_logical, only: MAG_X2C
      IMPLICIT NONE
************************************************************************
*     Objective: to compute the 'spin-dependent' part of the hyperfine *
*                from the magnetic integrals contributions             *
*     Output: BUFF                                                     *
************************************************************************
      CHARACTER(LEN=8), intent(in):: LABEL
      INTEGER(KIND=IWP), INTENT(IN):: IC,NBUFF,NSIZ
      INTEGER(KIND=IWP), INTENT(INOUT):: ISCHK
      REAL(KIND=WP), INTENT(INOUT):: BUFF(NBUFF)

      INTEGER(KIND=IWP) ICM,INBUFF,IOPT,IRC
      real(kind=wp) DA
      real(kind=wp), allocatable:: TA(:)

c Set MAG_X2C to avoid add_info in hfcts
      MAG_X2C=.True.
      IOPT=0
      CALL mma_allocate(TA,NBUFF,Label='TA')
c BUFF needs to be initialized
      do iNBUFF = 1, NBUFF
         BUFF(iNBUFF) = Zero
      enddo
c end of initialization
      DA = Two

      Select Case (IC)
      Case(1)
c EF2(1) = (2*MAG(1)-MAG(5)-MAG(9))*(2/3)
         ICM = 1
         DA = Four/Three
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo

         ICM = 5
         DA = -Two/Three
         Call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)

         ICM = 9
         DA = -Two/Three
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
      Case(2)
c EF2(2) = MAG(2)*2
         ICM = 2
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo
      Case(3)
c EF2(3) = MAG(3)*2
         ICM = 3
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo
      Case(4)
c EF2(4) = (2*MAG(5)-MAG(1)-MAG(9))*(2/3)
         ICM = 5
         DA = Four/Three
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo

         ICM = 1
         DA = -Two/Three
         Call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)

         ICM = 9
         DA = -Two/Three
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
      Case(5)
c EF2(5) = MAG(6)*2
         ICM = 6
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo
      Case(6)
c EF2(6) = (MAG(1)+MAG(5)+MAG(9))*2
         ICM = 1
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = TA(1+iNBUFF)
         enddo

         ICM = 5
         Call RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)

         ICM = 9
         CALL RDONE(IRC,IOPT,LABEL,ICM,TA,ISCHK)
         IF (IRC.NE.0) CALL ErrStop()
         CALL DAXPY_(NSIZ,DA,TA,1,BUFF,1)
      Case Default
         WRITE(u6,'(6X,A)')'*** ERROR IN SUBROUTINE HFCSD ***'
         Call ABend()
      End Select
      CALL mma_deallocate(TA)

      Contains
      Subroutine ErrStop()
        WRITE(u6,*)
        WRITE(u6,'(6X,A)')'*** ERROR IN SUBROUTINE HFCSD ***'
        WRITE(u6,'(6X,A)')'  FAILED IN READING FROM  ONEINT'
        WRITE(u6,'(6X,A)')' PLEASE MAKE SURE THE MAGNETIC'
        WRITE(u6,'(6X,A)')' HYPERFINE INTEGRALS ARE AVAILABLE'
        WRITE(u6,'(6X,A,A)')'  LABEL     = ',LABEL
        WRITE(u6,'(6X,A,I2)')'  COMPONENT = ',ICM
        WRITE(u6,*)
        CALL ABEND()
      END Subroutine ErrSTop

      END SUBROUTINE HFCSD

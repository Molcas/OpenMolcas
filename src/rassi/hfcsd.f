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
      IMPLICIT REAL*8 (A-H,O-Z)
************************************************************************
*     Objective: to compute the 'spin-dependent' part of the hyperfine *
*                from the magnetic integrals contributions             *
*     Output: BUFF                                                     *
************************************************************************
#include "WrkSpc.fh"
#include "hfc_logical.fh"
      CHARACTER*8 LABEL
      INTEGER IC,NBUFF,NSIZ,ISCHK
      REAL*8 BUFF(NBUFF)
      INTEGER ICM,INBUFF
      REAL*8 DA

c Set MAG_X2C to avoid add_info in hfcts
      MAG_X2C=.True.
      IOPT=0
      CALL GETMEM('MAG','Allo','Real',ITA,NBUFF)
c BUFF needs to be initialized
      do iNBUFF = 1, NBUFF
         BUFF(iNBUFF) = 0.0
      enddo
c end of initialization
      DA = 2.0D0
      If (IC.eq.1) then
c EF2(1) = (2*MAG(1)-MAG(5)-MAG(9))*(2/3)
         ICM = 1
         DA = 4.0D0/3.0D0
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
         ICM = 5
         DA = -2.0D0/3.0D0
         Call RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         ICM = 9
         DA = -2.0D0/3.0D0
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
      endif
      if (IC.eq.2) then
c EF2(2) = MAG(2)*2
         ICM = 2
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
      endif
      if (IC.eq.3) then
c EF2(3) = MAG(3)*2
         ICM = 3
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
      endif
      if (IC.eq.4) then
c EF2(4) = (2*MAG(5)-MAG(1)-MAG(9))*(2/3)
         ICM = 5
         DA = 4.0D0/3.0D0
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
         ICM = 1
         DA = -2.0D0/3.0D0
         Call RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         ICM = 9
         DA = -2.0D0/3.0D0
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
      endif
      if (IC.eq.5) then
c EF2(5) = MAG(6)*2
         ICM = 6
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
      endif
      if (IC.eq.6) then
c EF2(6) = (MAG(1)+MAG(5)+MAG(9))*2
         ICM = 1
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         do iNBUFF = NSIZ, NBUFF-1
             BUFF(iNBUFF+1) = WORK(ITA+iNBUFF)
         enddo
         ICM = 5
         Call RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
         ICM = 9
         CALL RDONE(IRC,IOPT,LABEL,ICM,WORK(ITA),ISCHK)
         IF (IRC.NE.0) GOTO 300
         CALL DAXPY_(NSIZ,DA,WORK(ITA),1,BUFF,1)
      endif
      CALL GETMEM('MAG','Free','Real',ITA,NBUFF)

300   CONTINUE
      IF (IRC.NE.0) THEN
        WRITE(6,*)
        WRITE(6,'(6X,A)')'*** ERROR IN SUBROUTINE HFCSD ***'
        WRITE(6,'(6X,A)')'  FAILED IN READING FROM  ONEINT'
        WRITE(6,'(6X,A)')' PLEASE MAKE SURE THE MAGNETIC'
        WRITE(6,'(6X,A)')' HYPERFINE INTEGRALS ARE AVAILABLE'
        WRITE(6,'(6X,A,A)')'  LABEL     = ',LABEL
        WRITE(6,'(6X,A,I2)')'  COMPONENT = ',ICM
        WRITE(6,*)
        CALL ABEND
       ENDIF

      END

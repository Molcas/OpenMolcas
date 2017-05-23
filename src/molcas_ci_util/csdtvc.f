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
* Copyright (C) 1999, Markus P. Fuelscher                              *
*               2012, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,
     *                  IREFSM,ICOPY)
C
C     PURPOSE: TRANDSORM FRON DETERMINANT TO CSF BASIS AND IVCE VERSA
C              IWAY = 1 : CSF TO DETERMINANT TRANSFORMATION
C              IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "ciinfo.fh"
#include "spinfo.fh"
C
      DIMENSION CSFVEC(*),DETVEC(*)
      DIMENSION DTOCMT(*),ICTSDT(*)

      Call qEnter('CSDTVC')
C
C     To avoid compiler complains
C
      IOFFCS = 0
      IOFFDT = 0
      IOFFCD = 0
C
      NTEST = 00
C
      NDET = NDTASM(IREFSM)
      NCSF = NCSASM(IREFSM)
      IF (NTEST .GE. 100) THEN
         WRITE(6,*) '======================================='
         WRITE(6,*) '         CSDTVC INFORMATIONS'
         IF (IWAY .EQ. 1) THEN
            WRITE(6,*) '   CSF TO DETERMINANT TRANSFORMATION'
         ELSE
            WRITE(6,*) '   DETERMINANT TO CSF TRANSFORMATION'
         END IF
         WRITE(6,*) '======================================='
         WRITE(6,*)
         WRITE(6,*) '  NDET                  = ',NDET
         WRITE(6,*) '  NCSF                  = ',NCSF
         WRITE(6,*)
      END IF
C
C     CSF ==> DET TRANSFORMATION
C
      IF( IWAY.EQ.1 ) THEN
        IF (NTEST .GE. 100) THEN
          WRITE(6,*) '   INPUT CSF VECTOR:'
          CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
          WRITE(6,*)
        END IF
        CALL DCOPY_(NDET,0.0d0,0,DETVEC,1)
        DO ITYP = 1,NTYP
          IDET = NDTFTP(ITYP)
          ICSF = NCSFTP(ITYP)
          ICNF = NCNFTP(ITYP,IREFSM)
          IF( ITYP.EQ.1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
            IOFFCD = IOFFCD + NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
          END IF
          IF( (IDET*ICNF*ICSF).GT.0 )
     &         CALL MATML4(DETVEC(IOFFDT),DTOCMT(IOFFCD),
     &         CSFVEC(IOFFCS), IDET,ICNF,IDET,ICSF,ICSF,ICNF,0)
        END DO
        Call Sort_Cdet(nDet,ICTSDT,DetVec)
        IF( ICOPY.NE.0 ) CALL DCOPY_(NDET,DETVEC,1,CSFVEC,1)
        IF (NTEST .GE. 100) THEN
          WRITE(6,*) '   OUTPUT DET VECTOR:'
          CALL WRTMAT(DETVEC,1,NDET,1,NDET)
          WRITE(6,*)
        END IF

      ELSE
C
C     DET ==> CSF TRANSFORMATION
C
        IF (NTEST .GE. 100) THEN
          WRITE(6,*) '   INPUT DET VECTOR:'
          CALL WRTMAT(DETVEC,1,NDET,1,NDET)
          WRITE(6,*)
        END IF
        CALL GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
        IF (NTEST .GE. 100) THEN
           WRITE(6,*) ' ICTSDT reorder array '
           CALL IWRTMA(ICTSDT,1,100,1,100)
        END IF
        CALL DCOPY_(NDET,CSFVEC,1,DETVEC,1)
        DO ITYP = 1,NTYP
          IDET = NDTFTP(ITYP)
          ICSF = NCSFTP(ITYP)
          ICNF = NCNFTP(ITYP,IREFSM)
          IF( ITYP.EQ. 1 ) THEN
            IOFFCS = 1
            IOFFDT = 1
            IOFFCD = 1
          ELSE
            IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
            IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
            IOFFCD = IOFFCD + NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
          END IF
          IF( (IDET*ICNF*ICSF).GT.0 )
     &         CALL MATML4(CSFVEC(IOFFCS),DTOCMT(IOFFCD),
     &         DETVEC(IOFFDT), ICSF,ICNF,IDET,ICSF,IDET,ICNF,1)
        END DO
        IF( ICOPY.NE.0 ) CALL DCOPY_(NCSF,CSFVEC,1,DETVEC,1)
        IF (NTEST .GE. 100) THEN
          WRITE(6,*) '   OUTPUT CSF VECTOR:'
          CALL WRTMAT(CSFVEC,1,NCSF,1,NCSF)
          WRITE(6,*)
        END IF
       END IF
C
      Call qExit('CSDTVC')
      RETURN
      END

      Subroutine Sort_Cdet(N,Index,X)
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     Sort CI-coefficients in determinant basis and change phase       *
*                                                                      *
*     calling arguments:                                               *
*     N       : integer                                                *
*               dimension of the CI-vector                             *
*     Index   : integer                                                *
*               reordering indices                                     *
*     X       : real*8                                                 *
*               CI-vector                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     G. Li Manni on 03 Feb 2012 on a desperate situation for saving   *
*     time!                                                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1999                                 *
*                                                                      *
* It was a great piece of code for those times when memory was a very  *
* big issue! Nowdays only time is important!                           "
************************************************************************
*      Do i = 1,N                                                      *
*        i_old = i                                                     *
*        i_new = ABS(Index(i_old))                                     *
*        Do while ( i_new.gt.i )                                       *
*          i_old = i_new                                               *
*          i_new = ABS(Index(i_old))                                   *
*        End Do                                                        *
*        If ( i_new.eq.i ) then                                        *
*          i_old = i                                                   *
*          X_old = X(i_old)                                            *
*          i_new = ABS(Index(i_old))                                   *
*          X_new = X(i_new)                                            *
*          alpha = DBLE(SIGN(1,Index(i_old)))                          *
*          Do while ( i_new.gt.i )                                     *
*            X(i_new) = alpha*X_old                                    *
*            i_old = i_new                                             *
*            X_old = X_new                                             *
*            i_new = ABS(Index(i_old))                                 *
*            X_new = X(i_new)                                          *
*            alpha = DBLE(SIGN(1,Index(i_old)))                        *
*          End Do                                                      *
*          X(i) = alpha*X_old                                          *
*        End If                                                        *
*      End Do                                                          *
************************************************************************

      Implicit Real*8 (A-H,O-Z)

      Dimension Index(N),X(N)
      Intrinsic SIGN
#include "WrkSpc.fh"
      call getmem('ReordSDs','Allo','Real',iReoSDs,N)
      Do i = 1,N
        i_new = ABS(Index(i))
        alpha = DBLE(SIGN(1,Index(i)))
        Work(iReoSDs+i_new -1) = alpha*X(i)
      End do
      call dcopy_(N,Work(iReoSDs),1,X,1)
      call getmem('ReordSDs','Free','Real',iReoSDs,N)
      Return
      End

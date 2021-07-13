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
      SUBROUTINE CSDIAG_CI_UTIL(CSFDIA,DETDIA,NCNFTP,NTYP,              &
     &                  ICTSDT,NDTFTP,NCSFTP,IPRINT)
!
!     PURPOSE: OBTAIN AVERAGE CI DIAGONAL ELEMENTS AND STORE IN CSFDIA
!
!     CALLS TO SUBROUTINES AND EXTERNAL FUNCTIONS:
!     DCOPY,WRTMAT
!
!     CALLING PARAMETERS.
!     CSFDIA  : CI DIAGONAL IN SCF BASIS
!     DETDIA  : CI DIAGONAL IN DETERMINENT BASIS
!     NCNFTP  :
!     NTYP    :
!     ICTSDT  :
!     NDTFTP  :
!     NCSFTP  :
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      DIMENSION CSFDIA(*),DETDIA(*)
      DIMENSION NCNFTP(NTYP),NDTFTP(NTYP),NCSFTP(NTYP)
      DIMENSION ICTSDT(*)
!
      ICSOFF = 1
      IDTOFF = 1
      JCNABS = 0
      DO 100 ITYP = 1, NTYP
        IDET = NDTFTP(ITYP)
        ICSF = NCSFTP(ITYP)
        ICNF = NCNFTP(ITYP)
        DO 80 JCNF = 1, ICNF
          JCNABS = JCNABS + 1
          EAVER = 0.0D0
          DO 70 JDET = 1, IDET
            EAVER = EAVER +DETDIA(ABS(ICTSDT(IDTOFF-1+JDET)) )
70        CONTINUE
          IF( IDET .NE. 0 )EAVER = EAVER/DBLE(IDET)
          CALL DCOPY_(ICSF,[EAVER],0,CSFDIA(ICSOFF),1)
          ICSOFF = ICSOFF + ICSF
          IDTOFF = IDTOFF + IDET
80      CONTINUE
100   CONTINUE
!
      IF( IPRINT.GE.40 ) THEN
        NCSTOT = ICSOFF-1
        NDTTOT = IDTOFF-1
        WRITE(6,*) ' '
        WRITE(6,*) ' CIDIAGONAL IN DET BASIS '
        CALL WRTMAT(DETDIA,1,NDTTOT,1,NDTTOT)
        WRITE(6,*) ' '
        WRITE(6,*) ' CIDIAGONAL IN CSF BASIS '
        CALL WRTMAT(CSFDIA,1,NCSTOT,1,NCSTOT)
      END IF
!
      RETURN
      END

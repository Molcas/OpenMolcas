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
FUNCTION fmm_second()

   USE fmm_global_paras
   IMPLICIT NONE
   REAL(REALK) :: fmm_second
   INTEGER :: count, count_rate, count_max, tmp
   INTEGER(INTK), SAVE :: parity=0
   INTEGER(INTK), SAVE :: last_count=0

! VV: order changed due to a bug in NAG compiler.
   CALL SYSTEM_CLOCK(count,count_rate,count_max)

   IF (parity == 0) parity = 1  ! first call to fmm_second()
   parity = -parity
   ! must check if count has "cycled"
   tmp = count
   IF ( parity == 1 ) THEN  ! second of a paired call to fmm_second()
      IF ( count < last_count ) tmp = count + (count_max - last_count)
   END IF
   fmm_second = REAL(tmp,REALK)/REAL(count_rate,REALK)
   last_count = count

END FUNCTION fmm_second

!-------------------------------------------------------------------------------

SUBROUTINE TIMTXT(TEXTIN,TIMUSD,IUNIT)

! TIMTXT based on TIMER by TUH //900709-hjaaj

   USE fmm_global_paras
   IMPLICIT NONE
   CHARACTER(LEN=*) :: TEXTIN
   CHARACTER :: AHOUR*6, ASEC*8, AMIN*8
   REAL(REALK) :: TIMUSD
   INTEGER :: ISECND, IUNIT, IHOURS, MINUTE
   CHARACTER(LEN=45) :: TEXT

   TEXT = TEXTIN

   ISECND = INT(TIMUSD)
   IF (ISECND .GE. 60) THEN
      MINUTE = ISECND/60
      IHOURS = MINUTE/60
      MINUTE = MINUTE - 60*IHOURS
      ISECND = ISECND - 3600*IHOURS - 60*MINUTE
      IF (IHOURS .EQ. 1) THEN
         AHOUR = ' hour '
      ELSE
         AHOUR = ' hours'
      END IF
      IF (MINUTE .EQ. 1) THEN
         AMIN = ' minute '
      ELSE
         AMIN = ' minutes'
      END IF
      IF (ISECND .EQ. 1) THEN
         ASEC = ' second '
      ELSE
         ASEC = ' seconds'
      END IF
      IF (IHOURS .GT. 0) THEN
         WRITE(IUNIT,100) TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
      ELSE
         WRITE(IUNIT,200) TEXT, MINUTE, AMIN, ISECND, ASEC
      END IF
   ELSE
      WRITE(IUNIT,300) TEXT,TIMUSD
   END IF
100 FORMAT(1X,A,I4,A,I3,A,I3,A)
200 FORMAT(1X,A,     I3,A,I3,A)
300 FORMAT(1X,A,F7.2,' seconds')
   RETURN
END SUBROUTINE TIMTXT

!-------------------------------------------------------------------------------

SUBROUTINE fmm_quit(msg)
   CHARACTER(LEN=*) msg
   write(6,*) msg
   WRITE(6,*) ">>> FATAL ERROR"
   CALL Abend()
END SUBROUTINE fmm_quit

!------------------------------------------------------------------------------

SUBROUTINE fmm_matrix_norm(label,matrix,ndim)

   USE fmm_global_paras

   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN) :: label
   INTEGER,      INTENT(IN) :: ndim
   REAL(REALK),  INTENT(IN) :: matrix(ndim)

   INTEGER :: i
   REAL(REALK) :: norm

   norm = 0d0
   DO i = 1, ndim
      norm = norm + matrix(i)*matrix(i)
   END DO
   WRITE(6,*) 'o fmm_matrix_norm: ', label,' = ', SQRT(norm)

END SUBROUTINE fmm_matrix_norm

!------------------------------------------------------------------------------

MODULE fmm_utils

   USE fmm_global_paras
   IMPLICIT NONE
   PRIVATE
   ! Public procedures
   PUBLIC :: fmm_distance_between_lines

CONTAINS

!-------------------------------------------------------------------------------

   FUNCTION fmm_cross_product(a,b)

      IMPLICIT NONE
      REAL(REALK), INTENT(IN)  :: a(3), b(3)
      REAL(REALK) :: fmm_cross_product(3)

      fmm_cross_product(1) =    a(2)*b(3) - a(3)*b(2)
      fmm_cross_product(2) = -( a(1)*b(3) - a(3)*b(1) )
      fmm_cross_product(3) =    a(1)*b(2) - a(2)*b(1)

   END FUNCTION fmm_cross_product

!-------------------------------------------------------------------------------

   FUNCTION fmm_distance_between_lines(p,q,ep,eq)

      IMPLICIT NONE
      REAL(REALK), INTENT(IN)  :: p(3), q(3), ep(3), eq(3)
      REAL(REALK) :: fmm_distance_between_lines

      REAL(REALK) :: d1, d2, n(3)

      n = fmm_cross_product(ep,eq)

      d1 = DOT_PRODUCT(p,n)
      d2 = DOT_PRODUCT(q,n)

      fmm_distance_between_lines = ABS(d1 - d2)

   END FUNCTION fmm_distance_between_lines

!-------------------------------------------------------------------------------

END MODULE fmm_utils


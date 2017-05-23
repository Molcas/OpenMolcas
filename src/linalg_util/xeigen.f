************************************************************************
* This file from Cascade:                                              *
*   http://www.netlib.org/ieeecss/cascade/                             *
*                                                                      *
* No warranty is made concerning the quality of the cascade software   *
* or its applicability.  The software available as a part of netlib    *
* has been placed in the public domain by the authors, as a part of    *
* the cascade development project, sponsored by the U.S.  Department   *
* of Energy.  Several authors  contributed to this project; they are   *
* cited in the manuals.                                                *
************************************************************************

C     SUBROUTINE EIGEN (NVEC,NA,N,A,EVR,EVI,VECS,SCR1,SCR2,IERR)
      SUBROUTINE XEIGEN (NVEC,NA,N,A,EVR,EVI,VECS,SCR1,SCR2,IERR)
      INTEGER NVEC,NA,N,IERR
      REAL*8 A(NA,N),EVR(N),EVI(N),VECS(NA,N),SCR1(N),SCR2(N)
C
C     ***** PURPOSE:
C     THIS SUBROUTINE COMPUTES THE EIGENVALUES AND EIGENVECTORS
C     (IF DESIRED) OF A REAL GENERAL MATRIX A BY THE DOUBLE FRANCIS
C     QR ALGORITHM AS IMPLEMENTED IN EISPACK.
C     REFERENCE:    SMITH, B.T., ET. AL., MATRIX EIGENSYSTEM ROUTINES--
C                   EISPACK GUIDE, SECOND EDITION, LECTURE NOTES IN
C                   COMPUTER SCIENCE, VOL. 6, SPRINGER-VERLAG, 1976.
C
C     ON ENTRY:
C
C        NVEC   INTEGER
C               SET = 0 IF NO EIGENVECTORS ARE DESIRED, I.E., TO
C               COMPUTE EIGENVALUES ONLY; OTHERWISE SET TO ANY
C               NONZERO INTEGER IF BOTH EIGENVALUES AND EIGENVECTORS
C               ARE DESIRED.
C
C        NA     INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING A AND VECS
C               AS DECLARED IN THE MAIN CALLING PROGRAM.
C
C        N      INTEGER
C               THE ORDER OF THE MATRIX A.
C
C        A      real*8(NA,N)
C               A REAL GENERAL MATRIX WHOSE EIGENVALUES AND EIGEN-
C               VECTORS (IF DESIRED) ARE TO BE COMPUTED.
C
C     ON RETURN:
C
C        EVR    real*8(N)
C               THE REAL PARTS OF THE EIGENVALUES OF A.
C
C        EVI    real*8(N)
C               THE CORRESPONDING IMAGINARY PARTS OF THE EIGENVALUES
C               OF A.  NOTE THAT COMPLEX CONJUGATE PAIRS OF EIGENVALUES
C               APPEAR CONSECUTIVELY WITH THE EIGENVALUE HAVING THE
C               POSITIVE IMAGINARY PART FIRST.
C
C        VECS   real*8(NA,N)
C               IF NVEC IS NONZERO, THIS ARRAY CONTAINS THE REAL AND
C               IMAGINARY PARTS OF THE EIGENVECTORS OF A.  IF THE J-TH
C               EIGENVALUE IS REAL, THE J-TH COLUMN OF VECS CONTAINS
C               THE CORRESPONDING EIGENVECTOR (NORMALIZED TO HAVE
C               EUCLIDEAN OR 2- NORM = 1 AND POSITIVE MAXIMUM COMP-
C               ONENT).  IF THE J-THE EIGENVALUE IS COMPLEX WITH
C               POSITIVE IMAGINARY PART, THE J-TH AND (J+1)-TH
C               COLUMNS OF VECS CONTAIN THE REAL AND IMAGINARY
C               PARTS OF THE CORRESPONDING COMPLEX EIGENVECTOR
C               (NORMALIZED TO HAVE COMPLEX EUCLIDEAN OR 2- NORM
C               =1 AND REAL, POSITIVE MAXIMUM COMPONENT).  THE CONJ-
C               UGATE OF THIS VECTOR IS THE EIGENVECTOR FOR THE
C               CONJUGATE EIGENVALUE.
C
C        SCR1   real*8(N)
C               THE I-TH COMPONENT OF THIS VECTOR CONTAINS THE
C               UNDAMPED NATURAL FREQUENCY (MODULUS) OF THE I-TH
C               EIGENVALUE;  SCR1 IS ALSO USED INTERNALLY
C               AS A SCRATCH VECTOR FOR THE EISPACK SUBROUTINE
C               BALANC.
C
C        SCR2   real*8(N)
C               THE I-TH COMPONENT OF THIS VECTOR CONTAINS THE
C               DAMPING RATIO OF THE I-TH EIGENVALUE;  SCR2 IS ALSO
C               USED INTERNALLY AS A SCRATCH VECTOR FOR THE
C               EISPACK SUBROUTINE ORTHES.
C
C        IERR   INTEGER
C               ERROR COMPLETION CODE RETURNED BY EISPACK SUBROUTINE
C               HQR OR HQR2.  NORMAL RETURN VALUE IS ZERO.  SEE THE
C               EISPACK GUIDE, P. 331, FOR A DISCUSSION OF NONZERO
C               VALUES OF IERR.
C
C     PROGRAM WRITTEN BY ALAN J. LAUB, DEPT. OF ELEC. AND COMP.ENGRG.,
C     UNIVERSITY OF CALIFORNIA, SANTA BARBARA, CA 93106,
C     PH.: (805) 961-3616.
C     JUNE 1981.
C     MOST RECENT MODIFICATION: JAN. 2, 1985
C
C     INTERNAL VARIABLES:
C
      INTEGER I,IGH,J,JM1,K,LOW
      REAL*8 ANORM,EI,EPS,EPSP1,ER,T,TIM,TRE,T1,T2
C
C     FORTRAN FUNCTIONS CALLED:
C
*     REAL*8 ABS,SQRT
C
C     SUBROUTINES AND FUNCTIONS CALLED:
C
C     BALANC,BALBAK,HQR,HQR2,ORTHES,ORTRAN (ALL FROM EISPACK)
C
C     ------------------------------------------------------------------
C
C     DETERMINE MACHINE PRECISION
C
      EPS = 1.0D0
   10 CONTINUE
      EPS = EPS/2.0D0
      EPSP1 = EPS+1.0D0
      IF (EPSP1 .GT. 1.0D0) GO TO 10
      EPS = 2.0D0*EPS
C
C     BALANCE A
C
      CALL BALANC (NA,N,A,LOW,IGH,SCR1)
C
C     COMPUTE 1-NORM OF THE BALANCED A
C
      ANORM = 0.0D0
      DO 30 J = 1,N
         T = 0.0D0
         DO 20 I = 1,N
            T = T+abs(A(I,J))
   20    CONTINUE
         IF (T .GT. ANORM) ANORM = T
   30 CONTINUE
C
C     REDUCE A TO UPPER HESSENBERG FORM
C
      CALL ORTHES (NA,N,LOW,IGH,A,SCR2)
      IF (NVEC .NE. 0) GO TO 40
C
C     COMPUTE EIGENVALUES USING QR ALGORITHM
C
      CALL HQR (NA,N,LOW,IGH,A,EVR,EVI,IERR)
      IF (IERR .NE. 0) RETURN
      GO TO 110
   40 CONTINUE
C
C     COMPUTE EIGENVALUES AND EIGENVECTORS USING QR ALGORITHM
C
      CALL ORTRAN (NA,N,LOW,IGH,A,SCR2,VECS)
      CALL HQR2 (NA,N,LOW,IGH,A,EVR,EVI,VECS,IERR)
      IF (IERR .NE. 0) RETURN
      CALL BALBAK (NA,N,LOW,IGH,SCR1,N,VECS)
C
C     NORMALIZE EIGENVECTORS TO HAVE EUCLIDEAN OR 2- NORM EQUAL TO 1
C
      DO 100 J = 1,N
         K = 0 ! dummy initialize
         IF (EVI(J) .NE. 0.0D0) GO TO 70
         T = 0.0D0
         T1 = 0.0D0
         DO 50 I = 1,N
            T2 = VECS(I,J)**2
            IF (T2 .LE. T1) GO TO 45
            K = I
            T1 = T2
   45       CONTINUE
            T = T+T2
   50    CONTINUE
         T = sign(sqrt(T),VECS(K,J))
         DO 60 I = 1,N
            VECS(I,J) = VECS(I,J)/T
   60    CONTINUE
         GO TO 100
   70    CONTINUE
         IF (EVI(J) .GT. 0.0D0) GO TO 100
         JM1 = J-1
         T = 0.0D0
         T1 = 0.0D0
         DO 80 I = 1,N
            T2 = VECS(I,JM1)**2 + VECS(I,J)**2
            IF (T2 .LE. T1) GO TO 75
            K = I
            T1 = T2
   75       CONTINUE
            T = T+T2
   80    CONTINUE
         T = sqrt(T)
         T1 = sqrt(T1)
         DO 90 I = 1,N
            TRE = VECS(I,JM1)*VECS(K,JM1) + VECS(I,J)*VECS(K,J)
            TIM = VECS(I,J)*VECS(K,JM1) - VECS(I,JM1)*VECS(K,J)
            VECS(I,JM1) = (TRE/T1)/T
            VECS(I,J) = (TIM/T1)/T
   90    CONTINUE
  100 CONTINUE
  110 CONTINUE
C
C     COMPUTE NATURAL FREQUENCIES AND DAMPING RATIOS.  SET
C     EIGENVALUES WITH NORM LESS THAN EPS*ANORM TO (0.0D0,0.0D0)
C
      EPS = EPS*ANORM
      DO 130 I = 1,N
         T = abs(EVR(I))+abs(EVI(I))
         IF (T .GT. EPS) GO TO 120
         EVR(I) = 0.0D0
         EVI(I) = 0.0D0
         SCR1(I) = 0.0D0
         SCR2(I) = 1.0D0
         GO TO 130
  120    CONTINUE
         ER = EVR(I)/T
         EI = EVI(I)/T
         SCR1(I) = sqrt(ER**2 + EI**2)
         SCR2(I) = -ER/SCR1(I)
         SCR1(I) = T*SCR1(I)
         IF (abs(EVI(I)) .LT. EPS) EVI(I) = 0.0D0
  130 CONTINUE
      RETURN
      END

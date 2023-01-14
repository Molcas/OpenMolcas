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
!
c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
      SUBROUTINE AF3X3LEV(RDIST,DELTAE,C3val,C6val,C8val,De,ULR)
c=======================================================================
c*** Simplified version of AF3x3potRet which does not return derivatives
      REAL*8  H(3,3),DM1(3,3),DM3(3,3),DM5(3,3),DR(3,3),
     1              DDe(3,3),Q(3,3)
      REAL*8  DEIGM1(1,1),DEIGM3(1,1),DEIGM5(1,1),DEIGR(1,1),
     1        DEIGDe(1,1), EIGVEC(3,1), W(3)
      REAL*8  RDIST,RDIST2,RDIST3,DELTAE,C3val,C6val,C8val,De,ULR,
     1   RET,RETSig,RETPi,Modulus,M1,M3,M5,Z
      INTEGER          I,J,L
      M1= C3val
      M3= C6val
      M5= C8val
      RET= 9.36423830d-4*RDIST
      RETSig= DCOS(RET) + (RET)*DSIN(RET)
      RETPi= RETSig - RET**2 *DCOS(RET)
      RDIST2= RDIST**2
      RDIST3= RDIST*RDIST2
*      WRITE(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)" '
*      WRITE(25,*) 'zone T = "U(r)"'
c  Initialize interaction matrix to 0.d0
      DO  I= 1,3
          H(I,I)=0.0D0
          ENDDO
ccccc Prepare interation matrix  H
      H(1,1)= -(M1*RETSig+ M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
      H(1,2)= -(DSQRT(2.D0))*H(1,1)
      H(2,1)= H(1,2)
      H(1,3)= M1*RETPi/(DSQRT(6.D0)*RDIST3)
      H(3,1)= H(1,3)
      H(2,2)= 2*H(1,1) + DELTAE
      H(2,3)= H(1,3)/DSQRT(2.d0)
      H(3,2)= H(2,3)
      H(3,3)= DELTAE
cccccc Prepare radial derivative of interaction matrix (? is it needed ?)
      DR(1,1)= (3.d0*M1*RETSig + 6.d0*M3/RDIST3
     1                  + 8.D0*M5/(RDIST3*RDIST2))/(3.d0*RDIST3*RDIST)
      DR(1,2)= -DSQRT(2.d0)*DR(1,1)
      DR(2,1)= DR(1,2)
      DR(2,2)= 2.d0*DR(1,1)
      DR(1,3)= -3.d0*H(1,3)/RDIST
      DR(3,1)= DR(1,3)
      DR(2,3)= -3.d0*H(2,3)/RDIST
      DR(3,2)= DR(2,3)
      DR(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C3
      DM1(1,1)= -(RETSig + M1/(2.d0*De*RDIST3))/(3.d0*RDIST3)
      DM1(1,2)= -DSQRT(2.d0)*DM1(1,1)
      DM1(2,1)= DM1(1,2)
      DM1(2,2)= 2.d0*DM1(1,1)
      DM1(1,3)= RETPi/(DSQRT(6.d0)*RDIST3)
      DM1(3,1)= DM1(1,3)
      DM1(2,3)= DM1(1,3)/DSQRT(2.d0)
      DM1(3,2)= DM1(2,3)
      DM1(3,3)= 0.d0
cccccc Partial derivative of interaction matric  H  w.r.t.  C6
      DM3(1,1)= -1.d0/(3.d0*RDIST3**2)
      DM3(1,2)= -SQRT(2.d0)*DM3(1,1)
      DM3(1,3)= 0.D0
      DM3(2,1)= DM3(1,2)
      DM3(2,2)= 2.d0*DM3(1,1)
      DM3(2,3)= 0.D0
      DM3(3,1)= DM3(1,3)
      DM3(3,2)= DM3(2,3)
      DM3(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  C8
      DM5(1,1)= DM3(1,1)/(RDIST2)
      DM5(1,2)= DM3(1,2)/(RDIST2)
      DM5(1,3)= 0.D0
      DM5(2,1)= DM3(1,2)
      DM5(2,2)= DM3(2,2)/(RDIST2)
      DM5(2,3)= 0.D0
      DM5(3,1)= DM5(1,3)
      DM5(3,2)= DM5(2,3)
      DM5(3,3)= 0.D0
cccccc Partial derivative of interaction matric  H  w.r.t.  De
      DDe(1,1)= M1**2/(12.D0*(RDIST3*De)**2)
      DDe(1,2)= -SQRT(2.D0)*DDe(1,1)
      DDe(1,3)= 0.D0
      DDe(2,1)= DDe(1,2)
      DDe(2,2)= 2.D0*DDe(1,1)
      DDe(2,3)= 0.d0
      DDe(3,1)= DDe(1,3)
      DDe(3,2)= DDe(2,3)
      DDe(3,3)= 0.D0
cccccc Call subroutine to prepare and invert interaction matrix  H
      CALL ZHEEVJ3(H,Q,W)
      L=1
ccc Nor - identify the lowest eigenvalue of  H  and label it  L
      DO J=2,3
          IF (W(J) .LT. W(L)) THEN
              L=J
              ENDIF
          ENDDO
      ULR= -W(L)
      DO I=1,3
          EIGVEC(I,1) = Q(I,L)
          ENDDO
      WRITE(6,*) EIGVEC
! Make sure the following variables are "referenced":
      DEIGM1=0.d0
      DEIGM3=DEIGM1*0.d0
      DEIGM5=DEIGM3*0.d0
      DEIGR =DEIGM5*0.d0
      DEIGDe=DEIGR*0.d0
      DEIGM1=DEIGDe*0.d0
c     WRITE(25,600) RDIST ,ULR
c 600 FORMAT(2D16.7)
c     WRITE(26,601) RDIST , DEIGM1, DEIGR ,DEIGDe
c 601 FORMAT(4D16.7)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Modulus = SQRABS(Z)
!     Modulus =  REAL(Z)**2
      Modulus = 0.d0
      Z = Modulus ! Make sure it's "referenced"
      Modulus = Z ! Make sure it's "referenecd"
      RETURN
      CONTAINS
*=======================================================================
      SUBROUTINE ZHEEVJ3(H,Q,W)
*=======================================================================
c** Subroutine to setup and invert the matrix  H  and return
c   eigenvalues W and eigenvector matric  Q
      INTEGER   N, I, X, Y, R
      PARAMETER (N=3)
      REAL*8    H(3,3),Q(3,3), W(3)
      REAL*8    SD,SO,S,T,C,G,B,Z,THRESH
c     DOUBLE PRECISION FUNCTION SQRABS
* Initialize Q to the identitity matrix
* --- This loop can be omitted if only the eigenvalues are desired ---
      DO  X = 1, N
          Q(X,X) = 1.0D0
          DO  Y = 1, X-1
              Q(X, Y) = 0.0D0
              Q(Y, X) = 0.0D0
              ENDDO
          ENDDO
* Initialize W to diag(A)
      DO  X = 1, N
          W(X) =  REAL(H(X, X))
          ENDDO
* Calculate SQR(tr(A))
      SD= 0.0D0
      DO  X = 1, N
          SD= SD + ABS(W(X))
          ENDDO
      SD = SD**2
* Main iteration loop
      DO  I = 1, 50
* Test for convergence
          SO = 0.0D0
          DO  X = 1, N
              DO  Y = X+1, N
                  SO = SO + ABS(REAL(H(X, Y)))
                  ENDDO
              ENDDO
          IF(SO.EQ.0.0D0) RETURN
          IF (I .LT. 4) THEN
              THRESH = 0.2D0 * SO / N**2
            ELSE
              THRESH = 0.0D0
            END IF
* Do sweep
          DO  X= 1, N
              DO  Y= X+1, N
                  G= 100.0D0*(ABS(REAL(H(X, Y))) )
                  IF((I.GT.4).AND.((ABS(W(X))+G).EQ.ABS(W(X)))
     $                         .AND.((ABS(W(Y))+G).EQ.ABS(W(Y)))) THEN
                      H(X, Y)= 0.0D0
                    ELSEIF(ABS(REAL(H(X, Y))).GT.THRESH) THEN
* Calculate Jacobi transformation
                      B= W(Y) - W(X)
                      IF((ABS(B)+G).EQ.ABS(B)) THEN
                          T= H(X, Y) / B
                        ELSE
                          IF(B .LE. 0.0D0) THEN
                              T= -2.0D0 * H(X, Y)
c    $                       /(SQRT(B**2 + 4.0D0*SQRABS(H(X, Y))) - B)
     $                       /(SQRT(B**2 + 4.0D0*REAL(H(X, Y))**2) - B)
                            ELSE IF (B .EQ. 0.0D0) THEN
                              T= H(X, Y) * (1.0D0 / ABS(H(X, Y)))
                            ELSE
                              T= 2.0D0 * H(X, Y)
c    $                       /(SQRT(B**2 + 4.0D0*SQRABS(H(X, Y))) + B)
     $                       /(SQRT(B**2 + 4.0D0*REAL(H(X, Y))**2) + B)
                            ENDIF
                        ENDIF
c                     C= 1.0D0 / SQRT( 1.0D0 + SQRABS(T) )
                      C= 1.0D0 / SQRT( 1.0D0 + REAL(T)**2 )
                      S= T * C
                      Z= REAL(T * (H(X, Y)))
* Apply Jacobi transformation
                      H(X, Y) = 0.0D0
                      W(X)    = W(X) - Z
                      W(Y)    = W(Y) + Z
                      DO  R = 1, X-1
                          T       = H(R, X)
                          H(R, X) = C * T - (S) * H(R, Y)
                          H(R, Y) = S * T + C * H(R, Y)
                          ENDDO
                      DO  R = X+1, Y-1
                          T       = H(X, R)
                          H(X, R) = C * T - S * (H(R, Y))
                          H(R, Y) = S * (T) + C * H(R, Y)
                          ENDDO
                      DO  R = Y+1, N
                          T       = H(X, R)
                          H(X, R) = C * T - S * H(Y, R)
                          H(Y, R) = (S) * T + C * H(Y, R)
                          ENDDO
* eigenvectors
* This loop can be omitted if only the eigenvalues are desired ---
                      DO  R = 1, N
                          T       = Q(R, X)
                          Q(R, X) = C * T - (S) * Q(R, Y)
                          Q(R, Y) = S * T + C * Q(R, Y)
                          ENDDO
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      WRITE(6,*) 'ZHEEVJ3: No convergence.'
      END SUBROUTINE ZHEEVJ3

c     CONTAINS
*=======================================================================
c      DOUBLE PRECISION FUNCTION SQRABS(Z)
c      SUBROUTINE SQRABS(Z)
*=======================================================================
* Calculates the squared absolute value of a complex number Z
* ----------------------------------------------------------------------
*  Parameters ..
c     REAL*8 Z
c     SQRABS = DREAL(Z)**2
c     RETURN
c     END SUBROUTINE SQRABS
c     END FUNCTION SQRABS
*
      END SUBROUTINE AF3X3LEV

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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine AF3X3LEV(RDIST,DELTAE,C3val,C6val,C8val,De,ULR)
!=======================================================================
!*** Simplified version of AF3x3potRet which does not return derivatives

real*8 H(3,3), DM1(3,3), DM3(3,3), DM5(3,3), DR(3,3), DDe(3,3), Q(3,3)
real*8 DEIGM1(1,1), DEIGM3(1,1), DEIGM5(1,1), DEIGR(1,1), DEIGDe(1,1), EIGVEC(3,1), W(3)
real*8 RDIST, RDIST2, RDIST3, DELTAE, C3val, C6val, C8val, De, ULR, RET, RETSig, RETPi, Modulus, M1, M3, M5, Z
integer I, J, L

M1 = C3val
M3 = C6val
M5 = C8val
RET = 9.36423830d-4*RDIST
RETSig = dcos(RET)+(RET)*dsin(RET)
RETPi = RETSig-RET**2*dcos(RET)
RDIST2 = RDIST**2
RDIST3 = RDIST*RDIST2
!write(25,*) 'Variables = "r", "U(r)","U(r)-U(r)^2/(4De)"'
!write(25,*) 'zone T = "U(r)"'
! Initialize interaction matrix to 0.d0
do I=1,3
  H(I,I) = 0.0d0
end do
! Prepare interation matrix  H
H(1,1) = -(M1*RETSig+M3/(RDIST3)+M5/(RDIST3*RDIST2))/(3.d0*RDIST3)
H(1,2) = -(dsqrt(2.d0))*H(1,1)
H(2,1) = H(1,2)
H(1,3) = M1*RETPi/(dsqrt(6.d0)*RDIST3)
H(3,1) = H(1,3)
H(2,2) = 2*H(1,1)+DELTAE
H(2,3) = H(1,3)/dsqrt(2.d0)
H(3,2) = H(2,3)
H(3,3) = DELTAE
! Prepare radial derivative of interaction matrix (? is it needed ?)
DR(1,1) = (3.d0*M1*RETSig+6.d0*M3/RDIST3+8.d0*M5/(RDIST3*RDIST2))/(3.d0*RDIST3*RDIST)
DR(1,2) = -dsqrt(2.d0)*DR(1,1)
DR(2,1) = DR(1,2)
DR(2,2) = 2.d0*DR(1,1)
DR(1,3) = -3.d0*H(1,3)/RDIST
DR(3,1) = DR(1,3)
DR(2,3) = -3.d0*H(2,3)/RDIST
DR(3,2) = DR(2,3)
DR(3,3) = 0.d0
! Partial derivative of interaction matric  H  w.r.t.  C3
DM1(1,1) = -(RETSig+M1/(2.d0*De*RDIST3))/(3.d0*RDIST3)
DM1(1,2) = -dsqrt(2.d0)*DM1(1,1)
DM1(2,1) = DM1(1,2)
DM1(2,2) = 2.d0*DM1(1,1)
DM1(1,3) = RETPi/(dsqrt(6.d0)*RDIST3)
DM1(3,1) = DM1(1,3)
DM1(2,3) = DM1(1,3)/dsqrt(2.d0)
DM1(3,2) = DM1(2,3)
DM1(3,3) = 0.d0
! Partial derivative of interaction matric  H  w.r.t.  C6
DM3(1,1) = -1.d0/(3.d0*RDIST3**2)
DM3(1,2) = -sqrt(2.d0)*DM3(1,1)
DM3(1,3) = 0.d0
DM3(2,1) = DM3(1,2)
DM3(2,2) = 2.d0*DM3(1,1)
DM3(2,3) = 0.d0
DM3(3,1) = DM3(1,3)
DM3(3,2) = DM3(2,3)
DM3(3,3) = 0.d0
! Partial derivative of interaction matric  H  w.r.t.  C8
DM5(1,1) = DM3(1,1)/(RDIST2)
DM5(1,2) = DM3(1,2)/(RDIST2)
DM5(1,3) = 0.d0
DM5(2,1) = DM3(1,2)
DM5(2,2) = DM3(2,2)/(RDIST2)
DM5(2,3) = 0.d0
DM5(3,1) = DM5(1,3)
DM5(3,2) = DM5(2,3)
DM5(3,3) = 0.d0
! Partial derivative of interaction matric  H  w.r.t.  De
DDe(1,1) = M1**2/(12.d0*(RDIST3*De)**2)
DDe(1,2) = -sqrt(2.d0)*DDe(1,1)
DDe(1,3) = 0.d0
DDe(2,1) = DDe(1,2)
DDe(2,2) = 2.d0*DDe(1,1)
DDe(2,3) = 0.d0
DDe(3,1) = DDe(1,3)
DDe(3,2) = DDe(2,3)
DDe(3,3) = 0.d0
! Call subroutine to prepare and invert interaction matrix  H
call ZHEEVJ3(H,Q,W)
L = 1
! Nor - identify the lowest eigenvalue of  H  and label it  L
do J=2,3
  if (W(J) < W(L)) then
    L = J
  end if
end do
ULR = -W(L)
do I=1,3
  EIGVEC(I,1) = Q(I,L)
end do
write(6,*) EIGVEC
! Make sure the following variables are "referenced":
DEIGM1 = 0.d0
DEIGM3 = DEIGM1*0.d0
DEIGM5 = DEIGM3*0.d0
DEIGR = DEIGM5*0.d0
DEIGDe = DEIGR*0.d0
DEIGM1 = DEIGDe*0.d0
!write(25,600) RDIST,ULR
!600 format(2D16.7)
!write(26,601) RDIST,DEIGM1,DEIGR,DEIGDe
!601 format(4D16.7)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Modulus = SQRABS(Z)
!Modulus = REAL(Z)**2
Modulus = 0.d0
Z = Modulus ! Make sure it's "referenced"
Modulus = Z ! Make sure it's "referenced"

return

contains

subroutine ZHEEVJ3(H,Q,W)
  !=======================================================================
  !** Subroutine to setup and invert the matrix  H  and return
  !   eigenvalues W and eigenvector matric  Q
  integer N, I, X, Y, R
  parameter(N=3)
  real*8 H(3,3), Q(3,3), W(3)
  real*8 SD, SO, S, T, C, G, B, Z, THRESH
  ! Initialize Q to the identitity matrix
  ! --- This loop can be omitted if only the eigenvalues are desired ---
  do X=1,N
    Q(X,X) = 1.0d0
    do Y=1,X-1
      Q(X,Y) = 0.0d0
      Q(Y,X) = 0.0d0
    end do
  end do
  ! Initialize W to diag(A)
  do X=1,N
    W(X) = real(H(X,X))
  end do
  ! Calculate SQR(tr(A))
  SD = 0.0d0
  do X=1,N
    SD = SD+abs(W(X))
  end do
  SD = SD**2
  ! Main iteration loop
  do I=1,50
    ! Test for convergence
    SO = 0.0d0
    do X=1,N
      do Y=X+1,N
        SO = SO+abs(real(H(X,Y)))
      end do
    end do
    if (SO == 0.0d0) return
    if (I < 4) then
      THRESH = 0.2d0*SO/N**2
    else
      THRESH = 0.0d0
    end if
    ! Do sweep
    do X=1,N
      do Y=X+1,N
        G = 100.0d0*(abs(real(H(X,Y))))
        if ((I > 4) .and. (abs(W(X))+G == abs(W(X))) .and. (abs(W(Y))+G == abs(W(Y)))) then
          H(X,Y) = 0.0d0
        else if (abs(real(H(X,Y))) > THRESH) then
          ! Calculate Jacobi transformation
          B = W(Y)-W(X)
          if ((abs(B)+G) == abs(B)) then
            T = H(X,Y)/B
          else
            if (B <= 0.0d0) then
              T = -2.0d0*H(X,Y)/(sqrt(B**2+4.0d0*real(H(X,Y))**2)-B)
              !                /(sqrt(B**2+4.0D0*SQRABS(H(X,Y)))-B)
            else if (B == 0.0d0) then
              T = H(X,Y)*(1.0d0/abs(H(X,Y)))
            else
              T = 2.0d0*H(X,Y)/(sqrt(B**2+4.0d0*real(H(X,Y))**2)+B)
              !               /(sqrt(B**2+4.0D0*SQRABS(H(X,Y)))+B)
            end if
          end if
          !C = 1.0D0/sqrt(1.0D0+SQRABS(T))
          C = 1.0d0/sqrt(1.0d0+real(T)**2)
          S = T*C
          Z = real(T*(H(X,Y)))
          ! Apply Jacobi transformation
          H(X,Y) = 0.0d0
          W(X) = W(X)-Z
          W(Y) = W(Y)+Z
          do R=1,X-1
            T = H(R,X)
            H(R,X) = C*T-(S)*H(R,Y)
            H(R,Y) = S*T+C*H(R,Y)
          end do
          do R=X+1,Y-1
            T = H(X,R)
            H(X,R) = C*T-S*(H(R,Y))
            H(R,Y) = S*(T)+C*H(R,Y)
          end do
          do R=Y+1,N
            T = H(X,R)
            H(X,R) = C*T-S*H(Y,R)
            H(Y,R) = (S)*T+C*H(Y,R)
          end do
          ! eigenvectors
          ! This loop can be omitted if only the eigenvalues are desired ---
          do R=1,N
            T = Q(R,X)
            Q(R,X) = C*T-(S)*Q(R,Y)
            Q(R,Y) = S*T+C*Q(R,Y)
          end do
        end if
      end do
    end do
  end do
  write(6,*) 'ZHEEVJ3: No convergence.'
end subroutine ZHEEVJ3

!real*8 function SQRABS(Z)
!subroutine SQRABS(Z)
!  !=====================================================================
!  ! Calculates the squared absolute value of a complex number Z
!  ! --------------------------------------------------------------------
!  !  Parameters ..
!  real*8 Z
!  SQRABS = dreal(Z)**2
!  return
!end subroutine SQRABS
!end function SQRABS

end subroutine AF3X3LEV

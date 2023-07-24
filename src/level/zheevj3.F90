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
subroutine ZHEEVJ3(H,Q,W)
!=======================================================================
! Subroutine to setup and invert the matrix  H  and return
! eigenvalues W and eigenvector matric  Q

use Constants, only: Zero, One, Two, Four
use Definitions, only: wp, iwp, u6

integer(kind=iwp), parameter :: N = 3
real(kind=wp), intent(inout) :: H(N,N)
real(kind=wp), intent(out) :: Q(N,N), W(N)
integer(kind=iwp) :: I, R, X, Y
real(kind=wp) :: B, C, G, S, SO, T, THRESH, Z

! Initialize Q to the identity matrix
! --- This loop can be omitted if only the eigenvalues are desired ---
call unitmat(Q,N)
! Initialize W to diag(A)
do X=1,N
  W(X) = H(X,X)
end do
! Main iteration loop
do I=1,50
  ! Test for convergence
  SO = Zero
  do X=1,N
    do Y=X+1,N
      SO = SO+abs(H(X,Y))
    end do
  end do
  if (SO == Zero) return
  if (I < 4) then
    THRESH = 0.2_wp*SO/N**2
  else
    THRESH = Zero
  end if
  ! Do sweep
  do X=1,N
    do Y=X+1,N
      G = 100.0_wp*(abs(H(X,Y)))
      if ((I > 4) .and. (abs(W(X))+G == abs(W(X))) .and. (abs(W(Y))+G == abs(W(Y)))) then
        H(X,Y) = Zero
      else if (abs(H(X,Y)) > THRESH) then
        ! Calculate Jacobi transformation
        B = W(Y)-W(X)
        if ((abs(B)+G) == abs(B)) then
          T = H(X,Y)/B
        else
          if (B <= Zero) then
            T = -Two*H(X,Y)/(sqrt(B**2+Four*H(X,Y)**2)-B)
            !              /(sqrt(B**2+Four*SQRABS(H(X,Y)))-B)
          else if (B == Zero) then
            T = H(X,Y)/abs(H(X,Y))
          else
            T = Two*H(X,Y)/(sqrt(B**2+Four*H(X,Y)**2)+B)
            !             /(sqrt(B**2+Four*SQRABS(H(X,Y)))+B)
          end if
        end if
        !C = One/sqrt(One+SQRABS(T))
        C = One/sqrt(One+T**2)
        S = T*C
        Z = T*H(X,Y)
        ! Apply Jacobi transformation
        H(X,Y) = Zero
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
write(u6,*) 'ZHEEVJ3: No convergence.'

end subroutine ZHEEVJ3

!function SQRABS(Z)
!!=====================================================================
!! Calculates the squared absolute value of a complex number Z
!! --------------------------------------------------------------------
!
!use Definitions, only: wp
!
!real(kind=wp) :: SQRABS
!real(kind=wp), intent(in) :: Z
!SQRABS = real(Z)**2
!
!return
!
!end function SQRABS

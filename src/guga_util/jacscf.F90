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

subroutine JACSCF(A,B,C,NAA,NQQ,EPSLON)
! VERSION 4 AUGUST 1971
! SUBROUTINE TO FIND ALL THE EIGENVALUES AND EIGENVECTORS OF A
! BY THE JACOBI METHOD
! THIS PROGRAM TREATS THE ORTHOGONAL CASE (S=1)
! NAA=DIMENSION OF A,B,C
! A TRIANGULAR, B MATRIX OF VECTORS, C EIGENVALUES
! B CLEARED TO UNIT MATRIX IF NQ=-1, NOT CLEARED IF NQ= NO. ROWS B
! NQ MUST NOT BE LESS THAN NAA
! EPSLON IS THE CONVERGENCE CRITERIA FOR OFF DIAGONAL ELEMENTS

use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: A(*)
real(kind=wp), intent(_OUT_) :: B(*), C(*)
integer(kind=iwp), intent(in) :: NAA, NQQ
real(kind=wp), intent(in) :: EPSLON
integer(kind=iwp) :: I, ID, J, JD, JJ, K, KA, KB, KC, KD, KDX, L, LL, LOOPC, N, NA, NN, NQ
real(kind=wp) :: ALPHA, AMAX, BETA, CC, CSQ, ROOT, RSUM, S, SSQ, TEMPA, TERM, THRESH, THRSHG, TWOSC

! POW: Unnecessary but warning stopping initialization
TERM = huge(TERM)

NQ = NQQ
LOOPC = 0
NA = NAA
NN = (NA*(NA+1))/2
if (NQ <= 0) then
  K = 1
  NQ = NA
  do I=1,NA
    do J=1,NA
      if (I /= J) then
        B(K) = Zero
      else
        B(K) = One
      end if
      K = K+1
    end do
  end do
end if
RSUM = Zero
if (NA < 1) return
if (NA > 1) then
  K = 1
  AMAX = Zero
  do I=1,NA
    do J=1,I
      if ((I /= J) .and. (abs(A(K)) > AMAX)) AMAX = abs(A(K))
      TERM = A(K)*A(K)
      RSUM = RSUM+TERM+TERM
      K = K+1
    end do
    RSUM = RSUM-TERM
  end do
  RSUM = sqrt(RSUM)
  THRESH = RSUM/sqrt(real(NA,kind=wp))
  THRSHG = THRESH*EPSLON
  if (THRSHG < AMAX) then
    THRESH = AMAX/Three
    if (THRESH < THRSHG) THRESH = THRSHG
    do
      K = 2
      N = 0
      JD = 1
      KDX = 0
      do J=2,NA
        ID = 0
        JD = JD+J
        JJ = J-1
        KC = 0
        KDX = KDX+NQ
        do I=1,JJ
          ID = ID+I
          if (abs(A(K)) <= THRESH) then
            KC = KC+NQ
          else
            N = N+1
            ALPHA = (A(JD)-A(ID))/(Two*A(K))
            BETA = One/(One+ALPHA*ALPHA)
            ROOT = One+abs(ALPHA)*sqrt(BETA)
            SSQ = BETA/(2*ROOT)
            CSQ = ROOT/2
            CC = sqrt(CSQ)
            S = -sqrt(SSQ)*sign(One,ALPHA)
            TWOSC = 2*CC*S
            TEMPA = CSQ*A(ID)+TWOSC*A(K)+SSQ*A(JD)
            A(JD) = CSQ*A(JD)-TWOSC*A(K)+SSQ*A(ID)
            A(ID) = TEMPA
            A(K) = Zero
            KA = JD-J
            KB = ID-I
            KD = KDX
            do L=1,NQ
              KC = KC+1
              KD = KD+1
              TEMPA = CC*B(KC)+S*B(KD)
              B(KD) = -S*B(KC)+CC*B(KD)
              B(KC) = TEMPA
              if (L <= NA) then
                if (I == L) then
                  KB = KB+1
                  KA = KA+1
                else if (I < L) then
                  KB = KB+L-1
                  if (J == L) then
                    KA = KA+1
                  else
                    if (J > L) then
                      KA = KA+1
                    else
                      KA = KA+L-1
                    end if
                    TEMPA = CC*A(KB)+S*A(KA)
                    A(KA) = -S*A(KB)+CC*A(KA)
                    A(KB) = TEMPA
                  end if
                else
                  KB = KB+1
                  KA = KA+1
                  TEMPA = CC*A(KB)+S*A(KA)
                  A(KA) = -S*A(KB)+CC*A(KA)
                  A(KB) = TEMPA
                end if
              end if
            end do
          end if
          K = K+1
        end do
        K = K+1
      end do
      LOOPC = LOOPC+1
      if (LOOPC >= 50) return
      if (N <= NN/8) then
        if (THRESH /= THRSHG) then
          THRESH = THRESH/Three
          if (THRESH < THRSHG) THRESH = THRSHG
        else if (N == 0) then
          exit
        end if
      end if
    end do
  end if
end if
LL = 0
do L=1,NA
  LL = LL+L
  C(L) = A(LL)
end do

return

end subroutine JACSCF

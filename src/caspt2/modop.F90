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

subroutine MODOP(OP1,NOP2,OP2,NOP3,OP3)
! Purpose: Modify the coefficients in OP1 and OP2, using the
! input values of OP2 and OP3, so that the operators can be
! calculated using products of elementary excitation
! operators rather than normal-ordered products.

use Index_Functions, only: iTri, nTri_Elem, nTri3_Elem
use caspt2_module, only: NACTEL, NASHT
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NOP2, NOP3
real(kind=wp), intent(inout) :: OP1(NASHT,NASHT), OP2(NOP2)
real(kind=wp), intent(in) :: OP3(NOP3)
integer(kind=iwp) :: I, I_N, IJ, IJKL, IJKLMN, IL, IND, J, K, KL, KN, L, M, MN, N
real(kind=wp) :: X

if (NACTEL > 2) then

  do I=1,NASHT
    do J=1,NASHT
      IJ = I+NASHT*(J-1)
      do K=1,NASHT
        do L=1,NASHT
          KL = K+NASHT*(L-1)
          if (KL > IJ) cycle

          do M=1,NASHT
            do N=1,NASHT
              MN = M+NASHT*(N-1)
              if (MN > KL) cycle
              IJKLMN = nTri3_Elem(IJ-1)+nTri_Elem(KL-1)+MN
              X = OP3(IJKLMN)
              if (abs(X) < 1.0e-15_wp) cycle

              if (K == J) then
                IL = I+NASHT*(L-1)
                IND = iTri(IL,MN)
                OP2(IND) = OP2(IND)-X
                if (M == L) OP1(I,N) = OP1(I,N)-X
              end if

              if (M == J) then
                I_N = I+NASHT*(N-1)
                IND = iTri(I_N,KL)
                OP2(IND) = OP2(IND)-X
              end if

              if (M == L) then
                KN = K+NASHT*(N-1)
                IND = iTri(KN,IJ)
                OP2(IND) = OP2(IND)-X
              end if

            end do
          end do

        end do
      end do
    end do
  end do

end if

if (NACTEL > 1) then

  do I=1,NASHT
    do J=1,NASHT
      IJ = I+NASHT*(J-1)
      do K=1,NASHT
        do L=1,NASHT
          KL = K+NASHT*(L-1)
          if (KL > IJ) cycle

          if (J == K) then
            IJKL = iTri(IJ,KL)
            OP1(I,L) = OP1(I,L)-OP2(IJKL)
          end if

        end do
      end do
    end do
  end do

end if

end subroutine MODOP

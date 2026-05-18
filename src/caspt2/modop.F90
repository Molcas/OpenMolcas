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

use caspt2_module, only: NASHT, NACTEL
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: NOP2, NOP3
real(kind=wp), intent(inout) :: OP1(NASHT,NASHT), OP2(NOP2)
real(kind=wp), intent(in) :: OP3(NOP3)
integer(kind=iwp) I, J, K, L, IJ, KL, M, N, MN, IJKLMN, IL, IN, KN, IND, IJKL
real(kind=wp) X

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
              IJKLMN = ((IJ+1)*IJ*(IJ-1))/6+(KL*(KL-1))/2+MN
              X = OP3(IJKLMN)
              if (abs(X) < 1.0e-15_wp) cycle

              if (K == J) then
                IL = I+NASHT*(L-1)
                if (IL >= MN) then
                  IND = (IL*(IL-1))/2+MN
                else
                  IND = (MN*(MN-1))/2+IL
                end if
                OP2(IND) = OP2(IND)-X
                if (M == L) then
                  OP1(I,N) = OP1(I,N)-X
                end if
              end if

              if (M == J) then
                IN = I+NASHT*(N-1)
                if (IN >= KL) then
                  IND = (IN*(IN-1))/2+KL
                else
                  IND = (KL*(KL-1))/2+IN
                end if
                OP2(IND) = OP2(IND)-X
              end if

              if (M == L) then
                KN = K+NASHT*(N-1)
                if (KN >= IJ) then
                  IND = (KN*(KN-1))/2+IJ
                else
                  IND = (IJ*(IJ-1))/2+KN
                end if
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
            IJKL = (IJ*(IJ-1))/2+KL
            OP1(I,L) = OP1(I,L)-OP2(IJKL)
          end if

        end do
      end do
    end do
  end do

end if

end subroutine MODOP

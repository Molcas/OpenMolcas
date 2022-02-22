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

subroutine PRMAT(IUOUT,R,N,M,HEAD)
! SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
! TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
! ((N+1)*N)/2 WHEN M IS ZERO

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IUOUT, N, M
real(kind=wp), intent(in) :: R(*)
character(len=*), intent(in) :: HEAD
integer(kind=iwp) :: I, I1, I2, IBL, IJ, IR, J, J1, J2, K, K1, K1S, K2, KD, KK, NKPB

write(IUOUT,1001) HEAD
NKPB = 4
if (M <= 0) then
  IBL = N/NKPB
  IR = N-IBL*NKPB
  J1 = 1
  K1S = 1
  KD = 0
  if (IBL /= 0) then
    J2 = NKPB
    do I=1,IBL
      write(IUOUT,1002) (J,J=J1,J2)
      K1 = K1S
      K2 = K1
      KK = 0
      do J=J1,J2
        write(IUOUT,1003) J,(R(K),K=K1,K2)
        KK = KK+1
        K1 = K1+KD+KK
        K2 = K1+KK
      end do
      J1 = J1+NKPB
      if (J1 > N) return
      J2 = J2+NKPB
      K2 = K1-1
      K1 = K2+1
      K2 = K1+(NKPB-1)
      K1S = K2+1
      KK = KD+NKPB
      do J=J1,N
        write(IUOUT,1003) J,(R(K),K=K1,K2)
        KK = KK+1
        K1 = K1+KK
        K2 = K2+KK
      end do
      KD = KD+NKPB
    end do
  end if
  if (IR /= 0) then
    K1 = K1S
    J2 = J1+IR-1
    KK = 0
    K2 = K1
    write(IUOUT,1002) (J,J=J1,J2)
    write(IUOUT,1003)
    do J=J1,J2
      write(IUOUT,1003) J,(R(K),K=K1,K2)
      KK = KK+1
      K1 = K1+KD+KK
      K2 = K1+KK
    end do
  end if
else
  IBL = M/NKPB
  IR = M-IBL*NKPB
  I2 = 0
  K2 = 0
  if (IBL /= 0) then
    do I=1,IBL
      I1 = (I-1)*N*NKPB+1
      I2 = I1+(NKPB-1)*N
      K1 = K2+1
      K2 = K1+(NKPB-1)
      write(IUOUT,1002) (K,K=K1,K2)
      do J=1,N
        write(IUOUT,1003) J,(R(IJ),IJ=I1,I2,N)
        I1 = I1+1
        I2 = I1+(NKPB-1)*N
      end do
    end do
  end if
  if (IR /= 0) then
    I1 = IBL*N*NKPB+1
    I2 = I1+(IR-1)*N
    K1 = K2+1
    K2 = M
    write(IUOUT,1002) (K,K=K1,K2)
    write(IUOUT,1003)
    do J=1,N
      write(IUOUT,1003) J,(R(IJ),IJ=I1,I2,N)
      I1 = I1+1
      I2 = I1+(IR-1)*N
    end do
  end if
  write(IUOUT,1003)
end if

return

1001 format(' MATRIX PRINTED:',2X,A)
1002 format(' ',4X,4(6X,I4,6X),/)
1003 format(' ',I4,4d16.8)

end subroutine PRMAT

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

subroutine HZ(ARR)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
parameter(IX1F=1,IX2F=2,IRR=3,IX1R=4,IX2R=5,IX1X1=6,IX2X1=7,IX2X2=8,IFDF=9,IFDR=10,IRDR=11)
dimension TMP(MXVEC,MXVEC)
dimension ARR(NRROOT,NRROOT,11)

! THIS SUBROUTINE FORMS THE OVERLAP AND H-ZERO MATRIX ELEMENTS
! IN THE BASIS OF PSI, RHO, XI1, AND XI2 FUNCTIONS.
!
!write(6,*)
!write(6,*)' CHECK PRINTS IN HZ.'
!write(6,*)' X1F ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX1F),J=1,NRROOT)
!end do
!write(6,*)' X2F ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX2F),J=1,NRROOT)
!end do
!write(6,*)' RR ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IRR ),J=1,NRROOT)
!end do
!write(6,*)' X1R ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX1R),J=1,NRROOT)
!end do
!write(6,*)' X2R ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX2R),J=1,NRROOT)
!end do
!write(6,*)' X1X1 ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX1X1),J=1,NRROOT)
!end do
!write(6,*)' X2X1 ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX2X1),J=1,NRROOT)
!end do
!write(6,*)' X2X2 ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IX2X2),J=1,NRROOT)
!end do
!write(6,*)' FDF ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IFDF),J=1,NRROOT)
!end do
!write(6,*)' FDR ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IFDR),J=1,NRROOT)
!end do
!write(6,*)' RDR ARRAY:'
!do I=1,NRROOT
!  write(6,'(1X,5F15.6)') (ARR(I,J,IRDR),J=1,NRROOT)
!end do
! FIRST, CREATE OVERLAP MATRIX, AND INITIALIZE HZERO MATRIX WITH ALL
! TERMS THAT DO NOT REQUIRE MATRIX MULTIPLIES:
do I1=1,NRROOT
  I2 = I1+NRROOT
  I3 = I1+2*NRROOT
  I4 = I1+3*NRROOT
  do J1=1,NRROOT
    J2 = J1+NRROOT
    J3 = J1+2*NRROOT
    J4 = J1+3*NRROOT
    SZERO(I1,J1) = 0.0d00
    SZERO(I2,J1) = 0.0d00
    SZERO(I3,J1) = ARR(I1,J1,IX1F)
    SZERO(I4,J1) = ARR(I1,J1,IX2F)
    SZERO(I2,J2) = ARR(I1,J1,IRR)
    SZERO(I3,J2) = ARR(I1,J1,IX1R)
    SZERO(I4,J2) = ARR(I1,J1,IX2R)
    SZERO(I3,J3) = ARR(I1,J1,IX1X1)
    SZERO(I4,J3) = ARR(I1,J1,IX2X1)
    SZERO(I4,J4) = ARR(I1,J1,IX2X2)
    HZERO(I1,J1) = 0.0d00
    if (I1 == J1) then
      SZERO(I1,J1) = 1.0d00
      HZERO(I1,J1) = ESMALL(I1)
    end if
    HZERO(I2,J1) = ARR(I1,J1,IRR)
    HZERO(I3,J1) = ARR(I1,J1,IX1F)*ESMALL(J1)+ARR(I1,J1,IX1R)
    HZERO(I4,J1) = ARR(I1,J1,IX2F)*ESMALL(J1)+ARR(I1,J1,IX2R)
    HZERO(I2,J2) = ARR(I1,J1,IRDR)
    HZERO(I3,J2) = ESMALL(I1)*ARR(I1,J1,IX1R)
    HZERO(I4,J2) = ESMALL(I1)*ARR(I1,J1,IX2R)+ARR(I1,J1,IRR)
    HZERO(I3,J3) = ARR(I1,J1,IX1X1)*ESMALL(J1)-ARR(J1,I1,IX1F)
    HZERO(I4,J3) = ARR(I1,J1,IX2X1)*ESMALL(J1)
    HZERO(I4,J4) = ARR(I1,J1,IX2X2)*ESMALL(J1)+ARR(I1,J1,IX2R)
  end do
end do
IO2 = NRROOT
IO3 = 2*NRROOT
IO4 = 3*NRROOT
do I=1,NRROOT
  do J=1,NRROOT
    SUM1 = HZERO(IO3+I,IO2+J)
    SUM2 = HZERO(IO4+I,IO2+J)
    do K=1,NRROOT
      SUM1 = SUM1+ARR(I,K,IX1F)*(ARR(K,J,IRR)-ARR(K,J,IFDR))
      SUM2 = SUM2+ARR(I,K,IX2F)*(ARR(K,J,IRR)-ARR(K,J,IFDR))
    end do
    HZERO(IO3+I,IO2+J) = SUM1
    HZERO(IO4+I,IO2+J) = SUM2
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM1 = 0.0d00
    SUM2 = 0.0d00
    do K=1,NRROOT
      SUM1 = SUM1+ARR(I,K,IX1F)*ARR(J,K,IX1F)
      SUM2 = SUM2+ARR(I,K,IX2F)*ARR(J,K,IX1F)
    end do
    HZERO(IO3+I,IO3+J) = HZERO(IO3+I,IO3+J)-SUM1*ESMALL(J)
    HZERO(IO3+J,IO3+I) = HZERO(IO3+J,IO3+I)-SUM1*ESMALL(J)
    HZERO(IO4+I,IO3+J) = HZERO(IO4+I,IO3+J)-SUM2*ESMALL(J)
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM1 = 0.0d00
    SUM2 = 0.0d00
    do K=1,NRROOT
      SUM1 = SUM1+ARR(I,K,IX2F)*ARR(J,K,IX1F)
      SUM2 = SUM2+ARR(I,K,IX2F)*ARR(J,K,IX2F)
    end do
    HZERO(IO4+I,IO3+J) = HZERO(IO4+I,IO3+J)-SUM1*ESMALL(I)
    HZERO(IO4+I,IO4+J) = HZERO(IO4+I,IO4+J)-SUM2*ESMALL(I)
    HZERO(IO4+J,IO4+I) = HZERO(IO4+J,IO4+I)-SUM2*ESMALL(I)
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM = 0.0d00
    do K=1,NRROOT
      SUM = SUM+ARR(I,K,IFDF)*ARR(J,K,IX1F)
    end do
    TMP(I,J) = SUM
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM1 = HZERO(IO3+I,IO3+J)
    SUM2 = HZERO(IO4+I,IO3+J)
    do K=1,NRROOT
      SUM1 = SUM1+ARR(I,K,IX1F)*TMP(K,J)
      SUM2 = SUM2+ARR(I,K,IX2F)*TMP(K,J)
    end do
    HZERO(IO3+I,IO3+J) = SUM1
    HZERO(IO4+I,IO3+J) = SUM2
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM = 0.0d00
    do K=1,NRROOT
      SUM = SUM+ARR(I,K,IFDF)*ARR(J,K,IX2F)
    end do
    TMP(I,J) = SUM
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM2 = HZERO(IO4+I,IO4+J)
    do K=1,NRROOT
      SUM2 = SUM2+ARR(I,K,IX2F)*TMP(K,J)
    end do
    HZERO(IO4+I,IO4+J) = SUM2
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    TMP(I,J) = ESMALL(I)*ARR(J,I,IX1F)+ARR(J,I,IX1R)
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM1 = HZERO(IO3+I,IO3+J)
    SUM2 = HZERO(IO4+I,IO3+J)
    do K=1,NRROOT
      SUM1 = SUM1+ARR(I,K,IX1F)*TMP(K,J)+ARR(I,K,IX1R)*ARR(J,K,IX1F)
      SUM2 = SUM2+ARR(I,K,IX2F)*TMP(K,J)+ARR(I,K,IX2R)*ARR(J,K,IX1F)
    end do
    HZERO(IO3+I,IO3+J) = SUM1
    HZERO(IO4+I,IO3+J) = SUM2
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    TMP(I,J) = ESMALL(I)*ARR(J,I,IX2F)+ARR(J,I,IX2R)
  end do
end do
do I=1,NRROOT
  do J=1,NRROOT
    SUM2 = HZERO(IO4+I,IO4+J)
    do K=1,NRROOT
      SUM2 = SUM2+ARR(I,K,IX2F)*TMP(K,J)+ARR(I,K,IX2R)*ARR(J,K,IX2F)
    end do
    HZERO(IO4+I,IO4+J) = SUM2
  end do
end do
do I1=1,NRROOT
  I2 = I1+NRROOT
  I3 = I1+2*NRROOT
  I4 = I1+3*NRROOT
  do J1=1,NRROOT
    J2 = J1+NRROOT
    J3 = J1+2*NRROOT
    J4 = J1+3*NRROOT
    HZERO(I1,J2) = HZERO(J2,I1)
    HZERO(I1,J3) = HZERO(J3,I1)
    HZERO(I1,J4) = HZERO(J4,I1)
    HZERO(I2,J3) = HZERO(J3,I2)
    HZERO(I2,J4) = HZERO(J4,I2)
    HZERO(I3,J4) = HZERO(J4,I3)
    SZERO(I1,J2) = SZERO(J2,I1)
    SZERO(I1,J3) = SZERO(J3,I1)
    SZERO(I1,J4) = SZERO(J4,I1)
    SZERO(I2,J3) = SZERO(J3,I2)
    SZERO(I2,J4) = SZERO(J4,I2)
    SZERO(I3,J4) = SZERO(J4,I3)
  end do
end do

return

end subroutine HZ

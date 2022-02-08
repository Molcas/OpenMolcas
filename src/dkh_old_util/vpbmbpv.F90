!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Bernd Artur Hess                                 *
!***********************************************************************

subroutine VPBMBPV(EIGA,EIGB,REVTA,SINVA,SINVB,NA,NB,ISIZEA,AAA,AAB,RRA,RRB,PVA,VPA,ISYMA,ISYMB,BU2,G2,AUX2,CMM1,BU4,CMM2,PVAA, &
                   VPAA,SCPV,SCVP)
! $Id: vpbmbpv.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NA, NB, ISIZEA, ISYMA, ISYMB
real(kind=wp), intent(in) :: EIGA(NA,NA), EIGB(NB,NB), REVTA(NA,NA), SINVA(NA,NA), SINVB(NB,NB), AAA(NA), AAB(NB), RRA(NA), &
                             RRB(NB), PVAA(ISIZEA), VPAA(ISIZEA), SCPV(NB,NA), SCVP(NB,NA)
real(kind=wp), intent(inout) :: PVA(NA,NB), VPA(NA,NB), CMM2(NA,NB)
real(kind=wp), intent(out) :: BU2(NA,NB), G2(NA,NB), AUX2(NA,NB), CMM1(NA,NB), BU4(NA,NB)
integer(kind=iwp) :: I, IJ, J, K

if (iSyma == iSymb) then
  IJ = 0
  do I=1,NA
    do J=1,I
      IJ = IJ+1
      PVA(I,J) = PVAA(IJ)
      PVA(J,I) = -VPAA(IJ)
      VPA(I,J) = VPAA(IJ)
      VPA(J,I) = -PVAA(IJ)
    end do
  end do
else if (iSyma < iSymb) then
  do I=1,NA
    do J=1,NB
      PVA(I,J) = SCPV(J,I)
      VPA(I,J) = SCVP(J,I)
    end do
  end do
  CMM2(:,:) = Zero
end if

! TRANSFORM pV TO T-BASIS

call TrSmrN(PVA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
call TrSmrN(G2,EIGA,EIGB,BU2,NA,NB,AUX2,CMM1)

G2(:,:) = Zero
Aux2(:,:) = Zero

! TRANSFORM Vp TO T-BASIS

call TrSmrN(VPA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
call TrSmrN(G2,EIGA,EIGB,BU4,NA,NB,AUX2,CMM1)

G2(:,:) = Zero
Aux2(:,:) = Zero
Cmm1(:,:) = Zero

! Multiply

do I=1,NA
  do J=1,NB
    G2(I,J) = BU4(I,J)*RRB(J)
    CMM1(I,J) = -BU2(I,J)*RRA(I)
  end do
end do

!write(u6,*)
!write(u6,*) 'pV part of commutator MATRIX'
!write(u6,*)  G2
!write(u6,*)

!write(u6,*)
!write(u6,*) '-Vp part commutator MATRIX'
!write(u6,*)  CMM1
!write(u6,*)

! Calculate the commutator <iSymA|[V,pb]|iSymB> and put into CMM1

CMM1(:,:) = CMM1+G2

! Multiply BY A MATRIX

do I=1,NA
  do J=1,NB
    CMM1(I,J) = CMM1(I,J)*AAA(I)*AAB(J)
  end do
end do

! Multiply BY REVTA  MATRIX

do I=1,NA
  do J=1,NB
    do K=1,NA
      CMM2(I,J) = CMM2(I,J)+REVTA(I,K)*CMM1(K,J)
    end do
  end do
end do

Cmm1(:,:) = Zero

!write(u6,*)
!write(u6,*) 'CMM2  MATRX FINAL'
!write(u6,*)  CMM2
!write(u6,*)

!write(u6,*) 'END OF VPBMBPV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

return

end subroutine VPBMBPV

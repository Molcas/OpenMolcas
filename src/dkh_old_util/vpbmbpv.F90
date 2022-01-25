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
subroutine VPBMBPV(idbg,epsilon,EIGA,EIGB,REVTA,SINVA,SINVB,NA,NB,ISIZEA,ISIZEB,VELIT,AAA,AAB,RRA,RRB,PVA,VPA,ISYMA,ISYMB,BU2,G2, &
                   AUX2,CMM1,BU4,CMM2,PVAA,VPAA,SCPV,SCVP)
! $Id: vpbmbpv.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de

implicit real*8(A-H,O-Z)
#include "real.fh"
dimension PVAA(ISIZEA), VPAA(ISIZEA)
dimension EIGA(NA,NA), EIGB(NB,NB), SINVA(NA,NA), SINVB(NB,NB), REVTA(NA,NA), AAA(NA), RRA(NA), AAB(NB), RRB(NB), PVA(NA,NB), &
          VPA(NA,NB), BU2(NA,NB), G2(NA,NB), CMM1(NA,NB), BU4(NA,NB)
dimension AUX2(NA,NB), CMM2(NA,NB), SCPV(NB,NA), SCVP(NB,NA)

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
end if

if (iSyma < iSymb) then
  do I=1,NA
    do J=1,NB
      CMM1(I,J) = SCPV(J,I)
      CMM2(I,J) = SCVP(J,I)
    end do
  end do
  do I=1,NA
    do J=1,NB
      PVA(I,J) = CMM1(I,J)
      VPA(I,J) = CMM2(I,J)
    end do
  end do
  call DCOPY_(NA*NB,[ZERO],0,CMM1,1)
  call DCOPY_(NA*NB,[ZERO],0,CMM2,1)

end if

! TRANSFORM pV TO T-BASIS

call TrSmrN(PVA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
call TrSmrN(G2,EIGA,EIGB,BU2,NA,NB,AUX2,CMM1)

call dcopy_(na*nb,[Zero],0,G2,1)
call dcopy_(na*nb,[Zero],0,Aux2,1)

! TRANSFORM Vp TO T-BASIS

call TrSmrN(VPA,SINVA,SINVB,G2,NA,NB,AUX2,CMM1)
call TrSmrN(G2,EIGA,EIGB,BU4,NA,NB,AUX2,CMM1)

call dcopy_(na*nb,[Zero],0,G2,1)
call dcopy_(na*nb,[Zero],0,Aux2,1)
call dcopy_(na*nb,[Zero],0,Cmm1,1)

! Multiply

do I=1,NA
  do J=1,NB
    G2(I,J) = BU4(I,J)*RRB(J)
    CMM1(I,J) = -BU2(I,J)*RRA(I)
  end do
end do

!write(6,*)
!write(6,*) 'pV part of commutator MATRIX'
!write(6,*)  G2
!write(6,*)

!write(6,*)
!write(6,*) '-Vp part commutator MATRIX'
!write(6,*)  CMM1
!write(6,*)

! Calculate the commutator <iSymA|[V,pb]|iSymB> and put into CMM1

call AddMar(NA*NB,G2,CMM1)

! MultiplY BY A MATRIX

do I=1,NA
  do J=1,NB
    CMM1(I,J) = CMM1(I,J)*AAA(I)*AAB(J)
  end do
end do

! MultiplY BY REVTA  MATRIX

do I=1,NA
  do J=1,NB
    do K=1,NA
      CMM2(I,J) = CMM2(I,J)+REVTA(I,K)*CMM1(K,J)
    end do
  end do
end do

call dcopy_(na*nb,[Zero],0,Cmm1,1)

!write(6,*)
!write(6,*) 'CMM2  MATRX FINAL'
!write(6,*)  CMM2
!write(6,*)

!write(6,*) 'END OF VPBMBPV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(idbg)
  call Unused_real(epsilon)
  call Unused_integer(ISIZEB)
  call Unused_real(VELIT)
end if

end subroutine VPBMBPV

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

subroutine SCFCLI4(idbg,eps,S,H,REVTA,SINVA,NA,NB,ISIZEA,ISIZEB,VELIT,AAA,AAB,ISYMA,ISYMB,CMM1,CMM2,EV4,BU2,BU6,EIG4,EW4,P)
! $Id: relsewc.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: idbg, NA, NB, ISIZEA, ISIZEB, ISYMA, ISYMB
real(kind=wp) :: eps, S(ISIZEA), H(ISIZEA), REVTA(NA,NA), SINVA(NA,NA), VELIT, AAA(NA), AAB(NB), CMM1(NB,NA), CMM2(NA,NB), &
                 EV4(ISIZEA), BU2(NA,NB), BU6(NA,NA), EIG4(NA,NA), EW4(NA), P(ISIZEA)
integer(kind=iwp) :: I, IJ, J, L

do I=1,NB
  do J=1,NA
    CMM1(I,J) = -CMM2(J,I)
  end do
end do

! Now we changed sign due to the imaginary character

IJ = 0
do I=1,NA
  do J=1,I
    IJ = IJ+1
    EV4(IJ) = Zero
    do L=1,NB
      EV4(IJ) = EV4(IJ)-CMM2(I,L)*CMM1(L,J)
    end do
  end do
end do

!write(u6,*)
!write(u6,*) 'Final BSS matrix no Hess part alpha**2 yet'
!write(u6,*) EV4

do IJ=1,ISIZEA
  EV4(IJ) = Half*(One/(VELIT*VELIT))*EV4(IJ)
end do

call AddMar(ISIZEA,EV4,H)

!ulf
if (idbg > 0) call PRMAT(IDBG,h,na,nb,'h   oper')
call Sogr(idbg,NA,S,SINVA,P,BU6,EW4)

call Diagr(H,NA,EIG4,EW4,SINVA,BU6,EV4)

!if (idbg > 0) write(idbg,*) '--- EIGENVALUES OF H MATRIX ---'
!if (idbg > 0) write(idbg,'(4D20.12)') EW4

!write(u6,*) 'END OF SCFCLI4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(eps)
  call Unused_real_array(REVTA)
  call Unused_integer(ISIZEB)
  call Unused_real_array(AAA)
  call Unused_real_array(AAB)
  call Unused_integer(ISYMA)
  call Unused_integer(ISYMB)
  call Unused_real_array(BU2)
end if

end subroutine SCFCLI4

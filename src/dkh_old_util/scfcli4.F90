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

subroutine SCFCLI4(idbg,S,H,SINVA,NA,NB,ISIZEA,VELIT,CMM1,CMM2,EV4,BU6,EIG4,EW4,P)
! $Id: relsewc.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de

use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: idbg, NA, NB, ISIZEA
real(kind=wp), intent(in) :: S(ISIZEA), VELIT, CMM2(NA,NB)
real(kind=wp), intent(inout) :: H(ISIZEA)
real(kind=wp), intent(out) :: SINVA(NA,NA), CMM1(NB,NA), EV4(ISIZEA), BU6(NA,NA), EIG4(NA,NA), EW4(NA), P(ISIZEA)
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

H(:) = H+EV4

!ulf
if (idbg > 0) call PRMAT(IDBG,h,na,nb,'h   oper')
call Sogr(idbg,NA,S,SINVA,P,BU6,EW4)

call Diagr(H,NA,EIG4,EW4,SINVA,BU6,EV4)

!if (idbg > 0) write(idbg,*) '--- EIGENVALUES OF H MATRIX ---'
!if (idbg > 0) write(idbg,'(4D20.12)') EW4

!write(u6,*) 'END OF SCFCLI4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

return

end subroutine SCFCLI4

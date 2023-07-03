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

subroutine Strtch(xyz,nCent,Avst,B,lWrite,Label,dB,ldB)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 B(3,nCent), xyz(3,nCent), dB(3,nCent,3,nCent), R(3)
logical lWrite, ldB
character*8 Label
#include "angstr.fh"

R(1) = xyz(1,2)-xyz(1,1)
R(2) = xyz(2,2)-xyz(2,1)
R(3) = xyz(3,2)-xyz(3,1)
R2 = R(1)**2+R(2)**2+R(3)**2
RR = sqrt(R2)
Avst = RR

aRR = RR*angstr
if (lWrite) write(6,'(1X,A,A,2(F10.6,A))') Label,' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
if (aRR < 1d-6) then
  call WarningMessage(2,'Abend in Strtch')
  write(6,*) '***************** ERROR **********************'
  write(6,*) ' Short (or negative) distance for coordinate: '
  write(6,'(1X,A,A,2(F10.6,A))') Label,' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
  write(6,*) '**********************************************'
  write(6,*)
  call Quit_OnUserError()
end if

! Compute the WDC B-matrix.

B(1,1) = -R(1)/RR
B(2,1) = -R(2)/RR
B(3,1) = -R(3)/RR

! Renormalize

xRR = sqrt(B(1,1)**2+B(2,1)**2+B(3,1)**2)
B(1,1) = B(1,1)/xRR
B(2,1) = B(2,1)/xRR
B(3,1) = B(3,1)/xRR

! Utilize translational invariance.
B(1,2) = -B(1,1)
B(2,2) = -B(2,1)
B(3,2) = -B(3,1)

! Compute the cartesian derivative of the B-matrix.

if (ldB) then

  do i=1,3
    do j=1,i
      if (i == j) then
        dB(i,1,j,1) = (One-B(j,1)*B(i,1))/RR
      else
        dB(i,1,j,1) = (-B(j,1)*B(i,1))/RR
      end if
      dB(j,1,i,1) = dB(i,1,j,1)

      dB(i,2,j,1) = -dB(i,1,j,1)
      dB(j,1,i,2) = dB(i,2,j,1)

      dB(i,1,j,2) = -dB(i,1,j,1)
      dB(j,2,i,1) = dB(i,1,j,2)

      dB(i,2,j,2) = -dB(i,2,j,1)
      dB(j,2,i,2) = dB(i,2,j,2)
    end do
  end do

end if

return

end subroutine Strtch

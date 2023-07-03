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

subroutine Bend(xyz,nCent,Fir,Bf,lWrite,lWarn,Label,dBf,ldB)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Bf(3,nCent), xyz(3,nCent), dBf(3,nCent,3,nCent), BRij(3,2), dBRij(3,2,3,2), BRjk(3,2), dBRjk(3,2,3,2)
logical lWrite, ldB, lWarn
character*8 Label

!                                                                      *
!***********************************************************************
!                                                                      *
!define _TIME_
!                                                                      *
!***********************************************************************
!                                                                      *
mCent = 2
call Strtch(xyz(1,1),mCent,Rij1,BRij,.false.,Label,dBRij,ldB)
call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.false.,Label,dBRjk,ldB)
Co = Zero
Crap = Zero
! BRij and BRjk should be normalized
do i=1,3
  Co = Co+BRij(i,1)*BRjk(i,2)
end do
do i=1,3
  Crap = Crap+(BRjk(i,2)-sign(One,Co)*BRij(i,1))**2
end do
Crap = sqrt(Crap)

! Special care for cases close to linearity

if (Crap < 1.0D-4) then
  Si = Crap
  if (Co < Zero) then
    Fir = Pi-ArSin(Si)
  else
    Fir = ArSin(Si)
  end if
else
  if (abs(Co) > One) Co = sign(One,Co)
  Fir = ArCos(Co)
  Si = sqrt(One-Co**2)
end if
!
if (Fir < 1.0d-13) then
  Fir = Zero
  return
else if (abs(Fir-Pi) < 1.0d-13) then
  Fir = Pi
  return
end if
dFir = 180.0d0*Fir/Pi
if (((abs(dFir) > 177.5) .or. (abs(dFir) < 2.5)) .and. lWarn) write(6,*) ' Valence angle close to end in range of definition'
if (lWrite) write(6,'(1X,A,A,F10.4,A,F10.6,A)') Label,' : Angle=    ',dFir,'   / Degree  ',Fir,' / rad'

! Compute the WDC B-matrix

if (Si == Zero) then
  ! Dummy assignment for a linear system!
  call dcopy_(3*nCent,[0.0d0],0,Bf,1)
else
  do i=1,3
    Bf(i,1) = (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
    Bf(i,3) = (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
    ! Utilize translational invariance.
    Bf(i,2) = -(Bf(i,1)+Bf(i,3))
  end do
end if
!call RecPrt('Bf',' ',Bf,9,1)

! Compute the cartesian derivative of the B-Matrix.

if (ldB) then

  !dBf = -11.11111
  if (Si == Zero) then
    call WarningMessage(2,'Bend: Si == 0.0D')
    call Abend()
  end if
  do i=1,3
    do j=1,i
      dBf(i,1,j,1) = (-Si*Bf(i,1)*BRij(j,1)+Co*dBRij(i,1,j,1)-Bf(j,1)*(Co*Bf(i,1)*Rij1+Si*BRij(i,1)))/(Si*Rij1)
      dBf(i,1,j,3) = (-Si*Bf(i,1)*BRjk(j,2)+dBRij(i,1,j,2)-Bf(j,3)*Co*Bf(i,1)*Rjk1)/(Si*Rjk1)
      !write(6,*) '13',dBf(i,1,j,3), i, j
      dBf(i,3,j,1) = (-Si*Bf(i,3)*BRij(j,1)+dBRjk(i,2,j,1)-Bf(j,1)*Co*Bf(i,3)*Rij1)/(Si*Rij1)
      dBf(i,3,j,3) = (-Si*Bf(i,3)*BRjk(j,2)+Co*dBRjk(i,2,j,2)-Bf(j,3)*(Co*Bf(i,3)*Rjk1+Si*BRjk(i,2)))/(Si*Rjk1)

      dBf(j,1,i,1) = dBf(i,1,j,1)
      dBf(j,3,i,1) = dBf(i,1,j,3)
      dBf(j,1,i,3) = dBf(i,3,j,1)
      dBf(j,3,i,3) = dBf(i,3,j,3)

      dBf(i,1,j,2) = -(dBf(i,1,j,1)+dBf(i,1,j,3))
      dBf(j,2,i,1) = dBf(i,1,j,2)
      dBf(j,1,i,2) = -(dBf(j,1,i,1)+dBf(j,1,i,3))
      dBf(i,2,j,1) = dBf(j,1,i,2)
      dBf(i,3,j,2) = -(dBf(i,3,j,1)+dBf(i,3,j,3))
      dBf(j,2,i,3) = dBf(i,3,j,2)
      dBf(j,3,i,2) = -(dBf(j,3,i,1)+dBf(j,3,i,3))
      dBf(i,2,j,3) = dBf(j,3,i,2)

      dBf(i,2,j,2) = -(dBf(i,2,j,1)+dBf(i,2,j,3))
      dBf(j,2,i,2) = dBf(i,2,j,2)

    end do
  end do
  !call RecPrt('dBf','(9F9.1)',dBf,9,9)

end if

return

end subroutine Bend

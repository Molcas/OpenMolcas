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

use Constants, only: Zero, One, Pi, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(in) :: xyz(3,nCent)
real(kind=wp), intent(out) :: Fir, Bf(3,nCent)
logical(kind=iwp), intent(in) :: lWrite, lWarn, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBf(3,nCent,3,nCent)
integer(kind=iwp) :: i, j, mCent
real(kind=wp) :: BRij(3,2), BRjk(3,2), Co, Crap, dBRij(3,2,3,2), dBRjk(3,2,3,2), dFir, Rij1, Rjk1, Si
real(kind=wp), external :: ArCos, ArSin

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

if (Crap < 1.0e-4_wp) then
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

if (Fir < 1.0e-13_wp) then
  Fir = Zero
  return
else if (abs(Fir-Pi) < 1.0e-13_wp) then
  Fir = Pi
  return
end if
dFir = Fir/deg2rad
if (((abs(dFir) > 177.5_wp) .or. (abs(dFir) < 2.5_wp)) .and. lWarn) write(u6,*) ' Valence angle close to end in range of definition'
if (lWrite) write(u6,'(1X,A,A,F10.4,A,F10.6,A)') Label,' : Angle=    ',dFir,'   / Degree  ',Fir,' / rad'

! Compute the WDC B-matrix

if (Si == Zero) then
  ! Dummy assignment for a linear system!
  Bf(:,:) = Zero
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

  !dBf = -11.11111_wp
  if (Si == Zero) then
    call WarningMessage(2,'Bend: Si == 0.0')
    call Abend()
  end if
  do i=1,3
    do j=1,i
      dBf(i,1,j,1) = (-Si*Bf(i,1)*BRij(j,1)+Co*dBRij(i,1,j,1)-Bf(j,1)*(Co*Bf(i,1)*Rij1+Si*BRij(i,1)))/(Si*Rij1)
      dBf(i,1,j,3) = (-Si*Bf(i,1)*BRjk(j,2)+dBRij(i,1,j,2)-Bf(j,3)*Co*Bf(i,1)*Rjk1)/(Si*Rjk1)
      !write(u6,*) '13',dBf(i,1,j,3), i, j
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

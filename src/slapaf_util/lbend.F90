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

subroutine LBend(Cent,nCent,Fir,Bf,lWrite,Label,dBf,ldB,Axis,Perp_Axis1,Force)

use Constants, only: Zero, One, Two, Pi, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Cent(3,3), Axis(3), Perp_Axis1(3)
integer(kind=iwp), intent(in) :: nCent
real(kind=wp), intent(out) :: Fir, Bf(3,nCent), dBf(3,nCent,3,nCent)
logical(kind=iwp), intent(in) :: lWrite, ldB, Force
character(len=8), intent(in) :: Label
#include "print.fh"
integer(kind=iwp) :: i, iPrint, iRout, j, Lu, mCent, Middle
real(kind=wp) :: Bfi1, Bfi3, Bfj1, Bfj3, BRij(3,2), BRjk(3,2), Co, Crap, dBRij(3,2,3,2), dBRjk(3,2,3,2), dFir, R1, R2, R3, Rij1, &
                 Rjk1, Scr1(3,3), Scr2(3,3), Si, uMtrx(3,3), uVec(3,3), xxx(3,3)
logical(kind=iwp) :: Linear
real(kind=wp), external :: ArCos, ArSin

iRout = 220
iPrint = nPrint(iRout)

Lu = u6

if (iPrint >= 99) then
  write(u6,*) 'LBend: Force ',Force
  call RecPrt('LBend: Axis',' ',Axis,3,1)
  call RecPrt('LBend: Perp_Axis1',' ',Perp_Axis1,3,1)
end if

uVec(:,1) = Axis(:)
uVec(:,2) = Perp_Axis1(:)
uVec(:,3) = Zero

! Project the coordinates to the plane

call DGEMM_('T','N',3,3,3,One,uVec,3,Cent,3,Zero,xxx,3)
xxx(3,1) = Zero
xxx(3,2) = Zero
xxx(3,3) = Zero
if (iPrint >= 99) then
  call RecPrt('Original coordinates','(3F24.12)',Cent,3,3)
  call RecPrt('uVec',' ',uVec,3,3)
  call RecPrt('Projected coordinates','(3F24.12)',xxx,3,3)
end if

! Swap atoms to ensure the complementary angle is always Pi

Middle = 2
if (Force) then
  R1 = (xxx(1,1)-xxx(1,2))**2+(xxx(2,1)-xxx(2,2))**2
  R2 = (xxx(1,2)-xxx(1,3))**2+(xxx(2,2)-xxx(2,3))**2
  R3 = (xxx(1,3)-xxx(1,1))**2+(xxx(2,3)-xxx(2,1))**2
  if ((R1 >= R3) .and. (R1 >= R2)) then
    Middle = 3
  else if (R2 >= R3) then
    Middle = 1
  end if
end if
if (Middle /= 2) then
  call DSwap_(3,xxx(1,2),1,xxx(1,Middle),1)
  if (iPrint >= 99) call RecPrt('Swapped coordinates','(3F24.12)',xxx,3,3)
end if

mCent = 2
call Strtch(xxx(1,1),mCent,Rij1,BRij,.false.,Label,dBRij,ldB)
call Strtch(xxx(1,2),mCent,Rjk1,BRjk,.false.,Label,dBRjk,ldB)

! We better be very careful here in order not to lose accuracy!

Co = Zero
Crap = Zero
do i=1,3
  Co = Co+BRij(i,1)*BRjk(i,2)
end do
do i=1,3
  Crap = Crap+(BRjk(i,2)-sign(One,Co)*BRij(i,1))**2
end do
Crap = sqrt(Crap)
Linear = .true.
if (iPrint >= 99) then
  call RecPrt('BRij','(3F24.12)',BRij,3,2)
  call RecPrt('BRjk','(3F24.12)',BRjk,3,2)
  write(u6,*) ' Rij1=',Rij1
  write(u6,*) ' Rjk1=',Rjk1
  write(u6,*) ' Diff=',abs(ArCos(Co)-Pi)
  write(u6,'(A,F24.16)') ' Co=',Co
  write(u6,'(A,F24.16)') ' Crap=',Crap
end if

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

!if (abs(Fir-Pi) > 1.0e-13_wp) then
if (abs(Si) > 1.0e-13_wp) then
  if (iPrint >= 99) write(Lu,*) ' LBend: Use nonlinear formulae'
  Linear = .false.
else
  if (iPrint >= 99) write(Lu,*) ' LBend: Use linear formulae'
end if

dFir = Fir/deg2rad
if (lWrite) write(Lu,'(1X,A,A,F10.6,A,F12.8,A)') Label,' : Projected Angle=',dFir,'/degree, ',Fir,'/rad'

uMtrx(:,:) = Zero
if (Linear) then
  uMtrx(2,1) = -sign(One,Co)/Rij1
else
  do i=1,3
    uMtrx(i,1) = (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
  end do
end if
call DGEMM_('N','N',3,1,3,One,uVec,3,uMtrx,3,Zero,Scr2,3)
Bf(1,1) = Scr2(1,1)
Bf(2,1) = Scr2(2,1)
Bf(3,1) = Scr2(3,1)

uMtrx(:,:) = Zero
if (Linear) then
  uMtrx(2,1) = One/Rjk1
else
  uMtrx(:,1) = (Co*BRjk(:,2)-BRij(:,1))/(Si*Rjk1)
end if
call DGEMM_('N','N',3,1,3,One,uVec,3,uMtrx,3,Zero,Scr2,3)
Bf(1,3) = Scr2(1,1)
Bf(2,3) = Scr2(2,1)
Bf(3,3) = Scr2(3,1)

! Utilize translational invariance.
Bf(:,2) = -(Bf(:,1)+Bf(:,3))

! Compute the cartesian derivative of the B-Matrix.

if (ldB) then

  ! 1,1 Block

  uMtrx(:,:) = Zero
  if (Linear) then
    if (Co > Zero) then
      uMtrx(1,2) = One/Rij1**2
      if (Rij1 < Rjk1) uMtrx(1,2) = Two*uMtrx(1,2)
    else
      uMtrx(1,2) = Two/Rij1**2-One/(Rij1**2+Rij1*Rjk1)
    end if
    uMtrx(2,1) = uMtrx(1,2)
  else
    do i=1,2
      Bfi1 = (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
      do j=1,2
        Bfj1 = (Co*BRij(j,1)-BRjk(j,2))/(Si*Rij1)
        uMtrx(i,j) = (-Si*Bfi1*BRij(j,1)+Co*dBRij(i,1,j,1)-Bfj1*(Co*Bfi1*Rij1+Si*BRij(i,1)))/(Si*Rij1)
      end do
    end do
  end if
  call DGEMM_('N','T',3,3,3,One,uMtrx,3,uVec,3,Zero,Scr1,3)
  call DGEMM_('N','N',3,3,3,One,uVec,3,Scr1,3,Zero,Scr2,3)
  dBf(:,1,:,1) = Scr2(:,:)

  ! 1,3 Block

  uMtrx(:,:) = Zero
  if (Linear) then
    if (Co > Zero) then
      if (Rij1 < Rjk1) then
        uMtrx(1,2) = -One/(Rij1*Rjk1)
        uMtrx(2,1) = Zero
      else
        uMtrx(1,2) = Zero
        uMtrx(2,1) = One/(Rij1*Rjk1)
      end if
    else
      uMtrx(1,2) = One/(Rij1**2+Rij1*Rjk1)
      uMtrx(2,1) = -One/(Rjk1**2+Rjk1*Rij1)
    end if
  else
    do i=1,2
      Bfi1 = (Co*BRij(i,1)-BRjk(i,2))/(Si*Rij1)
      do j=1,2
        Bfj3 = (Co*BRjk(j,2)-BRij(j,1))/(Si*Rjk1)
        uMtrx(i,j) = (-Si*Bfi1*BRjk(j,2)+dBRij(i,1,j,2)-Bfj3*Co*Bfi1*Rjk1)/(Si*Rjk1)
      end do
    end do
  end if
  call DGEMM_('N','T',3,3,3,One,uMtrx,3,uVec,3,Zero,Scr1,3)
  call DGEMM_('N','N',3,3,3,One,uVec,3,Scr1,3,Zero,Scr2,3)
  dBf(:,1,:,3) = Scr2(:,:)

  ! 3,1 Block

  uMtrx(:,:) = Zero
  if (Linear) then
    if (Co > Zero) then
      if (Rjk1 < Rij1) then
        uMtrx(1,2) = One/(Rjk1*Rij1)
        uMtrx(2,1) = Zero
      else
        uMtrx(1,2) = Zero
        uMtrx(2,1) = -One/(Rjk1*Rij1)
      end if
    else
      uMtrx(1,2) = -One/(Rjk1**2+Rjk1*Rij1)
      uMtrx(2,1) = One/(Rij1**2+Rij1*Rjk1)
    end if
  else
    do i=1,2
      Bfi3 = (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
      do j=1,2
        Bfj1 = (Co*BRij(j,1)-BRjk(j,2))/(Si*Rij1)
        uMtrx(i,j) = (-Si*Bfi3*BRij(j,1)+dBRjk(i,2,j,1)-Bfj1*Co*Bfi3*Rij1)/(Si*Rij1)
      end do
    end do
  end if
  call DGEMM_('N','T',3,3,3,One,uMtrx,3,uVec,3,Zero,Scr1,3)
  call DGEMM_('N','N',3,3,3,One,uVec,3,Scr1,3,Zero,Scr2,3)
  dBf(:,3,:,1) = Scr2(:,:)

  ! 3,3 Block

  uMtrx(:,:) = Zero
  if (Linear) then
    if (Co > Zero) then
      uMtrx(1,2) = -One/Rjk1**2
      if (Rjk1 < Rij1) uMtrx(1,2) = Two*uMtrx(1,2)
    else
      uMtrx(1,2) = -Two/Rjk1**2+One/(Rjk1**2+Rjk1*Rij1)
    end if
    uMtrx(2,1) = uMtrx(1,2)
  else
    do i=1,2
      Bfi3 = (Co*BRjk(i,2)-BRij(i,1))/(Si*Rjk1)
      do j=1,2
        Bfj3 = (Co*BRjk(j,2)-BRij(j,1))/(Si*Rjk1)
        uMtrx(i,j) = (-Si*Bfi3*BRjk(j,2)+Co*dBRjk(i,2,j,2)-Bfj3*(Co*Bfi3*Rjk1+Si*BRjk(i,2)))/(Si*Rjk1)
      end do
    end do
  end if
  call DGEMM_('N','T',3,3,3,One,uMtrx,3,uVec,3,Zero,Scr1,3)
  call DGEMM_('N','N',3,3,3,One,uVec,3,Scr1,3,Zero,Scr2,3)
  dBf(:,3,:,3) = Scr2(:,:)

  do i=1,3
    do j=1,i

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

end if

! Swap atoms back

if (Middle /= 2) then
  call DSwap_(3,Bf(:,2),1,Bf(:,Middle),1)
  if (ldB) then
    call DSwap_(3*nCent*3,dBf(:,:,:,2),1,dBf(:,:,:,Middle),1)
    call DSwap_(3*nCent,dBf(1,2,1,1),3*nCent,dBf(1,Middle,1,1),3*nCent)
    call DSwap_(3*nCent,dBf(2,2,1,1),3*nCent,dBf(2,Middle,1,1),3*nCent)
    call DSwap_(3*nCent,dBf(3,2,1,1),3*nCent,dBf(3,Middle,1,1),3*nCent)
  end if
end if

if (iPrint >= 99) then
  call RecPrt('Bf',' ',Bf,3,nCent)
  if (ldB) call RecPrt('dBf',' ',dBf,3*nCent,3*nCent)
end if

return

end subroutine LBend

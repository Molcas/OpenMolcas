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

subroutine CoSys(Cent,R,xyz)

use Constants, only: Zero, One
#ifdef _DEBUGPRINT_
use Constants, only: Pi
#endif
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Cent(3,3)
real(kind=wp), intent(out) :: r(3), xyz(3,2)
integer(kind=iwp) :: i, iComp(3), j, k, nComp
real(kind=wp) :: Co, Crap, r11, r12, r2, R21j, R21k, R23j, R23k, RR, RR1, RR2, Si
logical(kind=iwp) :: Linear, Retry
real(kind=wp), parameter :: ThrAcos = 1.0e-6_wp
#ifdef _DEBUGPRINT_
real(kind=wp) :: Fi
real(kind=wp), external :: ArCos, ArSin

call RecPrt('CoSys: Cent',' ',Cent,3,3)
#endif

! Check if linear

Co = Zero
Crap = Zero
RR1 = Zero
RR2 = Zero
do i=1,3
  Co = Co+(Cent(i,1)-Cent(i,2))*(Cent(i,3)-Cent(i,2))
  RR1 = RR1+(Cent(i,1)-Cent(i,2))**2
  RR2 = RR2+(Cent(i,3)-Cent(i,2))**2
end do
RR1 = sqrt(RR1)
RR2 = sqrt(RR2)
Co = Co/(RR1*RR2)
do i=1,3
  Crap = Crap+((Cent(i,3)-Cent(i,2))/RR2-sign(One,Co)*(Cent(i,1)-Cent(i,2))/RR1)**2
end do
Crap = sqrt(Crap)
#ifdef _DEBUGPRINT_
write(u6,*) 'Co=',Co
write(u6,*) 'Crap=',Crap
#endif
if (Crap < 1.0e-6_wp) then
  Si = Crap
# ifdef _DEBUGPRINT_
  if (Co < Zero) then
    Fi = Pi-ArSin(Si)
  else
    Fi = ArSin(Si)
  end if
# endif
else
  if ((Co > One) .and. (Co < One+ThrAcos)) Co = One
  if ((Co < -One) .and. (Co > -One-ThrAcos)) Co = -One
  if ((Co > One+ThrAcos) .or. (Co < -One-ThrAcos)) then
    call WarningMessage(2,'Error in CoSys')
    write(u6,*) 'Error in cosys: arcos(',Co,')'
    call Abend()
  end if
  Si = sqrt(One-Co**2)
# ifdef _DEBUGPRINT_
  Fi = ArCos(Co)
# endif
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'Fi,Pi=',Fi,Pi
write(u6,*) 'Pi-Fi=',Pi-Fi
#endif

Linear = abs(Si) < 1.0e-13_wp

! Form reference axis

R(:) = Cent(:,3)-Cent(:,1)
RR = sqrt(R(1)**2+R(2)**2+R(3)**2)
if ((RR1 >= RR2) .and. (RR1 >= RR)) then
  R(:) = Cent(:,1)-Cent(:,2)
  RR = RR1
else if (RR2 >= RR) then
  R(:) = Cent(:,3)-Cent(:,2)
  RR = RR2
end if
R(:) = R(:)/RR

#ifdef _DEBUGPRINT_
write(u6,*) 'Linear=',Linear
write(u6,*) 'RR=',RR
call RecPrt('R',' ',R,3,1)
#endif

Retry = .true.
do while (Retry)
  Retry = .false.
  if (Linear) then

    nComp = 0
    do i=1,3
      if (R(i) == Zero) then
        nComp = nComp+1
        iComp(nComp) = i
      end if
    end do

    ! Compute the WDC B-matrix

    xyz(:,:) = Zero
    if (nComp == 0) then
      !write(u6,*) ' Case nComp == 0'

      xyz(1,1) = R(1)
      xyz(2,1) = R(2)
      xyz(3,1) = -R(3)
      r12 = R(1)**2+R(2)**2-R(3)**2
      r11 = R(1)**2+R(2)**2+R(3)**2
      xyz(1,1) = xyz(1,1)-(r12/r11)*R(1)
      xyz(2,1) = xyz(2,1)-(r12/r11)*R(2)
      xyz(3,1) = xyz(3,1)-(r12/r11)*R(3)
      r2 = sqrt(xyz(1,1)**2+xyz(2,1)**2+xyz(3,1)**2)
      xyz(1,1) = xyz(1,1)/r2
      xyz(2,1) = xyz(2,1)/r2
      xyz(3,1) = xyz(3,1)/r2
      xyz(1,2) = R(2)*xyz(3,1)-R(3)*xyz(2,1)
      xyz(2,2) = R(3)*xyz(1,1)-R(1)*xyz(3,1)
      xyz(3,2) = R(1)*xyz(2,1)-R(2)*xyz(1,1)
      r2 = sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
      xyz(1,2) = xyz(1,2)/r2
      xyz(2,2) = xyz(2,2)/r2
      xyz(3,2) = xyz(3,2)/r2

    else if (nComp == 1) then
      !write(u6,*) ' Case nComp == 1'

      i = iComp(1)
      xyz(i,1) = One
      xyz(1,2) = R(2)*xyz(3,1)-R(3)*xyz(2,1)
      xyz(2,2) = R(3)*xyz(1,1)-R(1)*xyz(3,1)
      xyz(3,2) = R(1)*xyz(2,1)-R(2)*xyz(1,1)
      r2 = sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
      xyz(1,2) = xyz(1,2)/r2
      xyz(2,2) = xyz(2,2)/r2
      xyz(3,2) = xyz(3,2)/r2

    else if (nComp == 2) then
      !write(u6,*) ' Case nComp == 2'

      i = iComp(1)
      xyz(i,1) = One
      i = iComp(2)
      xyz(i,2) = One

    else
      call WarningMessage(2,'Error in CoSys')
      write(u6,*) ' CoSys: nComp == 3'
      call Abend()
    end if

  else     ! Non-linear

    ! Form the cross product R12xR32
    RR = Zero
    do i=1,3
      j = mod(i,3)+1
      k = mod(i+1,3)+1
      R21j = Cent(j,1)-Cent(j,2)
      R21k = Cent(k,1)-Cent(k,2)
      R23j = Cent(j,3)-Cent(j,2)
      R23k = Cent(k,3)-Cent(k,2)
      xyz(i,2) = R21j*R23k-R21k*R23j
      RR = RR+xyz(i,2)**2
    end do
    if (RR == Zero) then
      Linear = .true.
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Linear=',Linear
#     endif
      Retry = .true.
    else
      xyz(:,2) = xyz(:,2)/sqrt(RR)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'RR=',RR
#     endif

      RR = Zero
      do i=1,3
        j = mod(i,3)+1
        k = mod(i+1,3)+1
        xyz(i,1) = xyz(j,2)*R(k)-xyz(k,2)*R(j)
        RR = RR+xyz(i,1)**2
      end do
      xyz(:,1) = xyz(:,1)/sqrt(RR)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'RR=',RR
      call RecPrt('xyz',' ',xyz,3,2)
#     endif
    end if

  end if
end do

#ifdef _DEBUGPRINT_
call RecPrt(' Reference Axis',' ',R,3,1)
call RecPrt(' Perpendicular Axes',' ',xyz,3,2)
#endif

end subroutine CoSys

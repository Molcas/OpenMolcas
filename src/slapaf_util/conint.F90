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

subroutine ConInt(xyz,nCent,dE,Bf,lWrite_,Label,dBf,ldB,lIter)

use Slapaf_Info, only: ApproxNADC, Energy, Energy0, Gx0, NADC
use Constants, only: Zero, One, Two, auTokJmol
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nCent, lIter
real(kind=wp), intent(in) :: xyz(3,nCent)
real(kind=wp), intent(out) :: dE, Bf(3,nCent)
logical(kind=iwp), intent(in) :: lWrite_, ldB
character(len=8), intent(in) :: Label
real(kind=wp), intent(inout) :: dBf(3*nCent,3*nCent)
integer(kind=iwp) :: i, iCar, iCent, iOpt, ix, iy, j, jCent
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nAtoms
#endif
real(kind=wp) :: E0, E1, Fact, XX
integer(kind=iwp), external :: iDeg
real(kind=wp), external :: DDot_

E1 = Energy(lIter)
E0 = Energy0(lIter)
#ifdef _DEBUGPRINT_
nAtoms = size(Gx0,2)
write(u6,*) 'ConInt: lIter=',lIter
write(u6,*) 'ConInt: E1, E0=',E1,E0
call RecPrt('ConInt: Gx0',' ',Gx0(:,:,lIter),3,nAtoms)
#endif

! iOpt=1 -> Linear
! iOpt=2 -> Quadratic
! iOpt=3 -> Absolute value
if (NADC) then
  if (ApproxNADC) then
    iOpt = 2
  else
    iOpt = 3
  end if
else
  iOpt = 1
end if

! Observe that the program is storing the forces rather than the
! gradients!
!
! The average energy is stored in E1 and the energy difference in
! E0. Ditto for the gradients (see process_gradients).

dE = Zero
if (iOpt == 1) then
  ! Linear
  dE = E0
else if (iOpt == 2) then
  ! Quadratic
  dE = E0**2
else if (iOpt == 3) then
  ! Absolute value
  dE = abs(E0)
end if
if (lWrite_) then
  write(u6,'(2A,F18.8,A,F18.8,A)') Label,' : Energy difference = ',E0,' hartree, ',E0*auTokJmol,' kJ/mol'
  write(u6,'( A,F18.8,A)') '           Average energy    = ',E1,' hartree'
# ifdef _DEBUGPRINT_
  select case (iOpt)
    case (1)
      write(u6,*) 'Option: Linear'
    case (2)
      write(u6,*) 'Option: Quadratic'
    case (3)
      write(u6,*) 'Option: Absolute value'
  end select
# endif
end if

! Compute the WDC B-matrix

Bf(:,:) = Zero
do iCent=1,nCent
  Fact = real(iDeg(xyz(1,iCent)),kind=wp)
  !write(u6,*) 'Fact=',Fact
  do iCar=1,3
    if (iOpt == 1) then
      ! Linear
      Bf(iCar,iCent) = -Gx0(iCar,iCent,lIter)
    else if (iOpt == 2) then
      ! Quadratic

      ! When the energy difference becomes small, the true derivative vanishes.
      ! In such case use simply a scaled-down energy difference gradient

      if (abs(E0) > 1.0e-5_wp) then
        Bf(iCar,iCent) = -Two*E0*Gx0(iCar,iCent,lIter)
      else
        Bf(iCar,iCent) = -Two*1.0e-5_wp*Gx0(iCar,iCent,lIter)
      end if
    else if (iOpt == 3) then
      ! Absolute value
      Bf(iCar,iCent) = -sign(One,E0)*Gx0(iCar,iCent,lIter)
    end if

    Bf(iCar,iCent) = Fact*Bf(iCar,iCent)
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('ConInt: Bf',' ',Bf,3,nCent)
#endif
if (lWrite_ .and. (iOpt == 1)) then
  XX = sqrt(DDot_(3*nCent,Bf,1,Bf,1))
  if (XX <= 1.0e-3_wp) then
    write(u6,*)
    write(u6,*) '    Warning: PESs might be parallel!'
    write(u6,*)
  end if
end if

! Compute the cartesian derivative of the B-Matrix.

if (ldB) then
  dBf(:,:) = Zero
  if (iOpt == 1) then
    ! Linear

  else if (iOpt == 2) then
    ! Quadratic

    ix = 0
    do iCent=1,nCent
      do i=1,3
        ix = ix+1

        iy = 0
        do jCent=1,nCent
          do j=1,3
            iy = iy+1
            dBf(ix,iy) = -Two*Gx0(i,iCent,lIter)*Gx0(j,jCent,lIter)
          end do
        end do
      end do
    end do

  else if (iOpt == 3) then
    ! Absolute value

  end if
  !call RecPrt('dBf','(9F9.1)',dBf,3*nCent,3*nCent)

end if

return

end subroutine ConInt

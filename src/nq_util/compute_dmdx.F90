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

subroutine Compute_dMdx(ZA,RA,nAtoms,T,iAtom,iCar,dTdRAi,dMdx)

#ifdef _DEBUGPRINT_
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, iAtom, iCar
real(kind=wp), intent(in) :: ZA(nAtoms), RA(3,nAtoms), T(3), dTdRAi
real(kind=wp), intent(out) :: dMdx(3,3)
integer(kind=iwp) :: i, j, jAtom
real(kind=wp) :: RTx, RTy, RTz, tmp, ZB
#ifdef _DEBUGPRINT_
real(kind=wp) :: delta, M(3,3)
real(kind=wp), allocatable :: dRA(:,:)
#endif
real(kind=wp), parameter :: Thrs = 1.0e-14_wp

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
delta = 1.0e-4_wp
call mma_allocate(dRA,3,nAtoms,label='dRA')
dRA(:,:) = RA

dRA(iCar,iAtom) = RA(iCar,iAtom)+Delta
call Compute_M(ZA,nAtoms,dRA,T,M)

dRA(iCar,iAtom) = RA(iCar,iAtom)-Delta
call Compute_M(ZA,nAtoms,dRA,T,dMdx)

dRA(iCar,iAtom) = RA(iCar,iAtom)

dMdx(:,:) = (M-dMdx)/(Two*Delta)
call RecPrt('dMdx(Numerical)',' ',dMdx,3,3)
#endif

dMdx(:,:) = Zero
do jAtom=1,nAtoms
  ZB = ZA(jAtom)
  if (iAtom == jAtom) then
    tmp = (One-dTdRAi)*ZB
  else
    tmp = (-dTdRAi)*ZB
  end if

  RTx = RA(1,jAtom)-T(1)
  RTy = RA(2,jAtom)-T(2)
  RTz = RA(3,jAtom)-T(3)
  select case (iCar)
    case (1)
      dMdx(2,2) = dMdx(2,2)+Two*tmp*RTx
      dMdx(3,3) = dMdx(3,3)+Two*tmp*RTx
      dMdx(1,2) = dMdx(1,2)-tmp*RTy
      dMdx(2,1) = dMdx(2,1)-RTy*tmp
      dMdx(1,3) = dMdx(1,3)-tmp*RTz
      dMdx(3,1) = dMdx(3,1)-RTz*tmp
    case (2)
      dMdx(1,1) = dMdx(1,1)+Two*tmp*RTy
      dMdx(3,3) = dMdx(3,3)+Two*tmp*RTy
      dMdx(1,2) = dMdx(1,2)-RTx*tmp
      dMdx(2,1) = dMdx(2,1)-tmp*RTx
      dMdx(2,3) = dMdx(2,3)-tmp*RTz
      dMdx(3,2) = dMdx(3,2)-RTz*tmp
    case (3)
      dMdx(1,1) = dMdx(1,1)+Two*tmp*RTz
      dMdx(2,2) = dMdx(2,2)+Two*tmp*RTz
      dMdx(1,3) = dMdx(1,3)-RTx*tmp
      dMdx(3,1) = dMdx(3,1)-tmp*RTx
      dMdx(2,3) = dMdx(2,3)-RTy*tmp
      dMdx(3,2) = dMdx(3,2)-tmp*RTy
  end select
end do

! Remove noise

do i=1,3
  do j=1,3
    if (abs(dMdx(i,j)) < Thrs) dMdx(i,j) = Zero
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('dMdx',' ',dMdx,3,3)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Compute_dMdx

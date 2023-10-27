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
! Copyright (C) 2003, Roland Lindh                                     *
!***********************************************************************

subroutine RotGrd(RA,ZA,O,dOdx,d2Odx2,nAtoms,Do_Grad,Do_Hess)
!***********************************************************************
!                                                                      *
!     Object: Compute the principal axes system and to optionally      *
!             evaluate the gradient of the principal axes system.      *
!                                                                      *
!             See: B. G. Johnson et al., CPL, 220, 377 (1994).         *
!                                                                      *
!     Author: R. Lindh, Dept. of Chem. Phys., Univ. of Lund, Sweden.   *
!                                                                      *
!             Created on board M/S Polarlys on voyage from Tromso to   *
!             Trondheim, Sept. 2003.                                   *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(in) :: RA(3,nAtoms), ZA(nAtoms)
real(kind=wp), intent(out) :: O(3,3), dOdx(3,3,nAtoms,3)
real(kind=wp), intent(inout) :: d2Odx2(3,3,nAtoms,3,nAtoms,3)
logical(kind=iwp), intent(in) :: Do_Grad, Do_Hess
integer(kind=iwp) :: iAtom, iCar, jAtom, jCar, jCar_Max
real(kind=wp) :: dMdx(3,3), dMdy(3,3), dTdRAi, dTdRAj, EVal(3), Px(3,3), Py(3,3), T(3), Z_Tot
logical(kind=iwp) :: Rot_Corr
real(kind=wp), parameter :: Thrs = 1.0e-3_wp

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('RotGrd: RA',' ',RA,3,nAtoms)
call RecPrt('RotGrd: ZA',' ',ZA,1,nAtoms)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the total charge

Z_Tot = sum(ZA)

! Form the center of nuclear charge

call Compute_T(Z_Tot,T,ZA,RA,nAtoms)

! Form O

call Compute_O(ZA,RA,nAtoms,T,O,EVal)
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Do_Grad) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Turn off rotational correction to the gradient if the eigenvectors
! are close to degeneracy!

Rot_Corr = .true.
if (abs(Eval(1)-Eval(2))/(EVal(1)+EVal(2)) < Thrs) then
  write(u6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
  Rot_Corr = .false.
end if
if (abs(Eval(1)-Eval(3))/(EVal(1)+EVal(3)) < Thrs) then
  write(u6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
  Rot_Corr = .false.
end if
if (abs(Eval(2)-Eval(3))/(EVal(2)+EVal(3)) < Thrs) then
  write(u6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
  Rot_Corr = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the gradient

do iAtom=1,nAtoms
  dTdRAi = ZA(iAtom)/Z_Tot
  do iCar=1,3

    ! Form dO/dx

    call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx(:,:,iAtom,iCar),Px)

  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Do_Hess) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the Hessian

do iAtom=1,nAtoms
  dTdRAi = ZA(iAtom)/Z_Tot
  do iCar=1,3
    call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx(:,:,iAtom,iCar),Px)

    do jAtom=1,iAtom
      dTdRAj = ZA(jAtom)/Z_Tot
      jCar_Max = 3
      if (iAtom == jAtom) jCar_Max = iCar
      do jCar=1,jCar_Max
        call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,jAtom,jCar,dTdRAj,dMdy,dOdx(:,:,jAtom,jCar),Py)

        ! Form d2O/dx2

        call Compute_d2Odx2(ZA,nAtoms,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,Px,jAtom,jCar,dMdy,Py, &
                            d2Odx2(:,:,iAtom,iCar,jAtom,jCar))

        if ((iAtom /= jAtom) .or. (iCar /= jCar)) d2Odx2(:,:,jAtom,jCar,iAtom,iCar) = d2Odx2(:,:,iAtom,iCar,jAtom,jCar)

      end do ! jCar
    end do   ! iCar

  end do     ! jAtom
end do       ! iAtom
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine RotGrd

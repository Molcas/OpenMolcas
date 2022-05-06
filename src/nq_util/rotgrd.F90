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

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 RA(3,nAtoms), ZA(nAtoms), O(3,3), dOdx(3,3,nAtoms,3), d2Odx2(3,3,nAtoms,3,nAtoms,3), dMdx(3,3), dMdy(3,3), EVal(3), T(3), &
       Px(3,3), Py(3,3)
logical Do_Grad, Rot_Corr, Do_Hess

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('RotGrd: RA',' ',RA,3,nAtoms)
call RecPrt('RotGrd: ZA',' ',ZA,1,nAtoms)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the total charge

Z_Tot = DDot_(nAtoms,[One],0,ZA,1)

! Form the center of nuclear charge

call Compute_T(Z_Tot,T,ZA,RA,nAtoms)

! Form O

call Compute_O(ZA,RA,nAtoms,Z_Tot,T,O,EVal)
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
if (abs(Eval(1)-Eval(2))/(EVal(1)+EVal(2)) < 0.001d0) then
  write(6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
  Rot_Corr = .false.
end if
if (abs(Eval(1)-Eval(3))/(EVal(1)+EVal(3)) < 0.001d0) then
  write(6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
  Rot_Corr = .false.
end if
if (abs(Eval(2)-Eval(3))/(EVal(2)+EVal(3)) < 0.001d0) then
  write(6,*) 'Rotational correction to the DFT gradient is turned off due to close-to-degeneracy problems!'
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

    call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx(1,1,iAtom,iCar),Px)

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
    call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx(1,1,iAtom,iCar),Px)

    do jAtom=1,iAtom
      dTdRAj = ZA(jAtom)/Z_Tot
      jCar_Max = 3
      if (iAtom == jAtom) jCar_Max = iCar
      do jCar=1,jCar_Max
        call Compute_dOdx(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,jAtom,jCar,dTdRAj,dMdy,dOdx(1,1,jAtom,jCar),Py)

        ! Form d2O/dx2

        call Compute_d2Odx2(ZA,RA,nAtoms,T,O,EVal,Rot_Corr,iAtom,iCar,dTdRAi,dMdx,dOdx(1,1,iAtom,iCar),Px,jAtom,jCar,dTdRAj,dMdy, &
                            dOdx(1,1,jAtom,jCar),Py,d2Odx2(1,1,iAtom,iCar,jAtom,jCar))

        if ((iAtom /= jAtom) .or. ((iAtom == jAtom) .and. (iCar /= jCar))) then
          call dcopy_(9,d2Odx2(1,1,iAtom,iCar,jAtom,jCar),1,d2Odx2(1,1,jAtom,jCar,iAtom,iCar),1)
        end if

      end do ! jCar
    end do   ! iCar

  end do     ! jAtom
end do       ! iAtom
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine RotGrd

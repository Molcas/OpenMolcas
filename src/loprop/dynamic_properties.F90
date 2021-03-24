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

subroutine Dynamic_Properties(Temp,nAtoms,rMP,nij,nPert,nElem,Delta,EC,Polar,iANr,Bond_Threshold,ChPol,ChPolBB)

use Constants, only: Zero, One, Two, Three, Eight, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nij, nPert, nElem, iANr(nAtoms)
real(kind=wp), intent(out) :: Temp(nij), Polar(6,nij), ChPol(6,nij), ChPolBB(6,nij)
real(kind=wp), intent(in) :: rMP(nij,0:nElem-1,0:nPert-1), Delta, EC(3,nij), Bond_Threshold
integer(kind=iwp) :: iAtom, iCar, ii, ij, iPert, iPert_, iPol, jAtom, jCar, jj, jPert, jPert_
real(kind=wp) :: A(3), B(3), Pol1, Pol1a, Pol1b, Pol2, Rij_iCar

!                                                                      *
!***********************************************************************
!                                                                      *
!call RecPrt('rMP',' ',rMP,nij*nElem,nPert)
write(u6,*)
write(u6,*) ' D y n a m i c  P r o p e r t i e s'
write(u6,*)
write(u6,*) ' Properties computed with FFPT'
write(u6,*)

do iPol=1,6
  do iAtom=1,nAtoms
    do jAtom=1,iAtom
      ij = iAtom*(iAtom-1)/2+jAtom
      ChPol(iPol,ij) = Zero
      ChPolBB(iPol,ij) = Zero
    end do
  end do
end do

! iPol: index vector for polarizability
!       (1,2,3,4,5,6)=(xx,yx,yy,zx,zy,zz)

do iPol=1,6
  Temp(:) = Zero
  !write (u6,*)
  do iAtom=1,nAtoms
    ii = iAtom*(iAtom+1)/2
    call dcopy_(3,EC(1,ii),1,A,1)
    do jAtom=1,iAtom
      jj = jAtom*(jAtom+1)/2
      call dcopy_(3,EC(1,jj),1,B,1)

      ij = iAtom*(iAtom-1)/2+jAtom
      !                                                                *
      !*****************************************************************
      !                                                                *
      !        Polarizabilities: alpha(iAtom,jAtom,iCar,jCar)          *
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! iCar, jCar: index of cartesian for each perturbation
      !           1=x, 2=y, 3=z
      iCar = int((One+sqrt(Eight*real(iPol,kind=wp)-Three))/Two)
      jCar = iPol-iCar*(iCar-1)/2

      ! iPert, jPert: index vector to actuall perturbation
      !       (1,2,3,4,5,6)=(+dx,-dx,+dy,-dy,+dz,-dz)
      iPert = (jCar-1)*2+1
      jPert = iPert+1
      !
      ! Contribution due to change of localized dipole moment
      !
      !------- mu(iAtom,jAtom,iCar,+dF(jCar))
      !------- mu(iAtom,jAtom,iCar,-dF(jCar))
      !
      Pol1a = (rMP(ij,iCar,iPert)-rMP(ij,iCar,jPert))/(Two*Delta)
      !
      iPert_ = (iCar-1)*2+1
      jPert_ = iPert_+1
      Pol1b = (rMP(ij,jCar,iPert_)-rMP(ij,jCar,jPert_))/(Two*Delta)
      Pol1 = Half*(Pol1a+Pol1b)
      !write(u6,*) 'Pol1',ij,iCar,iPert,jPert
      !write(u6,*) rMP(ij,iCar,iPert),rMP(ij,iCar,jPert)
      !
      ! Contribution due to change of localized charges
      !
      if (iAtom /= jAtom) then
        Rij_iCar = B(iCar)-A(iCar)
        !write(u6,*) rMP(ij,0,iPert),rMP(ij,0,jPert)
        Pol2 = (rMP(ij,0,iPert)-rMP(ij,0,jPert))*Rij_iCar/(Two*Delta)
      else
        Pol2 = Zero
      end if

      !write(u6,*) Pol1, Pol2

      Temp(ij) = Temp(ij)+Pol1+Pol2
      Polar(iPol,ij) = Temp(ij)
      !                                                                *
      !*****************************************************************
      !        Compute Charge contrib                                  *
      !*****************************************************************
      !                                                                *
      ChPol(iPol,ij) = ChPol(iPol,ij)+Pol2
      ChPolBB(iPol,ij) = ChPolBB(iPol,ij)+Pol2
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do   ! jAtom
  end do     ! iAtom
!                                                                      *
!***********************************************************************
!                                                                      *
end do   ! iPol
!                                                                      *
!***********************************************************************
!                                                                      *
!            Move the polarizabilities to the atoms if needed          *
!                                                                      *
!***********************************************************************
!                                                                      *
call Move_Polar(Polar,EC,nAtoms,nij,iANr,Bond_Threshold)
call Move_Polar(ChPol,EC,nAtoms,nij,iANr,Bond_Threshold)

return

end subroutine Dynamic_Properties

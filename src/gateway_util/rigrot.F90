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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine RigRot(CoorIn,rM,nAtm)
!***********************************************************************
!                                                                      *
! Object: to compute the rotational constants and the rotational       *
!         spectrum within the rigid-rotor model.                       *
!                                                                      *
!    Reference: P.W. Atkins, in "Molecular Quantum Mechanics",         *
!               Oxford University Press, Oxford/New York,              *
!               ch. 11.2, pp. 289.                                     *
!               Michael Tinkham, in "Group Theory and Quantum          *
!               Mechanics", McGraw-Hill Book Company, New York,        *
!               ch. 7-14, pp. 250.                                     *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Sizes_of_Seward, only: S
use Gateway_Info, only: CoM, PAX, Prin, rMI, TMass
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Eight, Half, auTocm, auToHz, uToau
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtm
real(kind=wp), intent(in) :: CoorIn(3,nAtm), rM(nAtm)
#include "print.fh"
integer(kind=iwp) :: i, iAtom, iCar, iEn, ii, iPrint, iRout, j, jCar, k, k1, k2, kappa, kk, kk2, mDim, nEn, nHess, nTri
real(kind=wp) :: A, B, C, keep, rKappa, XI(3)
logical(kind=iwp) :: Linear, RR_Show
character(len=80) :: Label
real(kind=wp), allocatable :: Coor(:,:), En(:), Hess(:), Vec(:)
integer(kind=iwp), external :: iprintlevel

!#define _DEBUGPRINT_
iRout = 117
iPrint = nPrint(iRout)
RR_Show = iPrint >= 6
if (iprintlevel(-1) < 3) RR_Show = .false.
#ifdef _DEBUGPRINT_
iPrint = 99
RR_Show = .true.
#endif

if (RR_Show) then
  write(u6,*)
  call CollapseOutput(1,'   Rigid rotor info:')
  write(u6,'(3X,A)') '   -----------------'
  write(u6,*)
end if

PAx(:,:) = reshape([One,Zero,Zero,Zero,One,Zero,Zero,Zero,One],[3,3])
if (TMass == Zero) then
  call FinishUp()
  return
end if
Linear = .false.
if (iPrint >= 99) then
  call RecPrt(' In RigRot: CoorIn',' ',CoorIn,3,nAtm)
  call RecPrt(' In RigRot: Mass',' ',rM,1,nAtm)
end if
if (RR_Show) then
  write(u6,*)
  write(u6,'(19X,A,F10.5)') ' Total mass (a) :',TMass/uToau
  write(u6,*)
  write(u6,'(19X,A)') ' Center of mass '
  write(u6,'(19X,3A11)') '       X   ','       Y          Z   '
  write(u6,'(19X,3F11.5)') (CoM(i),i=1,3)
  write(u6,*)
end if

! Translate coordinate system to center of mass.

call mma_allocate(Coor,3,nAtm)
do iCar=1,3
  do iAtom=1,nAtm
    Coor(iCar,iAtom) = CoorIn(iCar,iAtom)-CoM(iCar)
  end do
end do

if (RR_Show) then
  write(u6,'(19X,A)') ' Reference system based on center of mass'
  write(u6,'(19X,A)') ' Coordinates and Masses of Atoms, in au and A'
  write(u6,'(19X,A)') '       X          Y          Z        Mass'
  write(u6,'(19X,4F11.5)') ((Coor(j,i),j=1,3),rM(i)/utoau,i=1,nAtm)
  write(u6,*)
end if
if (nAtm == 1) then
  call mma_deallocate(Coor)
  call FinishUp()
  return
end if

! Construct the moment of inertia tensor

ii = 0
do iCar=1,3
  do jCar=1,iCar-1
    ii = ii+1
    rMI(ii) = Zero
    do iAtom=1,nAtm
      rMI(ii) = rMI(ii)-Coor(iCar,iAtom)*Coor(jCar,iAtom)*rM(iAtom)
    end do
  end do
  ii = ii+1
  rMI(ii) = Zero
  do iAtom=1,nAtm
    rMI(ii) = rMI(ii)+rM(iAtom)*(Coor(1,iAtom)**2+Coor(2,iAtom)**2+Coor(3,iAtom)**2-Coor(iCar,iAtom)**2)
  end do
end do
call mma_deallocate(Coor)
if (RR_Show) then
  write(u6,'(19X,A)') ' The Moment of Inertia Tensor / au'
  write(u6,'(19X,14X,3A)') '    X     ','     Y    ','    Z     '
  write(u6,'(19X,A,12X,3(E11.4))') ' X',rMI(1)
  write(u6,'(19X,A,12X,3(E11.4))') ' Y',rMI(2),rMI(3)
  write(u6,'(19X,A,12X,3(E11.4))') ' Z',rMI(4),rMI(5),rMI(6)
  write(u6,*)
  call RecPrt('Pax',' ',Pax,3,3)
end if

! Diagonalize and find principle axis

call mma_Allocate(Hess,6)
Hess(:) = rMI(:)
call Jacob(Hess,Pax,3,3)
Prin(1) = Hess(1)
Prin(2) = Hess(3)
Prin(3) = Hess(6)

! Sort the prinipal axis such that z' is the one with the lowest eigenvalue.

do i=1,2
  do j=i+1,3
    if (Prin(i) < Prin(j)) then
      keep = Prin(i)
      Prin(i) = Prin(j)
      Prin(j) = keep
      call DSwap_(3,PAx(:,i),1,PAx(:,j),1)
    end if
  end do
end do
call mma_deallocate(Hess)
if (RR_Show) then
  write(u6,'(19X,A)') ' The Principal Axes and Moments of Inertia (au)'
  write(u6,'(19X,A,3(E11.4))') ' Eigenvalues :',(Prin(i),i=1,3)
  write(u6,'(19X,14X,3A)') '    X''    ','     Y''   ','    Z''    '
  write(u6,'(19X,A)') ' Eigenvectors:'
  write(u6,'(19X,A,3(E11.4))') ' X            ',Pax(1,:)
  write(u6,'(19X,A,3(E11.4))') ' Y            ',Pax(2,:)
  write(u6,'(19X,A,3(E11.4))') ' Z            ',Pax(3,:)
  write(u6,*)
  !call Put_dArray('PAX',Pax,9)
  write(u6,'(19X,A)') ' The Rotational Constants'
  write(u6,'(19X,A)') '         (cm-1)            (GHz)'
  if (Prin(1) >= 1.0e-3_wp) write(u6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(1),1.0e-9_wp*auToHz*Half/Prin(1)
  if (Prin(2) >= 1.0e-3_wp) write(u6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(2),1.0e-9_wp*auToHz*Half/Prin(2)
  if (Prin(3) >= 1.0e-3_wp) write(u6,'(19X,F16.3,1X,F16.3)') auTocm*Half/Prin(3),1.0e-9_wp*auToHz*Half/Prin(3)
end if
if ((Prin(1) == Zero) .and. (Prin(2) == Zero) .and. (Prin(3) == Zero)) then
  call FinishUp()
  return
end if
if (RR_Show) then
  write(u6,*)
  write(u6,*)
  write(u6,'(19X,A)') ' *******************************************'
  write(u6,'(19X,A)') ' *                                         *'
  write(u6,'(19X,A)') ' * R I G I D - R O T O R   A N A L Y S I S *'
  write(u6,'(19X,A)') ' *                                         *'
  write(u6,'(19X,A)') ' *******************************************'
  write(u6,*)
  write(u6,'(19X,A,I3)') ' j(Max):',S%jMax
  write(u6,*)
end if

! Order the three principal moments of inertia, Ia<=Ib<=Ic

XI(:) = Prin
call Order_Axis(XI,3)

! Asymmetry parameter, Tinkham formula 7-96

B = Half/XI(2)
C = Half/XI(3)
if (XI(1) <= 1.0e-3_wp) then
  Linear = .true.
  rKappa = -One
else
  A = Half/XI(1)
  if (abs(A-C) <= 1.0e-10_wp) then
    rKappa = Zero
  else
    rKappa = (Two*B-A-C)/(A-C)
  end if
end if

! Change order if molecule is oblate, Ia=Ib<=Ic

if ((abs(XI(1)-XI(2)) < 1.0e-6_wp) .and. (abs(XI(1)-XI(3)) > 1.0e-6_wp)) then
  keep = XI(1)
  XI(1) = XI(3)
  XI(3) = keep
end if

! Construct Hessian Matrix
! Iz=XI(1), Ix=XI(2), and Iy=XI(3)

! Set up constants, see Tinkham fomula 7-89a and 7-89b.

if (abs(XI(1)) > 1.0e-3_wp) then
  A = Half/XI(1)
else
  A = Zero
end if
B = (One/XI(2)+One/XI(3))/Four
C = (One/XI(2)-One/XI(3))/Eight
do i=1,80
  Label(i:i) = ' '
end do
! Iz=0 linear rotor
if (abs(A) < 1.0e-3_wp) then
  Label(1:13) = ' Linear Rotor'
! Ix=/=Iy asymmetric top
else if (C > 1.0e-3_wp) then
  Label(1:15) = ' Asymmetric Top'
! Iz=Ix=Iy
else if (abs(A-B) < 1.0e-10_wp) then
  Label(1:14) = ' Spherical Top'
else if (A > B) then
  Label(1:25) = ' Symmetric Top, Prolate'
else
  Label(1:24) = ' Symmetric Top, Oblate'
end if
if (RR_Show) then
  write(u6,'(19X,A,A)') ' Rotor Type:',trim(Label)
  write(u6,'(19X,A,F7.3)') ' Asymmetry parameter:',rKappa
  write(u6,'(19X,A)') ' Prolate = -1'
  write(u6,'(19X,A)') ' Oblate  =  1'
  write(u6,*)
end if

nEn = (S%jMax+1)*(S%jMax+2)*(S%jMax+3)/6
call mma_Allocate(En,nEn)
iEn = 1
nHess = (2*S%jMax+1)*(2*S%jMax+2)/2
call mma_allocate(Hess,nHess)
call mma_Allocate(Vec,(2*S%jMax+1)**2)
do j=0,S%jMax
  mDim = 2*j+1
  nTri = mDim*(mDim+1)/2
  Hess(1:nTri) = Zero
  call unitmat(Vec,mDim)
  if (iPrint >= 99) call RecPrt(' Vec',' ',Vec,mDim,mDim)
  k1 = 1
  do k=-j,j
    kk = k1*(k1+1)/2

    ! Formula 7-89a, (j,k|H|j,k), diagonal term

    Hess(kk) = B*real(j*(j+1),kind=wp)+(A-B)*real(k**2,kind=wp)
    if (Linear) Hess(kk) = B*real(j*(j+1),kind=wp)
    if ((k+2 <= j) .and. (.not. Linear)) then
      k2 = k1+2
      kk2 = k2*(k2-1)/2+k1

      ! Formula 7-89b, (j,k+2|H|j,k), the only off diagonal term.

      Hess(kk2) = C*sqrt(real((j*(j+1)-k*(k+1))*(j*(j+1)-(k+1)*(k+2)),kind=wp))
    end if
    k1 = k1+1
  end do
  if (iPrint >= 99) call TriPrt(' Hessian',' ',Hess,mDim)
  call Jacob(Hess,Vec,mDim,mDim)
  if (iPrint >= 99) call TriPrt(' Hessian',' ',Hess,mDim)
  do i=1,mDim
    En(iEn+I-1) = Hess(i*(i+1)/2)*auTocm
  end do
  call Order_Axis(En(iEn),mDim)
  iEn = iEn+mDim
end do
call mma_deallocate(Vec)
call mma_deallocate(Hess)

! Output

if (RR_Show) then
  write(u6,*)
  write(u6,'(19X,A)') ' Rotational energies / cm-1'
  iEn = 1
  do j=0,S%jMax
    write(u6,*)
    if (Linear) then
      write(u6,'(19X,A,I2,A,F8.3)') ' E(J=',J,') = ',En(iEn)
      iEn = iEn+(2*j+1)
    else
      do kappa=-j,j
        write(u6,'(19X,A,I2,A,I2,A,F12.3)') ' E(J=',J,',kappa=',kappa,') = ',En(iEn)
        iEn = iEn+1
      end do
    end if
  end do
end if
call mma_deallocate(En)

call FinishUp()

return

contains

subroutine FinishUp()

  if (RR_Show) then
    call CollapseOutput(0,'   Rigid rotor info:')
    write(u6,*)
  end if

end subroutine FinishUp

end subroutine RigRot

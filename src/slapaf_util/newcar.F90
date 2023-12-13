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

!#define _DEBUGPRINT_
subroutine NewCar(Iter,nAtom,Coor,mTtAtm,Error)
!***********************************************************************
!                                                                      *
!     Object: To compute the new symm. distinct Cartesian coordinates  *
!             from the suggested shift of the internal coordinates.    *
!                                                                      *
!***********************************************************************

use Symmetry_Info, only: VarR, VarT
use Slapaf_Info, only: AtomLbl, BMx, BSet, Curvilinear, Cx, Degen, HSet, Lbl, lOld, qInt, RefGeo, Shift, User_Def, &
                       WeightedConstraints
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: Iter, nAtom, mTtAtm
real(kind=wp), intent(inout) :: Coor(3,nAtom)
logical(kind=iwp), intent(inout) :: Error
#include "print.fh"
#include "warnings.h"
logical(kind=iwp) :: BSet_Save, Converged, HSet_Save, lOld_Save
integer(kind=iwp) :: i, iAtom, iInter, iMax, iPrint, iRout, iterMx, jter, Lu, M, N, nQQ, nWndw
real(kind=wp) :: denom, dx2, dx_RMS, rMax
real(kind=wp), allocatable :: DFC(:,:), dss(:), rInt(:)
integer(kind=iwp), parameter :: NRHS = 1

!                                                                      *
!***********************************************************************
!                                                                      *
nQQ = size(qInt,1)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('NewCar: q',' ',qInt,nQQ,Iter+1)
!call RecPrt('NewCar: Shift',' ',Shift,nQQ,iter)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(DFC,3,nAtom,Label='DFC')
call mma_allocate(dss,nQQ,Label='dss')
call mma_allocate(rInt,nQQ,Label='rInt')
!                                                                      *
!***********************************************************************
!                                                                      *
rInt(:) = qInt(:,iter)
dss(:) = Shift(:,iter)
rInt(:) = rInt(:)+dss(:)
!                                                                      *
!***********************************************************************
!                                                                      *
Lu = u6
iRout = 33
iPrint = nPrint(iRout)
#ifdef _DEBUGPRINT_
iPrint = 99
#endif
if (iPrint >= 11) then
  write(Lu,*)
  write(Lu,*) ' *** Transforming internal coordinates to Cartesian ***'
  write(Lu,*)
  write(Lu,*) ' Iter  Internal  Error'
end if

if (iPrint >= 99) then
  write(Lu,*)
  write(Lu,*) ' In NewCar: Shifts'
  write(Lu,*)
  write(Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter),iInter=1,nQQ)
  call RecPrt(' In NewCar: qInt',' ',qInt,nQQ,Iter+1)
end if

! Compute the final internal coordinates, plus sign due to the use
! of forces and not gradients.

rMax = Zero
iMax = 0
jter = 0
do i=1,nQQ
  if (abs(dss(i)) > abs(rMax)) then
    rMax = dss(i)
    iMax = i
  end if
end do
if (iPrint >= 11) then
  if (iMax /= 0) then
    write(Lu,300) jter,Lbl(iMax),rMax
  else
    write(Lu,300) jter,'N/A     ',rMax
  end if
end if

if (iPrint >= 19) then
  write(Lu,*)
  write(Lu,*) ' Internal coordinates of the next macro iteration'
  write(Lu,*)
  write(Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),rInt(iInter),iInter=1,nQQ)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the new Cartesian coordinates.

iterMx = 50
Converged = .false.
do jter=1,iterMx

  ! Compute the Cartesian shift, solve dq = B^T dx!

  M = 3*nAtom
  N = nQQ
  call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,dSS,DFC)
  call mma_deallocate(BMx)

  if (iPrint >= 99) call PrList('Symmetry Distinct Nuclear Displacements',AtomLbl,nAtom,DFC,3,nAtom)

  ! Compute the RMS in Cartesian Coordinates.

  dx2 = Zero
  denom = Zero
  do iAtom=1,nAtom
    do i=1,3
      dx2 = dx2+Degen(i,iAtom)*DFC(i,iAtom)**2
      denom = denom+Degen(i,iAtom)
    end do
  end do
  dx_RMS = sqrt(dx2/denom)

  ! Update the symmetry distinct Cartesian coordinates.

  Coor(:,:) = Coor(:,:)+DFC(:,:)

  ! Dirty fix of zeros

  do iAtom=1,nAtom
    if ((Cx(1,iAtom,Iter) == Zero) .and. (abs(Coor(1,iAtom)) < 1.0e-13_wp)) Coor(1,iAtom) = Zero
    if ((Cx(2,iAtom,Iter) == Zero) .and. (abs(Coor(2,iAtom)) < 1.0e-13_wp)) Coor(2,iAtom) = Zero
    if ((Cx(3,iAtom,Iter) == Zero) .and. (abs(Coor(3,iAtom)) < 1.0e-13_wp)) Coor(3,iAtom) = Zero
  end do

  Cx(:,:,Iter+1) = Coor(:,:)
  if (iPrint >= 99) call PrList('Symmetry Distinct Nuclear Coordinates / Bohr',AtomLbl,nAtom,Coor,3,nAtom)

  ! Compute new values q and the Wilson B-matrix for the new
  ! geometry with the current new set of Cartesian coordinates.

  nWndw = 1
  BSet_Save = BSet
  HSet_Save = HSet
  lOld_Save = lOld
  BSet = .false.
  HSet = .false.
  lOld = .false.
  call BMtrx(nAtom,Coor,iter+1,mTtAtm,nWndw)
  BSet = BSet_Save
  HSet = HSet_Save
  lOld = lOld_Save

  ! Check if the final structure is reached and get the
  ! difference between the present structure and the final.

  iMax = 1
  rMax = Zero
  do i=1,nQQ
    dSS(i) = rInt(i)-qInt(i,Iter+1)
    if (abs(dSS(i)) > abs(rMax)) then
      rMax = dSS(i)
      iMax = i
    end if
  end do

  ! Convergence based on the RMS of the Cartesian displacements.

  if (dx_RMS < 1.0e-6_wp) then
    Converged = .true.
    exit
  end if

  if (iPrint >= 99) then
    write(Lu,*)
    write(Lu,*) ' Displacement of internal coordinates'
    write(Lu,*)
    write(Lu,'(1X,A,2X,F10.4)') (Lbl(iInter),dss(iInter),iInter=1,nQQ)
  end if

end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. Converged) then
  ! On input, Error specifies whether an error should be signalled
  ! (.True.) or the calculation aborted (.False.)
  ! In the former case, the output value indicates if an error has occurred

  if (.not. User_Def) call RecPrt('NewCar: rInt  ','(10F15.10)',rInt,nQQ,1)
  if (.not. User_Def) call RecPrt('NewCar: qInt','(10F15.10)',qInt(:,Iter+1),nQQ,1)
  if (.not. Error) then
    call WarningMessage(2,'Error in NewCar')
    write(Lu,*)
    write(Lu,*) '***********************************************'
    write(Lu,*) ' ERROR: No convergence in NewCar !             '
    write(Lu,*) ' Strong linear dependency among Coordinates.   '
    write(Lu,*) ' Hint: Try to change the Internal Coordinates. '
    write(Lu,*) '***********************************************'
    if (.not. User_Def) call RecPrt('NewCar: rInt  ','(10F15.10)',rInt,nQQ,1)
    if (.not. User_Def) call RecPrt('NewCar: qInt','(10F15.10)',qInt(:,Iter+1),nQQ,1)
    write(Lu,*)
    call Quit(_RC_NOT_CONVERGED_)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *

if (iPrint >= 6) then
  write(Lu,*)
  write(Lu,'(A,i2,A)') ' New Cartesian coordinates were found in',jter,' Newton-Raphson iterations.'
  write(Lu,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Finally, just to be safe align the new Cartesian structure with
! the reference structure (see init2.f)

if (WeightedConstraints .and. (.not. (VarR .or. VarT))) call Align(Cx(:,:,iter+1),RefGeo,nAtom)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(rInt)
call mma_deallocate(dss)
call mma_deallocate(DFC)
!                                                                      *
!***********************************************************************
!                                                                      *
Error = .false.

return

300 format(1X,'Iter:',I5,2X,A,1X,ES11.4)

end subroutine NewCar

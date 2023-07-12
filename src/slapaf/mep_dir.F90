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
! Copyright (C) 2013-2015, Ignacio Fdez. Galvan                        *
!***********************************************************************
!  MEP_Dir
!
!> @brief
!>   Compute the new reference structure and initial coordinates for the next MEP point
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Compute the new reference structure (and direction) and initial geometry
!> for the next MEP point optimization (using the Gonz&aacute;lez--Schlegel or
!> M&uacute;ller--Brown method).
!> Calculate some properties of the path (length and curvature) between the previous
!> and current MEP points.
!> All calculations are done in weighted coordinates (mass-weighted by default).
!>
!> @param[in,out] Cx            Cartesian coordinates in all iterations
!> @param[in]     Gx            Cartesian gradient in all iterations
!> @param[in]     nAtom         Number of symmetry-unique atoms
!> @param[in]     iMEP          Number of this MEP point
!> @param[in]     iOff_iter     Iteration of the previous MEP point
!> @param[in]     iPrint        Print level
!> @param[in]     IRCRestart    Flag to mark the start of a backward IRC search
!> @param[out]    ResGrad       Residual gradient
!> @param[out]    BadConstraint Flag to signal constraint problems
!***********************************************************************

subroutine MEP_Dir(Cx,Gx,nAtom,iMEP,iOff_iter,iPrint,IRCRestart,ResGrad,BadConstraint)

use Symmetry_Info, only: nIrrep
use Slapaf_Info, only: dMEPStep, IRC, iter, MEP, MEP_Algo, MEP_Type, MEPNum, MF, nLambda, nMEP, RefGeo, rMEP, Weights
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, iMEP, iOff_iter, iPrint
real(kind=wp), intent(inout) :: Cx(3*nAtom,iter+1)
real(kind=wp), intent(in) :: Gx(3*nAtom,iter+1)
logical(kind=iwp), intent(in) :: IRCRestart
real(kind=wp), intent(out) :: ResGrad
logical(kind=iwp), intent(out) :: BadConstraint
integer(kind=iwp) :: iAd, iAtom, iDum(1), iLambda, iOff, iPrev_iter, ixyz, LudRdX, nCoor, nCoor_, nLambda_
real(kind=wp) :: ConstraintAngle, Curvature, dd, dDir, dDisp, dGrad, dPostDir, dPostDirGrad, dPrevDir, dPrevDirDisp, dPrevDirGrad, &
                 dPrevDirPostDir, drd, Fact, PathAngle, PathLength, TWeight, xWeight
real(kind=wp), allocatable :: Cen(:,:), Cur(:), Dir(:,:), Disp(:,:), drdx(:,:), Grad(:,:), PLn(:), PostDir(:,:), PrevDir(:,:)
integer(kind=iwp), external :: iDeg
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
nCoor = 3*nAtom
iPrev_iter = max(iOff_iter,1)
call mma_allocate(PrevDir,3,nAtom,Label='PrevDir')
call mma_allocate(PostDir,3,nAtom,Label='PostDir')
call mma_allocate(Disp,3,nAtom,Label='Disp')
call mma_allocate(Grad,3,nAtom,Label='Grad')
call mma_allocate(Dir,3,nAtom,Label='Dir')
call mma_allocate(Cen,3,nAtom,Label='Cen')

! Obtain some useful vectors:
! PrevDir: difference between ref. structure and previous MEP point
! PostDir: difference between current MEP point and ref. structure
! Disp:    difference between current and previous MEP points
! Grad:    gradient at current MEP point

if (iter > 1) then
  PrevDir(:,:) = RefGeo(:,:)-reshape(Cx(:,iPrev_iter),[3,nAtom])
  PostDir(:,:) = reshape(Cx(:,iter),[3,nAtom])-RefGeo(:,:)
  Disp(:,:) = PrevDir(:,:)+PostDir(:,:)
else
  PrevDir(:,:) = Zero
  PostDir(:,:) = Zero
  Disp(:,:) = Zero
end if
Grad(:,:) = reshape(Gx(:,iter),[3,nAtom])

! Normalize the vectors in weighted coordinates
! and compute some angles that provide information on the path
! shape and quality.
! Note that gradient and coordinates transform differently

TWeight = Zero
dPrevDir = Zero
dPostDir = Zero
dDisp = Zero
dGrad = Zero
dPostDirGrad = Zero
dPrevDirGrad = Zero
dPrevDirDisp = Zero
dPrevDirPostDir = Zero
iOff = 0
do iAtom=1,nAtom
  Fact = real(iDeg(Cx(1+iOff,iter)),kind=wp)
  xWeight = Weights(iAtom)
  TWeight = TWeight+Fact*xWeight
  do ixyz=1,3
    dPrevDir = dPrevDir+Fact*xWeight*PrevDir(ixyz,iAtom)**2
    dPostDir = dPostDir+Fact*xWeight*PostDir(ixyz,iAtom)**2
    dDisp = dDisp+Fact*xWeight*Disp(ixyz,iAtom)**2
    dGrad = dGrad+Fact*Grad(ixyz,iAtom)**2/xWeight
    dPostDirGrad = dPostDirGrad+Fact*PostDir(ixyz,iAtom)*Grad(ixyz,iAtom)
    dPrevDirGrad = dPrevDirGrad+Fact*PrevDir(ixyz,iAtom)*Grad(ixyz,iAtom)
    dPrevDirDisp = dPrevDirDisp+Fact*xWeight*PrevDir(ixyz,iAtom)*Disp(ixyz,iAtom)
    dPrevDirPostDir = dPrevDirPostDir+Fact*xWeight*PrevDir(ixyz,iAtom)*PostDir(ixyz,iAtom)
    iOff = iOff+1
  end do
end do
dPrevDir = sqrt(dPrevDir)
dPostDir = sqrt(dPostDir)
dDisp = sqrt(dDisp)
dGrad = sqrt(dGrad)
if (dPrevDir > Zero) PrevDir(:,:) = PrevDir(:,:)/dPrevDir
if (dPostDir > Zero) PostDir(:,:) = PostDir(:,:)/dPostDir
if (dDisp > Zero) Disp(:,:) = Disp(:,:)/dDisp
if (dGrad > Zero) Grad(:,:) = Grad(:,:)/dGrad
! Any zero vector is assumed to be parallel to any other
if (dPostDir*dGrad > Zero) then
  dPostDirGrad = dPostDirGrad/(dPostDir*dGrad)
else
  dPostDirGrad = One
end if
if (dPrevDir*dGrad > Zero) then
  dPrevDirGrad = dPrevDirGrad/(dPrevDir*dGrad)
else
  dPrevDirGrad = One
end if
if (dPrevDir*dDisp > Zero) then
  dPrevDirDisp = dPrevDirDisp/(dPrevDir*dDisp)
else
  dPrevDirDisp = One
end if
if (dPrevDir*dPostDir > Zero) then
  dPrevDirPostDir = dPrevDirPostDir/(dPrevDir*dPostDir)
else
  dPrevDirPostDir = One
end if
! A negative curvature means there is no appropriate value
Curvature = -1.0e-12_wp
PathLength = dDisp/sqrt(TWeight)
if (MEP .and. (MEP_Type == 'SPHERE    ')) then
  ! The curvature is the inverse radius of the circle tangent to
  ! both the current and previous MEP points
  ! The path length is the arc length between these points
  if (One-dPrevDirPostDir > Zero) Curvature = (One-dPrevDirPostDir)/sqrt(One-dPrevDirPostDir**2)/sqrt(dPrevDir*dPostDir/TWeight)
  if (Curvature > Zero) PathLength = acos(dPrevDirPostDir)/Curvature
end if

! Store the length and curvature values, and print results

call mma_allocate(PLn,nMEP+1,Label='PLn')
call mma_allocate(Cur,nMEP+1,Label='Cur')
if (iMEP >= 1) then
  call Get_dArray('MEP-Lengths   ',PLn,nMEP+1)
  call Get_dArray('MEP-Curvatures',Cur,nMEP+1)
  if (IRC == -1) then
    PLn(1+iMEP) = -PathLength
  else
    PLn(1+iMEP) = PathLength
  end if
  Cur(1+iMEP) = Curvature
  call Put_dArray('MEP-Lengths   ',PLn,nMEP+1)
  call Put_dArray('MEP-Curvatures',Cur,nMEP+1)
  if ((iMEP >= 1) .and. (iPrint >= 5)) then
    if (MEP_Type == 'TRANSVERSE') then
      ConstraintAngle = acos(dPrevDirGrad)/deg2rad
    else
      ConstraintAngle = acos(dPostDirGrad)/deg2rad
    end if
    if (MEP) then
      PathAngle = acos(dPrevDirPostDir)/deg2rad
    else
      PathAngle = acos(-dPrevDirDisp)/deg2rad
    end if
    ResGrad = dGrad
    write(u6,*)
    write(u6,'(a)') ' Last IRC/MEP step (in weighted coordinates / sqrt(total weight))'
    write(u6,'(a)') ' --------------------------------------------------------'
    write(u6,100) 'Residual gradient size:',ResGrad,'hartree/bohr'
    write(u6,100) 'Angle with constraint surface:',ConstraintAngle,'degrees'
    write(u6,100) 'Path angle:',PathAngle,'degrees'
    if (Curvature >= Zero) write(u6,100) 'Path curvature:',Curvature,'bohr^(-1)'
    write(u6,100) 'Path length:',PathLength,'bohr'
  end if
else
  PLn(:) = Zero
  Cur(:) = Zero
  call Put_dArray('MEP-Lengths   ',PLn,nMEP+1)
  call Put_dArray('MEP-Curvatures',Cur,nMEP+1)
end if
call mma_deallocate(PLn)
call mma_deallocate(Cur)

! Do not mess with the geometry or reference if the next iteration
! will be the start of a reverse IRC search.

if (.not. IRCRestart) then

  ! The new direction for the MEP should be the gradient in weighted
  ! coordinates, but this may break additional constraints if they
  ! exist, and the gradient may be close to zero.
  ! Instead, we will use the direction that, on convergence, should
  ! be parallel to the gradient when there are no other constraints.
  ! Note that we could be following the gradient uphill

  if (MEP_Type == 'TRANSVERSE') then
    ! In the TRANSVERSE case, PrevDir is the vector parallel to the
    ! gradient, but following this will never change the hyperplane
    ! orientation.
    ! We will try using a linear combination of PrevDir and Disp
    !   dp  = Disp.PrevDir
    !   Dir = dp*Disp + a*(PrevDir-dp*Disp)
    !       = a*PrevDir + (1-a)*dp*Disp
    ! Try different values for Fact (the "a" above).
    ! I believe the true MEP direction should be more slanted from
    ! the plane normal (PrevDir) than the Disp vector, therefore
    ! a negative value is probably better.
    Fact = -Half
    Dir(:,:) = Fact*PrevDir(:,:)+(One-Fact)*dPrevDirDisp*Disp(:,:)
  else
    ! In the SPHERE case, PostDir is the vector to use
    Dir(:,:) = PostDir(:,:)
  end if

  ! Special cases

  if (iMEP == 0) then
    if (IRC == 0) then
      ! In the initial iteration of a MEP, use the initial direction
      call Get_dArray('Transverse',Dir(:,:),nCoor)
    else
      ! In the initial iteration of an IRC branch, use the reaction vector
      Dir(:,:) = MF(:,:)
    end if
  end if

  ! Project any additional constraints out of the direction vector
  ! The constraint vectors are read from the dRdX file

  LudRdX = 30
  call DaName(LudRdX,'dRdX')
  iAd = 0
  call iDaFile(LudRdX,2,iDum,1,iAd)
  nLambda_ = iDum(1)
  call iDaFile(LudRdX,2,iDum,1,iAd)
  nCoor_ = iDum(1)
  call mma_allocate(drdx,nCoor_,nLambda_,Label='drdx')
  call dDaFile(LudRdX,2,drdx,nLambda_*nCoor_,iAd)
  call DaClos(LudRdX)
  do iLambda=1,nLambda
    if (iLambda /= MEPnum) then
      dd = dDot_(nCoor,drdx(:,iLambda),1,drdx(:,iLambda),1)
      drd = dDot_(nCoor,drdx(:,iLambda),1,Dir(:,:),1)
      Dir(:,:) = Dir(:,:)-drd/dd*reshape(drdx(:,iLambda),[3,nAtom])
    end if
  end do
  call mma_deallocate(drdx)

  ! Compute the length of the direction vector in weighted coordinates

  dDir = Zero
  iOff = 0
  do iAtom=1,nAtom
    Fact = real(iDeg(Cx(1+iOff,iter)),kind=wp)
    xWeight = Weights(iAtom)
    do ixyz=1,3
      dDir = dDir+Fact*xWeight*Dir(ixyz,iAtom)**2
      iOff = iOff+1
    end do
  end do
  dDir = sqrt(dDir)

  ! According to the Gonzalez-Schlegel method, the reference point
  ! is half-step away in the search direction.
  ! First set the new reference point (except for rMEP) and then
  ! compute the new starting structure at the full step distance
  ! For an IRC first step, keep the initial structure as reference

  Fact = dMEPStep*sqrt(TWeight)/dDir
  Cen(:,:) = reshape(Cx(:,iter),[3,nAtom])
  if (MEP_Algo == 'GS') then
    if ((IRC == 0) .or. (iMEP /= 0)) call Find_Distance(Cx(:,iter),Cen(:,:),Dir(:,:),Half*Fact,Half*dMEPStep,nAtom,BadConstraint)
    if (.not. rMEP) call Put_dArray('Ref_Geom',Cen(:,:),nCoor)
    call Find_Distance(Cen(:,:),Cx(:,iter+1),Dir(:,:),Half*Fact,Half*dMEPStep,nAtom,BadConstraint)
  else if (MEP_Algo == 'MB') then
    if (.not. rMEP) call Put_dArray('Ref_Geom',Cx(1,iter),nCoor)
    call Find_Distance(Cx(:,iter),Cx(:,iter+1),Dir(:,:),Fact,dMEPStep,nAtom,BadConstraint)
  end if

  ! Randomly displace the new geometry 1/20th of the MEP distance
  ! to try to break symmetry

  if (nIrrep == 1) then
    call Random_Vector(nCoor,Disp(:,:),.true.)
    dDir = Zero
    iOff = 0
    do iAtom=1,nAtom
      xWeight = Weights(iAtom)
      do ixyz=1,3
        dDir = dDir+xWeight*Disp(ixyz,iAtom)**2
        iOff = iOff+1
      end do
    end do
    Fact = dMEPStep*sqrt(TWeight/dDir)
    Cx(:,iter+1) = Cx(:,iter+1)+0.05_wp*Fact*reshape(Disp(:,:),[nCoor])
  end if
  call Put_dArray('Transverse',Dir(:,:),nCoor)
  if (iter == 1) BadConstraint = .false.
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(PrevDir)
call mma_deallocate(PostDir)
call mma_deallocate(Disp)
call mma_deallocate(Grad)
call mma_deallocate(Dir)
call mma_deallocate(Cen)

return

100 format(1X,A30,1X,F12.6,1X,A)

end subroutine MEP_Dir

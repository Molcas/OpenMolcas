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
! Copyright (C) 1999, Roland Lindh                                     *
!***********************************************************************

subroutine Subblock(iNQ,x_NQ,y_NQ,z_NQ,InBox,x_min_,x_max_,y_min_,y_max_,z_min_,z_max_,list_p,nlist_p,Grid,Weights,mGrid,Process, &
                    number_of_grid_points,R_box_min,R_box_max,ilist_p,xyz0,iAngular_Grid,nR_Eff)
!***********************************************************************
!                                                                      *
! Object:                                                              *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept of Chemical Physics,                                *
!             University of Lund, Sweden                               *
!             August 1999                                              *
!***********************************************************************

use NQ_structure, only: Info_Ang, NQ_Data
use nq_Info, only: Angular_Pruning, iOpt_Angular, ntotgp, On
use Grid_On_Disk, only: ExpandBatchInfo, iBatchInfo, iDisk_Grid, Lu_Grid, nBatch
use Constants, only: One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iNQ, nlist_p, list_p(nlist_p), mGrid, ilist_p, nR_Eff, iAngular_Grid(nR_Eff)
real(kind=wp), intent(in) :: x_NQ, y_NQ, z_NQ, x_min_, x_max_, y_min_, y_max_, z_min_, z_max_, R_box_min, R_box_max, xyz0(3,2)
logical(kind=iwp), intent(in) :: InBox, Process
real(kind=wp), intent(inout) :: Grid(3,mGrid), Weights(mGrid)
integer(kind=iwp), intent(inout) :: number_of_grid_points
integer(kind=iwp) :: iEnd_R, iPoint, iR, iR_End, iR_Start, iSet, iStart_R, iStrt, kSet, mGrid_, nGrid, nRemoved
real(kind=wp) :: Fact, R_Value, Radius, w_g, weight, x, xpt, y, ypt, z, zpt

!                                                                      *
!***********************************************************************
!                                                                      *
nGrid = (9*mGrid)/10
iStrt = number_of_grid_points+1
!                                                                      *
!***********************************************************************
!                                                                      *
! Start loop over the atomic grid

#ifdef _DEBUGPRINT_
write(u6,*) ' x_NQ=',x_NQ
write(u6,*) ' y_NQ=',y_NQ
write(u6,*) ' z_NQ=',z_NQ
write(u6,*) ' Process=',Process
write(u6,*) ' number_of_grid_points=',number_of_grid_points
write(u6,*) 'x:',x_min_,x_max_
write(u6,*) 'y:',y_min_,y_max_
write(u6,*) 'z:',z_min_,z_max_
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
iStart_R = nR_Eff
iEnd_R = 1
!write(u6,*)
!write(u6,*) 'Start range:',iEnd_R,iStart_R
if (.not. btest(iOpt_Angular,1)) then
  !write(u6,*) 'Find R subrange!'

  ! Compute valid subrange for R

  iR_End = iEnd_R
  do iR=iEnd_R,iStart_R
    R_Value = NQ_Data(iNQ)%R_Quad(1,iR)
    if (R_Value > R_box_Min) exit
    iR_End = iR
  end do

  iR_Start = iStart_R
  do iR=iStart_R,iR_End,-1
    R_Value = NQ_Data(iNQ)%R_Quad(1,iR)
    if (R_Value < R_Box_Max) exit
    iR_Start = iR
  end do

else
  !write(u6,*) 'Do whole R range!'

  ! Scan the whole range

  iR_Start = iStart_R
  iR_End = iEnd_R

end if

! Reset iStart_R and iEnd_R, these are not to be modified again!

iStart_R = iR_Start
iEnd_R = iR_End
iSet = -1
!write(u6,*) 'Actual range:',iEnd_R,iStart_R
!                                                                      *
!***********************************************************************
!                                                                      *
! Outer loop over angular grids

do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Determine a range (iStart_R,iEnd_R) where we will use a
  ! specific angular grid. As we get closer to the nuclei we will
  ! reduce the order of the angular grid.

  ! Start loop at the outermost point and iterate towards the
  ! nuclei.

  iR_End = iR_Start+1
  do iR=iR_Start,iEnd_R,-1
    kSet = iAngular_Grid(iR)

    ! Save new iSet for the first point of this subrange.

    if (iR == iR_Start) iSet = kSet

    ! Branch out if we hit on a range where we can reduce the
    ! angular grid further.

    if (kSet /= iSet) exit

    ! Update inner index

    iR_End = iR
  end do

  ! Branch out if subrange is outside the box

  if ((NQ_Data(iNQ)%R_Quad(1,iR_End) <= R_box_Max) .and. (NQ_Data(iNQ)%R_Quad(1,iR_Start) >= R_box_Min)) then

    !write(u6,*) 'Selected range:',iR_End,iR_Start
    !write(u6,*) 'l_max=',Info_Ang(iSet)%L_Eff
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Angular loop

    do iPoint=1,Info_Ang(iSet)%nPoints
#     ifdef _DEBUGPRINT_
      write(u6,*) 'X,Y,Z*',Info_Ang(iSet)%R(1,iPoint),Info_Ang(iSet)%R(2,iPoint),Info_Ang(iSet)%R(3,iPoint)
#     endif

      if ((.not. InBox) .and. (.not. btest(iOpt_Angular,1))) then
        !write(u6,*) 'Select angular points!'
        if ((Info_Ang(iSet)%R(1,iPoint) < xyz0(1,1)) .or. (Info_Ang(iSet)%R(1,iPoint) > xyz0(1,2)) .or. &
            (Info_Ang(iSet)%R(2,iPoint) < xyz0(2,1)) .or. (Info_Ang(iSet)%R(2,iPoint) > xyz0(2,2)) .or. &
            (Info_Ang(iSet)%R(3,iPoint) < xyz0(3,1)) .or. (Info_Ang(iSet)%R(3,iPoint) > xyz0(3,2))) cycle
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Radial loop over the reduced range

      do iR=iR_End,iR_Start
        Radius = NQ_Data(iNQ)%R_Quad(1,iR)
        ! In the atomic referential
        xpt = Radius*Info_Ang(iSet)%R(1,iPoint)
        ypt = Radius*Info_Ang(iSet)%R(2,iPoint)
        zpt = Radius*Info_Ang(iSet)%R(3,iPoint)
        ! In the system referential
        x = xpt+x_NQ
        y = ypt+y_NQ
        z = zpt+z_NQ
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Radius=',Radius
        write(u6,*) ' x,y,z:',x,y,z
        write(u6,*) x_NQ,y_NQ,z_NQ
#       endif

        ! Check if the point is inside the box.
        ! Points on the border of boxes are shared.

        if ((x >= x_min_) .and. (x <= x_max_) .and. (y >= y_min_) .and. (y <= y_max_) .and. (z >= z_min_) .and. (z <= z_max_)) then

          ! For shared points modify the weight.

          Fact = One
          if (x == x_min_) Fact = Fact*Half
          if (y == y_min_) Fact = Fact*Half
          if (z == z_min_) Fact = Fact*Half
          if (x == x_max_) Fact = Fact*Half
          if (y == y_max_) Fact = Fact*Half
          if (z == z_max_) Fact = Fact*Half
#         ifdef _DEBUGPRINT_
          write(u6,*) Info_Ang(iSet)%R(1,iPoint),Info_Ang(iSet)%R(2,iPoint),Info_Ang(iSet)%R(3,iPoint)
          write(u6,*) 'x:',xyz0(1,1),xyz0(1,2)
          write(u6,*) 'y:',xyz0(2,1),xyz0(2,2)
          write(u6,*) 'z:',xyz0(3,1),xyz0(3,2)
          write(u6,*) ' Inside:',x,y,z
#         endif

          ! Radial weight
          weight = NQ_Data(iNQ)%R_Quad(2,iR)
          ! Combine the radial and angular weight
          w_g = weight*Info_Ang(iSet)%R(4,iPoint)
          if (w_g*Fact >= 1.0e-15_wp) then
            number_of_grid_points = number_of_grid_points+1
            Grid(1,number_of_grid_points) = x
            Grid(2,number_of_grid_points) = y
            Grid(3,number_of_grid_points) = z
            ! Compute the partitioning weight
            Weights(number_of_grid_points) = w_g*Fact
          end if
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        ! The call of Do_Batch is done if the buffer is full:
        ! (number_of_grid_points > nGrid)
        !                                                              *
        !***************************************************************
        !                                                              *
        if (number_of_grid_points > mGrid) then
          call WarningMessage(2,'Subblock: Buffer overflowed!;Try a larger buffer size!')
          call Abend()
        end if
        if (number_of_grid_points > nGrid) then

          ! Dump grid information to disk
          nBatch = nBatch+1
          if (nBatch > size(iBatchInfo,2)) call ExpandBatchInfo()

          ! Generate weights

          mGrid_ = number_of_grid_points-iStrt+1
          call W(Grid(1,iStrt),ilist_p,Weights(iStrt),list_p,nList_p,mGrid_,nRemoved)
          number_of_grid_points = number_of_grid_points-nRemoved

          iBatchInfo(1,nBatch) = iDisk_Grid
          iBatchInfo(3,nBatch) = iNQ
          iBatchInfo(2,nBatch) = number_of_grid_points

          call dDaFile(Lu_Grid,1,Grid,3*number_of_grid_points,iDisk_Grid)
          call dDaFile(Lu_Grid,1,Weights,number_of_grid_points,iDisk_Grid)

          ntotgp = ntotgp+number_of_grid_points
          !write(u6,*) 'ntotgp=',ntotgp
          number_of_grid_points = 0
          iStrt = number_of_grid_points+1
        end if

      end do ! iR, Radial loop
    end do   ! iPoint
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  iR_Start = iR_End-1
  if ((Angular_Pruning /= On) .or. (iR_End == iEnd_R)) exit
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate weights
!
!write(u6,*) 'number_of_grid_points=',number_of_grid_points
!write(u6,*) 'iStrt=',iStrt
if (number_of_grid_points-iStrt+1 > 0) then
  mGrid_ = number_of_grid_points-iStrt+1
  call W(Grid(1,iStrt),ilist_p,Weights(iStrt),list_p,nList_p,mGrid_,nRemoved)
  number_of_grid_points = number_of_grid_points-nRemoved
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Process batch if not processed yet.

if (Process .and. (number_of_grid_points > 0)) then

  ! Dump grid information to disk
  nBatch = nBatch+1
  if (nBatch > size(iBatchInfo,2)) call ExpandBatchInfo()
  iBatchInfo(1,nBatch) = iDisk_Grid
  iBatchInfo(2,nBatch) = number_of_grid_points
  iBatchInfo(3,nBatch) = iNQ
  call dDaFile(Lu_Grid,1,Grid,3*number_of_grid_points,iDisk_Grid)
  call dDaFile(Lu_Grid,1,Weights,number_of_grid_points,iDisk_Grid)

  ntotgp = ntotgp+number_of_grid_points
  !write(u6,*) 'ntotgp=',ntotgp
  number_of_grid_points = 0
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Subblock

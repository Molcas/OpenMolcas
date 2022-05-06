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

use NQ_structure, only: NQ_Data, Info_Ang
use Grid_On_Disk
use nq_Info

implicit real*8(A-H,O-Z)
#include "itmax.fh"
#include "real.fh"
#include "setup.fh"
#include "nsd.fh"
#include "debug.fh"
integer list_p(nlist_p)
real*8 Grid(3,mGrid), Weights(mGrid), xyz0(3,2)
logical Process, InBox, Check
integer iAngular_Grid(nR_Eff)
! Statement functions
Check(i,j) = iand(i,2**(j-1)) /= 0
x_a(i,iSet) = Info_Ang(iSet)%R(1,i)
y_a(i,iSet) = Info_Ang(iSet)%R(2,i)
z_a(i,iSet) = Info_Ang(iSet)%R(3,i)
w_a(i,iSet) = Info_Ang(iSet)%R(4,i)

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
if (Debug) then
  write(6,*) ' x_NQ=',x_NQ
  write(6,*) ' y_NQ=',y_NQ
  write(6,*) ' z_NQ=',z_NQ
  write(6,*) ' Process=',Process
  write(6,*) ' number_of_grid_points=',number_of_grid_points
  write(6,*) 'x:',x_min_,x_max_
  write(6,*) 'y:',y_min_,y_max_
  write(6,*) 'z:',z_min_,z_max_
  write(6,*)
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
iStart_R = nR_Eff
iEnd_R = 1
!write(6,*)
!write(6,*) 'Start range:',iEnd_R,iStart_R
if (.not. Check(iOpt_Angular,2)) then
  !write(6,*) 'Find R subrange!'

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
  !write(6,*) 'Do whole R range!'

  ! Scan the whole range

  iR_Start = iStart_R
  iR_End = iEnd_R

end if

! Reset iStart_R and iEnd_R, these are not to be modified again!

iStart_R = iR_Start
iEnd_R = iR_End
iSet = -1
!write(6,*) 'Actual range:',iEnd_R,iStart_R
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

    !write(6,*) 'Selected range:',iR_End,iR_Start
    !write(6,*) 'l_max=',Info_Ang(iSet)%L_Eff
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Angular loop

    do iPoint=1,Info_Ang(iSet)%nPoints
#     ifdef _DEBUGPRINT_
      if (Debug) then
        write(6,*) 'X,Y,Z*',x_a(iPoint,iSet),y_a(iPoint,iSet),z_a(iPoint,iSet)
      end if
#     endif

      if ((.not. InBox) .and. (.not. Check(iOpt_Angular,2))) then
        !write(6,*) 'Select angular points!'
        if ((x_a(iPoint,iset) < xyz0(1,1)) .or. (x_a(iPoint,iset) > xyz0(1,2)) .or. (y_a(iPoint,iset) < xyz0(2,1)) .or. &
            (y_a(iPoint,iset) > xyz0(2,2)) .or. (z_a(iPoint,iset) < xyz0(3,1)) .or. (z_a(iPoint,iset) > xyz0(3,2))) cycle
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Radial loop over the reduced range

      do iR=iR_End,iR_Start
        Radius = NQ_Data(iNQ)%R_Quad(1,iR)
        ! In the atomic referential
        xpt = Radius*x_a(iPoint,iSet)
        ypt = Radius*y_a(iPoint,iSet)
        zpt = Radius*z_a(iPoint,iSet)
        ! In the system referential
        x = xpt+x_NQ
        y = ypt+y_NQ
        z = zpt+z_NQ
#       ifdef _DEBUGPRINT_
        if (Debug) then
          write(6,*) 'Radius=',Radius
          write(6,*) ' x,y,z:',x,y,z
          write(6,*) x_NQ,y_NQ,z_NQ
        end if
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
          if (Debug) then
            write(6,*) x_a(iPoint,iSet),y_a(iPoint,iSet),z_a(iPoint,iSet)
            write(6,*) 'x:',xyz0(1,1),xyz0(1,2)
            write(6,*) 'y:',xyz0(2,1),xyz0(2,2)
            write(6,*) 'z:',xyz0(3,1),xyz0(3,2)
            write(6,*) ' Inside:',x,y,z
          end if
#         endif

          ! Radial weight
          weight = NQ_Data(iNQ)%R_Quad(2,iR)
          ! Combine the radial and angular weight
          w_g = weight*w_a(iPoint,iSet)
          if (w_g*Fact >= 1.0D-15) then
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
          if (nBatch > nBatch_Max) then
            call WarningMessage(2,'Subblock: nBatch > nBatch_Max')
            call Abend()
          end if

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
          !write(6,*) 'ntotgp=',ntotgp
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
  if ((Angular_Prunning /= On) .or. (iR_End == iEnd_R)) exit
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate weights
!
!write(6,*) 'number_of_grid_points=',number_of_grid_points
!write(6,*) 'iStrt=',iStrt
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
  if (nBatch > nBatch_Max) then
    call WarningMessage(2,'Subblock: nBatch > nBatch_Max')
    call Abend()
  end if
  iBatchInfo(1,nBatch) = iDisk_Grid
  iBatchInfo(2,nBatch) = number_of_grid_points
  iBatchInfo(3,nBatch) = iNQ
  call dDaFile(Lu_Grid,1,Grid,3*number_of_grid_points,iDisk_Grid)
  call dDaFile(Lu_Grid,1,Weights,number_of_grid_points,iDisk_Grid)

  ntotgp = ntotgp+number_of_grid_points
  !write(6,*) 'ntotgp=',ntotgp
  number_of_grid_points = 0
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Subblock

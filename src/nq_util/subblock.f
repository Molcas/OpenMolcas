************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      SubRoutine Subblock(iNQ,x_NQ,y_NQ,z_NQ,InBox,
     &                    x_min_,x_max_,
     &                    y_min_,y_max_,
     &                    z_min_,z_max_,
     &                    list_p,nlist_p,
     &                    Grid,Weights,mGrid,
     &                    Process,
     &                    number_of_grid_points,
     &                    R_box_min,R_box_max,
     &                    ilist_p,xyz0,iAngular_Grid,nR_Eff)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "setup.fh"
#include "nq_info.fh"
#include "grid_on_disk.fh"
#include "nsd.fh"
#include "debug.fh"
      Integer list_p(nlist_p)
      Real*8 Grid(3,mGrid), Weights(mGrid), xyz0(3,2)
      Logical Process,InBox, Check
      Integer iAngular_Grid(nR_Eff)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions                                              *
*                                                                      *
#include "nq_structure.fh"
      declare_ip_r_quad
      iROff(i,ir)=2*(ir-1)+i-1
      Check(i,j)=iAnd(i,2**(j-1)).ne.0
      x_a(i,iSet)=Work(Info_Ang(3,iSet)+(i-1)*4  )
      y_a(i,iSet)=Work(Info_Ang(3,iSet)+(i-1)*4+1)
      z_a(i,iSet)=Work(Info_Ang(3,iSet)+(i-1)*4+2)
      w_a(i,iSet)=Work(Info_Ang(3,iSet)+(i-1)*4+3)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call QEnter('Subblock')
#endif
      nGrid=(9*mGrid)/10
      iStrt=number_of_grid_points+1
*                                                                      *
************************************************************************
*                                                                      *
*---- Start loop over the atomic grid
*
      ip_iRx=ip_of_iWork(Work(ip_R_Quad(iNQ)))
      ip_Rx=iWork(ip_iRx)
#ifdef _DEBUG_
      If (Debug) Then
         Write (6,*) 'ip_Rx=',ip_Rx
         Write (6,*) ' x_NQ=',x_NQ
         Write (6,*) ' y_NQ=',y_NQ
         Write (6,*) ' z_NQ=',z_NQ
         Write (6,*) ' Process=',Process
         Write (6,*) ' number_of_grid_points=',number_of_grid_points
         Write (6,*) 'x:', x_min_,x_max_
         Write (6,*) 'y:', y_min_,y_max_
         Write (6,*) 'z:', z_min_,z_max_
         Write (6,*)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iStart_R= nR_Eff
      iEnd_R=1
c     Write (*,*)
c     Write (*,*) 'Start range:',iEnd_R, iStart_R
      If (.Not.Check(iOpt_Angular,2)) Then
c        Write (*,*) 'Find R subrange!'
*
*------- Compute valid subrange for R
*
         iR_End = iEnd_R
         Do iR = iEnd_R, iStart_R
            R_Value = Work(ip_Rx+iROff(1,iR))
            If (R_Value.le.R_box_Min) Then
               iR_End = iR
            Else
               Go To 8888
            End If
         End Do
 8888    Continue
*
         iR_Start = iStart_R
         Do iR = iStart_R, iR_End, -1
            R_Value = Work(ip_Rx+iROff(1,iR))
            If (R_Value.ge.R_Box_Max) Then
               iR_Start = iR
            Else
               Go To 8889
            End If
         End Do
 8889    Continue
*
      Else
c        Write (*,*) 'Do whole R range!'
*
*------- Scan the whole range
*
         iR_Start=iStart_R
         iR_End  =iEnd_R
*
      End If
*
*---- Reset iStart_R and iEnd_R, these are not to be modified again!
*
      iStart_R=iR_Start
      iEnd_R  =iR_End
      iSet=-1
c     Write (*,*) 'Actual range:',iEnd_R, iStart_R
*                                                                      *
************************************************************************
*                                                                      *
*     Outer loop over angular grids
*
 999  Continue
*                                                                      *
************************************************************************
*                                                                      *
*        Determine a range (iStart_R,iEnd_R) where we will use a
*        specific angular grid. As we get closer to the nuclei we will
*        reduce the order of the angular grid.
*
*
*------- Start loop at the outermost point and iterate towards the
*        nuclei.
*
         iR_End=iR_Start+1
         Do iR = iR_Start, iEnd_R, -1
            kSet = iAngular_Grid(iR)
*
*---------- Save new iSet for the first point of this
*           subrange.
*
            If (iR.eq.iR_Start) iSet  = kSet
*
*---------- Branch out if we hit on a range where we can reduce the
*           angular grid further.
*
            If (kSet.ne.iSet) Go To 888
*
*---------- Update inner index
*
            iR_End = iR
         End Do
*
  888    Continue
*
*------ Branch out if subrange is outside the box
*
         R_Value_Min = Work(ip_Rx+iROff(1,iR_End))
         If (R_Value_Min.gt.R_box_Max) Go To 8887
         R_Value_Max = Work(ip_Rx+iROff(1,iR_Start))
         If (R_Value_Max.lt.R_box_Min) Go To 8887
*
c        Write (*,*) 'Selected range:',iR_End, iR_Start
c        Write (*,*) 'l_max=',Info_Ang(1,iSet)
*                                                                      *
************************************************************************
*                                                                      *
*---- Angular loop
*
      Do iPoint=1,Info_Ang(2,iSet)
#ifdef _DEBUG_
         If (Debug) Then
            Write (6,*) 'X,Y,Z*',x_a(iPoint,iSet),
     &                           y_a(iPoint,iSet),
     &                           z_a(iPoint,iSet)
         End If
#endif
*
         If (.Not.InBox.and..Not.Check(iOpt_Angular,2)) Then
c           Write (*,*) 'Select angular points!'
            If (x_a(iPoint,iset).lt.xyz0(1,1) .or.
     &          x_a(iPoint,iset).gt.xyz0(1,2) .or.
     &          y_a(iPoint,iset).lt.xyz0(2,1) .or.
     &          y_a(iPoint,iset).gt.xyz0(2,2) .or.
     &          z_a(iPoint,iset).lt.xyz0(3,1) .or.
     &          z_a(iPoint,iset).gt.xyz0(3,2) ) Go To 7777
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Radial loop over the reduced range
*
         Do iR = iR_End,iR_Start
            Radius=Work(ip_Rx+iROff(1,iR))
*           In the atomic referential
            xpt=Radius*x_a(iPoint,iSet)
            ypt=Radius*y_a(iPoint,iSet)
            zpt=Radius*z_a(iPoint,iSet)
*           In the system referential
            x=xpt+x_NQ
            y=ypt+y_NQ
            z=zpt+z_NQ
#ifdef _DEBUG_
            If (Debug) Then
               Write (6,*) 'Radius=',Radius
               Write (6,*) ' x,y,z:',x,y,z
               Write (6,*) x_NQ, y_NQ, z_NQ
            End If
#endif
*
*---------- Check if the point is inside the box.
*           Points on the border of boxes are shared.
*
            If (
     &           (x.ge.x_min_).and.(x.le.x_max_) .and.
     &           (y.ge.y_min_).and.(y.le.y_max_) .and.
     &           (z.ge.z_min_).and.(z.le.z_max_)
     &         ) Then
*
*------------- For shared points modify the weight.
*
               Fact=One
               If (x.eq.x_min_) Fact=Fact*Half
               If (y.eq.y_min_) Fact=Fact*Half
               If (z.eq.z_min_) Fact=Fact*Half
               If (x.eq.x_max_) Fact=Fact*Half
               If (y.eq.y_max_) Fact=Fact*Half
               If (z.eq.z_max_) Fact=Fact*Half
#ifdef _DEBUG_
               If (Debug) Then
                  Write (6,*) x_a(iPoint,iSet),
     &                        y_a(iPoint,iSet),
     &                        z_a(iPoint,iSet)
                  Write (6,*) 'x:', xyz0(1,1),xyz0(1,2)
                  Write (6,*) 'y:', xyz0(2,1),xyz0(2,2)
                  Write (6,*) 'z:', xyz0(3,1),xyz0(3,2)
                  Write (6,*) ' Inside:',x,y,z
               End If
#endif
*
               number_of_grid_points=number_of_grid_points+1
*              Radial weight
               weight=Work(ip_Rx+iROff(2,iR))
*              Combine the radial and angular weight
               w_g=weight*w_a(iPoint,iSet)
               Grid(1,number_of_grid_points)= x
               Grid(2,number_of_grid_points)= y
               Grid(3,number_of_grid_points)= z
*              Compute the partitioning weight
               Weights(number_of_grid_points)=w_g*Fact
            End If
*                                                                      *
************************************************************************
*                                                                      *
*    The call of Do_Batch is done if the buffer is full :
*     (number_of_grid_points.gt.nGrid)
*                                                                      *
************************************************************************
*                                                                      *
         If (number_of_grid_points.gt.mGrid) Then
            Call WarningMessage(2,'Subblock: Buffer overflowed!;'//
     &                  'Try a larger buffer size!')
            Call Abend()
         End If
         If (number_of_grid_points.gt.nGrid) Then
*
*---------- Dump grid information to disk
            nBatch = nBatch + 1
            If (nBatch.gt.nBatch_Max) Then
               Call WarningMessage(2,'Subblock: nBatch.gt.nBatch_Max')
               Call Abend()
            End If
            iBatchInfo(1,nBatch)=iDisk_Grid
            iBatchInfo(2,nBatch)=number_of_grid_points
            iBatchInfo(3,nBatch)=iNQ
*
            Call dDaFile(Lu_Grid,1,Grid,3*number_of_grid_points,
     &                   iDisk_Grid)
*
*---------- Generate weights
*
            Call W(Grid(1,iStrt),ilist_p,Weights(iStrt),list_p,
     &             nList_p,number_of_grid_points-iStrt+1)
            Call dDaFile(Lu_Grid,1,Weights,number_of_grid_points,
     &                   iDisk_Grid)
*
            ntotgp=ntotgp+number_of_grid_points
C           Write (*,*) 'ntotgp=',ntotgp
            number_of_grid_points=0
            iStrt=number_of_grid_points+1
         End If
*
         End Do           ! iR, Radial loop
 7777    Continue
      End Do              ! iPoint
*                                                                      *
************************************************************************
*                                                                      *
 8887 Continue
      iR_Start=iR_End-1
      If (Angular_Prunning.eq.On.and.iR_End.ne.iEnd_R) Go To 999
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate weights
*
*     Write (6,*) 'number_of_grid_points=',number_of_grid_points
*     Write (6,*) 'iStrt=',iStrt
      If (number_of_grid_points-iStrt+1.gt.0)
     &   Call W(Grid(1,iStrt),ilist_p,Weights(iStrt),list_p,
     &          nList_p,number_of_grid_points-iStrt+1)
*                                                                      *
************************************************************************
*                                                                      *
*---- Process batch if not processed yet.
*
      If (Process.and.number_of_grid_points.gt.0) Then
*
*------- Dump grid information to disk
         nBatch = nBatch + 1
         If (nBatch.gt.nBatch_Max) Then
            Call WarningMessage(2,'Subblock: nBatch.gt.nBatch_Max')
            Call Abend()
         End If
         iBatchInfo(1,nBatch)=iDisk_Grid
         iBatchInfo(2,nBatch)=number_of_grid_points
         iBatchInfo(3,nBatch)=iNQ
         Call dDaFile(Lu_Grid,1,Grid,3*number_of_grid_points,
     &                iDisk_Grid)
         Call dDaFile(Lu_Grid,1,Weights,number_of_grid_points,
     &                iDisk_Grid)
*
         ntotgp=ntotgp+number_of_grid_points
C        Write (*,*) 'ntotgp=',ntotgp
         number_of_grid_points=0
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call QExit('Subblock')
#endif
      Return
      End

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
!               2024, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Get_Subblock(Kernel,Func,ixyz,Maps2p,list_s,list_exp,list_bas,nShell,nSym,list_p,nNQ,FckInt,nFckDim,nFckInt,nD,mGrid, &
                        nP2_ontop,Do_Mo,Do_Grad,Grad,nGrad,mAO,mdRho_dR,EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)
!***********************************************************************
!                                                                      *
! Object: to generate the list of the shell and exponent that have an  *
!         influence on a subblock                                      *
!                                                                      *
!     Called from: Drvnq_inner                                         *
!                                                                      *
!     Author: Roland Lindh,                                            *
!             Dept of Chemical Physics,                                *
!             University of Lund, Sweden                               *
!             August 1999                                              *
!***********************************************************************

use iSD_data, only: iSD
use Basis_Info, only: Shells
use Center_Info, only: dc
use nq_Grid, only: dRho_dR, dW_dR, Grid, IndGrd, iTab, kAO, List_G, nR_Eff, TabAO, TabAO_Pack, TabAO_Short, Weights
use nq_Info, only: On, Rotational_Invariance
use NQ_Structure, only: NQ_Data
use nq_MO, only: nMOs
use nq_Info, only: Block_Size, Grid_Type, Moving_Grid, nPot1, nTotGP, nx, ny, nz, Off, On, Threshold, x_min, y_min, z_min
use Grid_On_Disk, only: Grid_Status, GridInfo, iBatchInfo, iDisk_Grid, Lu_Grid, LuGridFile, nBatch, nBatch_Max, Use_Old, WriteGrid
use DFT_Functionals, only: DFT_FUNCTIONAL
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
!#define _ANALYSIS_
#if defined (_DEBUGPRINT_) || defined (_ANALYSIS_)
use Definitions, only: u6
#endif

implicit none
procedure(DFT_FUNCTIONAL) :: Kernel
integer(kind=iwp), intent(in) :: ixyz, nShell, nSym, Maps2p(nShell,0:nSym-1), nNQ, nFckDim, nFckInt, nD, mGrid, nP2_ontop, nGrad, &
                                 mAO, mdRho_dR, nTmpPUVX
real(kind=wp), intent(inout) :: Func, FckInt(nFckInt,nFckDim), Grad(nGrad), EG_OT(nTmpPUVX), PDFTPot1(nPot1), PDFTFocI(nPot1), &
                                PDFTFocA(nPot1)
integer(kind=iwp), intent(out) :: list_s(2,nSym*nShell), list_exp(nSym*nShell), list_bas(2,nSym*nShell), list_p(nNQ)
logical(kind=iwp), intent(in) :: Do_Mo, Do_Grad
integer(kind=iwp) :: i, iAng, iBatch, iCar, iCmp, iExp, iGrad, iIndex, ilist_p, ilist_s, iNQ, iPseudo, iShell, iShll, iSkal, iSym, &
                     ix, iy, iyz, iz, jDisk_Grid, jlist_s, jNQ, jShell, jSym, kNQ, mdci, nAOs, nAOs_Eff, nBfn, nDegi, nExpTmp, &
                     nGrad_Eff, nIndex, nlist_d, nlist_p, nlist_s, nogp, NrBas, NrBas_Eff, NrExp, nTabMO, nTabSO, &
                     number_of_grid_points, nx_Roots, ny_Roots, nz_Roots
real(kind=wp) :: r, R_Box_Max, R_Box_Min, RMax, Roots(3,3), ValExp, X, x_box_max, x_box_min, x_max_, x_min_, x_NQ, Xref, &
                 xyz0(3,2), y, y_box_max, y_box_min, y_max_, y_min_, y_NQ, z, z_box_max, z_box_min, z_max_, z_min_, z_NQ
logical(kind=iwp) :: Added, InBox, More_to_come
integer(kind=iwp), allocatable :: Indx(:), invlist(:), invlist_g(:,:)
real(kind=wp), allocatable :: dPB(:,:,:), dW_Temp(:,:), R2_Trial(:), TabMO(:), TabSO(:)
integer(kind=iwp), external :: nBas_Eff, NrOpr
real(kind=wp), external :: Eval_RMax

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter Get_Subblock'
#endif

! Resolve triplet index

iyz = 1+(ixyz-1)/nx
ix = ixyz-(iyz-1)*nx
iz = 1+(iyz-1)/ny
iy = iyz-(iz-1)*ny

! Get the extreme coordinates of the box.

x_min_ = x_min+real(ix-2,kind=wp)*Block_Size
x_max_ = x_min_+Block_Size
y_min_ = y_min+real(iy-2,kind=wp)*Block_Size
y_max_ = y_min_+Block_Size
z_min_ = z_min+real(iz-2,kind=wp)*Block_Size
z_max_ = z_min_+Block_Size
if (ix == 1) x_min_ = -1.0e99_wp
if (ix == nx) x_max_ = 1.0e99_wp
if (iy == 1) y_min_ = -1.0e99_wp
if (iy == ny) y_max_ = 1.0e99_wp
if (iz == 1) z_min_ = -1.0e99_wp
if (iz == nz) z_max_ = 1.0e99_wp

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Block_Size=',Block_Size
write(u6,*) 'ix,iy,iz=',ix,iy,iz
write(u6,*) 'x_min_,x_max_',x_min_,x_max_
write(u6,*) 'y_min_,y_max_',y_min_,y_max_
write(u6,*) 'z_min_,z_max_',z_min_,z_max_
write(u6,*) 'nNQ=',nNQ
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Generate list over atoms which contribute to the subblock.
!                                                                      *
!***********************************************************************
!                                                                      *
ilist_p = 0
call mma_allocate(R2_Trial,nNQ,Label='R2_Trial')
do iNQ=1,nNQ
  ! Get the coordinates of the partitioning
  x_NQ = NQ_Data(iNQ)%Coor(1)
  y_NQ = NQ_Data(iNQ)%Coor(2)
  z_NQ = NQ_Data(iNQ)%Coor(3)

  R2_Trial(iNQ) = max(Zero,x_NQ-x_max_,x_min_-x_NQ)**2+max(Zero,y_NQ-y_max_,y_min_-y_NQ)**2+max(Zero,z_NQ-z_max_,z_min_-z_NQ)**2

  ! 1) R2_Trial == 0      : center is in the box
  ! 2) R2_Trial < RMax**2 : atomic grid of this center extends inside the box.

  RMax = NQ_Data(iNQ)%R_max
  if (R2_Trial(iNQ) <= RMax**2) then
    ilist_p = ilist_p+1
    list_p(ilist_p) = iNQ
  end if
end do
! nlist_p is the number of centers with a grid extending into this subblock
! nlist_d is the number of centers that affect the density or weights on the grid points in this subbblock (for derivatives)
nlist_p = ilist_p
nlist_d = nlist_p
if (nlist_p == 0) then
  call mma_deallocate(R2_Trial)
  return
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'Get_Subblock: List_p:',List_p
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Generate list over shells which contribute to the subblock.          *
!                                                                      *
!***********************************************************************
!                                                                      *
ilist_s = 0
do iShell=1,nShell
# ifdef _DEBUGPRINT_
  write(u6,*) 'iShell,nShell=',iShell,nShell
# endif
  NrExp = iSD(5,iShell)
  iAng = iSD(1,iShell)
  iShll = iSD(0,iShell)
  NrBas = iSD(3,iShell)
  mdci = iSD(10,iShell)
  nDegi = nSym/dc(mdci)%nStab

  do jSym=0,nDegi-1
    iSym = dc(mdci)%iCoSet(jSym,0)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iSym,nDegi-1=',iSym,nDegi-1
#   endif

    iNQ = Maps2p(iShell,NrOpr(iSym))
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iNQ=',iNQ
    write(u6,*) 'InBox(iNQ)=',R2_Trial(iNQ) == Zero
#   endif

    Added = .false.

    ! 1) the center of this shell is inside the box

    if (R2_Trial(iNQ) == Zero) then
      ilist_s = ilist_s+1
      list_s(1,ilist_s) = iShell
      list_s(2,ilist_s) = iSym
      list_exp(ilist_s) = NrExp
      list_bas(1,ilist_s) = NrBas
      Added = .true.
#     ifdef _ANALYSIS_
      write(u6,*) ' Shell is in box, ilist_s: ',ilist_s
#     endif
    else
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Passed here!'
      write(u6,*) 'Threshold:',Threshold
#     endif

      ! 2) the Gaussian has a grid point that extends inside the box.
      !    The Gaussians are ordered from the most diffuse to the most
      !    contracted.

      nExpTmp = 0
      do iExp=1,NrExp
        ! Get the value of the exponent
        ValExp = Shells(iShll)%Exp(iExp)
        ! If the exponent has an influence then increase the
        ! number of active exponents for this shell, else
        ! there is no other active exponent (they are ordered)
        RMax = Eval_RMax(ValExp,iAng,Threshold)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'iShell,iNQ=',iShell,iNQ
        write(u6,*) 'ValExp,iExp=',ValExp,iExp
        write(u6,*) 'RMax_Exp=',Eval_RMax(ValExp,iAng,Threshold)
        write(u6,*) 'RMax=',RMax
        write(u6,*) 'R2_Trial(iNQ),RMax**2=',R2_Trial(iNQ),RMax**2
#       endif
        if (R2_Trial(iNQ) > RMax**2) exit
        nExpTmp = nExpTmp+1
      end do  ! iExp
      if (nExpTmp /= 0) then
        ilist_s = ilist_s+1
        list_s(1,ilist_s) = iShell
        list_s(2,ilist_s) = iSym
        list_exp(ilist_s) = nExpTmp
        Added = .true.

        ! Examine if contracted basis functions can be ignored.
        ! This will be the case for segmented basis sets.

        list_bas(1,ilist_s) = nBas_Eff(NrExp,NrBas,Shells(iShll)%pCff,list_exp(ilist_s))
#       ifdef _ANALYSIS_
        write(u6,*) ' Shell is included, ilist_s: ',ilist_s
        write(u6,*) ' nExpTmp=',nExpTmp
        write(u6,*) 'R2_Trial(iNQ),RMax**2=',R2_Trial(iNQ),RMax**2
#       endif
      end if
    end if
    ! We may need to increase the number of centers considered for derivatives
    if (Added .and. Do_Grad .and. (Grid_Type == Moving_Grid)) then
      do ilist_p=1,nlist_d
        if (list_p(ilist_p) == iNQ) exit
      end do
      if (ilist_p == nlist_d+1) then
        ! If this center is not found in the current list_p, add it
        nlist_d = ilist_p
        list_p(ilist_p) = iNQ
      end if
    end if
  end do ! iSym
end do   ! iShell
nlist_s = ilist_s
#ifdef _DEBUGPRINT_
write(u6,*) 'nList_s,nList_p,nList_d=',nList_s,nList_p,nList_d
#endif
if (nlist_s == 0) then
  call mma_deallocate(R2_Trial)
  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate index arrays to address the density matrix, which is in
! the full shell and the reduced shell over which in the basis
! functions will be evaluated.

nIndex = 0
do ilist_s=1,nlist_s
  iShell = list_s(1,ilist_s)
  NrBas_Eff = list_bas(1,ilist_s)
  iCmp = iSD(2,iShell)
  nIndex = nIndex+NrBas_Eff*iCmp
end do

call mma_allocate(Indx,nIndex,Label='Indx')

iIndex = 1
nAOs = 0
nAOs_Eff = 0
do ilist_s=1,nlist_s
  iShell = list_s(1,ilist_s)
  NrBas = iSD(3,iShell)
  NrBas_Eff = list_bas(1,ilist_s)
  iCmp = iSD(2,iShell)
  nAOs = nAOs+NrBas*iCmp
  nAOs_Eff = nAOs_Eff+NrBas_Eff*iCmp
  list_bas(2,ilist_s) = iIndex
  call Do_Index(Indx(iIndex),NrBas,NrBas_Eff,iCmp)
  iIndex = iIndex+NrBas_Eff*iCmp
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'Contribution to the subblock :'
write(u6,*) 'NQ :',(list_p(ilist_p),ilist_p=1,nlist_p)
if (nlist_d > nlist_p) then
  write(u6,*) '   :',(list_p(ilist_p),ilist_p=nlist_p+1,nlist_d)
end if
write(u6,*) 'Sh :',(list_s(1,ilist_s),ilist_s=1,nlist_s)
write(u6,*) '   :',(list_s(2,ilist_s),ilist_s=1,nlist_s)
write(u6,*) 'Exp:',(list_exp(ilist_s),ilist_s=1,nlist_s)
#endif

nBfn = 0
do ilist_s=1,nlist_s
  iSkal = list_s(1,ilist_s)
  NrBas_Eff = list_bas(1,ilist_s)
  iCmp = iSD(2,iSkal)
  nBfn = nBfn+NrBas_Eff*iCmp
end do

if (Do_MO) then
  nTabMO = mAO*nMOs*mGrid
  nTabSO = mAO*nMOs*mGrid
else
  nTabMO = 1
  nTabSO = 1
end if
call mma_allocate(TabMO,nTabMO,Label='TabMO')
call mma_allocate(TabSO,nTabSO,Label='TabSO')
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(invlist,nNQ,Label='invlist')
invlist(:) = -1
! With rotational invariance, *all* centers with nonzero charge affect the grid, because they affect the orientation
if (Rotational_Invariance == On) then
  do iNQ=1,nNQ
    if (NQ_Data(iNQ)%Atom_Nr /= 0) invlist(iNQ) = 0
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (Grid_Status /= Use_Old) then
  call mma_allocate(iBatchInfo,3,nBatch_Max,label='iBatchInfo')
  iBatchInfo(:,:) = 0
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! For each partition active in the subblock create the grid
  ! and perform the integration on it.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  number_of_grid_points = 0
  nBatch = 0
  do ilist_p=1,nlist_p
    iNQ = list_p(ilist_p)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'ilist_p=',ilist_p
    write(u6,*) 'Get_SubBlock: iNQ=',iNQ
#   endif

    ! Get the coordinates of the partition
    x_NQ = NQ_Data(iNQ)%Coor(1)
    y_NQ = NQ_Data(iNQ)%Coor(2)
    z_NQ = NQ_Data(iNQ)%Coor(3)
    ! Get the maximum radius on which we have to integrate for the partition
    RMax = NQ_Data(iNQ)%R_max

    call Box_On_Sphere(x_Min_-x_NQ,x_Max_-x_NQ,y_Min_-y_NQ,y_Max_-y_NQ,z_Min_-z_NQ,z_Max_-z_NQ,xyz0(1,1),xyz0(1,2),xyz0(2,1), &
                       xyz0(2,2),xyz0(3,1),xyz0(3,2))
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Establish R_Box_Max and R_Box_Min, the longest and the shortest
    ! distance from the origin of the atomic grid to a point in the box

    R_Box_Max = Zero
    R_Box_Min = RMax

    x_box_min = x_min_-x_NQ
    x_box_max = x_max_-x_NQ
    y_box_min = y_min_-y_NQ
    y_box_max = y_max_-y_NQ
    z_box_min = z_min_-z_NQ
    z_box_max = z_max_-z_NQ

    Roots(1,1) = x_box_min
    Roots(2,1) = x_box_max
    if (x_box_max*x_box_min < Zero) then
      nx_Roots = 3
      Roots(3,1) = Zero
    else
      nx_Roots = 2
    end if

    Roots(1,2) = y_box_min
    Roots(2,2) = y_box_max
    if (y_box_max*y_box_min < Zero) then
      ny_Roots = 3
      Roots(3,2) = Zero
    else
      ny_Roots = 2
    end if

    Roots(1,3) = z_box_min
    Roots(2,3) = z_box_max
    if (z_box_max*z_box_min < Zero) then
      nz_Roots = 3
      Roots(3,3) = Zero
    else
      nz_Roots = 2
    end if

    ! Check all stationary points

    do ix=1,nx_Roots
      x = Roots(ix,1)
      do iy=1,ny_Roots
        y = Roots(iy,2)
        do iz=1,nz_Roots
          z = Roots(iz,3)

          r = sqrt(x**2+y**2+z**2)

          R_Box_Max = max(R_Box_Max,r)
          R_Box_Min = min(R_Box_Min,r)

        end do
      end do
    end do

    if (abs(R_Box_Min) < 1.0e-12_wp) R_Box_Min = Zero
    R_Box_Max = R_Box_Max+1.0e-15_wp
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Get_Subblock ----> Subblock'
#   endif

    ! Note that in gradient calculations we process the grid points for
    ! each atomic grid separately in order to used the translational
    ! invariance on the atomic contributions to the gradient.

    InBox = R2_Trial(iNQ) == Zero
    call Subblock(iNQ,x_NQ,y_NQ,z_NQ,InBox,x_min_,x_max_,y_min_,y_max_,z_min_,z_max_,nNQ,invlist,Grid,Weights,mGrid,.true., &
                  number_of_grid_points,R_Box_Min,R_Box_Max,xyz0,NQ_Data(iNQ)%Angular,nR_Eff(iNQ))

#   ifdef _DEBUGPRINT_
    write(u6,*) 'Subblock ----> Get_Subblock'
#   endif
  end do
  GridInfo(1,ixyz) = iDisk_Grid
  GridInfo(2,ixyz) = nBatch
  call iDaFile(Lu_Grid,1,iBatchInfo,3*nBatch,iDisk_Grid)
  call iDaFile(Lu_Grid,1,invlist,nNQ,iDisk_Grid)
  call mma_deallocate(iBatchInfo)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
iDisk_Grid = GridInfo(1,ixyz)
nBatch = GridInfo(2,ixyz)
call mma_allocate(iBatchInfo,3,nBatch,label='iBatchInfo')
call iDaFile(Lu_Grid,2,iBatchInfo,3*nBatch,iDisk_Grid)
call iDaFile(Lu_Grid,2,invlist,nNQ,iDisk_Grid)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate indexation of which shells contributes to which centers
! and center index for each gradient contribution which is computed.

nGrad_Eff = 0
if (Do_Grad) then
  call mma_allocate(invlist_g,3,nNQ,label='invlist_g')
  invlist_g(:,:) = -1
  List_G(:,:) = 0
  do ilist_s=1,nlist_s
    iShell = list_s(1,ilist_s)
    iSym = list_s(2,ilist_s)
    mdci = iSD(10,iShell)
    iNQ = Maps2p(iShell,NrOpr(iSym))
    do iCar=1,3
      if (List_G(iCar,ilist_s) /= 0) cycle
      iPseudo = iSD(12,iShell)
      if ((iSD(15+iCar,iShell) /= 0) .or. (iPseudo == 1)) then
        nGrad_Eff = nGrad_Eff+1

        ! For pseudo centers note that there will not be a
        ! gradient computed for this center.

        if (iPseudo == 0) then
          IndGrd(nGrad_Eff) = iSD(15+iCar,iShell)
        else
          IndGrd(nGrad_Eff) = -1
        end if
        List_G(iCar,ilist_s) = nGrad_Eff
        iTab(1,nGrad_Eff) = iCar
        iTab(3,nGrad_Eff) = iNQ
        invlist_g(iCar,iNQ) = nGrad_Eff
        kNQ = Maps2p(iShell,0)
        Xref = NQ_Data(kNQ)%Coor(iCar)
        X = NQ_Data(iNQ)%Coor(iCar)
        if (X == Xref) then
          iTab(4,nGrad_Eff) = dc(mdci)%nStab
        else
          iTab(4,nGrad_Eff) = -dc(mdci)%nStab
        end if

        ! Find all other shells which contribute to the same gradient.

        do jlist_s=ilist_s+1,nlist_s
          jShell = list_s(1,jlist_s)
          if ((iSD(15+iCar,iShell) == iSD(15+iCar,jShell)) .and. (iSym == list_s(2,jlist_s))) List_G(iCar,jlist_s) = nGrad_Eff
        end do

      else if (iSD(15+iCar,iShell) == 0) then

        ! Include derivatives which will be used for
        ! the translational invariance equation but which do not
        ! contribute directly to a symmetry adapted gradient

        nGrad_Eff = nGrad_Eff+1
        IndGrd(nGrad_Eff) = -1
        List_G(iCar,ilist_s) = nGrad_Eff
        iTab(1,nGrad_Eff) = iCar
        iTab(3,nGrad_Eff) = iNQ
        invlist_g(iCar,iNQ) = nGrad_Eff
        iTab(4,nGrad_Eff) = dc(mdci)%nStab

        ! Find all other shells which contribute to the same gradient.

        do jlist_s=ilist_s+1,nlist_s
          jShell = list_s(1,jlist_s)
          jSym = list_s(2,jlist_s)
          jNQ = Maps2p(jShell,NrOpr(jSym))
          if (iNQ == jNQ) List_G(iCar,jlist_s) = nGrad_Eff
        end do
      end if
    end do
  end do

  if (Grid_Type == Moving_Grid) then
    ! Add the centers that contribute through modification of the weights
    ! (i.e. those for which invlist is 0)
    ! Unfortunately, with rotational invariance this probably includes all centers
    ! (This should be revised with symmetry)
    do ilist_p=1,nlist_d
      invlist(list_p(ilist_p)) = ilist_p
    end do
    do iNQ=1,nNQ
      if (invlist(iNQ) < 0) cycle
      if (invlist(iNQ) == 0) then
        nlist_d = nlist_d+1
        list_p(nlist_d) = iNQ
        invlist(iNQ) = nlist_d
      end if
      do iCar=1,3
        if (invlist_g(iCar,iNQ) > 0) cycle
        nGrad_Eff = nGrad_Eff+1
        if (NQ_Data(iNQ)%Grad_idx(iCar) /= 0) then
          IndGrd(nGrad_Eff) = NQ_Data(iNQ)%Grad_idx(iCar)
        else
          IndGrd(nGrad_Eff) = -1
        end if
        iTab(1,nGrad_Eff) = iCar
        iTab(3,nGrad_Eff) = iNQ
        iTab(4,nGrad_Eff) = NQ_Data(iNQ)%Fact
        invlist_g(iCar,iNQ) = nGrad_Eff
      end do
    end do

    call mma_allocate(dW_dR,nGrad_Eff,mGrid,Label='dW_dR')
    call mma_allocate(dW_Temp,3,nlist_d,Label='dW_Temp')
    call mma_allocate(dPB,3,nlist_d,nlist_d,Label='dPB')
  end if
  call mma_deallocate(invlist_g)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if ((.not. Do_Grad) .or. (nGrad_Eff /= 0)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Process grid points

  iBatch = 0
  nogp = 0
  outer: do
    iBatch = iBatch+1
    if (iBatch > nBatch) exit outer
    do
      jDisk_Grid = iBatchInfo(1,iBatch)
      number_of_grid_points = iBatchInfo(2,iBatch)

      iNQ = iBatchInfo(3,iBatch)
#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'iNQ=',iNQ
      write(u6,*)
#     endif

      if (nogp+number_of_grid_points <= mGrid) then
        call dDaFile(Lu_Grid,2,Grid(1,nogp+1),3*number_of_grid_points,jDisk_Grid)
        call dDaFile(Lu_Grid,2,Weights(nogp+1),number_of_grid_points,jDisk_Grid)
        nogp = nogp+number_of_grid_points

        ! If this is not a gradient evaluation read next buffer if the
        ! current one is not the last one.

        More_to_Come = .false.
        if ((.not. Do_Grad) .and. (iBatch /= nBatch)) cycle outer
      else
        More_to_Come = .true.
      end if

      ! Here if it is a gradient evaluation or we have a buffer to process.

      if (Do_Grad) then
        iTab(2,1:nGrad_Eff) = On
        if (Grid_Type == Moving_Grid) then
          do iGrad=1,nGrad_Eff
            jNQ = iTab(3,iGrad)
            if (iNQ == jNQ) iTab(2,iGrad) = Off
          end do

          ! Generate derivative with respect to the weights if needed.

          call dWdR(Grid,iNQ,Weights,list_p,nlist_d,invlist,dW_dR,nGrad_Eff,iTab,dW_Temp,dPB,number_of_grid_points)
        end if
      end if

      call mma_allocate(TabAO,mAO,nogp,nBfn,Label='TabAO')
      if (Do_Grad) call mma_allocate(TabAO_Short,kAO,nogp,nBfn,Label='TabAO_Short')
      TabAO_Pack(1:mAO*nogp*nBfn) => TabAO(:,:,:)
      if (Do_Grad) then
        call mma_allocate(dRho_dR,mdRho_dR,nogp,nGrad_eff,Label='dRho_dR')
      else
        call mma_allocate(dRho_dR,1,1,1,Label='dRho_dR')
      end if

      call Do_Batch(Kernel,Func,nogp,list_s,nlist_s,list_exp,list_bas,Indx,nIndex,FckInt,nFckDim,nFckInt,mAO,nD,nP2_ontop,Do_Mo, &
                    TabMO,TabSO,nMOs,Do_Grad,Grad,nGrad,mdRho_dR,nGrad_Eff,iNQ,EG_OT,nTmpPUVX,PDFTPot1,PDFTFocI,PDFTFocA)

      call mma_deallocate(dRho_dR,safe='*')
      call mma_deallocate(TabAO_Short,safe='*')
      nullify(TabAO_Pack)
      call mma_deallocate(TabAO)

      nTotGP = nTotGP+nogp
      if (WriteGrid) then
        ! update the "LuGridFile":
        do i=1,nogp
          write(LuGridFile,'(3ES24.14,1x,ES24.14)') Grid(:,i),Weights(i)
        end do
      end if
      nogp = 0
      if (.not. More_To_Come) exit
    end do
  end do outer
end if
call mma_deallocate(iBatchInfo)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(R2_Trial)
call mma_deallocate(invlist)
call mma_deallocate(Indx)
call mma_deallocate(TabMO,safe='*')
call mma_deallocate(TabSO,safe='*')
if (Do_Grad .and. (Grid_Type == Moving_Grid)) then
  call mma_deallocate(dPB)
  call mma_deallocate(dW_Temp)
  call mma_deallocate(dW_dR)
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Get_Subblock

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

module nq_Info

use Definitions, only: wp, iwp

implicit none
private

!     R_Max : Maximum radius associated with the i'th center for the
!             radial loop.
! Integer
!
!     L_Quad  : Value of the angular momentum for which the grid is
!               generated unless the quadrature is pruned.
!     nR      : Initial number of radial points for which the grid is
!               generated.
!     nAtoms  : Number of atoms in the system.
!     nTotGP  : Total number of grid points generated during the program.

integer(kind=iwp), parameter :: LMax_NQ = 62, &
                                Other_Type = 0, LDA_Type = 1, GGA_Type = 2, meta_GGA_Type1 = 3, meta_GGA_Type2 = 4, &
                                Fixed_Grid = 0, Moving_Grid = 1, &
                                On = 1, Off = 0
integer(kind=iwp) :: Angular_Pruning, Functional_Type, Grid_Type, IOff_Ash(0:7), IOff_Bas(0:7), IOff_BasAct(0:7), iOpt_Angular, &
                     L_Quad, L_Quad_save, mBas(0:7), mIrrep, mOrb(0:7), mRad, nAngularGrids, nAsh(0:7), NASHT, nAtoms, ndc, &
                     nFro(0:7), nISh(0:7), nOrbt, nPot1, nPot2, NQ_Direct, nR, nR_Save, nTotGP, number_of_subblocks, &
                     nUVX(0:7,0:7,0:7), nUVXt, nVX(0:7,0:7), nVXt, nx, ny, nz, OffBas(0:7), OffBas2(0:7), OffBasFro(0:7), &
                     OffOrb(0:7), OffOrb2(0:7), OffOrbTri(0:7), OffPUVX(0:7), OffUVX(0:7,0:7,0:7), OffVX(0:7,0:7), Packing, &
                     Rotational_Invariance
real(kind=wp) :: Block_Size, Crowding, Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Energy_integrated, Fade, &
                 Grad_I, R_Max(0:LMax_NQ), T_Y, Tau_I, ThrC, Threshold, Threshold_save, x_min, y_min, z_min
character(len=10) :: Quadrature
character(len=8) :: MBC

public :: Angular_Pruning, Block_Size, Crowding, Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Energy_integrated, &
          Fade, Fixed_Grid, Functional_Type, GGA_Type, Grad_I, Grid_Type, IOff_Ash, IOff_Bas, IOff_BasAct, iOpt_Angular, L_Quad, &
          L_Quad_save, LDA_Type, LMax_NQ, mBas, MBC, meta_GGA_Type1, meta_GGA_Type2, mIrrep, mOrb, Moving_Grid, mRad, &
          nAngularGrids, nAsh, NASHT, nAtoms, ndc, nFro, nISh, nOrbt, nPot1, nPot2, NQ_Direct, NQ_Info_Dmp, NQ_Info_Get, nR, &
          nR_Save, nTotGP, number_of_subblocks, nUVX, nUVXt, nVX, nVXt, nx, ny, nz, Off, OffBas, OffBas2, OffBasFro, OffOrb, &
          OffOrb2, OffOrbTri, OffPUVX, OffUVX, OffVX, On, Other_Type, Packing, Quadrature, R_Max, Rotational_Invariance, T_Y, &
          Tau_I, ThrC, Threshold, Threshold_save, x_min, y_min, z_min

contains

subroutine NQ_Info_Dmp()

  use fortran_strings, only: char_array, str
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: i, lcDmp
  integer(kind=iwp), allocatable :: iDmp(:)
  real(kind=wp), allocatable :: rDmp(:)
  character, allocatable :: cDmp(:)
  integer(kind=iwp), parameter :: liDmp = 25+5*8, lrDmp = 20+(LMax_NQ+1)

  ! Real Stuff

  call mma_allocate(rDmp,lrDmp,Label='rDmp')
  i = 1
  rDmp(i) = Threshold_save
  i = i+1
  rDmp(i) = Crowding
  i = i+1
  rDmp(i) = Threshold
  i = i+1
  rDmp(i:i+LMax_NQ) = R_Max
  i = i+LMax_NQ+1
  rDmp(i) = Energy_integrated
  i = i+1
  rDmp(i) = Dens_I
  i = i+1
  rDmp(i) = Grad_I
  i = i+1
  rDmp(i) = Tau_I
  i = i+1
  rDmp(i) = Dens_a1
  i = i+1
  rDmp(i) = Dens_b1
  i = i+1
  rDmp(i) = Dens_a2
  i = i+1
  rDmp(i) = Dens_b2
  i = i+1
  rDmp(i) = Dens_t1
  i = i+1
  rDmp(i) = Dens_t2
  i = i+1
  rDmp(i) = Block_Size
  i = i+1
  rDmp(i) = x_min
  i = i+1
  rDmp(i) = y_min
  i = i+1
  rDmp(i) = z_min
  i = i+1
  rDmp(i) = Fade
  i = i+1
  rDmp(i) = ThrC
  i = i+1
  rDmp(i) = T_Y
  i = i+1
  call Put_dArray('Quad_r',rDmp,lrDmp)
  call mma_deallocate(rDmp)

  ! Integer Stuff

  call mma_allocate(iDmp,liDmp,Label='iDmp')
  i = 1
  iDmp(i) = NASHT
  i = i+1
  iDmp(i) = nPot1
  i = i+1
  iDmp(i) = nOrbt
  i = i+1
  iDmp(i) = nPot2
  i = i+1
  iDmp(i) = ndc
  i = i+1
  iDmp(i) = nAngularGrids
  i = i+1
  iDmp(i) = L_Quad_save
  i = i+1
  iDmp(i) = nR_Save
  i = i+1
  iDmp(i) = Angular_Pruning
  i = i+1
  iDmp(i) = nx
  i = i+1
  iDmp(i) = ny
  i = i+1
  iDmp(i) = nz
  i = i+1
  iDmp(i) = number_of_subblocks
  i = i+1
  iDmp(i) = L_Quad
  i = i+1
  iDmp(i) = nR
  i = i+1
  iDmp(i) = nAtoms
  i = i+1
  iDmp(i) = nTotGP
  i = i+1
  iDmp(i) = iOpt_Angular
  i = i+1
  iDmp(i) = mIrrep
  i = i+1
  iDmp(i:i+7) = nISh
  i = i+8
  iDmp(i:i+7) = nAsh
  i = i+8
  iDmp(i:i+7) = mBas
  i = i+8
  iDmp(i) = Functional_type
  i = i+1
  iDmp(i:i+7) = mOrb
  i = i+8
  iDmp(i) = Grid_Type
  i = i+1
  iDmp(i) = Rotational_Invariance
  i = i+1
  iDmp(i) = mRad
  i = i+1
  iDmp(i) = NQ_Direct
  i = i+1
  iDmp(i) = Packing
  i = i+1
  iDmp(i:i+7) = OffPUVX
  i = i+8
  call Put_iArray('Quad_i',iDmp,liDmp)
  call mma_deallocate(iDmp)

  ! Character Stuff

  lcDmp = len(Quadrature)+len(MBC)
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  i = 0
  cDmp(i+1:i+len(Quadrature)) = char_array(Quadrature)
  i = i+len(Quadrature)
  cDmp(i+1:i+len(MBC)) = char_array(MBC)
  i = i+len(MBC)
  call Put_cArray('Quad_c',str(cDmp),lcDmp)
  call mma_deallocate(cDmp)

end subroutine NQ_Info_Dmp

subroutine NQ_Info_Get()

  use fortran_strings, only: str
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: i, lcDmp
  integer(kind=iwp), allocatable :: iDmp(:)
  real(kind=wp), allocatable :: rDmp(:)
  character, allocatable :: cDmp(:)
  integer(kind=iwp), parameter :: liDmp = 25+5*8, lrDmp = 20+(LMax_NQ+1)

  ! Real Stuff

  call mma_allocate(rDmp,lrDmp,Label='rDmp')
  call Get_dArray('Quad_r',rDmp,lrDmp)
  i = 1
  Threshold_save = rDmp(i)
  i = i+1
  Crowding = rDmp(i)
  i = i+1
  Threshold = rDmp(i)
  i = i+1
  R_Max = rDmp(i:i+LMax_NQ)
  i = i+LMax_NQ+1
  Energy_integrated = rDmp(i)
  i = i+1
  Dens_I = rDmp(i)
  i = i+1
  Grad_I = rDmp(i)
  i = i+1
  Tau_I = rDmp(i)
  i = i+1
  Dens_a1 = rDmp(i)
  i = i+1
  Dens_b1 = rDmp(i)
  i = i+1
  Dens_a2 = rDmp(i)
  i = i+1
  Dens_b2 = rDmp(i)
  i = i+1
  Dens_t1 = rDmp(i)
  i = i+1
  Dens_t2 = rDmp(i)
  i = i+1
  Block_Size = rDmp(i)
  i = i+1
  x_min = rDmp(i)
  i = i+1
  y_min = rDmp(i)
  i = i+1
  z_min = rDmp(i)
  i = i+1
  Fade = rDmp(i)
  i = i+1
  ThrC = rDmp(i)
  i = i+1
  T_Y = rDmp(i)
  i = i+1
  call mma_deallocate(rDmp)

  ! Integer Stuff

  call mma_allocate(iDmp,liDmp,Label='iDmp')
  call Get_iArray('Quad_i',iDmp,liDmp)
  i = 1
  NASHT = iDmp(i)
  i = i+1
  NPot1 = iDmp(i)
  i = i+1
  nOrbt = iDmp(i)
  i = i+1
  nPot2 = iDmp(i)
  i = i+1
  ndc = iDmp(i)
  i = i+1
  nAngularGrids = iDmp(i)
  i = i+1
  L_Quad_save = iDmp(i)
  i = i+1
  nR_Save = iDmp(i)
  i = i+1
  Angular_Pruning = iDmp(i)
  i = i+1
  nx = iDmp(i)
  i = i+1
  ny = iDmp(i)
  i = i+1
  nz = iDmp(i)
  i = i+1
  number_of_subblocks = iDmp(i)
  i = i+1
  L_Quad = iDmp(i)
  i = i+1
  nR = iDmp(i)
  i = i+1
  nAtoms = iDmp(i)
  i = i+1
  nTotGP = iDmp(i)
  i = i+1
  iOpt_Angular = iDmp(i)
  i = i+1
  mIrrep = iDmp(i)
  i = i+1
  nISh = iDmp(i:i+7)
  i = i+8
  nAsh = iDmp(i:i+7)
  i = i+8
  mBas = iDmp(i:i+7)
  i = i+8
  Functional_type = iDmp(i)
  i = i+1
  mOrb = iDmp(i:i+7)
  i = i+8
  Grid_Type = iDmp(i)
  i = i+1
  Rotational_Invariance = iDmp(i)
  i = i+1
  mRad = iDmp(i)
  i = i+1
  NQ_Direct = iDmp(i)
  i = i+1
  Packing = iDmp(i)
  i = i+1
  OffPUVX = iDmp(i:i+7)
  i = i+8
  call mma_deallocate(iDmp)

  ! Character Stuff

  lcDmp = len(Quadrature)+len(MBC)
  call mma_allocate(cDmp,lcDmp,Label='cDmp')
  call Get_cArray('Quad_c',cDmp,lcDmp)
  i = 0
  Quadrature = str(cDmp(i+1:i+len(Quadrature)))
  i = i+len(Quadrature)
  MBC = str(cDmp(i+1:i+len(MBC)))
  i = i+len(MBC)
  call mma_deallocate(cDmp)

end subroutine NQ_Info_Get

end module nq_Info

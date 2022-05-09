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
!     nMaxExp : Maximum number of exponents over all the shells.
!     NbrMxBas: Maximum number of basis functions over all the shells.
!     nTotGP  : Total number of grid points generated during the program.
!     iAngMax : Maximum angulart momentum for the system.
!
! Pointer
!
!     ip_ioffsh   : Pointer to the offset of the shells in the
!                   overlap matrix.

integer(kind=iwp), parameter :: LMax_NQ = 62, &
                                Other_Type = 0, LDA_Type = 1, GGA_Type = 2, meta_GGA_Type1 = 3, meta_GGA_Type2 = 4, &
                                Fixed_Grid = 0, Moving_Grid = 1, &
                                On = 1, Off = 0
integer(kind=iwp) :: Angular_Pruning, Functional_Type, Grid_Type, iAngMax, IOff_Ash(0:7), IOff_Bas(0:7), IOff_BasAct(0:7), &
                     iOpt_Angular, ip_ioffsh, ip_nR_Eff, ip_OrbDip(3), ip_R, ipMem, L_Quad, L_Quad_save, maxUVX, mBas(0:7), &
                     mIrrep, mOrb(0:7), mRad, mTmp, nAngularGrids, nAOMax, nAsh(0:7), NASHT, NASHT4, nAtoms, nbrmxbas, ndc, &
                     nFro(0:7), nISh(0:7), nMaxExp, nMem, nOrbt, NPOt1, nPot2, NQ_Direct, nR, nR_Save, nTotGP, &
                     number_of_subblocks, nUVX(0:7,0:7,0:7), nUVXt, nVX(0:7,0:7), nVXt, nx, ny, nz, OffBas(0:7), OffBas2(0:7), &
                     OffBasFro(0:7), OffOrb(0:7), OffOrb2(0:7), OffOrbTri(0:7), OffPUVX(0:7), OffUVX(0:7,0:7,0:7), OffVX(0:7,0:7), &
                     Packing, Rotational_Invariance
real(kind=wp) :: Block_Size, Crowding, Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Energy_integrated, Fade, &
                 Grad_I, R_Max(0:LMax_NQ), T_Y, Tau_I, ThrC, Threshold, Threshold_save, x_max, x_min, y_max, y_min, z_max, z_min
character(len=10) :: Quadrature
character(len=8) :: MBC

public :: Angular_Pruning, Block_Size, Crowding, Dens_a1, Dens_a2, Dens_b1, Dens_b2, Dens_I, Dens_t1, Dens_t2, Energy_integrated, &
          Fade, Fixed_Grid, Functional_Type, GGA_Type, Grad_I, Grid_Type, iAngMax, IOff_Ash, IOff_Bas, IOff_BasAct, iOpt_Angular, &
          ip_ioffsh, ip_nR_Eff, ip_OrbDip, ip_R, ipMem, L_Quad, L_Quad_save, LDA_Type, LMax_NQ, maxUVX, mBas, MBC, meta_GGA_Type1, &
          meta_GGA_Type2, mIrrep, mOrb, Moving_Grid, mRad, mTmp, nAngularGrids, nAOMax, nAsh, NASHT, NASHT4, nAtoms, nbrmxbas, &
          ndc, nFro, nISh, nMaxExp, nMem, nOrbt, NPOt1, nPot2, NQ_Direct, nR, nR_Save, nTotGP, number_of_subblocks, nUVX, nUVXt, &
          nVX, nVXt, nx, ny, nz, Off, OffBas, OffBas2, OffBasFro, OffOrb, OffOrb2, OffOrbTri, OffPUVX, OffUVX, OffVX, On, &
          Other_Type, Packing, Quadrature, R_Max, Rotational_Invariance, T_Y, Tau_I, ThrC, Threshold, Threshold_save, x_max, &
          x_min, y_max, y_min, z_max, z_min

!IFG
integer, public :: iQStrt, iQEnd
common/Quad_i/iQStrt,NASHT,NASHT4,NPot1,nOrbt,nPot2,maxUVX,ndc,nAngularGrids,L_Quad_save,nR_Save,Angular_Pruning,nx,ny,nz, &
              number_of_subblocks,ip_nR_Eff,ip_R,ipMem,nMem,L_Quad,nR,nAtoms,nMaxExp,nTotGP,nbrmxbas,iAngMax,ip_ioffsh, &
              iOpt_Angular,mIrrep,nISh,nAsh,mBas,Functional_type,mOrb,Grid_Type,Rotational_Invariance,mTmp,mRad,nAOMax,NQ_Direct, &
              Packing,OffPUVX,iQEnd,ip_OrbDip,ioff_ash,ioff_bas,ioff_basact,OffBas,OffOrb,OffOrb2,OffOrbTri,OffBas2,OffBasFro, &
              OffUVX,nUVX,nUVXt,OffVX,nVX,nVXt
common/Quad_ii/nFro
real*8, public :: rQStrt, rQEnd
common/Quad_r/rQStrt,Threshold_save,Crowding,Threshold,R_Max,Energy_integrated,Dens_I,Grad_I,Tau_I,Dens_a1,Dens_b1,Dens_a2, &
              Dens_b2,Dens_t1,Dens_t2,Block_Size,x_min,x_max,y_min,y_max,z_min,z_max,Fade,ThrC,T_Y,rQEnd
integer, public :: cQStrt, cQEnd
character Pad_*6
common/Quad_c/cQStrt,Quadrature,MBC,Pad_,cQEnd

end module nq_Info

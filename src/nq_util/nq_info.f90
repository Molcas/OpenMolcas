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
!
!     R_Max : Maximum radius associated with the i'th center for the
!              radial loop.
! Integer
!
!     L_Quad  : Value of the angular momentum for which the grid is
!                generated unless the quadrature is pruned.
!     nR      : Initial number of radial points for which the grid is
!                generated.
!     nAtoms  : Number of atoms in the system.
!     nMaxExp : Maximum number of exponents over all the shells.
!     NbrMxBas: Maximum number of basis functions over all the shells.
!     nTotGP  : Total number of grid points generated during the program.
!     iAngMax : Maximum angulart momentum for the system.
!
! Pointer
!
!     ip_ioffsh   : Pointer to the offset of the shells in the
!                    overlap matrix.
!
Module nq_Info
      Integer LMax_NQ
      Parameter(LMax_NQ=62)
      Integer mBas(0:7),nISh(0:7), nAsh(0:7), nFro(0:7),mOrb(0:7)
      Integer IOff_Ash(0:7),IOff_Bas(0:7),IOff_BasAct(0:7)
      Integer OffOrb(0:7),OffOrb2(0:7),OffOrbTri(0:7),OffBas2(0:7)
      Integer OffBas(0:7),OffBasFro(0:7),OffPUVX(0:7)
      INTEGER OffUVX(0:7,0:7,0:7),nUVX(0:7,0:7,0:7)
      INTEGER OffVX(0:7,0:7),nVX(0:7,0:7)
#include "functional_types.fh"
      Integer Grid_Type, Fixed_Grid, Moving_Grid
      Parameter(Fixed_Grid=0,Moving_Grid=1)
      Integer Angular_Prunning, On, Off, Rotational_Invariance,         &
     &        Functional_Type
      Integer Packing
      Parameter(On=1, Off=0)
      Integer          iQStrt,NASHT,NASHT4,NPOt1,nOrbt,nPot2,maxUVX,    &
     &                          ndc,nUVXt ,nVXt,                        &
     &                           nAngularGrids,                         &
     &                 L_Quad_save, nR_Save,                            &
     &                      nx,ny,nz,number_of_subblocks,               &
     &                 ip_nR_Eff,ip_R,ipMem,nMem,                       &
     &                 L_Quad,nR,nAtoms,nMaxExp,                        &
     &                 nTotGP,nbrmxbas,iAngMax,ip_ioffsh,               &
     &                 iOpt_Angular,     mIrrep,                        &
     &                            mTmp,mRad,nAOMax,                     &
     &                            NQ_Direct,                            &
     &                 iQEnd,ip_OrbDip(3)
      Common /Quad_i / iQStrt,NASHT,NASHT4,NPot1,nOrbt,nPot2,maxUVX,    &
     &                          ndc,                                    &
     &                           nAngularGrids,                         &
     &                 L_Quad_save, nR_Save, Angular_Prunning,          &
     &                      nx,ny,nz,number_of_subblocks,               &
     &                 ip_nR_Eff,ip_R,ipMem,nMem,                       &
     &                 L_Quad,nR,nAtoms,nMaxExp,                        &
     &                 nTotGP,nbrmxbas,iAngMax,ip_ioffsh,               &
     &                 iOpt_Angular,                                    &
     &                 mIrrep, nISh, nAsh, mBas,Functional_type,mOrb,   &
     &                 Grid_Type, Rotational_Invariance,                &
     &                            mTmp,mRad,nAOMax,                     &
     &                            NQ_Direct,Packing,OffPUVX,            &
     &                 iQEnd,ip_OrbDip,ioff_ash,ioff_bas,ioff_basact,   &
     &                 OffBas,OffOrb,OffOrb2,OffOrbTri,OffBas2,         &
     &                 OffBasFro,OffUVX,nUVX,nUVXt,OffVX,nVX,nVXt
      Common /Quad_ii / nFro
!
      Real*8 R_Max(0:LMax_NQ)
      Real*8 rQStrt,Threshold_save,Crowding,Threshold,Energy_integrated, &
     &       Dens_I,Grad_I,Tau_I,Dens_a1,Dens_b1,Dens_a2,Dens_b2,        &
     &       Dens_t1,Dens_t2,Block_Size,x_min,x_max,y_min,y_max,z_min,   &
     &       z_max,Fade,ThrC,T_Y,rQEnd
      Common /Quad_r /rQStrt,                                           &
     &                Threshold_save, Crowding,                         &
     &                Threshold,R_Max,Energy_integrated,                &
     &                Dens_I,Grad_I,Tau_I,                              &
     &                Dens_a1,Dens_b1,Dens_a2,Dens_b2,Dens_t1,Dens_t2,  &
     &                Block_Size,x_min,x_max,y_min,y_max,z_min,z_max,   &
     &                Fade,                ThrC, T_Y,                   &
     &                rQEnd
!
      Integer cQStrt, cQEnd
      Character Quadrature*10, Pad*6, MBC*8
!
      Common /Quad_c /cQStrt,                                           &
     &                Quadrature, Pad, MBC,                             &
     &                cQEnd
!
End Module nq_Info

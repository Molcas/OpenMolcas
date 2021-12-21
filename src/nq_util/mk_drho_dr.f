***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2000,2002,2021, Roland Lindh                           *
************************************************************************
      Subroutine Mk_dRho_dR(list_s,nlist_s,list_g,list_bas)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
*             Modified December 2021 in Uppsala by RL.                 *
************************************************************************
      use iSD_data
      use nq_Grid, only: Grid_AO, Ind_Grd, TabAO, dRho_dR
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "nq_info.fh"
      Integer list_s(2,nlist_s), list_g(3,nlist_s), list_bas(2,nlist_s)
      Integer, Parameter :: Index_d2(3,3)=
     &    Reshape([5,6,7, 6,8,9, 7,9,10],[3,3])
      Integer, Parameter :: Index_d3(3,3) =
     &    Reshape([11,14,16, 12,17,19, 13,18,19],[3,3])
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      nAO = SIZE(Grid_AO,3)
      nD  = SIZE(Grid_AO,4)
      mGrid = SIZE(TabAO,2)
*                                                                      *
************************************************************************
*                                                                      *
      iOff = 0
      Do ilist_s=1,nlist_s
         iS       =list_s(1,ilist_s)
         iBas_Eff = list_bas(1,ilist_s)
         iCmp     =iSD( 2,iS)

         n_iBas=iBas_Eff*iCmp
*
         Ind_Grd(1,iOff+1:iOff+n_iBas) = List_g(1,ilist_s)
         Ind_Grd(2,iOff+1:iOff+n_iBas) = List_g(2,ilist_s)
         Ind_Grd(3,iOff+1:iOff+n_iBas) = List_g(3,ilist_s)
*
         iOff = iOff + iBas_Eff*iCmp
      End Do                         ! ilist_s
*                                                                      *
************************************************************************
*                                                                      *
      If (Functional_Type.eq.LDA_Type) Then
#ifdef _REDUNDANT_
         dRho_dR(:,:,:)=Zero
         Do iD = 1, nD
            Do iAO = 1, nAO
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)
                  j = iCar + 1

                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of the density.
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                                = dRho_dR(iD,iGrid,Ind_xyz)
     &                                + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
                     End Do
                  End If

               End Do
            End Do
         End Do
#endif
      Else If (Functional_Type.eq.GGA_Type) Then
         dRho_dR(:,:,:)=Zero
         Do iD = 1, nD                      ! index of rho
            Do iAO = 1, nAO
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)
                     End Do
                  End If

               End Do
            End Do
         End Do
       Else If (Functional_Type.eq.meta_GGA_Type1) Then
         dRho_dR(:,:,:)=Zero
         Do iD = 1, nD                      ! index of rho
            Do iAO = 1, nAO
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  iT  = nD*4 + iD      ! index of tau component

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)

                        dRho_dR(iT,iGrid,Ind_xyz)
     &                             = dRho_dR(iT,iGrid,Ind_xyz)
     &                       + Four* TabAO(idjx,iGrid,iAO)
     &                              * Grid_AO(2,iGrid,iAO,iD)
     &                       + Four* TabAO(idjy,iGrid,iAO)
     &                              * Grid_AO(3,iGrid,iAO,iD)
     &                       + Four* TabAO(idjz,iGrid,iAO)
     &                              * Grid_AO(4,iGrid,iAO,iD)
                     End Do
                  End If

               End Do
            End Do
         End Do
       Else If (Functional_Type.eq.meta_GGA_Type2) Then
         dRho_dR(:,:,:)=Zero
         Do iD = 1, nD                      ! index of rho
            Do iAO = 1, nAO
*
*------------- Loop over cartesian components
*
               Do iCar = 1, 3

                  Ind_xyz=Ind_Grd(iCar,iAO)! index of  nuclear gradient

                  j = iCar + 1             ! index derivative of AO

                  iDx = nD + (iD-1)*3 + 1  ! index of grad rho component
                  iDy = iDx + 1
                  iDz = iDy + 1

                  iT  = nD*4 + iD      ! index of tau component

                  iL  = nD*5 + iD      ! index if laplacian component

                  idjx = Index_d2(1,iCar)
                  idjy = Index_d2(2,iCar)
                  idjz = Index_d2(3,iCar)

                  idjx2 = Index_d3(1,iCar)
                  idjy2 = Index_d3(2,iCar)
                  idjz2 = Index_d3(3,iCar)
                  idx2  = Index_d2(1,1)
                  idy2  = Index_d2(2,2)
                  idz2  = Index_d2(3,3)
                  If (Ind_xyz/=0) Then
                     Do iGrid = 1, mGrid
*
*                       Cartesian derivative of rho
*
                        dRho_dR(iD,iGrid,Ind_xyz)
     &                             = dRho_dR(iD,iGrid,Ind_xyz)
     &                        + Two * Grid_AO(1,iGrid,iAO,iD)
     &                                * TabAO(j,iGrid,iAO)
*
*                       Cartesian derivatives of grad rho
*
                        dRho_dR(iDx,iGrid,Ind_xyz)
     &                              = dRho_dR(iDx,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjx,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(2,iGrid,iAO,iD)
                        dRho_dR(iDy,iGrid,Ind_xyz)
     &                              = dRho_dR(iDy,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjy,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(3,iGrid,iAO,iD)
                        dRho_dR(iDz,iGrid,Ind_xyz)
     &                              = dRho_dR(iDz,iGrid,Ind_xyz)
     &                         + Two * TabAO(idjz,iGrid,iAO)
     &                                * Grid_AO(1,iGrid,iAO,iD)
     &                            + Two * TabAO(j,iGrid,iAO)
     &                                * Grid_AO(4,iGrid,iAO,iD)

                        dRho_dR(iT,iGrid,Ind_xyz)
     &                             = dRho_dR(iT,iGrid,Ind_xyz)
     &                       + Four* TabAO(idjx,iGrid,iAO)
     &                              * Grid_AO(2,iGrid,iAO,iD)
     &                       + Four* TabAO(idjy,iGrid,iAO)
     &                              * Grid_AO(3,iGrid,iAO,iD)
     &                       + Four* TabAO(idjz,iGrid,iAO)
     &                              * Grid_AO(4,iGrid,iAO,iD)

                        dRho_dR(iL,iGrid,Ind_xyz)
     &                             = dRho_dR(iL,iGrid,Ind_xyz)
     &                             + Two * Grid_AO(1,iGrid,iAO,iD)
     &                             *(TabAO(idjx2,iGrid,iAO)
     &                             + TabAO(idjy2,iGrid,iAO)
     &                             + TabAO(idjz2,iGrid,iAO))
     &                             + Two *(Grid_AO(idx2,iGrid,iAO,iD)
     &                                    +Grid_AO(idy2,iGrid,iAO,iD)
     &                                    +Grid_AO(idz2,iGrid,iAO,iD))
     &                             *       TabAO(j,iGrid,iAO)

                     End Do
                  End If

               End Do
            End Do
         End Do
      Else
         Call abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      nGrad_Eff = SIZE(dRho_dR,3)
      Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       SIZE(dRho_dR,1)*mGrid,nGrad_Eff)
#endif
      Return
      End

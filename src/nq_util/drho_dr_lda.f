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
      Subroutine dRho_dR_LDA(nD,dRho_dR,ndRho_dR,mGrid,
     &                       list_s,nlist_s,nGrad_Eff,list_g,list_bas)
************************************************************************
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN.  2000                                   *
*                                                                      *
*             Modified May-June 2002 in Tokyo for DFT gradient by RL.  *
************************************************************************
      use iSD_data
      use nq_Grid, only: Grid_AO, Ind_Grd, TabAO
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "disp.fh"
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer On, Off
      Parameter (On=1, Off=0)
      Integer list_s(2,nlist_s), list_g(3,nlist_s), list_bas(2,nlist_s)
      Real*8 dRho_dR(ndRho_dR,mGrid,nGrad_Eff)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      dRho_dR(:,:,:)=Zero
      nAO = SIZE(Grid_AO,3)
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
      Do iD = 1, nD
         Do iAO = 1, nAO
*
*---------- Loop over cartesian components
*
            Do iCar = 1, 3

               Ind_xyz=Ind_Grd(iCar,iAO)
               j = iCar + 1

               If (Ind_xyz/=0) Then
                  Do iGrid = 1, mGrid
                     dRho_dR(iD,iGrid,Ind_xyz)=dRho_dR(iD,iGrid,Ind_xyz)
     &                             + Two * Grid_AO(1,iGrid,iAO,iD)
     &                             * TabAO(j,iGrid,iAO)
                  End Do
               End If

            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('dRho_dR_LDA: dRho_dR',' ',dRho_dR,
     &                       ndRho_dR*mGrid,nGrad_Eff)
*
#endif
      Return
      End

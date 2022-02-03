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
* Copyright (C) 2000,2022, Roland Lindh                                *
*               Ajitha Devarajan                                       *
************************************************************************
      Subroutine DFT_Int(list_s,nlist_s,FckInt,nFckInt,nD,Fact,ndc)
************************************************************************
*                                                                      *
* Object: to compute contributions to                                  *
*                                                                      *
*         <m|dF/drho|n> ; integrals over the potential                 *
*                                                                      *
*         where                                                        *
*                                                                      *
*         F(r)=rho(r)*e(rho(r),grad[rho(r)])                           *
*                                                                      *
*      Author:Roland Lindh, Department of Chemical Physics, University *
*             of Lund, SWEDEN. November 2000                           *
*             D.Ajitha:Modifying for the new Kernel outputs            *
************************************************************************
      use iSD_data
      use Symmetry_Info, only: nIrrep
      use nq_Grid, only: TabAO, Grid_AO,
     &                   AOIntegrals => Dens_AO
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "stdalloc.fh"
      Real*8 Fact(ndc**2), FckInt(nFckInt,nD)
      Integer list_s(2,nlist_s)
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the desired AO integrand here from the AOs, accumulate
*     contributions to the SO integrals on the fly.
*
      mGrid=SIZE(TabAO,2)
      mAO = SIZE(TabAO,1)
*
      nBfn = Size(AOIntegrals,1)
      nFn  = Size(Grid_AO,1)

      Call Do_NInt_d()
      Call Do_NIntX(AOIntegrals,mGrid,TabAO,nBfn,nD,mAO)
*                                                                      *
************************************************************************
*                                                                      *
*     Distribute result on to the full integral matrix.
*
      Do iD = 1, nD
         If (nIrrep.eq.1) Then
            Call AOAdd_Full(AOIntegrals(:,:,iD),nBfn,FckInt(:,iD),
     &                      nFckInt)
         Else
            Call SymAdp_Full(AOIntegrals(:,:,iD),nBfn,FckInt(:,iD),
     &                       nFckInt,list_s,nlist_s,Fact,ndc)
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      End

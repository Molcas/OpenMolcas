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
* Copyright (C) 1992, Roland Lindh                                     *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine AlloK2_Funi(nr_of_Densities)
************************************************************************
*                                                                      *
*  Object: Allocate space for K2 entities.                             *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden. November '92                 *
*             Martin Schuetz, Dept. of Theoretical Chemistry,          *
*             University of Lund, Sweden. Jun '95                      *
************************************************************************
      use iSD_data
      use k2_arrays
      use IOBUF
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "real.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
*
*
      Call Nr_Shells(nSkal)
*
*     determine memory size nDeDe, MaxDe, and MaxDRC
      nDeDe_DFT = 0
      MaxDe     = 0
*
************************************************************************
*                                                                      *
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      Do iS = 1, nSkal
C        iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
C        iPrim  = iSD( 5,iS)
         iAO    = iSD( 7,iS)
C        mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
*
         Do jS = 1, iS
C           jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
C           jPrim  = iSD( 5,jS)
            jAO    = iSD( 7,jS)
C           mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
*
C           iDeSiz = 1 + iPrim*jPrim + (iBas*jBas+1)*iCmp*jCmp
            iDeSiz = iBas*jBas*iCmp*jCmp
            MaxDe = Max(MaxDe,iDeSiz)
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            If (nSO.gt.0) Then
               nDeDe_DFT = nDeDe_DFT
     &                   + nr_of_Densities*iDeSiz*nIrrep
            End If
*
         End Do
      End Do
*
      Return
      End

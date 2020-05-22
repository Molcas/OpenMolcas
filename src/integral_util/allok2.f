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
      SubRoutine AlloK2()
************************************************************************
*                                                                      *
*  Object: Allocate space for K2 entities.                             *
*                                                                      *
* Called from: ReadIn (client) / DPSCF main (server)                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              DCR                                                     *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden. November '92                 *
*             Martin Schuetz, Dept. of Theoretical Chemistry,          *
*             University of Lund, Sweden. Jun '95                      *
************************************************************************
      use k2_setup
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "ndarray.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "IOBuf.fh"
#include "nsd.fh"
#include "setup.fh"
#include "status.fh"
#include "k2.fh"
*
*     declaration of local vars...
      Logical Debug
      Data Debug/.False./
*
*---- Statement function
*
      nElem(i)=(i+1)*(i+2)/2
      nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6  - 1
*
*     Call QEnter('AlloK2')
*
      If (Debug) Then
         If (Allocated(Data_k2)) Then
            Write (6,*) 'Enter Allok2, k2_Status=Active'
         Else If (k2_Status.eq.Produced) Then
            Write (6,*) 'Enter Allok2, k2_Status=Produced'
         Else
            Write (6,*) 'Enter Allok2, k2_Status=InActive'
         End If
      End If
      If (Allocated(Data_k2) .or. k2_Status.eq.Produced   ) Return
*
      Call Nr_Shells(nSkal)
*
*     determine memory size for K2 entities
*     for this, run dummy K2 loop...
      nk2 = 0
      nDeDe = 0
      MaxDe = 0
      nr_of_Densities=1 ! Hardwired option
*                                                                      *
************************************************************************
*                                                                      *
*-----Double loop over shells. These loops decide the integral type
*
      Do iS = 1, nSkal
         iShll  = iSD( 0,iS)
         If (AuxShell(iShll) .and. iS.ne.nSkal) Go To 100
         iAng   = iSD( 1,iS)
         iCmp   = iSD( 2,iS)
         iBas   = iSD( 3,iS)
         iPrim  = iSD( 5,iS)
         mdci   = iSD(10,iS)
         iShell = iSD(11,iS)
*
         Do jS = 1, iS
            jShll  = iSD( 0,jS)
            If (AuxShell(jShll) .and. jS.eq.nSkal) Go To 200
            jAng   = iSD( 1,jS)
            jCmp   = iSD( 2,jS)
            jBas   = iSD( 3,jS)
            jPrim  = iSD( 5,jS)
            mdcj   = iSD(10,jS)
            jShell = iSD(11,jS)
*
            If (nIrrep.eq.1) Then
               iDeSiz = 1 + iPrim*jPrim +               iCmp*jCmp
            Else
               iDeSiz = 1 + iPrim*jPrim + (iBas*jBas+1)*iCmp*jCmp
            End If
            MaxDe = Max(MaxDe,iDeSiz)
            iSmLbl = 1
            nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell)
            If (nSO.gt.0) Then
               nDeDe = nDeDe + nr_of_Densities*iDeSiz*nIrrep
            End If
*
            nZeta=iPrim*jPrim
            ijCmp=nElem(iAng)*nElem(jAng)
            If (.Not.DoGrad_) ijCmp=0
            nHm=iCmp*jCmp*(nabSz(iAng+jAng)-nabSz(Max(iAng,jAng)-1))
            nHm=nHm*nIrrep
            nData=nZeta*(nDArray+2*ijCmp)+nDScalar+nHm
            nk2 = nk2 + nData*nIrrep
 200        Continue
         End Do
 100     Continue
      End Do
*     now ... allocate memory
      Call mma_allocate(Data_k2,nk2)
      Call FZero(Data_k2,nk2)
      nIndk2=nShlls*(nShlls+1)/2
      call mma_allocate(Indk2,2,nIndk2)
*
      Return
      End

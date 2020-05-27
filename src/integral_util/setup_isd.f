************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine SetUp_iSD()
      use iSD_data
      Implicit Real*8 (A-H,O-Z)
#include "nsd.fh"
#include "setup.fh"
#include "stdalloc.fh"
*
*                                                                      *
************************************************************************
*                                                                      *
      If (Allocated(iSD)) Call mma_deallocate(iSD)
      Call Nr_Shells(nSkal)
      mSkal=nSkal
      nSkal_iSD=nSkal+4  ! Add four slots for future use.
      call mma_allocate(iSD,[0,nSD],[1,nSkal_iSD],label='iSD')
      Call Def_Shells(iSD,nSD,nSkal)
*                                                                      *
************************************************************************
*                                                                      *
*.... Compute the size of and allocate auxiliary memory
*
      Call Get_iScalar('nSym',nIrrep)
      MxPrm = 0
      MxFT = 0
      MxDij = 0
      Do iS = 1, nSkal
         iCmp =iSD(2,iS)
         iBas =iSD(3,iS)
         iPrim=iSD(5,iS)
         MxPrm=Max(MxPrm,iPrim)
         If (nIrrep.eq.1) Then
            MxFT=1 ! Dummay assignment
            MxDij= Max(MxDij,iCmp**2+iPrim**2+1)
         Else
            MxFT = Max(MxFT,6*(iBas*iCmp)**2)
            MxDij= Max(MxDij,(iBas**2+1)*iCmp**2+iPrim**2+1)
         End If
      End Do
      MxDij = 6 * nIrrep * MxDij
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

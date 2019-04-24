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
      Subroutine Print_EigenValues(H,nH)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 H(nH*(nH+1)/2)
*
      Call QEnter('PEV')
      Lu=6
*
      Call GetMem('EVal','Allo','Real',ipEVal,nH*(nH+1)/2)
      Call GetMem('EVec','Allo','Real',ipEVec,nH*nH)
*
*---- Copy elements for H
*
      call dcopy_(nH*(nH+1)/2,H,1,Work(ipEVal),1)
*
*---- Set up a unit matrix
*
      call dcopy_(nH*nH,[Zero],0,Work(ipEVec),1)
      call dcopy_(nH,[One],0,Work(ipEVec),nH+1)
*
*---- Compute eigenvalues and eigenvectors
*
      Call Jacob (Work(ipEVal),Work(ipEVec),nH,nH)
      Call Jacord(Work(ipEVal),Work(ipEVec),nH,nH)
*
*---- Print out the result
*
      Write (Lu,*)
      Write (Lu,*) 'Eigenvalues of the matrix'
      Write (Lu,*)
      Write (Lu,'(10F15.8)') (Work(i*(i+1)/2+ipEVal-1),i=1,nH)
*
      Call GetMem('EVec','Free','Real',ipEVec,nH*nH)
      Call GetMem('EVal','Free','Real',ipEVal,nH*(nH+1)/2)
*
      Call QExit('PEV')
      Return
      End

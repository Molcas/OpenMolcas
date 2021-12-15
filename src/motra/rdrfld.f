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
      Subroutine RdRfld(ipHOne)
*
      Implicit Real*8 (A-H,O-Z)
#include "files_motra.fh"

#include "motra_global.fh"
#include "trafo_motra.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      Logical Found
*
*----------------------------------------------------------------------*
*     If this is a perturbative reaction field calculation then        *
*     modifiy the one-electron Hamiltonian by the reaction field and   *
*     the nuclear attraction by the cavity self-energy                 *
*----------------------------------------------------------------------*
      nTemp=0
      Do iSym=1,nSym
         nTemp=nTemp+nBas(iSym)*(nBas(iSym)+1)/2
      End Do
      Call GetMem('RFFLD','Allo','Real',lTemp,nTemp)
      Call f_Inquire('RUNOLD',Found)
      If (Found) Call NameRun('RUNOLD')
      Call get_dscalar('RF Self Energy',ERFself)
      PotNuc=PotNuc+ERFself
      Call get_darray('Reaction field',Work(lTemp),nTemp)
      If (Found) Call NameRun('RUNFILE')
      Call Daxpy_(nTemp,1.0D0,Work(lTemp),1,Work(ipHone),1)
      Call GetMem('RFFLD','Free','Real',lTemp,nTemp)
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      Return
      End

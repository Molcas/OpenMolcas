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
      Subroutine MMEtoRMO(nAObas,nMObas,ipAvRed,iMME)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6)

*
*--- First all multipoles are transformed to MO-basis...
*
      Call GetMem('Squared','Allo','Real',iSq,nAObas**2)
      Call GetMem('TEMP','Allo','Real',iTEMP,nAObas*nMObas)
      Call GetMem('Final','Allo','Real',iMmeMO,nMObas**2)
      nUniqueM=1+3+6
      Do 11, iMlt=1,nUniqueM
        Call Square(Work(iMME(iMlt)),Work(iSq),iONE,nAObas,nAObas)
        Call Dgemm_('T','N',nMObas,nAObas,nAObas,ONE,Work(ipAvRed)
     &            ,nAObas,Work(iSq),nAObas,ZERO,Work(iTEMP),nMObas)
        Call Dgemm_('N','N',nMObas,nMObas,nAObas,ONE,Work(iTEMP)
     &            ,nMObas,Work(ipAvRed),nAObas,ZERO,Work(iMmeMO),nMObas)
        Call SqToTri_Q(Work(iMmeMO),Work(iMME(iMlt)),nMObas)
11    Continue
      Call GetMem('Squared','Free','Real',iSq,nAObas**2)
      Call GetMem('TEMP','Free','Real',iTEMP,nAObas*nMObas)
      Call GetMem('Final','Free','Real',iMmeMO,nMObas**2)

      Return
      End

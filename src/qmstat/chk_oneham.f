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
      Subroutine Chk_OneHam(nBas)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension nBas(MxSym)
      Character Label_Read*8, Label_Pure*8

      Lu_One=49
      Lu_One=IsFreeUnit(Lu_One)
      Label_Read='OneHam  '
      Label_Pure='OneHam 0'
      nBT=nBas(1)*(nBas(1)+1)/2
      Call OpnOne(irc,0,'ONEINT',Lu_One)
      Call GetMem('Read','Allo','Real',iOneR,nBT+4)
      Call GetMem('Pure','Allo','Real',iOneP,nBT+4)

      irc=-1
      iopt=0
      iSmLbl=0
      Call RdOne(irc,iopt,Label_Read,1,Work(iOneR),iSmLbl)
      irc=-1
      iopt=0
      iSmLbl=0
      Call RdOne(irc,iopt,Label_Pure,1,Work(iOneP),iSmLbl)
      Call ClsOne(irc,Lu_One)

      Call DaxPy_(nBT,-1.0d0,Work(iOneR),1,Work(iOneP),1)

      dNorm=dnrm2_(nBT,Work(iOneP),1)

      If(dNorm.gt.1d-8) then
        Write(6,*)
        Write(6,*)
        Write(6,*)' WARNING!'
        Write(6,*)
        Write(6,*)'   Your one-electron hamiltonian is not purely'
     &//' vacuum. This means that the Hamiltonian'
        Write(6,*)'   in QmStat can be contaminated. Is this'
     &//' intentional? If not, then make sure that the ONEINT'
        Write(6,*)'   file comes directly from a Seward calculation'
     &//' without any calls from'
        Write(6,*)'   FFPT (or similar) in between.'
        Write(6,*)
        Write(6,*)
      Endif

      Call GetMem('Read','Free','Real',iOneR,nBT+4)
      Call GetMem('Pure','Free','Real',iOneP,nBT+4)

      Return
      End

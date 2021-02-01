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
      Subroutine Desymmetrize(SOInt,nSOInt,Scr,nScr,AOInt,nBas,nBas1,
     &                        SymInv,nSym,iSyLbl)
      Implicit Real*8 (a-h,o-z)
*
#include "real.fh"
      Real*8 SOInt(nSOInt), Scr(nScr), AOInt(nBas1,nBas1),
     &       SymInv(nBas1**2)
      Integer nBas(0:nSym-1)
*
*define _DEBUGPRINT_
*
      Call FZero(AOInt,nBas1**2)
*
      iOffSO=1
      iOffPi=1
      Do iSym = 0, nSym-1
         iOffPj=1
         Do jSym = 0, iSym
            ijSym=iEor(iSym,jSym)
            If (iAnd(iSyLbl,2**ijSym).eq.0) Go To 20
            If (nBas(iSym)*nBas(jSym).eq.0) Go To 30
            If (iSym.eq.jSym) Then
#ifdef _DEBUGPRINT_
               Write (6,*) 'Diagonal Block'
               Call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),
     &                                               nBas(iSym))
#endif
*
               Call DGEMM_('N','T',
     &                     nBas(iSym),nBas1,nBas(iSym),
     &                     1.0d0,SOInt(iOffSO),nBas(iSym),
     &                     SymInv(iOffPi),nBas1,
     &                     0.0d0,Scr,nBas(iSym))
*
               Call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),
     &                    One,SymInv(iOffPi),nBas1,
     &                        Scr,nBas(iSym),
     &                    One,AOInt,nBas1)
*
            Else
#ifdef _DEBUGPRINT_
               Write (6,*) 'Off-Diagonal Block',iSym,jSym
               Call RecPrt('SOInt',' ',SOInt(iOffSO),nBas(iSym),
     &                                               nBas(jSym))
               Call RecPrt('Pi',' ',SymInv(iOffPi),nbas1,nBas(iSym))
               Call RecPrt('Pj',' ',SymInv(iOffPj),nbas1,nBas(jSym))
#endif
*
               Call DGEMM_('N','T',
     &                     nBas(iSym),nBas1,nBas(jSym),
     &                     1.0d0,SOInt(iOffSO),nBas(iSym),
     &                     SymInv(iOffPj),nBas1,
     &                     0.0d0,Scr,nBas(iSym))
*
               Call DGEMM_('N','N',nBas1,nBas1,nBas(iSym),
     &                    One,SymInv(iOffPi),nBas1,
     &                        Scr,nBas(iSym),
     &                    One,AOInt,nBas1)
*
               Call DGEMM_('T','T',nBas1,nBas1,nBas(iSym),
     &                    One,Scr,nBas(iSym),
     &                        SymInv(iOffPi),nBas1,
     &                    One,AOInt,nBas1)
*
            End If
*
 30         Continue
            iOffSO= iOffSO + nBas(iSym)*nBas(jSym)
 20         Continue
            iOffPj = iOffPj + nBas(jSym)*nBas1
         End Do
         iOffPi = iOffPi + nBas(iSym)*nBas1
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('AOInt',' ',AOInt,nBas1,nBas1)
#endif
*
      Return
      End

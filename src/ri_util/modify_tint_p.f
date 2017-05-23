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
      Subroutine Modify_TInt_p(TInt,nTheta_All,List2,mData)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Real*8 TInt(nTheta_All,nTheta_All)
      Integer List2(mData,nTheta_All)
*
#ifdef _DEBUG_
      Call RecPrt('Modify_TInt_p: TInt',' ',TInt,nTheta_All,nTheta_All)
#endif
      Do iTheta_All= 1, nTheta_All
         iPrim=List2(5,iTheta_All)
         iShll=List2(7,iTheta_All)
         nConti=nBasis_Cntrct(iShll)
         nPrimi=nExp(iShll)
         iOff = ipCff_Cntrct(iShll) + iPrim -1
         Coeff_i = DDot_(nConti,Work(iOff),nPrimi,Work(iOff),nPrimi)
         Coeff_i = Sqrt(Coeff_i)
*
         jPrim=List2(6,iTheta_All)
         jShll=List2(8,iTheta_All)
         nContj=nBasis_Cntrct(jShll)
         nPrimj=nExp(jShll)
         jOff = ipCff_Cntrct(jShll) + jPrim -1
         Coeff_j = DDot_(nContj,Work(jOff),nPrimj,Work(jOff),nPrimj)
         Coeff_j = Sqrt(Coeff_j)
*
         Do jTheta_All = 1, nTheta_All
            kPrim=List2(5,jTheta_All)
            kShll=List2(7,jTheta_All)
            nContk=nBasis_Cntrct(kShll)
            nPrimk=nExp(kShll)
            kOff = ipCff_Cntrct(kShll) + kPrim -1
            Coeff_k = DDot_(nContk,Work(kOff),nPrimk,Work(kOff),nPrimk)
            Coeff_k = Sqrt(Coeff_k)
*
            lPrim=List2(6,jTheta_All)
            lShll=List2(8,jTheta_All)
            nContl=nBasis_Cntrct(lShll)
            nPriml=nExp(lShll)
            lOff = ipCff_Cntrct(lShll) + lPrim -1
            Coeff_l = DDot_(nContl,Work(lOff),nPriml,Work(lOff),nPriml)
            Coeff_l = Sqrt(Coeff_l)
*
            TInt(iTheta_All,jTheta_All) = TInt(iTheta_All,jTheta_All)
     &                                  * Coeff_i * Coeff_j
     &                                  * Coeff_k * Coeff_l
#ifdef _DEBUG_
            Write (6,*)
            Write (6,*) Coeff_i,Coeff_j,Coeff_k,Coeff_l
            Write (6,*) iPrim,  jPrim,  kPrim,  lPrim
            Write (6,*) nPrimi, nPrimj, nPrimk, nPriml
#endif
*
         End Do
*
      End Do
#ifdef _DEBUG_
      Call RecPrt('Modify_TInt_p: TInt',' ',TInt,nTheta_All,nTheta_All)
#endif
*
      Return
      End

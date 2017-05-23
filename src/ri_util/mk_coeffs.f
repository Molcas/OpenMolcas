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
      Subroutine Mk_Coeffs(CoeffA,nPrimA,nConA,CoeffB,nPrimB,nConB,
     &                     Coeff,nTheta_Full,nPhi,iD,NumCho,
     &                     List2,mData,nPhi_All,
     &                     Indkl,nkl,nk,nl,iAng,jAng,
     &                     CoeffAP,CoeffBP)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 CoeffA(nPrimA,nConA), CoeffB(nPrimB,nConB),
     &       Coeff(nTheta_Full,nPhi),
     &       CoeffAP(nPrimA,nPrimA), CoeffBP(nPrimB,nPrimB)
      Integer List2(mData,nPhi_All), iD(NumCho), Indkl(nkl)
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
*
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('CoeffA',' ',CoeffA,nPrimA,nConA)
      Call RecPrt('CoeffB',' ',CoeffB,nPrimB,nConB)
      Call iVcPrt('Indkl',' ',Indkl,nkl)
      Call iVcPrt('Mk_Coeffs: iD',' ',iD,NumCho)
      Write (6,*) 'iAng,jAng=',iAng,jAng
#endif
      Do iCho = 1, NumCho
         iPhi_All = iD(iCho)
         If (List2(1,iPhi_All).eq.iAng .and.
     &       List2(2,iPhi_All).eq.jAng) Then
            ik=List2(5,iPhi_All)
            il=List2(6,iPhi_All)
*
            If (iAng.eq.jAng) Then
               iPhi_Full = iTri(ik,il)
               iPhi=Indkl(iPhi_Full)
               If (iPhi.eq.0) Go To 101
               Do iPrimA = 1, nPrimA
                  Do iPrimB = 1, iPrimA
                     Cff =(CoeffA(iPrimA,ik)*CoeffB(iPrimB,il)
     &                   + CoeffA(iPrimA,il)*CoeffB(iPrimB,ik))/
     &                    (CoeffAP(iPrimA,iPrimA)*
     &                     CoeffBP(iPrimB,iPrimB))
                     If (iPrimA.eq.iPrimB) Cff = Half * Cff
                     iTheta_Full=iPrimA*(iPrimA-1)/2 + iPrimB
                     Coeff(iTheta_Full,iPhi) = Cff
                  End Do
               End Do
            Else
               iPhi_Full = (il-1)*nk + ik
               iPhi=Indkl(iPhi_Full)
               If (iPhi.eq.0) Go To 101
               Do iPrimA = 1, nPrimA
                  Do iPrimB = 1, nPrimB
                     Cff = CoeffA(iPrimA,ik)*CoeffB(iPrimB,il)/
     &                    (CoeffAP(iPrimA,iPrimA)*
     &                     CoeffBP(iPrimB,iPrimB))
                     iTheta_Full=(iPrimB-1)*nPrimA + iPrimA
                     Coeff(iTheta_Full,iPhi) = Cff
                  End Do
               End Do
            End If
 101        Continue
*
         End If
      End Do
#ifdef _DEBUG_
      Call RecPrt('Coeff',' ',Coeff,nTheta_Full,nPhi)
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nl)
      End

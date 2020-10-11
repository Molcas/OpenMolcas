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
      Subroutine GF_Mult(G,F,GF,mInter)
      Implicit Real*8 (a-h,o-z)
      Real*8 G(mInter**2),F(mInter**2),GF(mInter*(mInter+1)/2)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Form the GF-matrix (actually G^(1/2)FG^(1/2))
*
      Do iX = 1, mInter
         ii = (iX-1)*mInter + iX
         XMass_i = Sqrt(G(ii))
         Do jX = 1, iX
            jj = (jX-1)*mInter + jX
            XMass_j = Sqrt(G(jj))
            ij = (jX-1)*mInter + iX
            ji = (iX-1)*mInter + jX
            ijT= iX*(iX-1)/2+jX
            GF(ijT) = XMass_i*XMass_j*F(ij)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call TriPrt('G^(1/2)FG^(1/2)',' ',GF,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

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
      Subroutine mk_G(G,GInv,nX)
      use Slapaf_Parameters, only: nDimBC
#include "info_slapaf.fh"
      Integer nX
      Real*8 G(nX,nX), GInv(nX*nX)
*
      Call mk_G_Internal(G,GInv,nDimBC)
*
      Return
      End
      Subroutine mk_G_Internal(G,GInv,nDoF)
      use Slapaf_Info, only: dMass, Degen, Smmtrc
      use Slapaf_Parameters, only: Curvilinear, User_Def
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "constants2.fh"
*
      Real*8 G(nDoF,nDoF), GInv(nDoF,nDoF)
      Logical Auto

      Auto=.Not.User_Def
      nsAtom=SIZE(Smmtrc,2)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the mass tensor
*
      G(:,:)=Zero
      GInv(:,:)=Zero
      ii = 0
      Do i = 1, nsAtom
         Do ix = 1, 3
            If (Smmtrc(ix,i)) Then
               ii = ii + 1
               If (Auto.and..Not.Curvilinear) Then
                  G(ii,ii) = Degen(ix,i)/dMass(i)
               Else
                  G(ii,ii) = One/(Degen(ix,i)*dMass(i))
               End If
               GInv(ii,ii) = One/(G(ii,ii)*UTOAU)
            End If
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('G (cartesian)',' ',G,nDoF,nDoF)
      Call RecPrt('G-1 (cartesian)',' ',GInv,nDoF,nDoF)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End

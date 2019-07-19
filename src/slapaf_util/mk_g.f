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
      Subroutine mk_G(ipG,ipGInv)
#include "info_slapaf.fh"
#include "WrkSpc.fh"
*
      nX=3*nsAtom
      mInter=mInt + mTROld
*
      Call Allocate_Work(ipG,nX**2)
      Call Allocate_Work(ipGInv,nX**2)
*
      Call mk_G_(Work(ipG),Work(ipGInv),nX,mInter,nsAtom,.Not.User_Def,
     &           CurviLinear,Smmtrc,Degen,Work(ipCM))
*
      Return
      End
      Subroutine mk_G_(G,GInv,nX,mInter,nAtom,Auto,nrc,Smmtrc,Degen,
     &                 dMass)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "constants2.fh"
*
      Real*8 G(mInter,mInter), GInv(mInter**2), dMass(nAtom),
     &       Degen(3,nAtom)
      Logical Auto, nrc, Smmtrc(3,nAtom)
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the mass tensor
*
      Call FZero(G,mInter**2)
      Call FZero(GInv,mInter**2)
      ii = 0
      Do i = 1, nAtom
         Do ix = 1, 3
            If (Smmtrc(ix,i)) Then
               ii = ii + 1
               If (Auto.and..Not.nrc) Then
                  G(ii,ii) = Degen(ix,i)/dMass(i)
               Else
                  G(ii,ii) = One/(Degen(ix,i)*dMass(i))
               End If
               jj = (ii-1)*mInter + ii
               GInv(jj) = One/(G(ii,ii)*UTOAU)
            End If
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('G (cartesian)',' ',G,mInter,mInter)
      Call RecPrt('G-1 (cartesian)',' ',GInv,mInter,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nX)
      End

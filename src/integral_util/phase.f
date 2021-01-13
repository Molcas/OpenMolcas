************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      Subroutine Phase(iCmp, jCmp, kCmp, lCmp, iAng,
     &                 iShll, kOp, ijkl, AOInt)
************************************************************************
*                                                                      *
*  Object: To change the phase of the integrals in accordance with the *
*          swapping of the operators operating on the integrals.       *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             June '90                                                 *
************************************************************************
      use Basis_Info
      use Real_Spherical, only: iSphCr
      use Symmetry_Info, only: iChBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iAng(4), iShll(4)
*
*     Statement Function
*
      iOff(ixyz)  = ixyz*(ixyz+1)*(ixyz+2)/6
*
*     Call RecPrt(' In Phase: AOInt ',' ',AOInt,ijkl,ijCmp*ijCmp)
*
*     Change phase factor. This is only necessary if T=/=E.
*
      If (kOp.eq.0 .or. iCmp*jCmp*kCmp*lCmp.eq.0) Go To 14
      ii = iOff(iAng(1))
      jj = iOff(iAng(2))
      kk = iOff(iAng(3))
      ll = iOff(iAng(4))
      Do 10 i1 = 1, iCmp
       iChBs = iChBas(ii+i1)
       If (Shells(iShll(1))%Transf) iChBs = iChBas(iSphCr(ii+i1))
       pa1T = DBLE(iPrmt(kOp,iChBs))
       Do 11 i2 = 1, jCmp
        jChBs = iChBas(jj+i2)
        If (Shells(iShll(2))%Transf) jChBs = iChBas(iSphCr(jj+i2))
        pb1T = DBLE(iPrmt(kOp,jChBs))
*
        Do 12 i3 = 1, kCmp
         kChBs = iChBas(kk+i3)
         If (Shells(iShll(3))%Transf) kChBs = iChBas(iSphCr(kk+i3))
         pa2T = DBLE(iPrmt(kOp,kChBs))
         Do 13 i4 = 1, lCmp
          lChBs = iChBas(ll+i4)
          If (Shells(iShll(4))%Transf) lChBs = iChBas(iSphCr(ll+i4))
          pb2T = DBLE(iPrmt(kOp,lChBs))
          Factor=pa1T*pb1T*pa2T*pb2T
          If (Factor.ne.One) Call DScal_(ijkl,Factor,
     &                                  AOInt(1,i1,i2,i3,i4),1)
 13      Continue
 12     Continue
 11    Continue
 10   Continue
 14   Continue
*
*     Call RecPrt(' Exit Phase: AOInt ',' ',AOInt,ijkl,
*    &            iCmp*jCmp*kCmp*lCmp)
      Return
      End

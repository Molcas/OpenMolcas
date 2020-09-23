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
* Copyright (C) 2004, Par Soderhjelm                                   *
************************************************************************
      Subroutine XFMoment(lMax,Cavxyz,Tmom,nCavxyz_,Org)

************************************************************************
*                                                                      *
*     Object:  Calculate total moment of XF multipoles, add to Cavxyz  *
*                                                                      *
*     Authors: P. Soderhjelm                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              November 2004                                           *
************************************************************************
      use External_Centers
      use Phase_Info
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"

      Real*8 Cavxyz(nCavxyz_),Tmom(nCavxyz_),Org(3)

      Real*8 Tco(3),A(3)
      Integer iStb(0:7), jCoSet(0:7,0:7)

      If(nOrd_XF.lt.0) Return
      If(nOrd_XF.gt.lMax) Then
         Call WarningMessage(2,'nOrd_XF.gt.lMax')
         Call Abend()
      EndIf
      nInp=(nOrd_XF+1)*(nOrd_XF+2)*(nOrd_XF+3)/6
      Call FZero(Org,3)
      do i=1,nXF
*
*------------- Generate Stabilizer of C
*
! IFG: "A" was undefined, is this the right point?
         A(1:3)=XF(1:3,i)
         iChxyz=iChAtm(A)
         iDum=0
         Call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)
         Do j = 0, nIrrep/nStb - 1
            Call FZero(Tmom,nCavxyz_)
            call dcopy_(nInp,XF(4,i),1,Tmom,1)
            TCo(1:3)=XF(1:3,i)
            Tco(1)=Tco(1)*DBLE(iPhase(1,jCoSet(j,0)))
            Tco(2)=Tco(2)*DBLE(iPhase(2,jCoSet(j,0)))
            Tco(3)=Tco(3)*DBLE(iPhase(3,jCoSet(j,0)))
            If(nOrd_XF.gt.0) Then
               Tmom(2)=Tmom(2)*DBLE(iPhase(1,jCoSet(j,0)))   !Dx
               Tmom(3)=Tmom(3)*DBLE(iPhase(2,jCoSet(j,0)))   !Dy
               Tmom(4)=Tmom(4)*DBLE(iPhase(3,jCoSet(j,0)))   !Dz
               If(nOrd_XF.gt.1) Then
                  Tmom(6)=Tmom(6)*DBLE(
     &                 iPhase(1,jCoSet(j,0))*iPhase(2,jCoSet(j,0))) !Qxy
                  Tmom(7)=Tmom(7)*DBLE(
     &                 iPhase(1,jCoSet(j,0))*iPhase(3,jCoSet(j,0))) !Qxz
                  Tmom(9)=Tmom(9)*DBLE(
     &                 iPhase(2,jCoSet(j,0))*iPhase(3,jCoSet(j,0))) !Qyz
               EndIf
            EndIf
            Call ReExpand(Tmom,1,nCavxyz_,Tco,Org,1,lMax)
            Call DaXpY_(nCavxyz_,One,Tmom,1,Cavxyz,1)
         EndDo
      EndDo

      Return
      End

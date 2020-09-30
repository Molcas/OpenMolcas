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
      Subroutine NucAtt_EMB(mGrid,Rho,nRho,P2_ontop,
     &                      nP2_ontop,iSpin,F_xc,dF_dRho,
     &                      ndF_dRho,dF_dP2ontop,ndF_dP2ontop)

      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "debug.fh"
#include "nsd.fh"
#include "setup.fh"
#include "nq_info.fh"
      Real*8  Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &        dF_dP2ontop(ndF_dP2ontop,mGrid),
     &        P2_ontop(nP2_ontop,mGrid), F_xc(mGrid)
*
*
      Call Get_nAtoms_All(mCenter)
      Call Allocate_Work(ip_RA,3*mCenter)
      Call Get_Coord_All(Work(ip_RA),mCenter)
*
      Call Allocate_Work(ip_ZA,mCenter)
      Call Get_iScalar('Unique atoms',nCenter)
*
      Call Allocate_iWork(ip_nStab,nCenter)
      Call Get_iArray('nStab',iWork(ip_nStab),nCenter)
      Call Allocate_Work(ip_Eff,nCenter)
      Call Get_dArray('Effective Nuclear Charge',Work(ip_Eff),nCenter)
      Call Get_iScalar('nSym',nSym)
*
      iOff = ip_ZA
      Do i = 1, nCenter
         n=nSym/iWork(ip_nStab+i-1)
         call dcopy_(n,Work(ip_Eff+i-1),0,Work(iOff),1)
         iOff = iOff + n
      End Do
*
      Call Free_Work(ip_Eff)
      Call Free_iWork(ip_nStab)
*
      Call Do_NucAtt_EMB(mGrid,Rho,nRho,P2_ontop,nP2_ontop,
     &                   iSpin,F_xc,dF_dRho,ndF_dRho,
     &                   dF_dP2ontop,ndF_dP2ontop,Work(ip_Grid),
     &                   Work(ip_RA),Work(ip_ZA),mCenter,T_X)
*
      Call Free_Work(ip_ZA)
      Call Free_Work(ip_RA)
*
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Do_NucAtt_EMB(mGrid,Rho,nRho,P2_ontop,nP2_ontop,
     &                         iSpin,F_xc,dF_dRho,ndF_dRho,
     &                         dF_dP2ontop,ndF_dP2ontop,Grid,RA,ZA,
     &                         mCenter,T_X)

      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid), dF_dRho(ndF_dRho,mGrid),
     &        dF_dP2ontop(ndF_dP2ontop,mGrid),
     &        P2_ontop(nP2_ontop,mGrid), F_xc(mGrid)
      Real*8 Grid(3,mGrid),RA(3,mCenter),ZA(mCenter)
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
        Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         DTot=Two*d_alpha
         If (DTot.lt.T_X) Go To 100
*
*------- Accumulate contributions to the nuclear attraction potential
*
         Attr=Zero
         Do i = 1, mCenter
            x=Grid(1,iGrid)-RA(1,i)
            y=Grid(2,iGrid)-RA(2,i)
            z=Grid(3,iGrid)-RA(3,i)
            Fact=ZA(i)/Sqrt(x**2+y**2+z**2)
            Attr=Attr+Fact
         End Do
         F_xc(iGrid)=F_xc(iGrid)-Attr*DTot
*
         dF_dRho(ipR,iGrid)=-Attr
*
100      Continue
*
        End Do
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
        Do iGrid = 1, mGrid
*
         d_alpha=Rho(1,iGrid)
         d_beta =Rho(2,iGrid)
         DTot=d_alpha+d_beta
         If (DTot.lt.T_X) Go To 200
*
*------- Accumulate contributions to the nuclear attraction potential
*
         Attr=Zero
         Do i = 1, mCenter
            x=Grid(1,iGrid)-RA(1,i)
            y=Grid(2,iGrid)-RA(2,i)
            z=Grid(3,iGrid)-RA(3,i)
            Fact=ZA(i)/Sqrt(x**2+y**2+z**2)
            Attr=Attr+Fact
         End Do
         F_xc(iGrid)=F_xc(iGrid)-Attr*DTot
*
         dF_dRho(ipRa,iGrid)=-Attr
         dF_dRho(ipRb,iGrid)=-Attr
*
200      Continue
*
        End Do
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End

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
* Copyright (C) 2002, Sergey Gusarov                                   *
*               2002, Roland Lindh                                     *
*               2002, Per Ake Malmqvist                                *
************************************************************************
      Subroutine Do_CS(Rho,nRho,mGrid,nP2,P2_ontop,P2II,
     &                 Coeff,iSpin,F_xc,Type1,Type2,Type3)
************************************************************************
*                                                                      *
* Object: Calculate functional from Colle-Salvetty approximation,      *
*         for MC wave function (basis formuls was written by PAM,      *
*         other by maple & GS (=/= CS!!!))                             *
*                                                                      *
*     Authors: Sergey Gusarov, University of Lund, SWEDEN,             *
*              Roland Lindh,   University of Lund, SWEDEN,             *
*              P.-A.Malmqvist, University of Lund, SWEDEN,             *
*              Dec. 2002                                               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Rho(nRho,mGrid), P2_ontop(nP2,mGrid),P2II(mGrid)
      Real*8 gradZ(3),
     &       gradBeta(3),
     &       dgradZdRho(3),
     &       dgradBetadRho(3),
     &       graddZdBeta(3),
     &       graddZdRho(3)
      Real*8 Type1(mGrid)
      Real*8 Type2(mGrid)
      Real*8 Type3(mGrid)
      Real*8 F_xc(mGrid)
*
*     Some functions:
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Some constants which we need...
*
      qCS = 2.29d0
      SqrtPi = sqrt(Pi)
      Pi32 = sqrt(Pi)*Pi
      Pi2  = Pi*Pi
      Pi52 = sqrt(Pi)*Pi2
      SqrtTwo = sqrt(Two)
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
         rho_alpha =Rho(1,iGrid)
         DTot    = Two*rho_alpha
         If (DTot.lt.1.0D-50) Go To 100
*
*        Beta:
*
         Beta  = qCS*Rho(1,iGrid)**(1.0d0/3.0d0)
         Beta2 = Beta*Beta
         Beta3 = Beta2*Beta
         Beta4 = Beta2*Beta2
         Beta5 = Beta2*Beta3
         Beta6 = Beta4*Beta2
        dBetadRho = qCS/(3.0d0*Rho(1,iGrid)**(2.0d0/3.0d0))
        Do i=1,3
          gradBeta(i) = qCS*Rho(1+i,iGrid)/
     &                  (3.0d0*Rho(1+i,iGrid)**(2.0d0/3.0d0))
*
          dgradBetadRho(i) = -2.0d0*qCS*Rho(1+i,iGrid)/
     &                    (9.0d0*Rho(1,iGrid)**(5.0d0/3.0d0))
        End Do
*
*
*
        D1  = (SqrtPi*Beta+1.D0)
        D2  = D1**2
        D3  = D1**3
        D4  = D1**4
*
        Y1  = -Two*Pi2*SqrtTwo +
     &         8.0d0*Pi52*Beta +
     &         9.0d0*Pi2 -
     &         24.0d0*Pi -
     &         32.0d0*Pi32*Beta
        Y   = Y1/(8.0d0*Beta2*D2)
*
        Z1  =  24.0d0*Pi52*Beta +
     &         26.0d0*Pi2 -
     &         3.0d0*Pi2*SqrtTwo -
     &         56.0d0*Pi -
     &         64.0d0*Pi32*Beta
        Z   = Z1/(96.0d0*Beta4*D2)
*
        Z2  = 60.0d0*Pi2*Beta2 +
     &        114.0d0*Pi32*Beta -
     &        160.0d0*Pi*Beta2 -
     &        264.0d0*SqrtPi*Beta +
     &        52.0d0*Pi -
     &        9.0d0*Pi32*SqrtTwo*Beta -
     &        6.0d0*Pi*SqrtTwo -
     &        112.0d0
*
        Z3  = -2328.0d0*Pi*Beta2 -
     &         1952.0d0*SqrtPi*Beta +
     &         978.0d0*Pi2*Beta2 -
     &         30.0d0*Pi*SqrtTwo +
     &         260.0d0*Pi +
     &         360.0d0*Pi52*Beta3 -
     &         960.0d0*Pi32*Beta3 -
     &         63.0d0*Pi2*Beta2*SqrtTwo -
     &         84.0d0*Pi32*SqrtTwo*Beta -
     &         560.0d0 +
     &         872.0d0*Pi32*Beta
*
        dZdBeta   = - Z2*Pi/(48.0d0*Beta5*D3)
        d2ZdBeta2 = Pi*Z3/(48.0d0*Beta6*D4)
        dZdRho    =  dZdBeta*dBetadRho
*
        Do i=1,3
*
          gradZ(i) = dZdBeta*gradBeta(i)
*
          dgradZdRho(i) = gradBeta(i)*d2ZdBeta2*dBetadRho +
     &                    dZdBeta*dgradBetadRho(i)
        End Do
*
*
        Z4 = 63.0d0*Pi2*Beta2*SqrtTwo +
     &       84.0d0*Pi32*SqrtTwo*Beta +
     &       560.0d0 +
     &       2328.0d0*Pi*Beta2 +
     &       1952.0d0*SqrtPi*Beta -
     &       260.0d0*Pi -
     &       978.0d0*Pi2*Beta2 +
     &       30.D0*Pi*SqrtTwo -
     &       360.0d0*Pi52*Beta3 -
     &       872.0d0*Pi32*Beta +
     &       960.0d0*Pi32*Beta3
       Z5 =  - Z4/(48.0d0*Beta6*D4)
       Do i=1,3
          graddZdBeta(i) = gradBeta(i)*Z5
       End Do
*
       Do i=1,3
          graddZdRho(i) = graddZdBeta(i)*dBetadRho +
     &                    dZdBeta*dgradBetadRho(i)
       End Do
*
        Y2 = -12.0d0*Pi2*Beta2 -
     &        22.0d0*Pi32*Beta +
     &        48.0d0*Pi*Beta2 +
     &        64.0d0*SqrtPi*Beta +
     &        4.0d0*Pi32*SqrtTwo*Beta +
     &        2.0d0*Pi*SqrtTwo -
     &        9.0d0*Pi +
     &        24.0d0
        dYdRho = Y2*Pi*dBetadRho/(4.0d0*Beta3*D3)
*
        grad2Z = Zero
*
        ii=0
        Do i=1,3
           If(i.eq.1) ii=5
           If(i.eq.2) ii=8
           If(i.eq.3) ii=9
           Betaxx = - Two*qCS*Rho(1+i,iGrid)**2/
     &                (9.0d0*Rho(1,iGrid)**(5.0d0/3.0d0)) +
     &                qCS*Rho(ii,iGrid)/
     &                (3.0d0*Rho(1,iGrid)**(2.0d0/3.0d0))
*
           Betaxx1= gradBeta(i)*gradBeta(i)
*
           Zxx = 60.0d0*Betaxx*Pi52*Beta4 +
     &           174.0d0*Betaxx*Pi2**Beta3 +
     &           166.0d0*Betaxx*Pi32*Beta2 -
     &           160.0d0*Betaxx*Pi32*Beta4 -
     &           424.0d0*Betaxx*Pi*Beta3 -
     &           376.0d0*Betaxx*SqrtPi*Beta2 -
     &           6.0d0*Betaxx*Pi*Beta*SqrtTwo -
     &           360.0d0*Betaxx1*Pi52*Beta3 -
     &           978.0d0*Betaxx1*Pi2*Beta2 -
     &           872.0d0*Betaxx1*Pi32*Beta +
     &           960.0d0*Betaxx1*Pi32*Beta3 +
     &           2328.0d0*Betaxx1*Pi*Beta2 +
     &           1952.0d0*Betaxx1*SqrtPi*Beta +
     &           63.0d0*Betaxx1*Pi2*Beta2*SqrtTwo +
     &           52.0d0*Betaxx*Pi*Beta +
     &           30.0d0*Betaxx1*Pi*SqrtTwo -
     &           112.0d0*Betaxx*Beta +
     &           560.0d0*Betaxx1 -
     &           260.0d0*Betaxx1*Pi +
     &           84.0d0*Betaxx1*Pi32*Beta*SqrtTwo -
     &           15.0d0*Betaxx*Pi32*Beta2*SqrtTwo -
     &           9.0d0*Betaxx*Pi2*Beta3*SqrtTwo
*
           grad2Z = grad2Z - Pi*Zxx/(48.0d0*Beta6*D4)
        End Do
************************************************************************
*                                                                      *
*    We have Y, dY/dRho, Z, grad Z, dZ/dRho, d grad(Z)/dRho, grad^2 Z  *
*       Rho type
*       From d /dRho:
*
        Type1(iGrid) = 0.5d0*(dYdRho*P2_ontop(1,iGrid) -
     &                 0.25d0*ddot_(3,dgradZdRho,1,
     &                        P2_ontop(2,iGrid),mGrid) -
     &                 2.0d0*dZdRho*P2II(iGrid) )
*
*       P2 type
*       From (- grad d /d gradRho):
*
        Type1(iGrid) = Type1(iGrid) + 0.25d0*(
     &                 ddot_(3,graddZdRho,1,
     &                 P2_ontop(2,iGrid),mGrid) +
     &                 dZdRho*(
     &                 P2_ontop( 5,iGrid) +
     &                 P2_ontop( 8,iGrid) +
     &                 P2_ontop(10,iGrid) )  )
*
*    From d / dP2  and grad (d / d(gradP2))
*
        Type2(iGrid) = .5d0*(Y + 0.25d0*grad2Z)
*
*      Other type:
*
        Type3(iGrid) = - Z
*
*
*
        F_xc(iGrid) = 0.50d0*(Y*P2_ontop(1,iGrid) -
     &               0.250d0*ddot_(3,gradZ,1,P2_ontop(2,iGrid),mGrid) -
     &               2.0d0*Z*P2II(iGrid))
*
 100     Continue
      End Do    ! iGrid




************************************************************************
*                                                                      *
*     iSpin=/=1
*
      Else
         Call WarningMessage(2,
     &              'Do_CS: Not implemented yet for iSpin=/=1')
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
*
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(Coeff)
      End

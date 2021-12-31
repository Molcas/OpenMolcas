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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Dec. 08, 2021, created this file.               *
* ****************************************************************
      Subroutine Calc_Pot2(Pot2,mGrid,Pi,nPi)
      use nq_Grid, only: Weights
      use nq_pdft
#include "nq_info.fh"
#include "stdalloc.fh"
******Input
      INTEGER mGrid,nPi
      Real*8,DIMENSION(nPi,mGrid)::Pi
******Output
      Real*8,DIMENSION(nPot2)::Pot2
******Internal
      INTEGER iGrid,nGOrb
      Real*8 ThrsZ2,ThrsPi,ggaterm,ftggaterm,predEdPip

      ThrsZ2=1.0d-15
      ThrsPi=1.0d-30

      IF(lGGA.and.lft) THEN
       CALL FZero(dEdPix,mGrid)
       CALL FZero(dEdPiy,mGrid)
       CALL FZero(dEdPiz,mGrid)
       CALL FZero(GdEdPiMO,mGrid*nOrbt)
      END IF


      DO iGrid=1,mGrid
       IF(Pass1(iGrid).and.(Pi(1,iGrid).gt.ThrsPi)) THEN
        If(Pass2(iGrid).or.Pass3(iGrid)) Then
         if(lGGA) then
          ggaterm=GradRhodFdRho(iGrid)
          if(lft) then
           ftggaterm=(d2ZdR2(iGrid)*dRdPi(iGrid)*GradRdFdRho(iGrid)+
     &           d2RdRhodPi(iGrid)*dZdR(iGrid)*GradRhodFdRho(iGrid))*
     &           RhoAB(iGrid)
          predEdPip=RhoAB(iGrid)*dZdR(iGrid)*dRdPi(iGrid)*Weights(iGrid)
           dEdPix(iGrid)=predEdPip*dF_dRhoxamb(iGrid)
           dEdPiy(iGrid)=predEdPip*dF_dRhoyamb(iGrid)
           dEdPiz(iGrid)=predEdPip*dF_dRhozamb(iGrid)
          else
           ftggaterm=0.0d0
          end if
         else
          ggaterm=0.0d0
          ftggaterm=0.0d0
         end if
          dEdPi(iGrid)=Weights(iGrid)*(dZdR(iGrid)*dRdPi(iGrid)*
     &   (RhoAB(iGrid)*dF_dRhoamb(iGrid)+ggaterm)+ftggaterm)
        Else
         dEdPi(iGrid)=0.0d0
        End If
       ELSE
        dEdPi(iGrid)=0.0d0
       END IF
      END DO
      nGOrb=mGrid*nOrbt


      CALL DSCal_(mGrid,0.5d0,dEdPi,1)
      IF(lGGA.and.lft) THEN
       CALL DSCal_(mGrid,0.5d0,dEdPix,1)
       CALL DSCal_(mGrid,0.5d0,dEdPiy,1)
       CALL DSCal_(mGrid,0.5d0,dEdPiz,1)
      END IF


      CALL DCopy_(nGOrb,MOas,1,dEdPiMO,1)

      DO iGrid=1,mGrid
       CALL DScal_(nOrbt,dEdPi(iGrid),dEdPiMO(iGrid),mGrid)
      END DO

      IF(lft.and.lGGA) THEN
       DO iGrid=1,mGrid
        CALL DAXpY_(nOrbt,dEdPix(iGrid),MOax(iGrid),mGrid,
     &                              GdEdPiMO(iGrid),mGrid)
        CALL DAXpY_(nOrbt,dEdPiy(iGrid),MOay(iGrid),mGrid,
     &                              GdEdPiMO(iGrid),mGrid)
        CALL DAXpY_(nOrbt,dEdPiz(iGrid),MOaz(iGrid),mGrid,
     &                              GdEdPiMO(iGrid),mGrid)
       END DO
       CALL DAXpY_(nGOrb,1.0d0,GdEdPiMO,1,dEdPiMO,1)
      END IF

*     dEdPiMO is practically (Phi_p*dEdPi+Phi_p'*dEdPi')
*     The subroutine below calculates
*     (Phi_p*dEdPi+Phi_p'*dEdPi')*Phi_u*Phi_v*Phi_x
      CALL Calc_Pot2_Inner(Pot2,mGrid,dEdPiMO,MOas,MOas,MOas,.false.)

      IF(lft.and.lGGA) THEN
       CALL Calc_Pot2_Inner(Pot2,mGrid,MOas,MOas,MOas,GdEdPiMO,.true.)
      END IF

      RETURN
      End Subroutine



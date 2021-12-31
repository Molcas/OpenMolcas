*********************************************************************
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
* Jie J. Bao, on Dec. 18, 2021, created this file.               *
* ****************************************************************
      Subroutine TranslateDens(Pi,dRho_dr,dPi,l_tanhr,nRho,mGrid,
     &                         nPi,ndRho_dr,nEGrad,DoGrad)
      use nq_Grid, only: Rho, GradRho, nGradRho
      use nq_pdft, only: ThrsRho,ThrsOMR,ThrsFT,ThrsNT,fta,ftb,ftc,
     &                 Pass1,Pass2,Pass3,
     &                 OnePZ,OneMZ,RatioA,ZetaA,RhoAB,dZdR,
     *                 lft,lGGA
******Input
      INTEGER nRho,mGrid,nPi,ndRho_dr,nEGrad
      Real*8,DIMENSION(nPi,mGrid)::Pi
      Real*8,DIMENSION(nPi,nEGrad,mGrid)::dPi
      Logical DoGrad,l_tanhr
******Input & Output
      Real*8,DIMENSION(ndRho_dr,mGrid,nEGrad)::dRho_dr
******In-subroutine
      INTEGER iGrid,iEGrad,ngragri,iOff1,nGRho
      Real*8 TempR,RRatio,Diff1,Rd2ZdR2,RdRdRho,RdRdPi,Rd2RdRho2,
     &       Rd2RdRhodPi,Rd2ZdRdZ,XAdd,YAdd,ZAdd,
C     &       RatioX,RatioY,RatioZ,GraddZdR,
     &       GraddZdR,
     &       GradRatio,GradRatioX,GradRatioY,GradRatioZ,
     &       ZetaX,ZetaY,ZetaZ,GradZetaX,GradZetaY,GradZetaZ

      Real*8,DIMENSION(mGrid)::dRhodx,dRhody,dRhodz,
     &                         tanhrx,tanhry,tanhrz,ftx23,fty23,ftz23,
     &                         RatioX,RatioY,RatioZ
      Real*8,DIMENSION(mGrid*nEGrad)::GradRhoAB,GradRhoX,GradRhoY,
     &                                GradRhoZ,dRatio,dZeta
******PassX
*     Pass1. Total density is greater than thresRho
*     Pass2. Do translation
*     Pass3. Do full translation
*     if Pass1 is false, Pass2, 3, are both false.
*     if Pass1 is true, Pass2 and 3 cannot both be true (can both be
*     false).
*********************************************************************

      nGRho=nGradRho
*********************************************************************
*     calculating total density at each grid
*********************************************************************
      CALL DCopy_(mGrid,Rho(1,1),nRho,RhoAB,1)
      CALL DAXPY_(mGrid,1.0d0,Rho(2,1),nRho,RhoAB,1)


*********************************************************************
*     calculating x, y, z components of density gradient
*********************************************************************
      IF(lGGA) THEN
       CALL DCopy_(mGrid,GradRho(1,1),nGRho,dRhodx,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(4,1),nGRho,dRhodx,1)
       CALL DCopy_(mGrid,GradRho(2,1),nGRho,dRhody,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(5,1),nGRho,dRhody,1)
       CALL DCopy_(mGrid,GradRho(3,1),nGRho,dRhodz,1)
       CALL DAXPY_(mGrid,1.0d0,GradRho(6,1),nGRho,dRhodz,1)
      END IF


*********************************************************************
*    Ratio and Zeta at each grid point
*********************************************************************
      CALL FZero( ZetaA,mGrid)
      CALL FZero(RatioA,mGrid)
      CALL FZero(dZdR  ,mGrid)
      DO iGrid=1,mGrid
       Pass1(iGrid)=.false.
       Pass2(iGrid)=.false.
       Pass3(iGrid)=.false.
      END DO

      IF(.not.lft) THEN
       DO iGrid=1,mGrid
        If(RhoAB(iGrid).ge.ThrsRho) Then
         Pass1(iGrid)=.true.
         RRatio=4.0d0*Pi(1,iGrid)/(RhoAB(iGrid)**2)
         RatioA(iGrid)=Rratio
         if(l_tanhr) RRatio=tanh(RRatio)
         if((1.0d0-Rratio).gt.ThrsOMR) then
          ZetaA(iGrid)=sqrt(1.0d0-Rratio)
          Pass2(iGrid)=.true.
          dZdR(iGrid)=-0.5d0/ZetaA(iGrid)
         end if
        End If
       END DO
      ELSE
       DO iGrid=1,mGrid
        If(RhoAB(iGrid).ge.ThrsRho) Then
         Pass1(iGrid)=.true.
         RRatio=4.0d0*Pi(1,iGrid)/(RhoAB(iGrid)**2)
         RatioA(iGrid)=Rratio
         if(RRatio.lt.ThrsFT) then ! do t-translation
          ZetaA(iGrid)=sqrt(1.0d0-Rratio)
          Pass2(iGrid)=.true.
          dZdR(iGrid)=-0.5d0/ZetaA(iGrid)
         else if(RRatio.le.ThrsNT) then ! do ft-translation
          Diff1=RRatio-ThrsNT
          ZetaA(iGrid)=(fta*Diff1**2+ftb*Diff1+ftc)*Diff1**3
          Pass3(iGrid)=.true.
          dZdR(iGrid)=
     &    (5.0d0*fta*Diff1**2+4.0d0*ftb*Diff1+3.0d0*ftc)*Diff1**2
         end if
        End If
       END DO
      END IF

*********************************************************************
*    (1 + zeta)/2 and (1 - zeta)/2
*********************************************************************
      CALL DCopy_(mGrid,[0.5d0],0,OnePZ,1)
      CALL DCopy_(mGrid,[0.5d0],0,OneMZ,1)
      CALL DAXPY_(mGrid, 0.5d0,ZetaA,1,OnePZ,1)
      CALL DAXPY_(mGrid,-0.5d0,ZetaA,1,OneMZ,1)


*********************************************************************
*     translating rho_a and rho_b
*********************************************************************
      DO iGrid=1,mGrid
       IF(Pass1(iGrid)) THEN
        Rho(1,iGrid)=OnePZ(iGrid)*RhoAB(iGrid)
        Rho(2,iGrid)=OneMZ(iGrid)*RhoAB(iGrid)
       END IF
      END DO

*********************************************************************
*     translating gradient component of rho_a and rho_b
*********************************************************************
      IF(lGGA) THEN
       DO iGrid=1,mGrid
        If(Pass1(iGrid)) Then
         GradRho(1,iGrid)=OnePZ(iGrid)*dRhodX(iGrid)
         GradRho(2,iGrid)=OnePZ(iGrid)*dRhodY(iGrid)
         GradRho(3,iGrid)=OnePZ(iGrid)*dRhodZ(iGrid)
         GradRho(4,iGrid)=OneMZ(iGrid)*dRhodX(iGrid)
         GradRho(5,iGrid)=OneMZ(iGrid)*dRhodY(iGrid)
         GradRho(6,iGrid)=OneMZ(iGrid)*dRhodZ(iGrid)
        End If
       END DO

       If(lft) Then
        DO iGrid=1,mGrid
         if(Pass1(iGrid)) then
          RatioX(iGrid)=(4.0d0*Pi(2,iGrid)/RhoAB(iGrid)-
     &                   2.0d0*RatioA(iGrid)*dRhodX(iGrid))/RhoAB(iGrid)
          RatioY(iGrid)=(4.0d0*Pi(3,iGrid)/RhoAB(iGrid)-
     &                   2.0d0*RatioA(iGrid)*dRhodY(iGrid))/RhoAB(iGrid)
          RatioZ(iGrid)=(4.0d0*Pi(4,iGrid)/RhoAB(iGrid)-
     &                   2.0d0*RatioA(iGrid)*dRhodZ(iGrid))/RhoAB(iGrid)
         else
          RatioX(iGrid)=0.0d0
          RatioY(iGrid)=0.0d0
          RatioZ(iGrid)=0.0d0
         end if
         ftx23(iGrid)=0.5d0*RhoAB(iGrid)*dZdR(iGrid)*RatioX(iGrid)
         fty23(iGrid)=0.5d0*RhoAB(iGrid)*dZdR(iGrid)*RatioY(iGrid)
         ftz23(iGrid)=0.5d0*RhoAB(iGrid)*dZdR(iGrid)*RatioZ(iGrid)
        END DO
        CALL DaXpY_(mGrid, 1.0d0,ftx23,1,GradRho(1,1),6)
        CALL DaXpY_(mGrid, 1.0d0,fty23,1,GradRho(2,1),6)
        CALL DaXpY_(mGrid, 1.0d0,ftz23,1,GradRho(3,1),6)
        CALL DaXpY_(mGrid,-1.0d0,ftx23,1,GradRho(4,1),6)
        CALL DaXpY_(mGrid,-1.0d0,fty23,1,GradRho(5,1),6)
        CALL DaXpY_(mGrid,-1.0d0,ftz23,1,GradRho(6,1),6)
       End If
      END IF

*********************************************************************
*     Additional terms in the tanh translation
*********************************************************************
      IF(l_tanhr) THEN
       CALL FZero(tanhrx,mGrid)
       CALL FZero(tanhry,mGrid)
       CALL FZero(tanhrz,mGrid)
       DO iGrid=1,mGrid
        If(Pass1(iGrid)) Then
         RRatio=RatioA(iGrid)
         TempR=4.0d0*Pi(1,iGrid)/RhoAB(iGrid)
         TanhrX(iGrid)=(RRatio**2-1.0d0)*(Pi(2,iGrid)-
     &(dRhodX(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
         TanhrY(iGrid)=(RRatio**2-1.0d0)*(Pi(3,iGrid)-
     &(dRhodY(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
         TanhrZ(iGrid)=(RRatio**2-1.0d0)*(Pi(4,iGrid)-
     &(dRhodZ(iGrid)*TempR))/(RhoAB(iGrid)*ZetaA(iGrid))
        End If
       END DO
       CALL DAXPY_(mGrid, 1.0d0,TanhrX,1,GradRho(1,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrX,1,GradRho(4,1),nRho)
       CALL DAXPY_(mGrid, 1.0d0,TanhrY,1,GradRho(2,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrY,1,GradRho(5,1),nRho)
       CALL DAXPY_(mGrid, 1.0d0,TanhrZ,1,GradRho(3,1),nRho)
       CALL DAXPY_(mGrid,-1.0d0,TanhrZ,1,GradRho(6,1),nRho)
      END IF



*********************************************************************
*     calculating terms needed in gradient calculation
*********************************************************************
*     if not doing gradient, code ends here
      IF(.not.DoGrad) RETURN
*********************************************************************
*     calculating density gradient wrt geometrical changes
*********************************************************************
      ngragri=mGrid*nEGrad
      CALL DCopy_(ngragri,dRho_dr(1,1,1),ndRho_dr,GradRhoAB,1)
      CALL DAXPY_(ngragri,1.0d0,dRho_dr(2,1,1),ndRho_dr,GradRhoAB,1)

      IF(lGGA) Then
       CALL DCopy_(ngragri,dRho_dr(3,1,1),ndRho_dr,GradRhoX,1)
       CALL DAXPY_(ngragri,1.0d0,dRho_dr(6,1,1),ndRho_dr,GradRhoX,1)
       CALL DCopy_(ngragri,dRho_dr(4,1,1),ndRho_dr,GradRhoY,1)
       CALL DAXPY_(ngragri,1.0d0,dRho_dr(7,1,1),ndRho_dr,GradRhoY,1)
       CALL DCopy_(ngragri,dRho_dr(5,1,1),ndRho_dr,GradRhoZ,1)
       CALL DAXPY_(ngragri,1.0d0,dRho_dr(8,1,1),ndRho_dr,GradRhoZ,1)
      END IF

*********************************************************************
*    dRatio and dZeta at each grid point
*********************************************************************
*     Calculate dRatio
      CALL Fzero(dRatio,nGraGri)
      DO iGrid=1,mGrid
       IF(Pass1(iGrid)) THEN
        Do iEGrad=1,nEGrad
         IOff1=(iEGrad-1)*mGrid
         dRatio(IOff1+iGrid)=4.0d0*dPi(1,iEGrad,iGrid)/(RhoAB(iGrid)**2)
     &       -8.0d0*Pi(1,iGrid)*GradRhoAB(IOff1+iGrid)/(RhoAB(iGrid)**3)
        End Do
       END IF
      END DO
*     alculate dZeta
      CALL Fzero(dZeta,nGraGri)
      DO iGrid=1,mGrid
       CALL DAxpy_(nEGrad,dZdR(iGrid),dRatio(iGrid),mGrid,
     &                                dZeta(iGrid),mGrid)
      END DO

      DO iEGrad=1,nEGrad
       IOff1=(iEGrad-1)*mGrid
       Do iGrid=1,mGrid
        If(Pass1(iGrid)) Then
         dRho_dr(1,iGrid,iEGrad)=OnePZ(iGrid)*GradRhoAB(IOff1+iGrid)+
     &                          0.50d0*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
         dRho_dr(2,iGrid,iEGrad)=OneMZ(iGrid)*GradRhoAB(IOff1+iGrid)-
     &                          0.50d0*dZeta(IOFf1+iGrid)*RhoAB(iGrid)
        End If
       End Do
      END DO

      IF(lGGA) THEN
       DO iEGrad=1,nEGrad
        IOff1=(iEGrad-1)*mGrid
        Do iGrid=1,mGrid
         If(Pass1(iGrid)) Then
          dRho_dr(3,iGrid,iEGrad)=OnePZ(iGrid)*GradRhoX(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
          dRho_dr(6,iGrid,iEGrad)=OneMZ(iGrid)*GradRhoX(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodx(iGrid)
          dRho_dr(4,iGrid,iEGrad)=OnePZ(iGrid)*GradRhoY(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhody(iGrid)
          dRho_dr(7,iGrid,iEGrad)=OneMZ(iGrid)*GradRhoY(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhody(iGrid)
          dRho_dr(5,iGrid,iEGrad)=OnePZ(iGrid)*GradRhoZ(IOff1+iGrid)+
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
          dRho_dr(8,iGrid,iEGrad)=OneMZ(iGrid)*GradRhoZ(IOff1+iGrid)-
     &                           0.50d0*dZeta(IOFf1+iGrid)*dRhodz(iGrid)
         End If
        End Do
       END DO
       If(lft) Then
        DO iGrid=1,mGrid
         if(.not.Pass1(iGrid)) cycle
         if(.not.(Pass2(iGrid).or.Pass3(iGrid))) cycle
         ZetaX=dZdR(iGrid)*RatioX(iGrid)
         ZetaY=dZdR(iGrid)*RatioY(iGrid)
         ZetaZ=dZdR(iGrid)*RatioZ(iGrid)
         RdRdRho=-2.0d0*RatioA(iGrid)/RhoAB(iGrid)
         RdRdPi=4.0d0/RhoAB(iGrid)**2
         Rd2RdRho2=-3.0d0*RdRdRho/RhoAB(iGrid)
         Rd2RdRhodPi=-2.0d0*RdRdPi/RhoAB(iGrid)
         Rd2ZdRdZ=0.0d0
         Rd2ZdR2=0.0d0
         if(Pass2(iGrid)) Rd2ZdRdZ=0.5d0/ZetaA(iGrid)**2
         if(Pass3(iGrid)) then
          Diff1=RatioA(iGrid)-ThrsNT
          Rd2ZdR2=(2.0d1*fta*Diff1**2+1.2d1*ftb*Diff1+6.0d0*ftc)*Diff1
         end if
         Do iEGrad=1,nEGrad
          IOff1=(iEGrad-1)*mGrid
          GradRatio=dRatio(iOff1+iGrid)
          GradRatioX=(Rd2RdRho2*dRhodX(iGrid)+Rd2RdRhodPi*Pi(2,iGrid))*
     &               GradRhoAB(iOff1+iGrid)+
     &               Rd2RdRhodPi*dRhodX(iGrid)*dPi(1,iEGrad,iGrid)+
     &               RdRdRho*GradRhoX(iOff1+iGrid)+
     &               RdRdPi*dPi(2,iEGrad,iGrid)

          GradRatioY=(Rd2RdRho2*dRhodY(iGrid)+Rd2RdRhodPi*Pi(3,iGrid))*
     &               GradRhoAB(iOff1+iGrid)+
     &               Rd2RdRhodPi*dRhodY(iGrid)*dPi(1,iEGrad,iGrid)+
     &               RdRdRho*GradRhoY(iOff1+iGrid)+
     &               RdRdPi*dPi(3,iEGrad,iGrid)

          GradRatioZ=(Rd2RdRho2*dRhodZ(iGrid)+Rd2RdRhodPi*Pi(4,iGrid))*
     &               GradRhoAB(iOff1+iGrid)+
     &               Rd2RdRhodPi*dRhodZ(iGrid)*dPi(1,iEGrad,iGrid)+
     &               RdRdRho*GradRhoZ(iOff1+iGrid)+
     &               RdRdPi*dPi(4,iEGrad,iGrid)
          GraddZdR=0.0d0
          if(Pass2(iGrid)) then
           GraddZdR=Rd2ZdRdZ*dZeta(iOff1+iGrid)
          else if(Pass3(iGrid)) then
           GraddZdR=Rd2ZdR2*GradRatio
          end if

          GradZetax=GraddZdR*RatioX(iGrid)+dZdR(iGrid)*GradRatioX
          GradZetaY=GraddZdR*RatioY(iGrid)+dZdR(iGrid)*GradRatioY
          GradZetaZ=GraddZdR*RatioZ(iGrid)+dZdR(iGrid)*GradRatioZ

      XAdd=0.5d0*(RhoAB(iGrid)*GradZetaX+ZetaX*GradRhoAB(iOff1+iGrid))
      YAdd=0.5d0*(RhoAB(iGrid)*GradZetaY+ZetaY*GradRhoAB(iOff1+iGrid))
      ZAdd=0.5d0*(RhoAB(iGrid)*GradZetaZ+ZetaZ*GradRhoAB(iOff1+iGrid))

          dRho_dr(3,iGrid,iEGrad)=dRho_dr(3,iGrid,iEGrad)+XAdd
          dRho_dr(6,iGrid,iEGrad)=dRho_dr(6,iGrid,iEGrad)-XAdd
          dRho_dr(4,iGrid,iEGrad)=dRho_dr(4,iGrid,iEGrad)+YAdd
          dRho_dr(7,iGrid,iEGrad)=dRho_dr(7,iGrid,iEGrad)-YAdd
          dRho_dr(5,iGrid,iEGrad)=dRho_dr(5,iGrid,iEGrad)+ZAdd
          dRho_dr(8,iGrid,iEGrad)=dRho_dr(8,iGrid,iEGrad)-ZAdd
         End Do
        END DO
       End If
      END IF

      RETURN
      END SUBROUTINE

***********************************************************************

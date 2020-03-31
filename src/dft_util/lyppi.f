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
* Copyright (C) Per Ake Malmqvist                                      *
*               Ajitha Devarajan                                       *
************************************************************************
      Subroutine LYPPI(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                 Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object:  Lyp Functional(Formula taken from Molpro Manual)            *
*                                                                      *
* Called from:Do_batch                                                 *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Author: Per-AAke Malmquist,Department of Theoretical Chemistry  *
*              University of LUnd, SWEDEN                              *
*              D. Ajitha , Department of Theoretical Chemistry         *
*              University of Lund, SWEDEN                              *
*              Modify Per-AAke's code for open shell case              *
*              and adopt for closed shell case                         *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),Rho(nRho,mGrid),F_xc(mGrid)
*
      dimension ec(0:4)
      data Cfconst / 2.8712340001881918D0 /
      data aconst,bconst,cconst,dconst
     &        /0.04918d0,0.132d0,0.2533d0,0.349d0/
* Constants aconst,..,dconst are from ???
* Cfconst is = (3/10)*(3*Pi**2)**(2/3)
      Cfconst2=Cfconst*2.0D0**(11.0D0/3)
      Rho_min=T_X*1.0D-2
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
      aConst=0.049180d0
      bConst=0.1320d0
      dConst=0.3490d0
      cConst=0.25330d0
      Cf=2.8712340001881918D0
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
************************************************************************
      Call dScal_(nRho*iSpin*mGrid,2.0d0,Rho,1)
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
      Dns = Rho(ipR,iGrid)
      if(Dns.lt.T_X) Go To 100
      gradRho2=ddot_(3,Rho(ipdRx,iGrid),1,Rho(ipdRx,iGrid),1)
      gradRho=sqrt(gradRho2)
      grad2Rho = Rho(ipL,iGrid)

      t1 = aConst*bConst
      t2 = Dns**(1.D0/3.D0)
      t3 = 1/t2
      t5 = exp(-cConst*t3)
      t7 = 1+dConst*t3
      t8 = 1/t7
      t9 = t5*t8
      t10 = t2**2
      t13 = t9/t10/Dns
      t16 = Dns**2
      t17 = 1/t16
      t18 = t1*t17
      t19 = cConst*t5
      t23 = t7**2
      t24 = 1/t23
      t25 = t5*t24
      t29 = 7.D0/36.D0*t1*t13-7.D0/72.D0*t18*t19*t8-7.D0/72.D0*t18*t25*d
     &Const
      t30 = gradRho2
      t34 = t10*t16
      t38 = t1*(Cf*t34-17.D0/72.D0*t30)
      Fxc = -t29*t30-aConst*Dns*t8-t38*t13
      t40 = t1*cConst
      t41 = t16*Dns
      t42 = 1/t41
      t47 = t1*t5
      t52 = 1/t34
      t53 = t9*t52
      t58 = t1/t2/t41
      t59 = cConst**2
      t64 = t24*dConst
      t71 = dConst**2
      t76 = (7.D0/27.D0*t40*t42*t5*t8+7.D0/27.D0*t47*t24*t42*dConst-35.D
     &0/108.D0*t1*t53-7.D0/216.D0*t58*t59*t5*t8-7.D0/108.D0*t58*t19*t64-
     &7.D0/108.D0*t58*t5/t23/t7*t71)*t30
      dFdDns = -t76-aConst*t8-aConst*t3*t64/3-8.D0/3.D0*t1*Cf*t5*t8-t38*
     &cConst*t42*t9/3-t38*t25*t42*dConst/3+5.D0/3.D0*t38*t53
      t100 = t2+dConst
      t101 = 1/t100
      t117 = t100**2
      Grad = 2*t29*grad2Rho+2*t76-17.D0/36.D0*t47/t2/Dns*t101*grad2Rho-2
     &*(17.D0/216.D0*t40*t52*t5*t101-17.D0/54.D0*t1*t5/t2/t16*t101-17.D0
     &/216.D0*t1*t5*t17/t117)*t30

*
      F_xc(iGrid) = Fxc
      dF_dRho(ipR,iGrid)= dFdDns + Grad
*
100   Continue
      End Do

*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/= 1
      Else
*                                                                      *
************************************************************************
*                                                                      *
*
      Do iGrid = 1, mGrid
      rhoa=Max(Rho_min,Rho(ipRa,iGrid))
      rhob=Max(Rho_min,Rho(ipRb,iGrid))
      rho_tot=rhoa+rhob
      if(rho_tot.lt.T_X) Go To 200
      gxa=Rho(ipdRxa,iGrid)
      gya=Rho(ipdRya,iGrid)
      gza=Rho(ipdRza,iGrid)
      gxb=Rho(ipdRxb,iGrid)
      gyb=Rho(ipdRyb,iGrid)
      gzb=Rho(ipdRzb,iGrid)

      gx=gxa+gxb
      gy=gya+gyb
      gz=gza+gzb

      rho3=rho_tot**(-1.0D0/3)
      crho3=cconst*rho3
      R=1/(1+dconst*rho3)
      if(crho3.lt.709.0d0) then
        expcr=exp(-crho3)
      else
        expcr=0.0D0
      end if
      omega=expcr*R*rho_tot**(-11.0D0/3)
      delta=crho3+1-R
      dlogodr=(delta-11)/(3*rho_tot)
      ddeltadr=-(crho3+R*(1-R))/(3*rho_tot)
      p=aconst*bconst*omega

      ec1=-(4*aconst*R)*rhoa*(rhob/rho_tot)

      dec1dra=(ec1*(1-R)/(3*rho_tot))-(4*aconst*R)*(rhob/rho_tot)**2
      dec1drb=(ec1*(1-R)/(3*rho_tot))-(4*aconst*R)*(rhoa/rho_tot)**2

      tmp1=-Cfconst2*p*rhoa**(11.0D0/3)*rhob
      tmp2=-Cfconst2*p*rhob**(11.0D0/3)*rhoa
      ec2=tmp1+tmp2
      dec2dra=(ec2*dlogodr)+(11*tmp1+3*tmp2)/(3*rhoa)
      dec2drb=(ec2*dlogodr)+(3*tmp1+11*tmp2)/(3*rhob)

      sa=gxa**2+gya**2+gza**2
      sb=gxb**2+gyb**2+gzb**2
      s =gx**2+gy**2+gz**2
      pp=-p*rhoa*rhob/(18*rho_tot)
      dppdra=((dlogodr-1/rho_tot)+1/rhoa)*pp
      dppdrb=((dlogodr-1/rho_tot)+1/rhob)*pp
      dqq1ds =47*rho_tot
      dqq1dsa=22*rhoa-45*rho_tot
      dqq1dsb=22*rhob-45*rho_tot
      dqq2ds =-7*rho_tot
      dqq2dsa= rhob-rhoa
      dqq2dsb= rhoa-rhob
      dqqds =dqq1ds +delta*dqq2ds
      dqqdsa=dqq1dsa+delta*dqq2dsa
      dqqdsb=dqq1dsb+delta*dqq2dsb
      dqq1dra=47*s-23*sa-45*sb
      dqq1drb=47*s-45*sa-23*sb
      qq1=dqq1dra*rhoa+dqq1drb*rhob
      dqq2dra=-7*s -sa +sb
      dqq2drb=-7*s +sa -sb
      qq2=dqq2dra*rhoa+dqq2drb*rhob
      dqqdra=dqq1dra+dqq2dra*delta+qq2*ddeltadr
      dqqdrb=dqq1drb+dqq2drb*delta+qq2*ddeltadr
      qq=qq1+qq2*delta
      ec3=pp*qq
      dec3dra=dppdra*qq+pp*dqqdra
      dec3drb=dppdrb*qq+pp*dqqdrb
      dec3ds =pp*dqqds
      dec3dsa=pp*dqqdsa
      dec3dsb=pp*dqqdsb

      dqdra=((4.0D0/3)*rho_tot*(s-sa-sb))+2*rhoa*sb
      dqdrb=((4.0D0/3)*rho_tot*(s-sa-sb))+2*rhob*sa
      q=(dqdra*rhoa+dqdrb*rhob)/2
      ec4=p*q
      dec4dra=dlogodr*ec4+p*dqdra
      dec4drb=dlogodr*ec4+p*dqdrb
      dec4ds =p*(2.0D0/3)*(rho_tot**2)
      dec4dsa=p*(rhob**2)-dec4ds
      dec4dsb=p*(rhoa**2)-dec4ds

      dec34ds =dec3ds +dec4ds
      dec34dsa=dec3dsa+dec4dsa
      dec34dsb=dec3dsb+dec4dsb

      ec(0)=ec1+ec2+ec3+ec4
      ec(1)=ec1
      ec(2)=ec2
      ec(3)=ec3
      ec(4)=ec4

      Functional= ec1+ec2+ec3+ec4
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*                                                                      *
************************************************************************
*                                                                      *
*
      dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
      dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)
     &                   +Coeff*(dec1drb+dec2drb+dec3drb+dec4drb)

      dF_dRho(ipdRxa,iGrid)=dF_dRho(ipdRxa,iGrid)
     &                   +2*Coeff*(dec34ds*gx+dec34dsa*gxa)
      dF_dRho(ipdRya,iGrid)=dF_dRho(ipdRy,iGrid)
     &                   +2*Coeff*(dec34ds*gy+dec34dsa*gya)
      dF_dRho(ipdRza,iGrid)=dF_dRho(ipdRza,iGrid)
     &                   +2*Coeff*(dec34ds*gz+dec34dsa*gza)
      dF_dRho(ipdRxb,iGrid)=dF_dRho(ipdRxb,iGrid)
     &                   +2*Coeff*(dec34ds*gx+dec34dsb*gxb)
      dF_dRho(ipdRyb,iGrid)=dF_dRho(ipdRyb,iGrid)
     &                   +2*Coeff*(dec34ds*gy+dec34dsb*gyb)
      dF_dRho(ipdRzb,iGrid)=dF_dRho(ipdRzb,iGrid)
     &                   +2*Coeff*(dec34ds*gz+dec34dsb*gzb)

200   Continue
      End Do
      Endif
*
      Return
      End

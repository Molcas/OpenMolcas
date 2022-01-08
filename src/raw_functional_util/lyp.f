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
      Subroutine LYP(mGrid,Coeff,iSpin,F_xc)
************************************************************************
*                                                                      *
* Object:  Lyp Functional(Formula taken from Molpro Manual)            *
*                                                                      *
*      Author: Per-AAke Malmquist,Department of Theoretical Chemistry  *
*              University of LUnd, SWEDEN                              *
*              D. Ajitha , Department of Theoretical Chemistry         *
*              University of Lund, SWEDEN                              *
*              Modify Per-AAke's code for open shell case              *
*              and adopt for closed shell case                         *
************************************************************************
      use KSDFT_Info, only: tmpB
      use nq_Grid, only: Rho, Sigma, l_casdft
      use nq_Grid, only: vRho, vSigma
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
*
      data Cfconst / 2.8712340001881918D0 /
      data aconst,bconst,cconst,dconst
     &        /0.04918d0,0.132d0,0.2533d0,0.349d0/
* Constants aconst,..,dconst are from ???
* Cfconst is = (3/10)*(3*Pi**2)**(2/3)
      Cfconst2=Cfconst*2.0D0**(11.0D0/3.D0)
      Rho_Min=T_X*1.0D-2
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
************************************************************************
*      write(6,*) 'iSpin in LYP routine', iSpin
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, mGrid
*
      rhoa=Rho(1,iGrid)
      rhob=rhoa
      rho_tot=rhoa+rhob
      if(rho_tot.lt.T_X) Go To 101

      rho3=rho_tot**(-1.0D0/3.D0)
      crho3=cconst*rho3
      R=1.D0/(1.D0+dconst*rho3)
      if(crho3.lt.709.0d0) then
        expcr=exp(-crho3)
      else
        expcr=0.0D0
      end if
      omega=expcr*R*rho_tot**(-11.0D0/3.D0)
      delta=crho3+1.D0-R
      dlogodr=(delta-11.D0)/(3.D0*rho_tot)
      ddeltadr=-(crho3+R*(1.D0-R))/(3.D0*rho_tot)
      p=aconst*bconst*omega

      ec1=-(4.D0*aconst*R)*rhoa*(rhob/rho_tot)


      dec1dra=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhob/rho_tot)**2
      dec1drb=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhoa/rho_tot)**2
      tmp1=-Cfconst2*p*rhoa**(11.0D0/3.D0)*rhob
      tmp2=-Cfconst2*p*rhob**(11.0D0/3.D0)*rhoa
      ec2=tmp1+tmp2
      dec2dra=(ec2*dlogodr)+(11.D0*tmp1+3.D0*tmp2)/(3.D0*rhoa)
      dec2drb=(ec2*dlogodr)+(3.D0*tmp1+11.D0*tmp2)/(3.D0*rhob)
      sa=Sigma(1,iGrid)
      sb=sa
      s =4.D0*sa
      pp=-p*rhoa*rhob/(18.D0*rho_tot)
      dppdra=((dlogodr-1.D0/rho_tot)+1.D0/rhoa)*pp
      dppdrb=((dlogodr-1.D0/rho_tot)+1.D0/rhob)*pp
      dqq1ds =47.D0*rho_tot
      dqq1dsa=22.D0*rhoa-45.D0*rho_tot
      dqq1dsb=22.D0*rhob-45.D0*rho_tot
      dqq2ds =-7.D0*rho_tot
      dqq2dsa= rhob-rhoa
      dqq2dsb= rhoa-rhob
      dqqds =dqq1ds +delta*dqq2ds
      dqqdsa=dqq1dsa+delta*dqq2dsa
      dqqdsb=dqq1dsb+delta*dqq2dsb
      dqq1dra=47.D0*s-23.D0*sa-45.D0*sb
      dqq1drb=47.D0*s-45.D0*sa-23.D0*sb
      qq1=dqq1dra*rhoa+dqq1drb*rhob
      dqq2dra=-7.D0*s -sa +sb
      dqq2drb=-7.D0*s +sa -sb
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

      dqdra=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhoa*sb
      dqdrb=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhob*sa
      q=(dqdra*rhoa+dqdrb*rhob)/2.D0
      ec4=p*q
      dec4dra=dlogodr*ec4+p*dqdra
      dec4drb=dlogodr*ec4+p*dqdrb
      dec4ds =p*(2.0D0/3.D0)*(rho_tot**2)
      dec4dsa=p*(rhob**2)-dec4ds
      dec4dsb=p*(rhoa**2)-dec4ds

      dec34ds =dec3ds +dec4ds
      dec34dsa=dec3dsa+dec4dsa
      dec34dsb=dec3dsb+dec4dsb

      Functional= ec1+ec2+ec3+ec4
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*                                                                      *
************************************************************************
*                                                                      *
*      dF/dRho
       vRho(1,iGrid)=vRho(1,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
*      dF/dGaa
       vSigma(1,iGrid)=vSigma(1,iGrid)
     &                     +Coeff*(dec34dsa+dec34ds)
     &                     +Coeff*(dec34ds)
*      dF/dGab
101   Continue
      End Do

*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=/= 1
      Else
*                                                                      *
************************************************************************
*                                                                      *
*         write(6,*) 'mGrid',mGrid
*         write(6,*) 'Rho_min',Rho_min
      If (l_casdft) Then
      Do iGrid = 1, mGrid
      rhoa=Max(Rho_min,Rho(1,iGrid))
      rhob=Max(Rho_min,Rho(2,iGrid))
      rho_tot=rhoa+rhob
*         if(iGrid.le.10) then
*           write(6,*) rhoa, rhob, rho_tot
*         end if
      if(rho_tot.lt.T_X) Cycle

      rho3=rho_tot**(-1.0D0/3.D0)
      crho3=cconst*rho3
      R=1.D0/(1.D0+dconst*rho3)
      if(crho3.lt.709.0d0) then
        expcr=exp(-crho3)
      else
        expcr=0.0D0
      end if
      omega=expcr*R*rho_tot**(-11.0D0/3.D0)
      delta=crho3+1.D0-R
      dlogodr=(delta-11.D0)/(3.D0*rho_tot)
      ddeltadr=-(crho3+R*(1.D0-R))/(3.D0*rho_tot)
      p=aconst*bconst*omega

      ec1=-(4.D0*aconst*R)*rhoa*(rhob/rho_tot)

      dec1dra=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhob/rho_tot)**2
      dec1drb=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhoa/rho_tot)**2

      tmp1=-Cfconst2*p*rhoa**(11.0D0/3.D0)*rhob
      tmp2=-Cfconst2*p*rhob**(11.0D0/3.D0)*rhoa
      ec2=tmp1+tmp2
      dec2dra=(ec2*dlogodr)+(11.D0*tmp1+3.D0*tmp2)
     &       /(3.D0*Max(rhoa,0.5D-50))
      dec2drb=(ec2*dlogodr)+(3.D0*tmp1+11.D0*tmp2)
     &       /(3.D0*Max(rhob,0.5D-50))

      sa=Sigma(1,iGrid)
      sb=Sigma(3,iGrid)
      s = sa + sb + Two*Sigma(2,iGrid)

      pp=-p*rhoa*rhob/(18.D0*rho_tot)
      dppdra=((dlogodr-1.D0/rho_tot)+1.D0/Max(rhoa,0.5D-50))*pp
      dppdrb=((dlogodr-1.D0/rho_tot)+1.D0/Max(rhob,0.5D-50))*pp
      dqq1ds =47.D0*rho_tot
      dqq1dsa=22.D0*rhoa-45.D0*rho_tot
      dqq1dsb=22.D0*rhob-45.D0*rho_tot
      dqq2ds =-7.D0*rho_tot
      dqq2dsa= rhob-rhoa
      dqq2dsb= rhoa-rhob
      dqqds =dqq1ds +delta*dqq2ds
      dqqdsa=dqq1dsa+delta*dqq2dsa
      dqqdsb=dqq1dsb+delta*dqq2dsb
      dqq1dra=47.D0*s-23.D0*sa-45.D0*sb
      dqq1drb=47.D0*s-45.D0*sa-23.D0*sb
      qq1=dqq1dra*rhoa+dqq1drb*rhob
      dqq2dra=-7.D0*s -sa +sb
      dqq2drb=-7.D0*s +sa -sb
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

      dqdra=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhoa*sb
      dqdrb=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhob*sa
      q=(dqdra*rhoa+dqdrb*rhob)/2.D0
      ec4=p*q
      dec4dra=dlogodr*ec4+p*dqdra
      dec4drb=dlogodr*ec4+p*dqdrb
      dec4ds =p*(2.0D0/3.D0)*(rho_tot**2)
      dec4dsa=p*(rhob**2)-dec4ds
      dec4dsb=p*(rhoa**2)-dec4ds

      dec34ds =dec3ds +dec4ds
      dec34dsa=dec3dsa+dec4dsa
      dec34dsb=dec3dsb+dec4dsb

      Functional= ec1+ec2+ec3+ec4
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
      tmpB(iGrid)=F_xc(iGrid)-tmpB(iGrid)
*                                                                      *
************************************************************************
*                                                                      *
*       dF/dRhoa, dF/dRhob
        vRho(1,iGrid)=vRho(1,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
        vRho(2,iGrid)=vRho(2,iGrid)
     &                   +Coeff*(dec1drb+dec2drb+dec3drb+dec4drb)
*       dF/dGaa, dF/dGab, and dF/dGbb
        vSigma(1,iGrid)=vSigma(1,iGrid)
     &                      +Coeff*(dec34dsa+dec34ds)
        vSigma(3,iGrid)=vSigma(3,iGrid)
     &                      +Coeff*(dec34dsb+dec34ds)
        vSigma(2,iGrid)=vSigma(2,iGrid)+2.0d0*Coeff*dec34ds
      End Do
      Else
      Do iGrid = 1, mGrid
      rhoa=Max(Rho_min,Rho(1,iGrid))
      rhob=Max(Rho_min,Rho(2,iGrid))
      rho_tot=rhoa+rhob
*         if(iGrid.le.10) then
*           write(6,*) rhoa, rhob, rho_tot
*         end if
      if(rho_tot.lt.T_X) Cycle

      rho3=rho_tot**(-1.0D0/3.D0)
      crho3=cconst*rho3
      R=1.D0/(1.D0+dconst*rho3)
      if(crho3.lt.709.0d0) then
        expcr=exp(-crho3)
      else
        expcr=0.0D0
      end if
      omega=expcr*R*rho_tot**(-11.0D0/3.D0)
      delta=crho3+1.D0-R
      dlogodr=(delta-11.D0)/(3.D0*rho_tot)
      ddeltadr=-(crho3+R*(1.D0-R))/(3.D0*rho_tot)
      p=aconst*bconst*omega

      ec1=-(4.D0*aconst*R)*rhoa*(rhob/rho_tot)

      dec1dra=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhob/rho_tot)**2
      dec1drb=(ec1*(1.D0-R)/(3.D0*rho_tot))
     &       -(4.D0*aconst*R)*(rhoa/rho_tot)**2

      tmp1=-Cfconst2*p*rhoa**(11.0D0/3.D0)*rhob
      tmp2=-Cfconst2*p*rhob**(11.0D0/3.D0)*rhoa
      ec2=tmp1+tmp2
      dec2dra=(ec2*dlogodr)+(11.D0*tmp1+3.D0*tmp2)
     &       /(3.D0*Max(rhoa,0.5D-50))
      dec2drb=(ec2*dlogodr)+(3.D0*tmp1+11.D0*tmp2)
     &       /(3.D0*Max(rhob,0.5D-50))

      sa=Sigma(1,iGrid)
      sb=Sigma(3,iGrid)
      s = sa + sb + Two*Sigma(2,iGrid)

      pp=-p*rhoa*rhob/(18.D0*rho_tot)
      dppdra=((dlogodr-1.D0/rho_tot)+1.D0/Max(rhoa,0.5D-50))*pp
      dppdrb=((dlogodr-1.D0/rho_tot)+1.D0/Max(rhob,0.5D-50))*pp
      dqq1ds =47.D0*rho_tot
      dqq1dsa=22.D0*rhoa-45.D0*rho_tot
      dqq1dsb=22.D0*rhob-45.D0*rho_tot
      dqq2ds =-7.D0*rho_tot
      dqq2dsa= rhob-rhoa
      dqq2dsb= rhoa-rhob
      dqqds =dqq1ds +delta*dqq2ds
      dqqdsa=dqq1dsa+delta*dqq2dsa
      dqqdsb=dqq1dsb+delta*dqq2dsb
      dqq1dra=47.D0*s-23.D0*sa-45.D0*sb
      dqq1drb=47.D0*s-45.D0*sa-23.D0*sb
      qq1=dqq1dra*rhoa+dqq1drb*rhob
      dqq2dra=-7.D0*s -sa +sb
      dqq2drb=-7.D0*s +sa -sb
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

      dqdra=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhoa*sb
      dqdrb=((4.0D0/3.D0)*rho_tot*(s-sa-sb))+2.D0*rhob*sa
      q=(dqdra*rhoa+dqdrb*rhob)/2.D0
      ec4=p*q
      dec4dra=dlogodr*ec4+p*dqdra
      dec4drb=dlogodr*ec4+p*dqdrb
      dec4ds =p*(2.0D0/3.D0)*(rho_tot**2)
      dec4dsa=p*(rhob**2)-dec4ds
      dec4dsb=p*(rhoa**2)-dec4ds

      dec34ds =dec3ds +dec4ds
      dec34dsa=dec3dsa+dec4dsa
      dec34dsb=dec3dsb+dec4dsb

      Functional= ec1+ec2+ec3+ec4
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
*                                                                      *
************************************************************************
*                                                                      *
*       dF/dRhoa, dF/dRhob
        vRho(1,iGrid)=vRho(1,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
        vRho(2,iGrid)=vRho(2,iGrid)
     &                   +Coeff*(dec1drb+dec2drb+dec3drb+dec4drb)
*       dF/dGaa, dF/dGab, and dF/dGbb
        vSigma(1,iGrid)=vSigma(1,iGrid)
     &                      +Coeff*(dec34dsa+dec34ds)
        vSigma(3,iGrid)=vSigma(3,iGrid)
     &                      +Coeff*(dec34dsb+dec34ds)
        vSigma(2,iGrid)=vSigma(2,iGrid)+2.0d0*Coeff*dec34ds
      End Do
      EndIf
      Endif
*
      Return
      End

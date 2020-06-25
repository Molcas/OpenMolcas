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
      Subroutine LYP(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &               Coeff,iSpin,F_xc,T_X)
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
#include "WrkSpc.fh"
#include "ksdft.fh"
      Real*8 dF_dRho(ndF_dRho,mGrid),Rho(nRho,mGrid),F_xc(mGrid)
cGLM     &               F_xca(mGrid),F_xcb(mGrid),tmpB(mGrid)
*
      dimension ec(0:4)
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
C     Call QEnter('LYP')
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
      rhoa=Rho(ipR,iGrid)
      rhob=rhoa
      rho_tot=rhoa+rhob
      if(rho_tot.lt.T_X) Go To 101

      gxa=Rho(ipdRx,iGrid)
      gya=Rho(ipdRy,iGrid)
      gza=Rho(ipdRz,iGrid)
      gxb=gxa
      gyb=gya
      gzb=gza

* New LYP code.
      gx=gxa+gxb
      gy=gya+gyb
      gz=gza+gzb

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
      sa=gxa**2+gya**2+gza**2
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
*      dF/dRho
       dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
*      dF/dGaa
       dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)
     &                     +Coeff*(dec34dsa+dec34ds)
*      dF/dGab
       dF_dRho(ipGxy,iGrid)=dF_dRho(ipGxy,iGrid)+2.0d0*Coeff*(dec34ds)
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
      tmpC_tot = 0.0d0
      Do iGrid = 1, mGrid
      rhoa=Max(Rho_min,Rho(ipRa,iGrid))
      rhob=Max(Rho_min,Rho(ipRb,iGrid))
      rho_tot=rhoa+rhob
*         if(iGrid.le.10) then
*           write(6,*) rhoa, rhob, rho_tot
*         end if
      if(rho_tot.lt.T_X) Go To 201
      gxa=Rho(ipdRxa,iGrid)
      gya=Rho(ipdRya,iGrid)
      gza=Rho(ipdRza,iGrid)
      gxb=Rho(ipdRxb,iGrid)
      gyb=Rho(ipdRyb,iGrid)
      gzb=Rho(ipdRzb,iGrid)

      gx=gxa+gxb
      gy=gya+gyb
      gz=gza+gzb

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

      sa=gxa**2+gya**2+gza**2
      sb=gxb**2+gyb**2+gzb**2
      s =gx**2+gy**2+gz**2
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

      ec(0)=ec1+ec2+ec3+ec4
      ec(1)=ec1
      ec(2)=ec2
      ec(3)=ec3
      ec(4)=ec4

      Functional= ec1+ec2+ec3+ec4
      F_xc(iGrid)=F_xc(iGrid)+Coeff*functional
      Work(ip_tmpB+iGrid-1)=F_xc(iGrid)-Work(ip_tmpB+iGrid-1)
*                                                                      *
************************************************************************
*                                                                      *
*       dF/dRhoa, dF/dRhob
        dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)
     &                   +Coeff*(dec1dra+dec2dra+dec3dra+dec4dra)
        dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)
     &                   +Coeff*(dec1drb+dec2drb+dec3drb+dec4drb)
*       dF/dGaa, dF/dGab, and dF/dGbb
        dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)
     &                      +Coeff*(dec34dsa+dec34ds)
        dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)
     &                      +Coeff*(dec34dsb+dec34ds)
        dF_dRho(ipGab,iGrid)=dF_dRho(ipGab,iGrid)+2.0d0*Coeff*dec34ds
201   Continue
      End Do
      Endif
*
C     Call QExit('LYP')
      Return
      End

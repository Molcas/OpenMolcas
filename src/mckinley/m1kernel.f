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
      subroutine m1kernel(Final,Hess,nHess,DAO,nDAO,
     &                    iAng,nRys,nZeta,
     &                    Alpha,Beta,Zeta,ZInv,
     &                    rKappa,P,TC,Coor,CoorAc,
     &                    Array,nArray,
     &                    ifgrd,indgrd,ifhss,indhss,
     &                    ifg,tr,nop,iuvwx,
     &                    kCnttp,fact,loper,idcar)

      use Real_Spherical
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
      External TNAI1, Fake, Cff2D
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"

      Integer iAng(4),nop(4),iuvwx(4)
      Real*8 Alpha(nZeta),Beta(nZeta),Zeta(nZeta),
     &        ZInv(nZeta),rKappa(nZeta),P(nZeta,*),
     &        TC(3),Coor(3,4),Array(nArray),Final(*),
     &        CoorAC(3,2), coori(3,4),Hess(*),DAO(nzeta,*)
      Logical ifg(4),tr(4),ifgrd(3,4),ifhss(3,4,3,4),eq,
     &         jfgrd(3,4),jfg(4),lGrad,lHess,jfhss(3,4,3,4)
      Integer indgrd(3,4,0:7),index(3,4),indhss(3,4,3,4,8),
     &        jndgrd(3,4,0:7),jndhss(3,4,3,4,8)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

      iprint = 0
      lGrad = idcar.ne.0
      lHess = nHess.ne.0
      call dcopy_(12,coor,1,coori,1)
      If (.Not.EQ(coor(1,1),coor(1,2)) .or.
     &    .Not.EQ(coor(1,1),coor(1,3))) Then
         Coori(1,1) = Coori(1,1)+One
      End If

      ip=1
      ipK = ip
      ip = ip + nZeta
      ipZ = ip
      ip = ip + nZeta
      ipZI = ip
      ip = ip + nZeta
      ipPx = ip
      ip = ip + nZeta
      ipPy = ip
      ip = ip + nZeta
      ipPz = ip
      ip = ip + nZeta
      ipDAO = ip
      ip = ip + nDAO*nZeta
      if (ip.ge.narray) then
        write(6,*) 'Out of memory in m1kernel (',narray,',',ip,')'
        Call QTrace
        Call Abend()
      endif



      if (iprint.ge.49)
     &  Write(6,*)'nM1=',dbsc(kCnttp)%nM1,'kCnttp=',kCnttp

      Do 1011 iM1xp=1, dbsc(kCnttp)%nM1
           Gamma =   dbsc(kCnttp)%M1xp(iM1xp)
           FactECP = dbsc(kCnttp)%M1cf(iM1xp)* Fact


           if (iprint.ge.49) then
              write(6,*) 'Fact=',FactECP,Fact
              write(6,*) 'im1xp=',iM1xp
              write(6,*) 'Gamma=',Gamma
           endif
*
*-----------Modify the original basis. Observe that
*           simplification due to A=B are not valid for the
*           exponent index, eq. P-A=/=0.
*
            Do 1012 iZeta = 1, nZeta
                PTC2 = (P(iZeta,1)-TC(1))**2
     &               + (P(iZeta,2)-TC(2))**2
     &               + (P(iZeta,3)-TC(3))**2
                Tmp0 = Zeta(iZeta)+Gamma
                Tmp1 = Exp(-Zeta(iZeta)*Gamma*PTC2/Tmp0)
                Array(ipK +iZeta-1) = rKappa(iZeta) * Tmp1
                Array(ipZ +iZeta-1) = Tmp0
                Array(ipZI+iZeta-1) = One/Tmp0
                Array(ipPx+iZeta-1) =
     &               (Zeta(iZeta)*P(iZeta,1)+Gamma*TC(1))/Tmp0
                Array(ipPy+iZeta-1) =
     &              (Zeta(iZeta)*P(iZeta,2)+Gamma*TC(2))/Tmp0
                Array(ipPz+iZeta-1) =
     &                  (Zeta(iZeta)*P(iZeta,3)+Gamma*TC(3))/Tmp0
 1012       Continue

            Do iDAO = 1, nDAO
              Do iZeta = 1, nZeta
                 Fac = FactECP * Array(ipK+iZeta-1) *
     &                        Array(ipZI+iZeta-1)*  Two * Pi
                 ipDAOt = nZeta*(iDAO-1) + iZeta-1 + ipDAO
                 Array(ipDAOt)= Fac * DAO(iZeta,iDAO)
              End Do
            End Do

*
*-----------Compute integrals with the Rys quadrature.
*
            call lcopy(4,ifg,1,jfg,1)

            call lcopy(12,ifgrd,1,jfgrd,1)
            call lcopy(12*12,ifhss,1,jfhss,1)

            call icopy(12*nirrep,indgrd,1,jndgrd,1)
            call icopy(12*12*nirrep,indhss,1,jndhss,1)

            Call Rysg2(iAng,nRys,nZeta,
     &                 Alpha,Beta,[One],[One],
     &                 Array(ipZ),Array(ipZI),nZeta,[One],[One],1,
     &                 Array(ipPx),nZeta,TC,1,Coori,Coor,CoorAC,
     &                 Array(ip),nArray-ip+1,
     &                 TNAI1,Fake,Cff2D,
     &                 Array(ipDAO),nDAO,Hess,nhess,jfGrd,jndGrd,
     &                 jfHss,jndHss,nOp,iuvwx,jfg,
     &                 nGr,Index,lgrad,lhess,tr)
          if (lGrad)  Then
            nb = nzeta*nElem(iAng(1))*nElem(iAng(2))
            Do iElem = 1, nElem(iAng(1))*nElem(iAng(2))*ngr
              Do iZeta = 1, nZeta
                tfac = Two*PI*Array(ipK+iZeta-1)*Array(ipZI-1+iZeta)
                indi=(iElem-1)*nZeta+iZeta
                Array(ip+indi-1) = FactECP*tfac * Array(ip+indi-1)
              End Do
            End Do


            Call SmAdNa(Array(ip),nb,Final,
     &            nop(1:3),loper,jndGrd,iuvwx(1:3),jfGrd,Index,
     &            idcar,One,jFG,tr)
      End if

 1011 Continue
      return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(ZInv)
      end

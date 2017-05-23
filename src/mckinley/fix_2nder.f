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
      Subroutine Fix_2nder(Fa1,Fa2,Fb1,Fb2,Final,nalpha,nbeta,
     &                     ishll,la,lb,iang,jfhess,nfun,fact)
*----------------------------------------------------------------------*
C
C       Fa1 includes  (  grad   < a | c > , < a | c >  )
C       Fb1 includes  (  grad   < c | b > , < c | b >  )
C
C
C
*
      Implicit Real*8(a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "real.fh"
CBS   Logical jfhess(3,3,2,2)   ! now same dimensions as in prjhss.f
      Logical jfhess(3,3,4,4)
      Real*8 FA1((2*iAng+1)*((la+1)*(la+2)/2)*nAlpha*nFun*4),
     &       FA2((2*iAng+1)*((la+1)*(la+2)/2)*nAlpha*nFun*6),
     &       FB1((2*iAng+1)*((lb+1)*(lb+2)/2)*nBeta*nFun*4),
     &       FB2((2*iAng+1)*((lb+1)*(lb+2)/2)*nBeta*nFun*6),
CBS  &Final(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6)
     &Final(nAlpha*nbeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,21)
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)


*              OK NOW EVERYTHING SHOULD BE MERGED TOGETHER
*
*
CBS   write(*,*) 'nFun ',nfun
CBS   write(*,*) 'FA1 ',FA1
CBS   write(*,*) 'FA2 ',FA2
CBS   write(*,*) 'FB1 ',FB1
CBS   write(*,*) 'FB2 ',FB2
      Do iCar = 1, 3
        Do jCar = 1, 3
         mvec=itri(iCar+3,jcar)
         If (jfHess(iCar,jcar,2,1)) Then
           ipFA1a = 1+ iCar *
     &        nAlpha*nFun*nElem(la)*(2*iAng+1)
           ipFB1a = 1 + jcar *
     &             nFun*nBeta*(2*iAng+1)*nElem(lb)

           Do  ib = 1, nElem(lb)
            Do ia = 1, nElem(la)
*
             Do  iC = 1, (2*iAng+1)
               iaC = (iC-1)*nElem(la) + ia
               ipaC = (iaC-1)*nAlpha*nFun + ipFA1a
               iCb = (ib-1)*(2*iAng+1) + iC
               ipCb = (iCb-1)*nFun*nBeta  + ipFB1a

               Call DGEMM_('N','N',
     &                    nAlpha,nBeta,nFun,
     &                    Fact,FA1(ipaC),nAlpha,
     &                        FB1(ipCb),nFun,
     &                    One,Final(1,ia,ib,mVec),nAlpha)
*
*
             End do
            End do
           End do
         End If
        End do
       End do

CBS   goto 4711
      Do iCar = 1, 3
        Do jCar = 1, icar
         mvec=itri(iCar,jcar)
         If (jfHess(iCar,jcar,1,1)) Then
           ipFA1a = 1 + (itri(iCar,jCar)-1) *
     &        nAlpha*nFun*nElem(la)*(2*iAng+1)
*
           Do  ib = 1, nElem(lb)
            Do ia = 1, nElem(la)

             Do  iC = 1, (2*iAng+1)
               iaC = (iC-1)*nElem(la) + ia
               ipaC = (iaC-1)*nAlpha*nFun + ipFA1a
               iCb = (ib-1)*(2*iAng+1) + iC
               ipCb = (iCb-1)*nFun*nBeta  + 1

               Call DGEMM_('N','N',
     &                    nAlpha,nBeta,nFun,
     &                    Fact,FA2(ipaC),nAlpha,
     &                    FB1(ipCb),nFun,
     &                    One,Final(1,ia,ib,mVec),nAlpha)
*
*
             End do
            End do
           End do
         End If
        End do
       End do
*4711 continue

CBS   goto 4712      ! for testing the different contributions
      Do iCar = 1, 3
        Do jCar = 1, icar
         mvec=itri(3+icar,3+jcar)
         If (jfHess(iCar,jcar,2,2)) Then
           ipFB1a = 1 + (itri(icar,jcar)-1) *
CBS  &        nAlpha*nFun*nElem(la)*(2*iAng+1)
CBS           looks much better ....  (Oh Anders ...)

     &        nBeta*nFun*nElem(lb)*(2*iAng+1)

           Do  ib = 1, nElem(lb)
            Do ia = 1, nElem(la)
*
             Do  iC = 1, (2*iAng+1)
               iaC = (iC-1)*nElem(la) + ia
               ipaC = (iaC-1)*nAlpha*nFun + 1
               iCb = (ib-1)*(2*iAng+1) + iC
               ipCb = (iCb-1)*nFun*nBeta  + ipFB1a

                        Call DGEMM_('N','N',
     &                    nAlpha,nBeta,nFun,
     &                    Fact,FA1(ipaC),nAlpha,
     &                        FB2(ipCb),nFun,
     &                    One,Final(1,ia,ib,mVec),nAlpha)

*
             End do
            End do
           End do
         End If
        End do
       End do
*4712  continue

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ishll)
      End

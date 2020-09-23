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
      Subroutine coreB(iang,lb,ishll,nordop,TC,RB,Array,narr,
     &                 Beta,nbeta,fb1,fb2,jfgrad,jfhess,
     &                 ld,debug)
*
*
*  Calculates <core|B'> and <core|B">
*
* @parameter iang  Angular momenta for core
* @parameter lb    Angular momenta for ket
* @parameter ishll identification for core shell
* @parameter nordop order for operator
* @parameter TC   Cartesian coordinates for core
* @parameter  RB Cartesian coordinates for ket
* @parameter Array Scratch
* @parameter narr size for scratch
* @parameter Beta  Ket exponents
* @parameter nbeta  number of exponents
* @parameter FB1 First derivatives (out)
* @parameter FB2 2nd derivatives (out)
* @parameter jfgrad true for all 1-deriavtives that are needed
* @parameter jfhess true for all 2-deriavtives that are needed
* @parameter ld Order of derivatives
* @parameter debug guess
*
      Use Basis_Info
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
      Logical ABeq(3),jfgrad(3),jfhess(4,3,4,3),debug
      Real*8 TC(3),RB(3),Array(*),fb1(*),fb2(*),beta(*)

      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      nExpi=Shells(iShll)%nExp
      if (debug) then
        Write(6,*) 'Shell: ',ishll,' nBeta:',nbeta,' nExp:',nExpi,
     &             'Angular',lb,iang
      Endif
*
*
      ip=1
      ipP2 = ip
      ip = ip + 3 * nExpi*nBeta
      ipZ2 = ip
      ip = ip + nExpi*nBeta
      ipK2 = ip
      ip = ip + nExpi*nBeta
      ipZI2 = ip
      ip = ip + nExpi*nBeta
      If (ip-1.gt.nArr) Then
         Write (6,*) '  ip-1.gt.nArr*nZeta(2) in bcore (',ip,',',
     &    narr,')'
         Call Abend
      End If
*
*--------------Calculate Effective center and exponent for <core|B>
*
      Call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,
     &                   Shells(iShll)%Exp,Beta)
      Call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,
     &                    TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
*--------------Calculate Overlap <core|B> and <core|B'>
*
CBS   nHer = (iAng+(lb+1)+2)/2
CBS   now it looks comparable to acore.f
      nHer = (lb+1+iAng+1+ld)/2
      ipCxyz = ip
      ip = ip + nBeta*nExpi*3*nHer*(iAng+1)
      ipBxyz = ip
      ip = ip + nBeta*nExpi*3*nHer*(lb+1+ld)
      ipRxyz = ip
      ip = ip + nBeta*nExpi*3*nHer*(nOrdOp+1)
      ipQ1 = ip
      ip = ip +
     &      nBeta*nExpi*3*(iAng+1)*(lb+1+ld)*(nOrdOp+1)
      ipB = ip
      ip = ip + nBeta*nExpi
      If (ip-1.gt.nArr) Then
         Write (6,*) '  ip-1.gt.nArr*nZeta(2b) in PrjGrd'
         Call Abend
      End If
      ABeq(1) = TC(1).eq.RB(1)
      ABeq(2) = TC(2).eq.RB(2)
      ABeq(3) = TC(3).eq.RB(3)
      if (debug)
     &write(6,*) 'shll=',ishll,' nExp=',nExpi,
     &            ' nBeta=',nBeta
      Call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,
     &            TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      Call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,
     &            RB,Array(ipBxyz),lb+ld,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call CrtCmp(Array(ipZ2),Array(ipP2),nExpi*nBeta,
     &            TC,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      If (debug) Then
            Write (6,*) ' nbeta  = ',nbeta ,' nExp(',ishll,')=',
     &              nExpi,' nHer=',nHer,' lb=',lb,' iAng=',
     &              iAng,' nOrdOp=',nOrdOp

                  Write (6,*) ' Array(ipCxyz)=',
     &             DNrm2_(nBeta*nExpi*3*nHer*(iAng+1),
     &             Array(ipCxyz),1)
                  Write (6,*) ' Array(ipBxyz)=',
     &             DNrm2_(nBeta*nExpi*3*nHer*(lb+2),
     &             Array(ipBxyz),1)
                  Write (6,*) ' Array(ipRxyz)=',
     &             DNrm2_(nBeta*nExpi*3*nHer*(nOrdOp+1),
     &             Array(ipRxyz),1)
      End If

      Call Assmbl(Array(ipQ1),
     &            Array(ipCxyz),iAng,
     &            Array(ipRxyz),nOrdOp,
     &            Array(ipBxyz),lb+ld,
     &            nExpi*nBeta,HerW(iHerW(nHer)),nHer)
      iStrt = ipB
      Do  iGamma = 1, nExpi
         call dcopy_(nBeta,Beta,1,Array(iStrt),nExpi)
         iStrt = iStrt + 1
      Enddo
      If (debug) Then
                  Write (6,*) ' Array(ipB)=',
     &            DNrm2_(nExpi*nBeta,Array(ipB),1)
      End If

      Call rKappa_Zeta(Array(ipK2),Array(ipZ2),
     &                 nExpi*nBeta)
      Call CmbnCB(Array(ipQ1),nExpi*nBeta,iAng,lb,
     &            Array(ipK2),FB1,
     &            Array(ipB),JfGrad,ld,nVecCB)
      If (debug) Then
         Write (6,*) ' Array(ipQ1)=',
     &        DNrm2_(
     &        nExpi*nBeta*3*(lb+1+ld+2)*(iAng+1)*(nOrdOp+1),
     &        Array(ipQ1),1)
         Write (6,*) ' Array(ipB)=',
     &            DNrm2_(nExpi*nBeta,Array(ipB),1)
      End If
      If (ld.ge.2) Then
        Call CmbnS2b(Array(ipQ1),nBeta *nExpi,iang,lb,
     &            Array(ipK2),FB2,
     &            Array(ipB),jfHess,ld)
        If (debug) Then
           Do i=1,6
              ipV=1
              n=nBeta*nExpi*nElem(lb)*nElem(iAng)
              Write(6,*)n,nBeta,nExpi,nElem(lb),nElem(iAng)

              Write(6,*) 'CmbnB2(',n,')=',DNrm2_(n,FB2(ipV),1)
              ipV=ipV+n
           End do
        End If
      End If

      Return
      End

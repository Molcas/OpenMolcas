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
      Subroutine Acore(iang,la,ishll,nordop,TC,A,Array,narr,
     &                 Alpha,nalpha,fa1,fa2,jfgrad,jfhess,
     &                 ld,debug)
*
*
*  Calculates <A'|core> and <A"|core>
*
* @parameter iang  Angular momenta for core
* @parameter la    Angular momenta for bra
* @parameter ishll identification for core shell
* @parameter nordop order for operator
* @parameter TC   Cartesian coordinates for core
* @parameter  A Cartesian coordinates for bra
* @parameter Array Scratch
* @parameter narr size for scratch
* @parameter Alpha Bra exponents
* @parameter nalpha number of exponents
* @parameter FA1 First derivatives (out)
* @parameter FA2 2nd derivatives (out)
* @parameter jfgrad true for all 1-deriavtives that are needed
* @parameter jfhess true for all 2-deriavtives that are needed
* @parameter ld Order of derivatives
* @parameter debug guess
      Use Basis_Info
      use Her_RW
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"
      Logical ABeq(3),jfgrad(3),jfhess(4,3,4,3),
     &        debug
      Real*8 TC(3),A(3),Array(*),fa1(*),fa2(*),alpha(*)
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
      nExpi=Shells(iShll)%nExp
      ip = 1
      ipP1 = ip
      ip = ip + 3 * nAlpha*nExpi
      ipZ1 = ip
      ip = ip + nAlpha*nExpi
      ipK1 = ip
      ip = ip + nAlpha*nExpi
      ipZI1 = ip
      ip = ip + nAlpha*nExpi
      If (ip-1.gt.nArr) Then
         Write (6,*) ' ip-1.gt.nArr in acore  (',ip,
     &   ',',narr,')'
         Call Abend
      End If

C------Calculate Effective center and exponent for <A|core>

      Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,
     &          Alpha,Shells(iShll)%Exp)
      Call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,
     &            A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
*--------------Calculate Overlap <A|core> and derivative <A'|core>
*
      nHer = (la+1+iAng+1+ld)/2
      ipAxyz = ip
      ip = ip + nAlpha*nExpi*3*nHer*(la+1+ld)
      ipCxyz = ip
      ip = ip + nAlpha*nExpi*3*nHer*(iAng+1)
      ipRxyz = ip
      ip = ip + nAlpha*nExpi*3*nHer*(nOrdOp+1)
      ipQ1 = ip
      ip = ip +
     &      nAlpha*nExpi*3*(la+1+ld)*(iAng+1)*(nOrdOp+1)
      ipA = ip
      ip = ip + nAlpha*nExpi
      If (ip-1.gt.nArr) Then
         Write (6,*) '  ip-1.gt.nArr (1b) in acore (',
     &    ip,',',narr,')','Order',ld,Shells(ishll)%nExp,nalpha
         Call Abend
      End If
      ABeq(1) = A(1).eq.TC(1)
      ABeq(2) = A(2).eq.TC(2)
      ABeq(3) = A(3).eq.TC(3)
      Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,
     &            A,Array(ipAxyz),la+ld,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,
     &            TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      ABeq(1) = .False.
      ABeq(2) = .False.
      ABeq(3) = .False.
      Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*nExpi,
     &            A,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &            nHer,ABeq)
      If (debug) Then
        Write (6,*) ' nAlpha = ',nAlpha,' nExp(',ishll,')=',
     &              nExpi,' nHer=',nHer,' la=',la,' iAng=',
     &              iAng,' nOrdOp=',nOrdOp

        Write (6,*) ' Array(ipAxyz)=',
     &             DNrm2_(nAlpha*nExpi*3*nHer*(la+ld+1),
     &             Array(ipAxyz),1)
        Write (6,*) ' Array(ipCxyz)=',
     &             DNrm2_(nAlpha*nExpi*3*nHer*(iAng+1),
     &             Array(ipCxyz),1)
        Write (6,*) ' Array(ipRxyz)=',
     &             DNrm2_(nAlpha*nExpi*3*nHer*(nOrdOp+1),
     &             Array(ipRxyz),1)
      End If

      Call Assmbl(Array(ipQ1),
     &            Array(ipAxyz),la+ld,
     &            Array(ipRxyz),nOrdOp,
     &            Array(ipCxyz),iAng,
     &            nAlpha*nExpi,HerW(iHerW(nHer)),nHer)
      iStrt = ipA
      Do 20 iGamma = 1, nExpi
         call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
         iStrt = iStrt + nAlpha
 20   Continue
      If (debug) Then
                  Write (6,*) ' Array(ipA)=',
     &            DNrm2_(nAlpha*nExpi,Array(ipA),1)
      End If

      Call rKappa_Zeta(Array(ipK1),Array(ipZ1),nExpi*nAlpha)
      Call CmbnAC(Array(ipQ1),nAlpha*nExpi,la,iAng,
     &            Array(ipK1),FA1,
     &            Array(ipA),JfGrad,ld,nVecAC)
      If (debug) Then
      write(6,*) 'nVecAC',nvecac
                  Write (6,*) ' Array(ipQ1)=',
     &            DNrm2_(
     &            nAlpha*nExpi*3*(la+ld+1)*(iAng+1)*(nOrdOp+1),
     &            Array(ipQ1),1)
                  Write (6,*) ' Array(ipA)=',
     &            DNrm2_(nAlpha*nExpi,Array(ipA),1)
      Do i=1,nvecac
        ipV=1
        n=nAlpha*nExpi*nElem(la)*nElem(iAng)
        Write(6,*) 'Cmbn(',i,')=',DNrm2_(n,FA1(ipV),1)
        ipV=ipV+n
      End do
      End If

      If (ld.ge.2) Then
        Call CmbnS2a(Array(ipQ1),nAlpha*nExpi,la,iAng,
     &              Array(ipK1),FA2,
     &              Array(ipA),jfHess,ld)
        If (debug) Then
          Do i=1,6
            ipV=1
            n=nAlpha*nExpi*nElem(la)*nElem(iAng)
            Write(6,*) 'Cmbn2(',i,')=',DNrm2_(n,FA2(ipV),1)
            ipV=ipV+n
          End do
        End If
      End If


      Return
      End

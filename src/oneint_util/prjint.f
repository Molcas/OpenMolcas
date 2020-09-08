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
* Copyright (C) 1993, Roland Lindh                                     *
*               1993, Per Boussard                                     *
************************************************************************
      SubRoutine PrjInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of ECP integrals.         *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              ZXia                                                    *
*              SetUp1                                                  *
*              MltPrm                                                  *
*              DGeTMO  (ESSL)                                          *
*              DGEMM_  (ESSL)                                          *
*              DScal   (ESSL)                                          *
*              DGEMM_  (ESSL)                                          *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*      Alpha : exponents of bra gaussians                              *
*      nAlpha: number of primitives (exponents) of bra gaussians       *
*      Beta  : as Alpha but for ket gaussians                          *
*      nBeta : as nAlpha but for the ket gaussians                     *
*      Zeta  : sum of exponents (nAlpha x nBeta)                       *
*      ZInv  : inverse of Zeta                                         *
*      rKappa: gaussian prefactor for the products of bra and ket      *
*              gaussians.                                              *
*      P     : center of new gaussian from the products of bra and ket *
*              gaussians.                                              *
*      Final : array for computed integrals                            *
*      nZeta : nAlpha x nBeta                                          *
*      nComp : number of components in the operator (e.g. dipolmoment  *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by ECPMem                 *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, Sweden, and Per Boussard, Dept. of Theoretical  *
*             Physics, University of Stockholm, Sweden, October '93.   *
************************************************************************
      use Basis_Info
      use Center_Info
      use Real_Spherical
      use Symmetry_Info, only: iChTbl
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3), C(3), TC(3)
      Integer iStabM(0:nStabM-1), lOper(nComp), iDCRT(0:7),
     &          iChO(nComp), iTwoj(0:7)
!#define _DEBUG_
#ifdef _DEBUG_
      Character*80 Label
#endif
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
#ifdef _DEBUG_
      Call RecPrt(' In PrjInt: Zeta',' ',Zeta,1,nZeta)
      Call RecPrt(' In PrjInt: A',' ',A,1,3)
      Call RecPrt(' In PrjInt: RB',' ',RB,1,3)
      Call RecPrt(' In PrjInt: Ccoor',' ',Ccoor,1,3)
      Call RecPrt(' In PrjInt: P',' ',P,nZeta,3)
      Write (6,*) ' In PrjInt: la,lb=',' ',la,lb
#endif
*
*     call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
      Final(:,:,:,:)=Zero
*
      llOper = lOper(1)
      iComp = 1
      mdc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.dbsc(iCnttp)%ECP) Then
            mdc = mdc + dbsc(iCnttp)%nCntr
            Cycle
         End If
         Do iCnt = 1,dbsc(iCnttp)%nCntr
            C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               dc(mdc+iCnt)%iStab,dc(mdc+iCnt)%nStab,iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
         Do lDCRT = 0, nDCRT-1
            Call OA(iDCRT(lDCRT),C,TC)
            Do iAng = 0, dbsc(iCnttp)%nPrj-1
               iShll = dbsc(iCnttp)%iPrj + iAng
               nExpi=Shells(iShll)%nExp
               nBasisi=Shells(iShll)%nBasis
               If (nExpi.eq.0 .or. nBasisi.eq.0) Cycle
*
#ifdef _DEBUG_
               Call RecPrt('Cff',' ',Shells(iShll)%pCff,nExpi,
     &                     nBasisi)
#endif
               ip = 1
               ipF1 = ip
               nac = nElem(la)*nElem(iAng)
               ip = ip + nAlpha*nExpi*nac
               ipP1 = ip
               ip = ip + 3 * nAlpha*nExpi
               ipZ1 = ip
               ip = ip + nAlpha*nExpi
               ipK1 = ip
               ip = ip + nAlpha*nExpi
               ipZI1 = ip
               ip = ip + nAlpha*nExpi
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'PrjInt: ip-1.gt.nArr*nZeta(1)')
                  Call Abend()
               End If
               mArr = (nArr*nZeta-(ip-1))/nZeta
*
*--------------Calculate Effective center and exponent for <A|core>
*
               Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExpi,
     &                   Alpha,Shells(iShll)%Exp)
               Call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,
     &                     A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
*--------------Calculate Overlap <A|core>
*
               nHer = (la+iAng+2)/2
               Call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,nExpi,
     &                   Array(ipZ1),Array(ipZI1),
     &                   Array(ipK1),Array(ipP1),
     &                   Array(ipF1),nAlpha*nExpi,iComp,
     &                   la,iAng,A,TC,nHer,Array(ip),
     &                   mArr,CCoor,nOrdOp)
               ip = ip - 6 * nAlpha*nExpi
*
               ipF2 = ip
               ncb = nElem(iAng)*nElem(lb)
               ip = ip + nExpi*nBeta*ncb
               ipP2 = ip
               ip = ip + 3 * nExpi*nBeta
               ipZ2 = ip
               ip = ip + nExpi*nBeta
               ipK2 = ip
               ip = ip + nExpi*nBeta
               ipZI2 = ip
               ip = ip + nExpi*nBeta
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'PrjInt: ip-1.gt.nArr*nZeta(2)')
                  Call Abend()
               End If
               mArr = (nArr*nZeta-(ip-1))/nZeta
*
*--------------Calculate Effective center and exponent for <core|B>
*
               Call ZXia(Array(ipZ2),Array(ipZI2),nExpi,nBeta,
     &                   Shells(iShll)%Exp,Beta)
               Call SetUp1(Shells(iShll)%Exp,nExpi,Beta,nBeta,
     &                    TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
*--------------Calculate Overlap <core|B>
*
               nHer = (iAng+lb+2)/2
               Call MltPrm(Shells(iShll)%Exp,nExpi,Beta,nBeta,
     &                   Array(ipZ2),Array(ipZI2),
     &                   Array(ipK2),Array(ipP2),
     &                   Array(ipF2),nExpi*nBeta,iComp,
     &                   iAng,lb,TC,RB,nHer,Array(ip),
     &                   mArr,CCoor,nOrdOp)
               ip = ip - 6 * nExpi*nBeta
               ipTmp = ip
               ip = ip + Max(nAlpha*nExpi*nac,
     &                       nBeta*ncb*nBasisi)
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'PrjInt: ip-1.gt.nArr*nZeta(3)')
                  Call Abend()
               End If
*
*--------------Calculate Contraction over components of the core
*              orbitals of type <A|core>Bc<core|B> where we now have in
*              Array(ipF1) the cartesian components of <A|core>, and
*              similarily, in Array(ipF2), we have stored the cartesian
*              components of <core|B>. Observe that the core orbitals
*              orthonomal atomic orbitals. Hence, the transformation
*              to the spherical harmonics has to be for normilized
*              spherical harminics.
*
*--------------From the lefthandside overlap, form iKaC from ikac by
*              1) i,kac -> k,aci
*
               Call DgeTMo(Array(ipF1),nAlpha,nAlpha,
     &                     nExpi*nac,Array(ipTmp),nExpi*nac)
*
*--------------2) aciK =  k,aci * k,K
*
               Call DGEMM_('T','N',
     &                     nAlpha*nac,nBasisi,nExpi,
     &                     1.0d0,Array(ipTmp),nExpi,
     &                     Shells(iShll)%pCff,nExpi,
     &                     0.0d0,Array(ipF1),nAlpha*nac)
*
*--------------3) Mult by shiftoperators aci,K -> Bk(K) * aci,K
*
               Do iBk = 1, nBasisi
                  Bk = Shells(ishll)%Bk(iBk)
                  Call DScal_(nAlpha*nac,Bk,
     &                       Array(ipF1+(iBk-1)*nAlpha*nac),1)
               End Do ! iBk
*
*--------------4) a,ciK -> ciKa
*
               Call DgeTMo(Array(ipF1),nElem(la),nElem(la),
     &                     nElem(iAng)*nAlpha*nBasisi,
     &                     Array(ipTmp),
     &                     nElem(iAng)*nAlpha*nBasisi)
*
*--------------5) iKa,C = c,iKa * c,C
*
               Call DGEMM_('T','N',
     &                     nAlpha*nBasisi*nElem(la),
     &                     (2*iAng+1),nElem(iAng),
     &                     1.0d0,Array(ipTmp),nElem(iAng),
     &                     RSph(ipSph(iAng)),nElem(iAng),
     &                     0.0d0,Array(ipF1),
     &                     nAlpha*nBasisi*nElem(la))
*
*--------------And (almost) the same thing for the righthand side, form
*              KjCb from kjcb
*              1) jcb,K = k,jcb * k,K
*
               Call DGEMM_('T','N',
     &                     nBeta*ncb,nBasisi,nExpi,
     &                     1.0d0,Array(ipF2),nExpi,
     &                     Shells(iShll)%pCff,nExpi,
     &                     0.0d0,Array(ipTmp),nBeta*ncb)
*
*--------------2)  j,cbK -> cbK,j
*
               Call DgeTMo(Array(ipTmp),nBeta,nBeta,
     &                     nBasisi*ncb,Array(ipF2),
     &                     nBasisi*ncb)
*
*--------------3) bKj,C = c,bKj * c,C
*
               Call DGEMM_('T','N',
     &                     nElem(lb)*nBasisi*nBeta,
     &                     (2*iAng+1),nElem(iAng),
     &                     1.0d0,Array(ipF2),nElem(iAng),
     &                     RSph(ipSph(iAng)),nElem(iAng),
     &                     0.0d0,Array(ipTmp),
     &                     nElem(lb)*nBasisi*nBeta)
*
*--------------4) b,KjC -> KjC,b
*
               Call DgeTMo(Array(ipTmp),nElem(lb),nElem(lb),
     &                     nBasisi*nBeta*(2*iAng+1),Array(ipF2),
     &                     nBasisi*nBeta*(2*iAng+1))
*
*--------------Next Contract (iKaC)*(KjCb) over K and C, producing ijab,
*              by the following procedure:
*              Loop over a and b
*                Loop over C
*                  Contract iK(aC)*Kj(Cb), over K producing ij(aCb),
*                    accumulate to ij(ab)
*                End loop C
*              End Loop b and a
*
               Do ib = 1, nElem(lb)
                  Do ia = 1, nElem(la)
*
                     Do iC = 1, (2*iAng+1)
                        iaC = (iC-1)*nElem(la) + ia
                        ipaC = (iaC-1)*nAlpha*nBasisi + ipF1
                        iCb = (ib-1)*(2*iAng+1) + iC
                        ipCb = (iCb-1)*nBasisi*nBeta  + ipF2
*
                        iIC = 0
                        Do iIrrep = 0, nIrrep-1
                           If (iAnd(llOper,iTwoj(iIrrep)).eq.0) Cycle
                           iIC = iIC + 1
                           nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
                           Xg=DBLE(iChTbl(iIrrep,nOp         ))
                           Factor=Xg*Fact
                           Call DGEMM_('N','N',
     &                                nAlpha,nBeta,nBasisi,
     &                                Factor,Array(ipaC),nAlpha,
     &                                    Array(ipCb),nBasisi,
     &                                One,Final(1,ia,ib,iIC),nAlpha)
                        End Do ! iIrrep
*
                     End Do ! iC
                  End Do    ! ia
               End Do       ! ib
            End Do ! iAng
*
         End Do ! lDCRT
         End Do ! iCnt
         mdc = mdc + dbsc(iCnttp)%nCntr
      End Do ! iCnttp
*
#ifdef _DEBUG_
         Write (6,*) ' Result in PrjInt'
         Do 100 ia = 1, (la+1)*(la+2)/2
            Do 200 ib = 1, (lb+1)*(lb+2)/2
               Write (Label,'(A,I2,A,I2,A)')
     &               ' Final(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,1),nAlpha,nBeta)
 200        Continue
 100     Continue
#endif
*
*     Call QExit('PrjInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P)
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_integer(nRys)
         Call Unused_integer_array(iChO)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End

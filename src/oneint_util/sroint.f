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
* Copyright (C) 1994, Roland Lindh                                     *
*               1994, Luis Seijo                                       *
************************************************************************
      SubRoutine SROInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of MP integrals.          *
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
*             of Lund, Sweden, and Luis Seijo, Dept. of Applied Phys-  *
*             ical Chemistry, the Free University of Madrid, Spain,    *
*             September '94.                                           *
************************************************************************
      use Basis_Info
      use Real_Spherical
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3), C(3), TC(3)
      Integer iStabM(0:nStabM-1), lOper(nComp), iDCRT(0:7),
     &          iChO(nComp), iTwoj(0:7)
      Character*80 Label
      Logical EQ
      Data iTwoj/1,2,4,8,16,32,64,128/
*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
*
*     Call qEnter('SROInt')
      iRout = 191
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In SROInt: A',' ',A,1,3)
         Call RecPrt(' In SROInt: RB',' ',RB,1,3)
         Call RecPrt(' In SROInt: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In SROInt: P',' ',P,nZeta,3)
         Write (6,*) ' In SROInt: la,lb=',' ',la,lb
      End If
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
      llOper = lOper(1)
      iComp = 1
      mdc = 0
      Do iCnttp = 1, nCnttp
         If (.Not.ECP(iCnttp) .or.
    &        nSRO_Shells(iCnttp).le.0) Then
            mdc = mdc + dbsc(iCnttp)%nCntr
            Cycle
         End If
         Do iCnt = 1,dbsc(iCnttp)%nCntr
            C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
*
            Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &               jStab(0,mdc+iCnt),nStab(mdc+iCnt),iDCRT,nDCRT)
            Fact = DBLE(nStabM) / DBLE(LmbdT)
*
         Do lDCRT = 0, nDCRT-1
            TC(1) = DBLE(iPhase(1,iDCRT(lDCRT)))*C(1)
            TC(2) = DBLE(iPhase(2,iDCRT(lDCRT)))*C(2)
            TC(3) = DBLE(iPhase(3,iDCRT(lDCRT)))*C(3)
            Do iAng = 0, nSRO_Shells(iCnttp)-1
               iShll = ipSRO(iCnttp) + iAng
               If (nExp(iShll).eq.0) Cycle
*
*
               ip = 1
               ipC = ip
               ip = ip + nExp(iShll)**2
*
               If (iPrint.ge.49)
     &            Call RecPrt(' The Akl matrix',' ',
     &                        Shells(iShll)%Akl(1,1,1),
     &                        nExp(iShll),nExp(iShll))
               call dcopy_(nExp(iShll)**2,Shells(iShll)%Akl(1,1,1),1,
     &                                    Array(ipC),1)
               If (EQ(A,RB).and.EQ(A,TC).and.
     &            lNoPair.and.NoPairL(iCnttp)) Then
                  If (iPrint.ge.49)
     &               Call RecPrt(' The Adl matrix',' ',
     &                           Shells(iShll)%Akl(1,1,2),
     &                           nExp(iShll),nExp(iShll))
                  Call DaXpY_(nExp(iShll)**2,One,
     &                                       Shells(iShll)%Akl(1,1,2),1,
     &                                       Array(ipC),1)
               End If
*
               ipF1 = ip
               nac = nElem(la)*nElem(iAng)
               ip = ip + nAlpha*nExp(iShll)*nac
               ipP1 = ip
               ip = ip + 3 * nAlpha*nExp(iShll)
               ipZ1 = ip
               ip = ip + nAlpha*nExp(iShll)
               ipK1 = ip
               ip = ip + nAlpha*nExp(iShll)
               ipZI1 = ip
               ip = ip + nAlpha*nExp(iShll)
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'SROInt: ip-1.gt.nArr*nZeta(1)')
                  Write (6,*) ' nArr, nZeta=',nArr, nZeta
                  Write (6,*) ' nac, nAlpha=', nac, nAlpha
                  Write (6,*) ' nExp(iShll)=',nExp(iShll)
                  Call Abend()
               End If
               mArr = (nArr*nZeta-(ip-1))/nZeta
*
*--------------Calculate Effective center and exponent for <A|alm>
*
               Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,nExp(iShll),
     &                   Alpha,Shells(iShll)%Exp)
               Call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,nExp(iShll),
     &                     A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
*--------------Calculate Overlap <A|alm>
*
               nHer = (la+iAng+2)/2
               Call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,nExp(iShll),
     &                   Array(ipZ1),Array(ipZI1),
     &                   Array(ipK1),Array(ipP1),
     &                   Array(ipF1),nAlpha*nExp(iShll),iComp,
     &                   la,iAng,A,TC,nHer,Array(ip),
     &                   mArr,CCoor,nOrdOp)
               If (iPrint.ge.99) Call RecPrt('<a|srbs>',' ',
     &                   Array(ipF1),nAlpha*nExp(iShll),nac)
               ip = ip - 6 * nAlpha*nExp(iShll)
*
               ipF2 = ip
               ncb = nElem(iAng)*nElem(lb)
               ip = ip + nExp(iShll)*nBeta*ncb
               ipP2 = ip
               ip = ip + 3 * nExp(iShll)*nBeta
               ipZ2 = ip
               ip = ip + nExp(iShll)*nBeta
               ipK2 = ip
               ip = ip + nExp(iShll)*nBeta
               ipZI2 = ip
               ip = ip + nExp(iShll)*nBeta
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'SROInt: ip-1.gt.nArr*nZeta(2)')
                  Call Abend()
               End If
               mArr = (nArr*nZeta-(ip-1))/nZeta
*
*--------------Calculate Effective center and exponent for <blm|B>
*
               Call ZXia(Array(ipZ2),Array(ipZI2),nExp(iShll),nBeta,
     &                   Shells(iShll)%Exp,Beta)
               Call SetUp1(Shells(iShll)%Exp,nExp(iShll),Beta,nBeta,
     &                    TC,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
*--------------Calculate Overlap <blm|B>
*
               nHer = (iAng+lb+2)/2
               Call MltPrm(Shells(iShll)%Exp,nExp(iShll),Beta,nBeta,
     &                   Array(ipZ2),Array(ipZI2),
     &                   Array(ipK2),Array(ipP2),
     &                   Array(ipF2),nExp(iShll)*nBeta,iComp,
     &                   iAng,lb,TC,RB,nHer,Array(ip),
     &                   mArr,CCoor,nOrdOp)
               If (iPrint.ge.99) Call RecPrt('<srbs|b>',' ',
     &                   Array(ipF2),nExp(iShll)*nBeta,ncb)
               ip = ip - 6 * nExp(iShll)*nBeta
               ipTmp = ip
               ip = ip + Max(nAlpha*nExp(iShll)*nac,
     &                       nExp(iShll)*nBeta*ncb)
               If (ip-1.gt.nArr*nZeta) Then
                  Call WarningMessage(2,'SROInt: ip-1.gt.nArr*nZeta(3)')
                  Call Abend()
               End If
*
*--------------Calculate Contraction over the spectral resolvent basis
*              set of the type <A|alm>A(l;ab)<blm|B> where we now have in
*              Array(ipF1) the cartesian components of <A|alm>, and
*              similarily, in Array(ipF2), we have stored the cartesian
*              components of <alm|B>. Observe that as opposed to the
*              projection operator that this contraction is done in the
*              primitive basis.
*
*--------------From the lefthandside overlap, form ikaC from ikac by
*              1) ika,c -> c,ika
*
               Call DgeTMo(Array(ipF1),nAlpha*nExp(iShll)*nElem(la),
     &                     nAlpha*nExp(iShll)*nElem(la),
     &                     nElem(iAng),Array(ipTmp),nElem(iAng))
*
*--------------2) ika,C = c,ika * c,C
*
               Call DGEMM_('T','N',
     &                     nAlpha*nExp(iShll)*nElem(la),
     &                     (2*iAng+1),nElem(iAng),
     &                     1.0d0,Array(ipTmp),nElem(iAng),
     &                     RSph(ipSph(iAng)),nElem(iAng),
     &                     0.0d0,Array(ipF1),
     &                     nAlpha*nExp(iShll)*nElem(la))
               If (iPrint.ge.99) Call RecPrt('<A|srbs>',' ',
     &                   Array(ipF1),nAlpha*nExp(iShll),
     &                   nElem(la)*(2*iAng+1))
*
*--------------And (almost) the same thing for the righthand side, form
*              kjCb from kjcb
*-------------1) kj,cb -> cb,kj
*
               Call DgeTMo(Array(ipF2),
     &                     nBeta*nExp(iShll),nBeta*nExp(iShll),
     &                     nElem(iAng)*nElem(lb),Array(ipTmp),
     &                     nElem(iAng)*nElem(lb))
*
*--------------2) bkj,C = c,bkj * c,C
*
               Call DGEMM_('T','N',
     &                  nElem(lb)*nExp(iShll)*nBeta,
     &                  (2*iAng+1),nElem(iAng),
     &                  1.0d0,Array(ipTmp),nElem(iAng),
     &                  RSph(ipSph(iAng)),nElem(iAng),
     &                  0.0d0,Array(ipF2),nElem(lb)*nExp(iShll)*nBeta)
*
*--------------3) b,kjC -> kjC,b
*
               Call DgeTMo(Array(ipF2),nElem(lb),nElem(lb),
     &                     nExp(iShll)*nBeta*(2*iAng+1),Array(ipTmp),
     &                     nExp(iShll)*nBeta*(2*iAng+1))
               call dcopy_(nExp(iShll)*nBeta*(2*iAng+1)*nElem(lb),
     &                    Array(ipTmp),1,Array(ipF2),1)
               If (iPrint.ge.99) Call RecPrt('<srbs|B>',' ',
     &                   Array(ipF2),nExp(iShll)*nBeta,
     &                   (2*iAng+1)*nElem(lb))
*
*--------------Next Contract (ikaC)*(klC)*(ljCb) over k,l and C,
*              producing ijab,
*              by the following procedure:
*              Loop over a and b
*                Loop over C
*                  Contract ik(aC)*kl(C), over k producing il(aC),
*                  Contract il(aC)*lj(Cb), over l producing ij(aCb)
*                    accumulate to ij(ab)
*                End loop C
*              End Loop b and a
*
               Do ib = 1, nElem(lb)
                  Do ia = 1, nElem(la)
                  If (iPrint.ge.99) Write (6,*) ' ia,ib=',ia,ib
*
                     Do iC = 1, (2*iAng+1)
                        If (iPrint.ge.99) Write (6,*) ' iC,=',iC
                        iaC = (iC-1)*nElem(la) + ia
                        ipaC = (iaC-1)*nAlpha*nExp(iShll) + ipF1
                        iCb = (ib-1)*(2*iAng+1) + iC
                        ipCb = (iCb-1)*nExp(iShll)*nBeta  + ipF2
*
                        iIC = 0
                        If (iPrint.ge.99) Then
                           Call RecPrt('<ia|iC>',' ',Array(ipaC),
     &                                  nAlpha,nExp(iShll))
                           Call RecPrt('<iC|ib>',' ',Array(ipCb),
     &                                  nExp(iShll),nBeta)
                        End If
                        Do iIrrep = 0, nIrrep-1
                           If (iAnd(llOper,iTwoj(iIrrep)).eq.0) Cycle
                           If (iPrint.ge.99) Write (6,*) ' iIC=',iIC
                           iIC = iIC + 1
                           nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
                           Xg=rChTbl(iIrrep,nOp         )
                           Factor=Xg*Fact
                           Call DGEMM_('N','N',
     &                                nAlpha,nExp(iShll),nExp(iShll),
     &                                One,Array(ipaC),nAlpha,
     &                                  Array(ipC),nExp(iShll),
     &                                Zero,Array(ipTmp),nAlpha)
                           Call DGEMM_('N','N',
     &                                nAlpha,nBeta,nExp(iShll),
     &                                Factor,Array(ipTmp),nAlpha,
     &                                    Array(ipCb),nExp(iShll),
     &                                One,Final(1,ia,ib,iIC),nAlpha)
                        End Do
*
                     End Do
                  End Do
               End Do
*
            End Do
         End Do
         End Do
         mdc = mdc + dbsc(iCnttp)%nCntr
      End Do
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in SROInt'
         Do ia = 1, (la+1)*(la+2)/2
            Do ib = 1, (lb+1)*(lb+2)/2
               Write (Label,'(A,I2,A,I2,A)')
     &               ' Final(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(1,ia,ib,1),nAlpha,nBeta)
            End Do
         End Do
      End If
*
*     Call GetMem(' Exit SROInt','LIST','REAL',iDum,iDum)
*     Call QExit('SROInt')
      Return
c Avoid unused argument warnings
      If (.False.) Then
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

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
* Copyright (C) Ben Swerts                                             *
************************************************************************
      SubRoutine FragPGrd(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                    Final,nZeta,la,lb,A,RB,nRys,
     &                    Array,nArr,Ccoor,nOrdOp,Grad,nGrad,
     &                    IfGrad,IndGrd,DAO,mdc,ndc,kOp,lOper,nComp,
     &                    iStabM,nStabM)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of FAIEMP Projection      *
*         operator integrals.                                          *
*                                                                      *
* Called from: OneEl_g                                                 *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              DCopy   (ESSL)                                          *
*              ZXia                                                    *
*              SetUp1                                                  *
*              Mlt1                                                    *
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
*      nComp : number of components in the operator (e.g. dipolemoment *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by FragMMG                *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*                                                                      *
*     based on PrjGrd                                                  *
*                                                                      *
************************************************************************
      use Her_RW
      use Real_Spherical
      use iSD_data
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: iOper
      Implicit None
#include "Molcas.fh"
#include "real.fh"
#include "print.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
      Integer  nAlpha,nBeta,nZeta,la,lb,nRys,nArr,nOrdOp,nGrad,mdc,ndc,
     &         nComp,nStabM
      Real*8   Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &         Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &         rKappa(nZeta), P(nZeta,3), A(3), RB(3), Grad(nGrad),
     &         Array(nZeta*nArr), Ccoor(3), C(3), TC(3), B(3), TB(3),
     &         DAO(nZeta,(la+1)*(la+2)/2*(lb+1)*(lb+2)/2)
      Integer  iStabM(0:nStabM-1), iDCRT(0:7), lOper(nComp),
     &         iuvwx(4), kOp(2), lOp(4),
     &         IndGrd(3,2), JndGrd(3,4)
      Character*80 Label
      Logical  IfGrad(3,2), JfGrad(3,4), ABeq(3), EQ
      Logical  EnergyWeight
      Integer  i,j,iIrrep,iComp,nElem,ia,ib,iAng,iAO,iBas
      Integer  iRout,iPrint,nSkal,iCar
      Integer  iCent,iCmp,iCnttp,iCurCenter,iCurCnttp,iCurMdc
      Integer  iGamma,iLoc,ip,ipA,ipAxyz,ipB,ipBxyz,ipCxyz,ipF1,ipF2
      Integer  ipF1a,ipF2a,ipIJ,ipK1,ipK2,ipP1,ipP2,ipQ1,iPrim,ipRxyz
      Integer  ipTmp,ipZ1,ipZ2,ipZI1,ipZI2,iS,iSbasis,iSEnd,iShell,iShll
      Integer  iSize,iSlocal,iSstart,iStemp,iStrt,iVec,jAng,jBas
      Integer  jCmp,jCnttp,jPrim,jS,jSbasis,jShell,jShll,jSize
      Integer  jSlocal,ld,lDCRT,LmbdT,mdci,mdcj,mGrad,mVec,mVecAC
      Integer  mVecCB,nac,ncb,nDAO,nDCRT,nDisp,nHer,jAO,maxDensSize
      Integer  nVecAC,nVecCB,iTri,iCnt, jCnt
      Real*8   Fact,DNrm2_
      External DNrm2_
*
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      iRout = 202
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
         Call RecPrt(' In FragPGrd: Grad',' ',Grad,1,nGrad)
         Call RecPrt(' In FragPGrd: A',' ',A,1,3)
         Call RecPrt(' In FragPGrd: RB',' ',RB,1,3)
         Call RecPrt(' In FragPGrd: Ccoor',' ',Ccoor,1,3)
         Call RecPrt(' In FragPGrd: P',' ',P,nZeta,3)
         Call RecPrt(' In FragPGrd: Alpha',' ',Alpha,nAlpha,1)
         Call RecPrt(' In FragPGrd: Beta',' ',Beta,nBeta,1)
         Write (6,*) ' In FragPGrd: la,lb=',' ',la,lb
      End If
*                                                                      *
************************************************************************
*                                                                      *
* Setup the fragment shells
*
      Call Set_Basis_Mode('Fragments')
      Call SetUp_iSD
      Call Nr_Shells(nSkal)
      If(iPrint.ge.99) Then
        write(6,*) 'looping over ',nSkal,' shells'
        write(6,*) 'Shells()%Frag = ',(Shells(i)%Frag,i=1,10)
      End If

*                                                                      *
************************************************************************
*                                                                      *
* Reserve space for the largest possible fragment energy weighted
* density matrix
      maxDensSize = 0
      Do iCnttp = 1, nCnttp
        If(dbsc(iCnttp)%nFragType.gt.0) maxDensSize = Max(maxDensSize,
     &      dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
* Loop over all shells belonging to the fragments
*
      nDAO = nElem(la)*nElem(lb)
      iIrrep = 0
      iuvwx(1) = dc(mdc)%nStab
      iuvwx(2) = dc(ndc)%nStab
      lOp(1) = iOper(kOp(1))
      lOp(2) = iOper(kOp(2))
*
      iComp = 1
      iCurMdc = 0
c      ! The mdc value of the current fragment placeholder
      iCurCnttp = 0
c      ! The Cnttp of the fragment placeholder
      iCurCenter = 999999
c      ! The index of the fragment in the fragment placeholder list of centers
      iSstart = 0
c      ! The index into the full shells list for the first shell of a fragment
      iSbasis = 0
c      ! The basis function index relative to the start of the fragment
      iSEnd = -1
c      ! Dummy initialize
      Do 1965 iS = 1, nSkal
        iShll  = iSD( 0,iS)
        iAng   = iSD( 1,iS)
        iCmp   = iSD( 2,iS)
        iBas   = iSD( 3,iS)
        iPrim  = iSD( 5,iS)
        iAO    = iSD( 7,iS)
        mdci   = iSD(10,iS)
        iShell = iSD(11,iS)
        iCnttp = iSD(13,iS)
        iCnt   = iSD(14,iS)
        C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

        iSize = nElem(iAng)
        if(Shells(iShll)%Transf.and.
     &     Shells(iShll)%Prjct ) iSize = 2*iAng+1
        If(Abs(dbsc(iCnttp)%nFragCoor).ne.iCurMdc) Then
* update fragment related quantities
          iCurMdc = Abs(dbsc(iCnttp)%nFragCoor)
          iSstart = iS
          iSend = nSkal
          Do iStemp = iSstart + 1,nSkal
            If(Abs(dbsc(iSD(13,iStemp))%nFragCoor).ne.iCurMdc) Then
              iSend = iStemp - 1
              goto 101
            End If
          End Do
 101      Continue
          iSbasis= 1
          iCurCenter = iCurCenter + 1
          If(iCurCenter.gt.dbsc(iCurCnttp)%nCntr) Then
            iCurCenter = 1
            Do jCnttp = iCurCnttp+1, nCnttp
              If(dbsc(jCnttp)%nFragType.gt.0) Then
                iCurCnttp = jCnttp
                goto 102
              End If
            End Do
 102        Continue
* update the energy weighted density matrix of the current fragment
            EnergyWeight = .true.
            Call MakeDens(dbsc(iCurCnttp)%nFragDens,
     &                    dbsc(iCurCnttp)%nFragEner,
     &                    dbsc(iCurCnttp)%FragCoef,
     &                    dbsc(iCurCnttp)%FragEner,
     &                    EnergyWeight,Array)
           If(iPrint.ge.49) Call TriPrt('Energy weighted fragment dens',
     &        ' ',Array,dbsc(iCurCnttp)%nFragDens)
* include the minus sign of -2eta_i
           Call DScal_(dbsc(iCurCnttp)%nFragDens
     &               *(dbsc(iCurCnttp)%nFragDens+1)/2,
     &                -One,Array,1)
           If(maxDensSize.lt.dbsc(iCurCnttp)%nFragDens*
     &        (dbsc(iCurCnttp)%nFragDens+1)/2) Stop 'maxIJSize'
          End If
        End If
c        write(*,*) '  iShll,iAng,mdci,iShell,iCnttp,iCurMdc,iCurCnttp',
c     &              iShll,iAng,mdci,iShell,iCnttp,iCurMdc,iCurCnttp
c        write(*,*) '  iPrim,iBas =',iPrim,iBas
*
* extra derivative stuff
            iuvwx(3) = dc(mdci)%nStab
            iuvwx(4) = dc(mdci)%nStab
            Call ICopy(6,IndGrd,1,JndGrd,1)
            Do i = 1, 3
               Do j = 1, 2
                 JfGrad(i,j) = IfGrad(i,j)
               End Do
            End Do
*
            nDisp = IndDsp(mdci,iIrrep)
            Do iCar = 0, 2
              JfGrad(iCar+1,3) = .False.
              iCmp = 2**iCar
* always equivalent of pChrg's
              JndGrd(iCar+1,3) = 0
            End Do
            Call ICopy(3,[0],0,JndGrd(1,4),1)
            JfGrad(1,4) = .False.
            JfGrad(2,4) = .False.
            JfGrad(3,4) = .False.
            mGrad = 0
            Do iCar = 1, 3
               Do i = 1, 2
                  If (JfGrad(iCar,i)) mGrad = mGrad + 1
               End Do
            End Do
            If (iPrint.ge.99) Write (6,*) ' mGrad=',mGrad
            If (mGrad.eq.0) Go To 1965
*                                                                      *
************************************************************************
*                                                                      *
* Loop over all other shells belonging to the same fragment
        jSbasis = 1
        Do jS = iSstart, iSend
          jShll  = iSD( 0,jS)
          jAng   = iSD( 1,jS)
          jCmp   = iSD( 2,jS)
          jBas   = iSD( 3,jS)
          jPrim  = iSD( 5,jS)
          jAO    = iSD( 7,iS)
          mdcj   = iSD(10,jS)
          jShell = iSD(11,jS)
          jCnttp = iSD(13,jS)
          jCnt   = iSD(14,jS)
          jSize = nElem(jAng)
          If(Shells(jShll)%Transf.and.
     &       Shells(jShll)%Prjct ) jSize = 2*jAng+1
          B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
c         write(*,*) '    jShll,jAng,mdcj,jShell,jCnttp =',
c    &                    jShll,jAng,mdcj,jShell,jCnttp
c         write(*,*) '    jPrim,jBas =',jPrim,jBas
*                                                                      *
************************************************************************
*                                                                      *
* Create a rectangular matrix sized (iBas*nElem(iAng),jBas*nElem(jAng))
* from the energy weighted density matrix (desymmetrized)
* contains values from iSbasis to iSbasis + iBas*nElem(iAng) - 1
*             and from jSbasis to jSbasis + jBas*nElem(jAng) - 1
           ipIJ = 1 + maxDensSize
c          write(*,*) '    extracting values from',iSbasis,' to',iSbasis
c    &                + iBas*iSize - 1,', and from',jSbasis,' to',
c    &                jSbasis + jBas*jSize - 1
           Do iSlocal = iSbasis, iSbasis + iBas*iSize - 1
             Do jSlocal = jSbasis, jSbasis + jBas*jSize - 1
               iLoc = ipIJ + (jSlocal-jSbasis)*iBas*iSize + iSlocal
     &                - iSbasis
               Array(iLoc) = Array(iTri(iSlocal,jSlocal))
               If(iSlocal.ne.jSlocal) Array(iLoc) = Array(iLoc)/Two
c              write(*,*) 'Filling (',iSlocal-iSbasis+1,',',
c    &           jSlocal-jSbasis+1,') from (',iSlocal,',',jSlocal,')'
             End Do
           End Do
           If(iPrint.ge.99) Call RecPrt('W(KC,LD)',
     &        ' ',Array(ipIJ),iBas*iSize,jBas*jSize)
*                                                                      *
************************************************************************
*                                                                      *
* DCR stuff (iS and jS have always the same symmetry character)
*
          Call DCR(LmbdT,iStabM,nStabM,
     &                   dc(mdci)%iStab,dc(mdci)%nStab,iDCRT,nDCRT)
          Fact = DBLE(nStabM) / DBLE(LmbdT)
*                                                                      *
************************************************************************
*                                                                      *
* Loop over symmetry operations acting on the basis.
*
          Do 1967 lDCRT = 0, nDCRT-1
            lOp(3) = iDCRT(lDCRT)
            lOp(4) = lOp(3)
            Call OA(iDCRT(lDCRT),C,TC)
            Call OA(iDCRT(lDCRT),B,TB)
            If (EQ(A,RB).and.EQ(A,TC)) Go To 1967
*                                                                      *
************************************************************************
*                                                                      *
* Calculate the overlap integral < alpha | is > and its derivative
*
**** Storage
*
            ip = ipIJ + maxDensSize
            ipF1 = ip
            nac = nElem(la)*nElem(iAng)*4
            ip = ip + nAlpha*nac*iPrim
            ipP1 = ip
            ip = ip + 3 * nAlpha*iPrim
            ipZ1 = ip
            ip = ip + nAlpha*iPrim
            ipK1 = ip
            ip = ip + nAlpha*iPrim
            ipZI1 = ip
            ip = ip + nAlpha*iPrim
            If (ip-1.gt.nArr*nZeta) Then
               Write (6,*) '  ip-1.gt.nArr*nZeta(1) in FragPGrd'
               Call Abend()
            End If
*
**** Effective center and exponent
*
            Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,iPrim,
     &                Alpha,Shells(iShll)%Exp)
            Call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,
     &                  A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
**** Overlap and derivative
*
            nHer = ((la+1)+iAng+2)/2
            ipAxyz = ip
            ip = ip + nAlpha*iPrim*3*nHer*(la+2)
            ipCxyz = ip
            ip = ip + nAlpha*iPrim*3*nHer*(iAng+1)
            ipRxyz = ip
            ip = ip + nAlpha*iPrim*3*nHer*(nOrdOp+1)
            ipQ1 = ip
            ip = ip + nAlpha*iPrim*3*(la+2)*(iAng+1)*(nOrdOp+1)
            ipA = ip
            ip = ip + nAlpha*iPrim
            If (ip-1.gt.nArr*nZeta) Then
               Write (6,*) '  ip-1.gt.nArr*nZeta(1b) in FragPGrd'
               Call Abend()
            End If
            ABeq(1) = A(1).eq.TC(1)
            ABeq(2) = A(2).eq.TC(2)
            ABeq(3) = A(3).eq.TC(3)
            Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,
     &                  A,Array(ipAxyz),la+1,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,
     &                  TC,Array(ipCxyz),iAng,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            ABeq(1) = .False.
            ABeq(2) = .False.
            ABeq(3) = .False.
            Call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,
     &                  Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            If (iPrint.ge.49) Then
               Write (6,*) ' Array(ipAxyz)=',
     &           DNrm2_(nAlpha*iPrim*3*nHer*(la+2),Array(ipAxyz),1)
               Write (6,*) ' Array(ipCxyz)=',
     &           DNrm2_(nAlpha*iPrim*3*nHer*(iAng+1),Array(ipCxyz),1)
               Write (6,*) ' Array(ipRxyz)=',
     &           DNrm2_(nAlpha*iPrim*3*nHer*(nOrdOp+1),Array(ipRxyz),1)
            End If
            Call Assmbl(Array(ipQ1),
     &                  Array(ipAxyz),la+1,
     &                  Array(ipRxyz),nOrdOp,
     &                  Array(ipCxyz),iAng,
     &                  nAlpha*iPrim,HerW(iHerW(nHer)),nHer)
            iStrt = ipA
            Do iGamma = 1, iPrim
              call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
              iStrt = iStrt + nAlpha
            End Do
            If (iPrint.ge.49) Write (6,*) ' Array(ipA)=',
     &        DNrm2_(nAlpha*iPrim,Array(ipA),1)
            Call rKappa_Zeta(Array(ipK1),Array(ipZ1),iPrim*nAlpha)
            ld=1
            Call CmbnAC(Array(ipQ1),nAlpha*iPrim,la,iAng,
     &                  Array(ipK1),Array(ipF1),
     &                  Array(ipA),JfGrad(1,1),ld,nVecAC)
            If (iPrint.ge.49) Then
              Write (6,*) ' Array(ipQ1)=',
     &          DNrm2_(nAlpha*iPrim*3*(la+2)*(iAng+1)*(nOrdOp+1),
     &                Array(ipQ1),1)
              Write (6,*) ' Array(ipA)=',
     &          DNrm2_(nAlpha*iPrim,Array(ipA),1)
            End If
            ip = ip - nAlpha*iPrim
     &            * ( 6 + 3*nHer*(la+2) + 3*nHer*(iAng+1)
     &            + 3*nHer*(nOrdOp+1) + 3*(la+2)*(iAng+1)*(nOrdOp+1) +1)
*                                                                      *
************************************************************************
*                                                                      *
* Calculate the overlap integral < js | beta > and its derivative
*
**** Storage
*
            ipF2 = ip
            ncb = nElem(jAng)*nElem(lb)*4
            ip = ip + jPrim*nBeta*ncb
            ipP2 = ip
            ip = ip + 3 * jPrim*nBeta
            ipZ2 = ip
            ip = ip + jPrim*nBeta
            ipK2 = ip
            ip = ip + jPrim*nBeta
            ipZI2 = ip
            ip = ip + jPrim*nBeta
            If (ip-1.gt.nArr*nZeta) Then
              Write (6,*) '  ip-1.gt.nArr*nZeta(2) in FragPGrd'
              Call Abend()
            End If
*
**** Effective center and exponent
*
            Call ZXia(Array(ipZ2),Array(ipZI2),jPrim,nBeta,
     &                Shells(jShll)%Exp,Beta)
            Call SetUp1(Shells(jShll)%Exp,jPrim,Beta,nBeta,
     &                  TB,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
**** Overlap and derivative
*
            nHer = (jAng+(lb+1)+2)/2
            ipCxyz = ip
            ip = ip + nBeta*jPrim*3*nHer*(jAng+1)
            ipBxyz = ip
            ip = ip + nBeta*jPrim*3*nHer*(lb+2)
            ipRxyz = ip
            ip = ip + nBeta*jPrim*3*nHer*(nOrdOp+1)
            ipQ1 = ip
            ip = ip + nBeta*jPrim*3*(jAng+1)*(lb+2)*(nOrdOp+1)
            ipB = ip
            ip = ip + nBeta*jPrim
            If (ip-1.gt.nArr*nZeta) Then
              Write (6,*) '  ip-1.gt.nArr*nZeta(2b) in FragPGrd'
              Call Abend()
            End If
            ABeq(1) = TB(1).eq.RB(1)
            ABeq(2) = TB(2).eq.RB(2)
            ABeq(3) = TB(3).eq.RB(3)
            Call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,
     &                  TB,Array(ipCxyz),jAng,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            Call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,
     &                  RB,Array(ipBxyz),lb+1,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            ABeq(1) = .False.
            ABeq(2) = .False.
            ABeq(3) = .False.
            Call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,
     &                  Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(nHer)),
     &                  nHer,ABeq)
            If (iPrint.ge.49) Then
              Write (6,*) ' Array(ipCxyz)=',
     &          DNrm2_(nBeta*jPrim*3*nHer*(jAng+1),Array(ipCxyz),1)
              Write (6,*) ' Array(ipBxyz)=',
     &          DNrm2_(nBeta*jPrim*3*nHer*(lb+2),Array(ipBxyz),1)
              Write (6,*) ' Array(ipRxyz)=',
     &          DNrm2_(nBeta*jPrim*3*nHer*(nOrdOp+1),Array(ipRxyz),1)
            End If
            Call Assmbl(Array(ipQ1),
     &                  Array(ipCxyz),jAng,
     &                  Array(ipRxyz),nOrdOp,
     &                  Array(ipBxyz),lb+1,
     &                  jPrim*nBeta,HerW(iHerW(nHer)),nHer)
            iStrt = ipB
            Do iGamma = 1, jPrim
              call dcopy_(nBeta,Beta,1,Array(iStrt),jPrim)
              iStrt = iStrt + 1
            End Do
            If (iPrint.ge.49) Write (6,*) ' Array(ipB)=',
     &        DNrm2_(jPrim*nBeta,Array(ipB),1)
            Call rKappa_Zeta(Array(ipK2),Array(ipZ2),jPrim*nBeta)
            ld=1
            Call CmbnCB(Array(ipQ1),jPrim*nBeta,jAng,lb,
     &                  Array(ipK2),Array(ipF2),
     &                  Array(ipB),JfGrad(1,2),ld,nVecCB)
            If (iPrint.ge.49) Then
              Write (6,*) ' Array(ipQ1)=',
     &          DNrm2_(jPrim*nBeta*3*(la+2)*(jAng+1)*(nOrdOp+1),
     &                Array(ipQ1),1)
              Write (6,*)' Array(ipB)=',DNrm2_(JPrim*nBeta,Array(ipB),1)
            End If
            ip = ip - nBeta*jPrim
     &            * ( 6 + 3*nHer*(lb+2) + 3*nHer*(jAng+1)
     &            + 3*nHer*(nOrdOp+1) + 3*(lb+2)*(jAng+1)*(nOrdOp+1) +1)
            nac = nElem(la)*nElem(iAng)*nVecAC
            ncb = nElem(jAng)*nElem(lb)*nVecCB
            ipTmp = ip
            ip = ip + Max(nAlpha*Max(iPrim,jBas)*nac,nBeta*ncb*jBas)
            If (ip-1.gt.nArr*nZeta) Then
              Write (6,*) '  ip-1.gt.nArr*nZeta(3) in FragPGrd'
              Call Abend()
            End If
            nac = nElem(la)*nElem(iAng)
            ncb = nElem(jAng)*nElem(lb)
*                                                                      *
************************************************************************
*                                                                      *
* Assemble the calculated quantities and contract
*
* Calculate Contraction over components of the fragment
* orbitals of type <A|iS>coef<jS|B> where we now have in
* Array(ipF1) the cartesian components of <A|iS>, and
* similarily, in Array(ipF2), we have stored the cartesian
* components of <jS|B>. Observe that the fragment orbitals are
* orthonomal atomic orbitals. Hence, the transformation
* to the spherical harmonics has to be for normalized
* spherical harmonics.
*
* nAlpha = i               nElem(la) = a
* nBeta  = j               nElem(lb) = b
* iPrim = k (iBas = K)     nElem(iAng) = c (iSize = C)
* jPrim = l (jBas = l)     nElem(jAng) = d (jSize = D)
*
*-----------From the lefthandside overlap, form iKaC from ikac by
*           1) i,kac -> k,aci
*
            Call DgeTMo(Array(ipF1),nAlpha,nAlpha,
     &                  iPrim*nac*nVecAC,Array(ipTmp),
     &                  iPrim*nac*nVecAC)
*
*-----------2) aciK =  k,aci * k,K (Contract over core orbital)
*
            Call DGEMM_('T','N',
     &                  nac*nVecAC*nAlpha,iBas,iPrim,
     &                  1.0d0,Array(ipTmp),iPrim,
     &                  Shells(iShll)%pCff,iPrim,
     &                  0.0d0,Array(ipF1),nac*nVecAC*nAlpha)
*
*-----------3) a,ciK -> ciKa
*
            Call DgeTMo(Array(ipF1),nElem(la),nElem(la),
     &                  nElem(iAng)*nVecAC*nAlpha*iBas,
     &                  Array(ipTmp),
     &                  nElem(iAng)*nVecAC*nAlpha*iBas)
*
*-----------4) iKa,C = c,iKa * c,C
*
            If(Shells(iShll)%Transf.and.
     &         Shells(iShll)%Prjct ) Then
              Call DGEMM_('T','N',
     &                    nVecAC*nAlpha*iBas*nElem(la),iSize,
     &                    nElem(iAng),
     &                    1.0d0,Array(ipTmp),nElem(iAng),
     &                    RSph(ipSph(iAng)),nElem(iAng),
     &                    0.0d0,Array(ipF1),
     &                    nVecAC*nAlpha*iBas*nElem(la))
            Else
              Call DgeTMo(Array(ipTmp),nElem(iAng),nElem(iAng),
     &                    nVecAC*iBas*nElem(la)*nAlpha,
     &                    Array(ipF1),
     &                    nVecAC*iBas*nElem(la)*nAlpha)
            End If
C what does this do and is it needed? (from PrjGrd)
            Call DgeTMo(Array(ipF1),nVecAC,nVecAC,
     &                  nAlpha*iBas*nElem(la)*iSize,
     &                  Array(ipTmp),
     &                  nAlpha*iBas*nElem(la)*iSize)
            call dcopy_(nVecAC*nAlpha*iBas*nElem(la)*iSize,Array(ipTmp),
     &                 1,Array(ipF1),1)
*
*-----------And (almost) the same thing for the righthand side, form
*           LjDb from ljdb
*           1) jdb,L = l,jdb * l,L
*
            Call DGEMM_('T','N',
     &                  nBeta*ncb*nVecCB,jBas,jPrim,
     &                  1.0d0,Array(ipF2),jPrim,
     &                  Shells(jShll)%pCff,jPrim,
     &                  0.0d0,Array(ipTmp),nBeta*ncb*nVecCB)
*
*-----------2)  j,dbL -> dbL,j
*
            Call DgeTMo(Array(ipTmp),nBeta,nBeta,
     &                  ncb*nVecCB*jBas,Array(ipF2),
     &                  ncb*nVecCB*jBas)
*
*-----------3) bLj,D = d,bLj * d,D
*
            If(Shells(jShll)%Transf.and.
     &         Shells(jShll)%Prjct ) Then
              Call DGEMM_('T','N',
     &                    nElem(lb)*nVecCB*jBas*nBeta,jSize,nElem(jAng),
     &                    1.0d0,Array(ipF2),nElem(jAng),
     &                    RSph(ipSph(jAng)),nElem(jAng),
     &                   0.0d0,Array(ipTmp),nElem(lb)*nVecCB*jBas*nBeta)
            Else
              Call DgeTMo(Array(ipF2),nElem(jAng),nElem(jAng),
     &                    nVecCB*jBas*nElem(lb)*nBeta,
     &                    Array(ipTmp),
     &                    nVecCB*jBas*nElem(lb)*nBeta)
            End If
*
*-----------4) b,LjD -> LjD,b
*
            Call DgeTMo(Array(ipTmp),nElem(lb)*nVecCB,
     &                  nElem(lb)*nVecCB,
     &                  jBas*nBeta*jSize,Array(ipF2),
     &                  jBas*nBeta*jSize)
*
*-----------Next Contract (iKaC)*W(KLCD)*(LjDb) producing ijab
*
            call dcopy_(nZeta*nElem(la)*nElem(lb)*6,[Zero],0,Final,1)
*
            If(iPrint.ge.99) Then
              Call RecPrt('ipF1 (nVecAC x X)',' ',Array(ipF1),
     &          nVecAC,iBas*nAlpha*iSize)
              Call RecPrt('ipF2 (nVecCB x Y)',' ',Array(ipF2),
     &          nVecCB,jBas*nBeta*jSize)
            End If
*
            mVec = 0
            mVecAC = 1
            mVecCB = 1
            Do iCar = 1, 3
              Do iCent = 1, 2
c               write(*,*) 'iCar, iCent = ',iCar,iCent
                If (JfGrad(iCar,iCent)) Then
                  mVec = mVec + 1
                  If (iCent.eq.1) Then
                    mVecAC = mVecAC+1
                    ipF1a = ipF1 + (mVecAC-1) *
     &              nAlpha*jBas*nElem(la)*iSize
                    ipF2a = ipF2
                  Else
                    ipF1a = ipF1
                    mVecCB = mVecCB+1
                    ipF2a = ipF2 + (mVecCB-1) *
     &                      jBas*nBeta*jSize*nElem(lb)
                  End If
            If(iPrint.ge.99) Then
              write(6,*) 'mVecAC, mVecCB = ',mVecAC,mVecCB
              Call RecPrt('ipF1a (nAlpha*aAng x iBas*iSize)',' ',
     &          Array(ipF1a),nAlpha*nElem(la),iBas*iSize)
              Call RecPrt('ipF2a (nBeta*bAng x jBas*jSize)',' ',
     &          Array(ipF2a),nBeta*nElem(lb),jBas*jSize)
            End If
*
               Call FragPCont(Array(ipF1a), nAlpha,iBas,nElem(la),iSize,
     &                        Array(ipF2a), jBas,nBeta,jSize,nElem(lb),
c     &                        Array(ipIJ), Final(1,1,1,mVec), Fact*Half)
     &                        Array(ipIJ), Final(:,:,:,mVec), Fact*Half)
                End If
              End Do !iCent
            End Do !iCar
*
            If (iPrint.ge.49) Then
              Do iVec = 1, mVec
                Write (6,*) iVec,
     &            Sqrt(DNrm2_(nZeta*nElem(la)*nElem(lb),
     &            Final(1,1,1,iVec),1))
              End Do
            End If
            If (iPrint.ge.99) Then
              Write (6,*) ' Result in FragPGrd'
              Do ia = 1, nElem(la)
                Do ib = 1, nElem(lb)
                  Do iVec = 1, mVec
      Write (Label,'(A,I2,A,I2,A,I2,A)')' Final(',ia,',',ib,',',iVec,')'
                    Call RecPrt(Label,' ',Final(1,ia,ib,iVec),
     &                          nAlpha,nBeta)
                  End Do
                End Do
              End Do
            End If
*
*-----------Distribute contributions to the gradient
*
            Call Distg1X(Final,DAO,nZeta,nDAO,mVec,Grad,nGrad,
     &                   JfGrad,JndGrd,iuvwx,lOp)
*
 1967     Continue !lDCRT
          jSbasis = jSbasis + jBas * jSize
        End Do !jS
        iSbasis = iSbasis + iBas * iSize
 1965 Continue !iS
*
* Revert to the valence shells
*
      Call Free_iSD()
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Zeta)
        Call Unused_real_array(ZInv)
        Call Unused_real_array(rKappa)
        Call Unused_integer(nRys)
        Call Unused_integer_array(lOper)
      End If
      End

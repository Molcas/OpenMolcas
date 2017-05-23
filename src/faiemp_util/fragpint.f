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
*               2016, Liviu Ungur                                      *
************************************************************************
      SubRoutine FragPInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                    Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,
     &                    Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                    iStabM,nStabM,
     &                    PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of Fragment AIEMP         *
*         projection integrals                                         *
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
*      nComp : number of components in the operator (e.g. dipolemoment *
*              operator has three components)                          *
*      la    : total angular momentum of bra gaussian                  *
*      lb    : total angular momentum of ket gaussian                  *
*      A     : center of bra gaussian                                  *
*      B     : center of ket gaussian                                  *
*      nRys  : order of Rys- or Hermite-Gauss polynomial               *
*      Array : Auxiliary memory as requested by PrjMem                 *
*      nArr  : length of Array                                         *
*      Ccoor : coordinates of the operator, zero for symmetric oper.   *
*      NOrdOp: Order of the operator                                   *
*                                                                      *
*     Author: Ben Swerts                                               *
*     based on seward/prjint.f                                         *
*   Modified: Liviu Ungur                                              *
*                                                                      *
* additional info:                                                     *
*       mdci, mdci = show the type number (i.e. from 1 to nCnttp) of   *
*                    the corresponding "fragment" atom. Must never take*
*                    values of the atoms from the "central region"     *
*   iCnttp, jCnttp = index number of the the type of the fragment atom *
*  nFragType(mdci) = number of different basis sets (kinds of atoms)   *
*                    are in the "fragment"                             *
*  nFragCoor(mdci) = number of atoms of the fragment "mdci"            *
*                    must be 1, as we expanded atoms in one atom/group *
*                    check also the "FragExpand" subroutine            *
*  nFragDens(mdci) = nBas of the entire fragment type "mdci"           *
************************************************************************
      use Real_Spherical
      use iSD_data
      Implicit None
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
#include "para_info.fh"
      Integer nZeta,la,lb,nIC,nArr,nComp,nAlpha,nBeta,nRys,nOrdOp,
     &        nStabM
      Real*8  Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &        Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &        rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &        Array(nZeta*nArr), Ccoor(3), C(3), TC(3), B(3), TB(3)
      Integer iStabM(0:nStabM-1), lOper(nComp), iDCRT(0:7),iTwoj(0:7),
     &        iChO(nComp),iAddPot
      Logical EnergyWeight
c      Character*80 Label
      Character*24 Label
      Data    iTwoj/1,2,4,8,16,32,64,128/
*
      Integer i,j,ixyz,nElem,iTri,nGrid,
     &        iRout,iPrint,ia,iAng,ib,iBas,iAO,iCff,iCmp,iCnttp,iComp,
     &        iCurCenter,iCurCnttp,iCurMdc,iExp,iIC,iIrrep,iLoc,iPrim,
     &        ip,ipF1,ipF2,ipIJ,ipK1,ipK2,ipP1,ipP2,ipTmp,ipZ1,ipZ2,
     &        ipZI1,ipZI2,iS,iSbasis,iSend,iShell,iShll,iSize,iSlocal,
     &        iSstart,iStemp,jAng,jAO,jBas,jCff,jCmp,jCnttp,jExp,jPrim,
     &        jS,jShell,jShll,jSize,jSlocal,jxyz,lDCRT,llOper,LmbdT,
     &        mArr,maxDensSize,mdci,mdcj,nac,ncb,nDCRT,nHer,nOp,nSkal,
     &        jSbasis
      Real*8  Fact,Factor,PtChrg,Xg
* external functions:
      Integer NrOpr
      Real*8  DNRM2_
      External DNRM2_,NrOpr
*
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
      iRout = 202
      iPrint = nPrint(iRout)
*
      If (iPrint.ge.49) Then
c data for individual fragments:
         Write (6,*) ' In FragPInt:    nCnttp          = ',nCnttp
         Write (6,*) ' In FragPInt: mdciCnttp(nCnttp)  = ',
     &                                    (mdciCnttp(i),i=1,nCnttp)
         Write (6,*) ' In FragPInt:     nCntr(nCnttp)  = ',
     &                                        (nCntr(i),i=1,nCnttp)
         Write (6,*) ' In FragPInt: nFragType(nCnttp)  = ',
     &                                    (nFragType(i),i=1,nCnttp)
         Write (6,*) ' In FragPInt: nFragCoor(nCnttp)  = ',
     &                                    (nFragCoor(i),i=1,nCnttp)
         Write (6,*) ' In FragPInt: nFragDens(nCnttp)  = ',
     &                                    (nFragDens(i),i=1,nCnttp)
         Write (6,*) ' In FragPInt: nFragDens(nCnttp)  = ',
     &                                    (nFragDens(i),i=1,nCnttp)
c
         Write (6,*) ' In FragPInt: nAlpha,nBeta=',' ',nAlpha,nBeta
         Write (6,*) ' In FragPInt: Alpha=',' ',Alpha
         Write (6,*) ' In FragPInt: Beta=',' ',Beta
c
         Write (6,*) ' In FragPInt: nZeta=',' ',nZeta
         Write (6,*) ' In FragPInt:  nArr=',' ',nArr
         Write (6,*) ' In FragPInt:   nIC=',' ',nIC
         Write (6,*) ' In FragPInt: la,lb=',' ',la,lb
         Write (6,*) ' In FragPInt: nElem(la)=',' ',nElem(la)
         Write (6,*) ' In FragPInt: nElem(lb)=',' ',nElem(lb)
c         Call RecPrt(' In FragPInt: A     ',' ',A     ,1,3)
         Call RecPrt(' In FragPInt: RB    ',' ',RB    ,1,3)
         Call RecPrt(' In FragPInt: Ccoor ',' ',Ccoor ,1,3)
c         Call RecPrt(' In FragPInt: P     ',' ',P     ,nZeta,3)
c         Call RecPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
c         Call TrcPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
      End If
*
      Call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,Zero,0,Final,1)
*                                                                      *
************************************************************************
*                                                                      *
* Setup the fragment shells
*
      Call Free_iSD
      Call Set_Basis_Mode('Fragments')
      Call SetUp_iSD
      Call Nr_Shells(nSkal)
      If (iPrint.ge.49) print *, 'nSkal_Fragment,nAlpha,nBeta = ',
     &                              nSkal,nAlpha,nBeta
*                                                                      *
************************************************************************
*                                                                      *
* Reserve space for the largest possible fragment energy weighted
* density matrix and one combination of shells of it.
      maxDensSize = 0
      Do iCnttp = 1, nCnttp
        If(nFragType(iCnttp).gt.0) then
        maxDensSize = Max( maxDensSize,
     &                     nFragDens(iCnttp)*(nFragDens(iCnttp)+1)/2)
c        write(6,'(A,i2,A,i6)') 'nFragDens(',iCnttp,')=',
c     &                                              nFragDens(iCnttp)
        endif
      End Do
*                                                                      *
************************************************************************
*                                                                      *
* Loop over all shells belonging to the fragments
*
      llOper = lOper(1)
      iComp = 1
      iCurMdc = 0
c      ! The mdc value of the current fragment placeholder
      iCurCnttp = 1
c     ! The Cnttp of the fragment placeholder
      iCurCenter = 999999
c     ! The index of the fragment in the fragment placeholder list of centers
      iSstart = 0
c     ! The index into the full shells list for the first shell of a fragment
      isEnd = nSkal
c     ! The index into the full shells list for the last shell of a fragment
      iSbasis = 0
c     ! The basis function index relative to the start of the fragment
      Do iS = 1, nSkal
        iShll  = iSD( 0,iS)
        iAng   = iSD( 1,iS)
        iCmp   = iSD( 2,iS)
        iBas   = iSD( 3,iS)
        iCff   = iSD( 4,iS)
        iPrim  = iSD( 5,iS)
        iExp   = iSD( 6,iS)
        iAO    = iSD( 7,iS)
        ixyz   = iSD( 8,iS)
        mdci   = iSD(10,iS)
        iShell = iSD(11,iS)
        iCnttp = iSD(13,iS)
        iSize = nElem(iAng)
c some printouts:
c      If (iPrint.ge.49) Then
c        print *, 'In FragPInt: iS=',iS,' iShll =',iShll
c        print *, 'In FragPInt: iS=',iS,' iAng  =',iAng
c        print *, 'In FragPInt: iS=',iS,' iCmp  =',iCmp
c        print *, 'In FragPInt: iS=',iS,' iBas  =',iBas
c        print *, 'In FragPInt: iS=',iS,' iCff  =',iCff
c        print *, 'In FragPInt: iS=',iS,' iPrim =',iPrim
c        print *, 'In FragPInt: iS=',iS,' iExp  =',iExp
c        print *, 'In FragPInt: iS=',iS,' iAO   =',iAO
c        print *, 'In FragPInt: iS=',iS,' ixyz  =',ixyz
c        print *, 'In FragPInt: iS=',iS,' mdci  =',mdci
c        print *, 'In FragPInt: iS=',iS,' iShell=',iShell
c        print *, 'In FragPInt: iS=',iS,' iCnttp=',iCnttp
c        print *, 'In FragPInt: iS=',iS,' iSize =',iSize
c        print *, 'In FragPInt: iS=',iS,' iCurMdc =',iCurMdc
c        print *, 'In FragPInt: nFragCoor(',mdci,') =',nFragCoor(mdci)
c      End If

        If(Transf(iShll).and.Prjct(iShll)) iSize = 2*iAng+1
        If(nFragCoor(mdci).ne.iCurMdc) Then
* update fragment related quantities
          iCurMdc = nFragCoor(mdci)
          iSstart = iS
          iSend = nSkal
          Do iStemp = iSstart + 1,nSkal
            If(nFragCoor(iSD(10,iStemp)).ne.iCurMdc) Then
              iSend = iStemp - 1
              goto 101
            End If
          End Do
 101      Continue
          iSbasis= 1
          iCurCenter = iCurCenter + 1

      If (iPrint.ge.49) Then
         write(6,*) 'start of new fragment encountered'
         write(6,*) 'iSstart,iSend,iCurCnttp,iCurCenter =',
     &               iSstart,iSend,iCurCnttp,iCurCenter
      End If

          If(iCurCenter.gt.nCntr(iCurCnttp)) Then
            iCurCenter = 1
            Do jCnttp = iCurCnttp+1, nCnttp
              If(nFragType(jCnttp).gt.0) Then
                iCurCnttp = jCnttp
                goto 102
              End If
            End Do

 102        Continue
* update the energy weighted density matrix of the current fragment
            EnergyWeight = .true.
c            write(6,*) 'Calling MakeDens 1: iCurCnttp=',iCurCnttp
c            write(6,*) 'nFragDens(iCurCnttp)=',nFragDens(iCurCnttp)
c            write(6,*) 'nFragEner(iCurCnttp)=',nFragEner(iCurCnttp)
c           Call MakeDens(                nBas,                nOrb,
c                                     Cff,                      OrbEn,
c             EnergyWeight, Dens)
            Call MakeDens(nFragDens(iCurCnttp),nFragEner(iCurCnttp),
     &        Work(ipFragCoef(iCurCnttp)),Work(ipFragEner(iCurCnttp)),
     &        EnergyWeight,Array)
            If(iPrint.ge.49) Call TriPrt('Energy weighted fragment dens'
     &                                  ,' ',Array,nFragDens(iCurCnttp))
* include the minus sign of -2eta_i
            Call DScal_(nFragDens(iCurCnttp)*(nFragDens(iCurCnttp)+1)/2,
     &                 -One,Array,1)
            If(iPrint.ge.49)
     &                    Call TriPrt('-1 Energy weighted fragment dens'
     &                                  ,' ',Array,nFragDens(iCurCnttp))
            If(maxDensSize.lt.nFragDens(iCurCnttp)*
     &         (nFragDens(iCurCnttp)+1)/2) Call Abend !'maxIJSize'
          End If
        End If
        call dcopy_(3,Work(ixyz),1,C,1)
c       write(*,*) '  iShll,iAng,mdci,iShell,iCnttp,iCurMdc,iCurCnttp',
c    &              iShll,iAng,mdci,iShell,iCnttp,iCurMdc,iCurCnttp
c       write(*,*) '  iPrim,iBas =',iPrim,iBas
*                                                                      *
************************************************************************
*                                                                      *
* Loop over all other shells belonging to the same fragment
        jSbasis = 1
      If (iPrint.ge.49) Then
      write(6,'(3(A,i4))') 'iS = ',iS,' iSstart=',iSstart,' iSEnd=',
     &                      iSend
      End If
        Do jS = iSstart, iSend
          jShll  = iSD( 0,jS)
          jAng   = iSD( 1,jS)
          jCmp   = iSD( 2,jS)
          jBas   = iSD( 3,jS)
          jCff   = iSD( 4,jS)
          jPrim  = iSD( 5,jS)
          jExp   = iSD( 6,jS)
          jAO    = iSD( 7,iS)
          jxyz   = iSD( 8,jS)
          mdcj   = iSD(10,jS)
          jShell = iSD(11,jS)
          jCnttp = iSD(13,jS)
          jSize = nElem(jAng)
      If (iPrint.ge.49) Then
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jShll =',jShll
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jAng  =',jAng
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jCmp  =',jCmp
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jBas  =',jBas
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jCff  =',jCff
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jPrim =',jPrim
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jExp  =',jExp
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jAO   =',jAO
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jxyz  =',jxyz
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' mdcj  =',mdcj
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jShell=',jShell
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jCnttp=',jCnttp
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jSize =',jSize
       End If

          if(Transf(jShll).and.Prjct(jShll)) jSize = 2*jAng+1
          call dcopy_(3,Work(jxyz),1,B,1)
c         write(*,*) '    jShll,jAng,mdcj,jShell,jCnttp =',
c     &                    jShll,jAng,mdcj,jShell,jCnttp
c         write(*,*) '    jPrim,jBas =',jPrim,jBas
*                                                                      *
************************************************************************
*                                                                      *
* Create a rectangular matrix sized (iBas*nElem(iAng),jBas*nElem(jAng))
* from the energy weighted density matrix (desymmetrized)
* contains values from iSbasis to iSbasis + iBas*nElem(iAng) - 1
*             and from jSbasis to jSbasis + jBas*nElem(jAng) - 1
           ipIJ = 1 + maxDensSize
c          write(6,*) '    ipIJ=',ipIJ
c          write(6,*) '    extracting values from',iSbasis,' to',iSbasis
c     &                + iBas*iSize - 1,', and from',jSbasis,' to',
c     &                jSbasis + jBas*jSize - 1
           Do iSlocal = iSbasis, iSbasis + iBas*iSize - 1
             Do jSlocal = jSbasis, jSbasis + jBas*jSize - 1
               iLoc = ipIJ + (jSlocal-jSbasis)*iBas*iSize + iSlocal
     &                - iSbasis
c               write(6,'(A,i3,A,i3,A,i4,A,i8)')
c     &         'iTri(',iSlocal,',',jSlocal,')=',iTri(iSlocal,jSlocal),
c     &        ' iLoc=',iLoc
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
          Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,
     &             jStab(0,mdci),nStab(mdci),iDCRT,nDCRT)
          Fact = DBLE(nStabM) / DBLE(LmbdT)
*                                                                      *
************************************************************************
*                                                                      *
* Loop over symmetry operations acting on the basis.
*
          Do lDCRT = 0, nDCRT-1
            TC(1) = iPhase(1,iDCRT(lDCRT))*C(1)
            TC(2) = iPhase(2,iDCRT(lDCRT))*C(2)
            TC(3) = iPhase(3,iDCRT(lDCRT))*C(3)
            TB(1) = iPhase(1,iDCRT(lDCRT))*B(1)
            TB(2) = iPhase(2,iDCRT(lDCRT))*B(2)
            TB(3) = iPhase(3,iDCRT(lDCRT))*B(3)

*                                                                      *
************************************************************************
*                                                                      *
* Calculate the overlap integral < alpha | is >
*
            ip = ipIJ + maxDensSize
            ipF1 = ip
            nac = nElem(la)*nElem(iAng)
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
               Write (6,*) '  ip-1.gt.nArr*nZeta(1) in FragPInt'
               Call Abend
            End If
            mArr = (nArr*nZeta-(ip-1))/nZeta
*
            Call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,iPrim,
     &                Alpha,Work(iExp))
            Call SetUp1(Alpha,nAlpha,Work(iExp),iPrim,
     &                  A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
            nHer = (la+iAng+2)/2
            Call MltPrm(Alpha,nAlpha,Work(iExp),iPrim,
     &                Array(ipZ1),Array(ipZI1),
     &                Array(ipK1),Array(ipP1),
     &                Array(ipF1),nAlpha*iPrim,iComp,
     &              la,iAng,A,TC,nHer,Array(ip),
     &              mArr,CCoor,nOrdOp)
            If(iPrint.ge.99) Call RecPrt('<alpha|iS> (aBas x X)',
     &        ' ',Array(ipF1),nAlpha*iPrim,nac)
*                                                                      *
************************************************************************
*                                                                      *
* Calculate the overlap integral < jS | beta >
*
            ip = ip - 6 * nAlpha*iPrim
            ipF2 = ip
            ncb = nElem(jAng)*nElem(lb)
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
               Write (6,*) '  ip-1.gt.nArr*nZeta(2) in FragPInt'
               Call Abend
            End If
            mArr = (nArr*nZeta-(ip-1))/nZeta
*
            Call ZXia(Array(ipZ2),Array(ipZI2),jPrim,nBeta,
     &                Work(jExp),Beta)
            Call SetUp1(Work(jExp),jPrim,Beta,nBeta,
     &                 TB,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
            nHer = (jAng+lb+2)/2
            Call MltPrm(Work(jExp),jPrim,Beta,nBeta,
     &                Array(ipZ2),Array(ipZI2),
     &                Array(ipK2),Array(ipP2),
     &                Array(ipF2),jPrim*nBeta,iComp,
     &                jAng,lb,TB,RB,nHer,Array(ip),
     &                mArr,CCoor,nOrdOp)
            ip = ip - 6 * jPrim*nBeta
            ipTmp = ip
            ip = ip + Max(nAlpha*nac*Max(iPrim,jBas),
     &                    nBeta*ncb*jBas)
            If (ip-1.gt.nArr*nZeta) Then
               Write (6,*) '  ip-1.gt.nArr*nZeta(3) in FragPInt'
               Call Abend
            End If
            If(iPrint.ge.99) Call RecPrt('<jS|beta> (bBas x Y)',
     &        ' ',Array(ipF2),nBeta*jPrim,ncb)
            If((ipTmp-ipF2).lt.nBeta*jPrim*ncb) stop 'sizetest 1'
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
     &                  iPrim*nac,Array(ipTmp),iPrim*nac)
            If(iPrint.ge.99) Call RecPrt('<alpha|iS>^T (X x aBas)',
     &        ' ',Array(ipTmp),iPrim*nac,nAlpha)
            If((ip-ipTmp).lt.nAlpha*iPrim*nac) stop 'sizetest 2'
*
*-----------2) aciK =  k,aci * k,K
*
            Call DGEMM_('T','N', nac*nAlpha, iBas, iPrim,
     &                  1.0d0, Array(ipTmp),       iPrim,
     &                          Work(iCff ),       iPrim,
     &                  0.0d0, Array(ipF1 ),  nAlpha*nac  )
            If(iPrint.ge.99) Then
              Call RecPrt('<alpha|iS>(regrouped, X x iPrim)',' ',
     &                    Array(ipTmp),nAlpha*nac,iPrim)
              Call RecPrt('Coeffs of iS (iPrim x iBas)',' ',Work(iCff),
     &                    iPrim,iBas)
             Call RecPrt('<alpha|iS>(re) * Coeffs of iS (X x iBas)',' ',
     &                   Array(ipF1),nAlpha*nac,iBas)
            End If
            If((ipF2-ipF1).lt.nAlpha*iBas*nac) stop 'sizetest 3'
*
*-----------4) a,ciK -> ciKa
*
            Call DgeTMo(Array(ipF1),nElem(la),nElem(la),
     &                  nElem(iAng)*nAlpha*iBas,
     &                  Array(ipTmp),
     &                  nElem(iAng)*nAlpha*iBas)
            If(iPrint.ge.99) Then
              Call RecPrt('result (regrouped, nElem(la) x X)',' ',
     &                  Array(ipF1),nElem(la),nElem(iAng)*nAlpha*iBas)
              Call RecPrt('transpose of result (X x nElem(la))',' ',
     &                  Array(ipTmp),nElem(iAng)*nAlpha*iBas,nElem(la))
            End If
            If((ip-ipTmp).lt.nAlpha*iBas*nElem(iAng)*nElem(la))
     &        stop 'sizetest 6'
*
*-----------5) iKa,C = c,iKa * c,C
*
            If(Transf(iShll).and.Prjct(iShll)) Then
              Call DGEMM_('T','N',
     &                    iBas*nElem(la)*nAlpha,iSize,nElem(iAng),
     &                    1.0d0,Array(ipTmp),nElem(iAng),
     &                    RSph(ipSph(iAng)),nElem(iAng),
     &                    0.0d0,Array(ipF1),nAlpha*iBas*nElem(la))
              If(iPrint.ge.99) Then
                Call RecPrt('result (regrouped, X x nElem(iAng))',' ',
     &                   Array(ipTmp),nElem(la)*nAlpha*iBas,nElem(iAng))
                Call RecPrt('Spher of iS (nElem(iAng) x (2*iAng+1))',' '
     &                      ,RSph(ipSph(iAng)),nElem(iAng),(2*iAng+1))
                Call RecPrt('result in spherical gaussians (X x iSize',
     &                    ' ',Array(ipF1),nAlpha*iBas*nElem(la),iSize)
              End If
            Else
* in this case nElem(iAng) = iSize
              Call DgeTMo(Array(ipTmp),nElem(iAng),nElem(iAng),
     &                    iBas*nElem(la)*nAlpha,
     &                    Array(ipF1),
     &                    iBas*nElem(la)*nAlpha)
            End If
            If((ipF2-ipF1).lt.nAlpha*iBas*nElem(la)*iSize)
     &        stop 'sizetest 7'
*
*-----------And (almost) the same thing for the righthand side, form
*           LjDb from ljdb
*           1) jdb,L = l,jdb * l,L
            Call DGEMM_('T','N',
     &                  nBeta*ncb,jBas,jPrim,
     &                  1.0d0,Array(ipF2),jPrim,
     &                  Work(jCff),jPrim,
     &                  0.0d0,Array(ipTmp),nBeta*ncb)
            If(iPrint.ge.99) Then
              Call RecPrt('<jS|beta>(regrouped, X x jPrim)',' ',
     &                    Array(ipF2),nBeta*ncb,jPrim)
              Call RecPrt('Coeffs of jS (jPrim x jBas)',' ',Work(jCff),
     &                    jPrim,jBas)
              Call RecPrt('<jS|beta>(re) * Coeffs of jS (Y x jBas)',' ',
     &                   Array(ipTmp),nBeta*ncb,jBas)

            End If
            If((ip-ipTmp).lt.nBeta*jBas*ncb) stop 'sizetest 8'
*
*-----------2)  j,dbL -> dbL,j
*
            Call DgeTMo(Array(ipTmp),nBeta,nBeta,
     &                  jBas*ncb,Array(ipF2),
     &                  jBas*ncb)
            If(iPrint.ge.99) Call RecPrt('transposed right 1 (Y x bBas)'
     &        ,' ',Array(ipF2),jBas*ncb,nBeta)
            If((ipTmp-ipF2).lt.nBeta*jBas*ncb) stop 'sizetest 9'
*
*-----------3) bLj,D = d,bLj * d,D
*
            If(Transf(jShll).and.Prjct(jShll)) Then
              Call DGEMM_('T','N',
     &                    nElem(lb)*jBas*nBeta,jSize,nElem(jAng),
     &                    1.0d0,Array(ipF2),nElem(jAng),
     &                    RSph(ipSph(jAng)),nElem(jAng),
     &                    0.0d0,Array(ipTmp),nElem(lb)*jBas*nBeta)
              If(iPrint.ge.99)
     &          Call RecPrt('multiply right 2 (Y x jSize)',' ',
     &                      Array(ipTmp),nBeta*jBas*nElem(lb),jSize)
            Else
* in this case nElem(jAng) = jSize
              Call DgeTMo(Array(ipF2),nElem(jAng),nElem(jAng),
     &                    jBas*nElem(lb)*nBeta,
     &                    Array(ipTmp),
     &                    jBas*nElem(lb)*nBeta)
            End If
            If((ip-ipTmp).lt.nBeta*jBas*nElem(lb)*jSize)
     &        stop 'sizetest 10'
*
*-----------4) b,LjD -> LjD,b
*
            Call DgeTMo(Array(ipTmp),nElem(lb),nElem(lb),
     &                  jBas*nBeta*jSize,Array(ipF2),
     &                  jBas*nBeta*jSize)
            If(iPrint.ge.99) Then
              Call RecPrt('transposed right 2 (Y x nElem(lb)',' ',
     &                    Array(ipF2),jBas*nBeta*jSize,nElem(lb))
            End If

            If((ipTmp-ipF2).lt.nBeta*jBas*nElem(lb)*jSize)
     &        stop 'sizetest 11'

*
*-----------Next Contract (iKaC)*W(KLCD)*(LjDb) producing ijab,
*           by the following procedure:
*           Loop over a and b
*             Loop over C
*               Contract iK(aC)*Kj(Cb), over K producing ij(aCb),
*                 accumulate to ij(ab)
*             End loop C
*           End Loop b and a
*
* Total size of ipF1 = nAlpha*nElem(la) * iBas*iSize  -> ordered as (nAlpha, iBas,  nElem(la), iSize)
*               ipF2 = nBeta*nElem(lb)  * jBas*jSize                (jBas,   nBeta, jSize,     nElem(lb))
*                  W = iBas*iSize * jBas*jSize                      (iBas,   iSize, jBas,      jSize)
*
            If(iPrint.ge.100) Then
              Write (6,*) ' Current contents of Final():'
              Do ia = 1, nElem(la)
                Do ib = 1, nElem(lb)
                  Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
                  Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
                End Do
              End Do
            End If

            iIC = 0
            Do iIrrep = 0, nIrrep - 1
              If (iAnd(llOper,iTwoj(iIrrep)).ne.0) Then
                iIC = iIC + 1
                nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
                Xg=rChTbl(iIrrep,nOp)
* Half is needed because we do a complete loop over iS,jS
                Factor=Xg*Fact*Half
                ! write(6,'(A,i24)') 'FragPInt:  ipIJ=', ipIJ
                ! print *,'CALL FragPCont'
                Call xFlush(6)

                Call FragPCont(Array(ipF1), nAlpha,iBas,nElem(la),iSize,
     &                         Array(ipF2), jBas,nBeta,jSize,nElem(lb),
c     &                         Array(ipIJ), Final(1,1,1,iIC), Factor)
     &                         Array(ipIJ), Final(:,:,:,iIC), Factor)
              End If
            End Do
            If(iIC.ne.nIC) stop 'iIC.ne.nIC'

            If(iPrint.ge.99) Then
              Write (6,*) ' After contraction:'
              Do ia = 1, nElem(la)
                Do ib = 1, nElem(lb)
                  Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
                  Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
                End Do
              End Do
            End If
*
          End Do ! end loop over lDCR
*
          jSbasis = jSbasis + jBas * jSize
        End Do ! end loop over jS
        Call xFlush(6)

        iSbasis = iSbasis + iBas * iSize
      End Do ! end loop over iS
      Call xFlush(6)
*
      If (iPrint.ge.49) Then
         Write (6,*) ' Result in FragPInt'
         Do ia = 1, nElem(la)
            Do ib = 1, nElem(lb)
               Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
               Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
               Call xFlush(6)
c               Call Add_Info(Label,Final(1,ia,ib,1),nAlpha*nBeta,5)
            End Do
         End Do
      End If
c add some verification data
      ! this check is to ensure that Add_Info is called only on the Master node
      ! The reason being that the execution of this Kernel routine is split on
      ! diffeent nodes eariler in the code.
      ! OneEl -> OneEl_ -> OneEl_IJ -> this routine.
      ! Normally, the Add_Info must be called after the parallelization is finalized,
      ! i.e. in the OneEl function.
c      If (MyRank.eq.0) then
c         Do ia = 1, nElem(la)
c            Do ib = 1, nElem(lb)
c            dA=0.d0
c            dA=dnrm2_(nAlpha*nBeta,Final(1,ia,ib,1),1)
c               If (dA .gt. 1.d-6 ) then
c                  write(label,'(A,i2,A,i2)') 'Fragpint: ',ia,' ib ',ib
c                  Call Add_Info(label, dA, 1, 6)
c               End if
c               Call xFlush(6)
c            End Do
c         End Do
c      End If
      Call Free_iSD
      Call Set_Basis_Mode('Valence')
      Call SetUp_iSD
      Call Nr_Shells(nSkal)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Zeta)
        Call Unused_real_array(ZInv)
        Call Unused_real_array(rKappa)
        Call Unused_real_array(P)
        Call Unused_integer(nRys)
        Call Unused_integer_array(iCho)
        Call Unused_real(PtChrg)
        Call Unused_integer(nGrid)
        Call Unused_integer(iAddPot)
      End If
      End

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
      SubRoutine FragPInt(
#define _CALLING_
#include "int_interface.fh"
     &                    )
************************************************************************
*                                                                      *
* Object: kernel routine for the computation of Fragment AIEMP         *
*         projection integrals                                         *
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
      use Basis_Info
      use Center_Info
      use Symmetry_Info, only: nIrrep, iChTbl
      Implicit None
#include "real.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"

#include "int_interface.fh"

*     Local variables
      Real*8  C(3), TC(3), B(3), TB(3)
      Integer iDCRT(0:7),iTwoj(0:7)
      Logical EnergyWeight
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*24 Label
      Integer ia, ib
#endif
      Data    iTwoj/1,2,4,8,16,32,64,128/
*
      Integer i,j,ixyz,nElem,iTri,
     &        iAng,iBas,iCnttp,iComp,
     &        iCurCenter,iCurCnttp,iCurMdc,iIC,iIrrep,iLoc,iPrim,
     &        ip,ipF1,ipF2,ipIJ,ipK1,ipK2,ipP1,ipP2,ipTmp,ipZ1,ipZ2,
     &        ipZI1,ipZI2,iS,iSbasis,iSend,iShll,iSize,iSlocal,
     &        iSstart,iStemp,jAng,jBas,jCnttp,jPrim,
     &        jS,jShll,jSize,jSlocal,lDCRT,llOper,LmbdT,
     &        mArr,maxDensSize,mdci,nac,ncb,nDCRT,nOp,nSkal,
     &        jSbasis,iCnt,jCnt
      Real*8  Fact,Factor,Xg
* external functions:
      Integer NrOpr
      External NrOpr
c      Real*8  DNRM2_
c      External DNRM2_
*
*     Statement functions
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
*
#ifdef _DEBUGPRINT_
c     data for individual fragments:
      Write (6,*) ' In FragPInt:    nCnttp          = ',nCnttp
      Write (6,*) ' In FragPInt: dbsc(nCnttp)%mdci  = ',
     &                                 (dbsc(i)%mdci,i=1,nCnttp)
      Write (6,*) ' In FragPInt:     dbsc(nCnttp)%nCntr  = ',
     &                                     (dbsc(i)%nCntr,i=1,nCnttp)
      Write (6,*) ' In FragPInt: nFragType(nCnttp)  = ',
     &                                 (dbsc(i)%nFragType,i=1,nCnttp)
      Write (6,*) ' In FragPInt: nFragCoor(nCnttp)  = ',
     &                                 (dbsc(i)%nFragCoor,i=1,nCnttp)
      Write (6,*) ' In FragPInt: nFragDens(nCnttp)  = ',
     &                                 (dbsc(i)%nFragDens,i=1,nCnttp)
      Write (6,*) ' In FragPInt: nFragDens(nCnttp)  = ',
     &                                 (dbsc(i)%nFragDens,i=1,nCnttp)
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
      Call RecPrt(' In FragPInt: A     ',' ',A     ,1,3)
      Call RecPrt(' In FragPInt: RB    ',' ',RB    ,1,3)
      Call RecPrt(' In FragPInt: Ccoor ',' ',Ccoor ,1,3)
      Call RecPrt(' In FragPInt: P     ',' ',P     ,nZeta,3)
      Call RecPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
      Call TrcPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
#endif
*
      Call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*                                                                      *
************************************************************************
*                                                                      *
* Setup the fragment shells
*
      Call Free_iSD
      Call Set_Basis_Mode('Fragments')
      Call SetUp_iSD
      Call Nr_Shells(nSkal)
#ifdef _DEBUGPRINT_
      Write(6,*)'nSkal_Fragment,nAlpha,nBeta = ',nSkal,nAlpha,nBeta
#endif
*                                                                      *
************************************************************************
*                                                                      *
* Reserve space for the largest possible fragment energy weighted
* density matrix and one combination of shells of it.
      maxDensSize = 0
      Do iCnttp = 1, nCnttp
        If(dbsc(iCnttp)%nFragType.gt.0) then
        maxDensSize = Max( maxDensSize,dbsc(iCnttp)%nFragDens
     &                               *(dbsc(iCnttp)%nFragDens+1)/2)
#ifdef _DEBUGPRINT_
        write(6,'(A,i2,A,i6)') 'nFragDens(',iCnttp,')=',
     &                                           dbsc(iCnttp)%nFragDens
#endif
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
c     ! The mdc value of the current fragment placeholder
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
        iBas   = iSD( 3,iS)
        iPrim  = iSD( 5,iS)
        mdci   = iSD(10,iS)
        iCnttp = iSD(13,iS)
        iCnt   = iSD(14,iS)
        iSize = nElem(iAng)
        C(1:3)=dbsc(iCnttp)%Coor(1:3,iCnt)
c some printouts:
#ifdef _DEBUGPRINT_
        Write(6,*) 'In FragPInt: iS=',iS,' iShll =',iShll
        Write(6,*) 'In FragPInt: iS=',iS,' iAng  =',iAng
        Write(6,*) 'In FragPInt: iS=',iS,' iBas  =',iBas
        Write(6,*) 'In FragPInt: iS=',iS,' iPrim =',iPrim
        Write(6,*) 'In FragPInt: iS=',iS,' ixyz  =',ixyz
        Write(6,*) 'In FragPInt: iS=',iS,' mdci  =',mdci
        Write(6,*) 'In FragPInt: iS=',iS,' iCnttp=',iCnttp
        Write(6,*) 'In FragPInt: iS=',iS,' iSize =',iSize
        Write(6,*) 'In FragPInt: iS=',iS,' iCurMdc =',iCurMdc
        Write(6,*) 'In FragPInt: nFragCoor(',iCnttp,') =',
     &            dbsc(iCnttp)%nFragCoor
#endif

        If(Shells(iShll)%Transf.and.
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

#ifdef _DEBUGPRINT_
      write(6,*) 'start of new fragment encountered'
      write(6,*) 'iSstart,iSend,iCurCnttp,iCurCenter =',
     &            iSstart,iSend,iCurCnttp,iCurCenter
#endif

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
#ifdef _DEBUGPRINT_
            Call TriPrt('Energy weighted fragment dens',' ',
     &                   Array,dbsc(iCurCnttp)%nFragDens)
#endif
* include the minus sign of -2eta_i
            Call DScal_(dbsc(iCurCnttp)%nFragDens
     &                *(dbsc(iCurCnttp)%nFragDens+1)/2,-One,Array,1)
#ifdef _DEBUGPRINT_
            Call TriPrt('-1 Energy weighted fragment dens',' ',
     &                  Array,dbsc(iCurCnttp)%nFragDens)
#endif
            If (maxDensSize.lt.dbsc(iCurCnttp)%nFragDens*
     &                        (dbsc(iCurCnttp)%nFragDens+1)/2)
     &         Call Abend !'maxIJSize'
          End If
        End If
#ifdef _DEBUGPRINT_
       write(6,*) '  iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp',
     &              iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp
       write(6,*) '  iPrim,iBas =',iPrim,iBas
#endif
*                                                                      *
************************************************************************
*                                                                      *
* Loop over all other shells belonging to the same fragment
        jSbasis = 1
#ifdef _DEBUGPRINT_
      write(6,'(3(A,i4))') 'iS = ',iS,' iSstart=',iSstart,' iSEnd=',
     &                      iSend
#endif
        Do jS = iSstart, iSend
          jShll  = iSD( 0,jS)
          jAng   = iSD( 1,jS)
          jBas   = iSD( 3,jS)
          jPrim  = iSD( 5,jS)
          jCnttp = iSD(13,jS)
          jCnt   = iSD(14,jS)
          jSize = nElem(jAng)
          B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
#ifdef _DEBUGPRINT_
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jShll =',jShll
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jAng  =',jAng
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jBas  =',jBas
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jPrim =',jPrim
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jCnttp=',jCnttp
        write(6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jSize =',jSize
#endif

          if(Shells(jShll)%Transf.and.
     &       Shells(jShll)%Prjct ) jSize = 2*jAng+1
#ifdef _DEBUGPRINT_
         write(6,*) '    jShll,jAng,jCnttp =',
     &                    jShll,jAng,jCnttp
         write(6,*) '    jPrim,jBas =',jPrim,jBas
#endif
*                                                                      *
************************************************************************
*                                                                      *
* Create a rectangular matrix sized (iBas*nElem(iAng),jBas*nElem(jAng))
* from the energy weighted density matrix (desymmetrized)
* contains values from iSbasis to iSbasis + iBas*nElem(iAng) - 1
*             and from jSbasis to jSbasis + jBas*nElem(jAng) - 1
           ipIJ = 1 + maxDensSize
#ifdef _DEBUGPRINT_
           write(6,*) '    ipIJ=',ipIJ
           write(6,*) '    extracting values from',iSbasis,' to',iSbasis
     &                + iBas*iSize - 1,', and from',jSbasis,' to',
     &                jSbasis + jBas*jSize - 1
#endif
           Do iSlocal = iSbasis, iSbasis + iBas*iSize - 1
             Do jSlocal = jSbasis, jSbasis + jBas*jSize - 1
               iLoc = ipIJ + (jSlocal-jSbasis)*iBas*iSize + iSlocal
     &                - iSbasis
#ifdef _DEBUGPRINT_
                write(6,'(A,i3,A,i3,A,i4,A,i8)')
     &         'iTri(',iSlocal,',',jSlocal,')=',iTri(iSlocal,jSlocal),
     &        ' iLoc=',iLoc
#endif
               Array(iLoc) = Array(iTri(iSlocal,jSlocal))
               If(iSlocal.ne.jSlocal) Array(iLoc) = Array(iLoc)/Two
#ifdef _DEBUGPRINT_
               write(6,*) 'Filling (',iSlocal-iSbasis+1,',',
     &           jSlocal-jSbasis+1,') from (',iSlocal,',',jSlocal,')'
#endif
             End Do
           End Do
#ifdef _DEBUGPRINT_
           Call RecPrt('W(KC,LD)',' ',Array(ipIJ),iBas*iSize,jBas*jSize)
#endif
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
          Do lDCRT = 0, nDCRT-1
            Call OA(iDCRT(lDCRT),C,TC)
            Call OA(iDCRT(lDCRT),B,TB)

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
     &                Alpha,Shells(iShll)%Exp)
            Call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,
     &                  A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))
*
            nHer = (la+iAng+2)/2
            Call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,
     &                Array(ipZ1),Array(ipZI1),
     &                Array(ipK1),Array(ipP1),
     &                Array(ipF1),nAlpha*iPrim,iComp,
     &              la,iAng,A,TC,nHer,Array(ip),
     &              mArr,CCoor,nOrdOp)
#ifdef _DEBUGPRINT_
            Call RecPrt('<alpha|iS> (aBas x X)',
     &        ' ',Array(ipF1),nAlpha*iPrim,nac)
#endif
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
     &                Shells(jShll)%Exp,Beta)
            Call SetUp1(Shells(jShll)%Exp,jPrim,Beta,nBeta,
     &                 TB,RB,Array(ipK2),Array(ipP2),Array(ipZI2))
*
            nHer = (jAng+lb+2)/2
            Call MltPrm(Shells(jShll)%Exp,jPrim,Beta,nBeta,
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
#ifdef _DEBUGPRINT_
            Call RecPrt('<jS|beta> (bBas x Y)',
     &        ' ',Array(ipF2),nBeta*jPrim,ncb)
#endif
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
#ifdef _DEBUGPRINT_
            Call RecPrt('<alpha|iS>^T (X x aBas)',
     &        ' ',Array(ipTmp),iPrim*nac,nAlpha)
#endif
            If((ip-ipTmp).lt.nAlpha*iPrim*nac) stop 'sizetest 2'
*
*-----------2) aciK =  k,aci * k,K
*
            Call DGEMM_('T','N', nac*nAlpha, iBas, iPrim,
     &                  1.0d0, Array(ipTmp),       iPrim,
     &                         Shells(iShll)%pCff, iPrim,
     &                  0.0d0, Array(ipF1 ),  nAlpha*nac  )
#ifdef _DEBUGPRINT_
             Call RecPrt('<alpha|iS>(regrouped, X x iPrim)',' ',
     &                    Array(ipTmp),nAlpha*nac,iPrim)
             Call RecPrt('Coeffs of iS (iPrim x iBas)',' ',
     &                    Shells(iShll)%pCff,iPrim,iBas)
             Call RecPrt('<alpha|iS>(re) * Coeffs of iS (X x iBas)',' ',
     &                   Array(ipF1),nAlpha*nac,iBas)
#endif
            If((ipF2-ipF1).lt.nAlpha*iBas*nac) stop 'sizetest 3'
*
*-----------4) a,ciK -> ciKa
*
            Call DgeTMo(Array(ipF1),nElem(la),nElem(la),
     &                  nElem(iAng)*nAlpha*iBas,
     &                  Array(ipTmp),
     &                  nElem(iAng)*nAlpha*iBas)
#ifdef _DEBUGPRINT_
            Call RecPrt('result (regrouped, nElem(la) x X)',' ',
     &                Array(ipF1),nElem(la),nElem(iAng)*nAlpha*iBas)
            Call RecPrt('transpose of result (X x nElem(la))',' ',
     &                Array(ipTmp),nElem(iAng)*nAlpha*iBas,nElem(la))
#endif
            If((ip-ipTmp).lt.nAlpha*iBas*nElem(iAng)*nElem(la))
     &        stop 'sizetest 6'
*
*-----------5) iKa,C = c,iKa * c,C
*
            If(Shells(iShll)%Transf.and.
     &         Shells(iShll)%Prjct ) Then
              Call DGEMM_('T','N',
     &                    iBas*nElem(la)*nAlpha,iSize,nElem(iAng),
     &                    1.0d0,Array(ipTmp),nElem(iAng),
     &                    RSph(ipSph(iAng)),nElem(iAng),
     &                    0.0d0,Array(ipF1),nAlpha*iBas*nElem(la))
#ifdef _DEBUGPRINT_
              Call RecPrt('result (regrouped, X x nElem(iAng))',' ',
     &                 Array(ipTmp),nElem(la)*nAlpha*iBas,nElem(iAng))
              Call RecPrt('Spher of iS (nElem(iAng) x (2*iAng+1))',' '
     &                    ,RSph(ipSph(iAng)),nElem(iAng),(2*iAng+1))
              Call RecPrt('result in spherical gaussians (X x iSize',
     &                  ' ',Array(ipF1),nAlpha*iBas*nElem(la),iSize)
#endif
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
     &                        Shells(jShll)%pCff,jPrim,
     &                  0.0d0,Array(ipTmp),nBeta*ncb)
#ifdef _DEBUGPRINT_
            Call RecPrt('<jS|beta>(regrouped, X x jPrim)',' ',
     &                  Array(ipF2),nBeta*ncb,jPrim)
            Call RecPrt('Coeffs of jS (jPrim x jBas)',' ',
     &                   Shells(jShll)%pCff,jPrim,jBas)
            Call RecPrt('<jS|beta>(re) * Coeffs of jS (Y x jBas)',' ',
     &                 Array(ipTmp),nBeta*ncb,jBas)
#endif
            If((ip-ipTmp).lt.nBeta*jBas*ncb) stop 'sizetest 8'
*
*-----------2)  j,dbL -> dbL,j
*
            Call DgeTMo(Array(ipTmp),nBeta,nBeta,
     &                  jBas*ncb,Array(ipF2),
     &                  jBas*ncb)
#ifdef _DEBUGPRINT_
            Call RecPrt('transposed right 1 (Y x bBas)',
     &         ' ',Array(ipF2),jBas*ncb,nBeta)
#endif
            If((ipTmp-ipF2).lt.nBeta*jBas*ncb) stop 'sizetest 9'
*
*-----------3) bLj,D = d,bLj * d,D
*
            If(Shells(jShll)%Transf.and.
     &         Shells(jShll)%Prjct ) Then
              Call DGEMM_('T','N',
     &                    nElem(lb)*jBas*nBeta,jSize,nElem(jAng),
     &                    1.0d0,Array(ipF2),nElem(jAng),
     &                    RSph(ipSph(jAng)),nElem(jAng),
     &                    0.0d0,Array(ipTmp),nElem(lb)*jBas*nBeta)
#ifdef _DEBUGPRINT_
              Call RecPrt('multiply right 2 (Y x jSize)',' ',
     &                    Array(ipTmp),nBeta*jBas*nElem(lb),jSize)
#endif
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
#ifdef _DEBUGPRINT_
            Call RecPrt('transposed right 2 (Y x nElem(lb)',' ',
     &                  Array(ipF2),jBas*nBeta*jSize,nElem(lb))
#endif

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
#ifdef _DEBUGPRINT_
            Write (6,*) ' Current contents of Final():'
            Do ia = 1, nElem(la)
              Do ib = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
                Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
              End Do
            End Do
#endif

            iIC = 0
            Do iIrrep = 0, nIrrep - 1
              If (iAnd(llOper,iTwoj(iIrrep)).ne.0) Then
                iIC = iIC + 1
                nOp = NrOpr(iDCRT(lDCRT))
                Xg=DBLE(iChTbl(iIrrep,nOp))
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

#ifdef _DEBUGPRINT_
            Write (6,*) ' After contraction:'
            Do ia = 1, nElem(la)
              Do ib = 1, nElem(lb)
                Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
                Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
              End Do
            End Do
#endif
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
#ifdef _DEBUGPRINT_
      Write (6,*) ' Result in FragPInt'
      Do ia = 1, nElem(la)
         Do ib = 1, nElem(lb)
            Write (Label,'(A,I2,A,I2,A)') ' Final(',ia,',',ib,')'
            Call RecPrt(Label,' ',Final(:,ia,ib,:),nAlpha,nBeta)
         End Do
      End Do
#endif
c add some verification data
!      this check is to ensure that Add_Info is called only on the Master node
!      The reason being that the execution of this Kernel routine is split on
!      diffeent nodes eariler in the code.
!      OneEl -> OneEl_ -> OneEl_IJ -> this routine.
!      Normally, the Add_Info must be called after the parallelization is finalized,
!      i.e. in the OneEl function.
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
        Call Unused_integer(nHer)
        Call Unused_integer_array(iCho)
        Call Unused_real_array(PtChrg)
        Call Unused_integer(iAddPot)
      End If
      End

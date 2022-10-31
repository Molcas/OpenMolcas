!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Ben Swerts                                             *
!***********************************************************************

subroutine FragPGrd( &
#                   define _CALLING_
#                   include "grd_interface.fh"
                   )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of FAIEMP Projection      *
!         operator integrals.                                          *
!                                                                      *
!     Author: Ben Swerts                                               *
!                                                                      *
!     based on PrjGrd                                                  *
!***********************************************************************

use Her_RW, only: HerR, HerW, iHerR, iHerW
use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iOper
use Index_Functions, only: iTri, nTri_Elem1
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "grd_interface.fh"
#include "print.fh"
real(kind=wp) :: C(3), TC(3), B(3), TB(3), Fact
integer(kind=iwp) :: i, ia, iAng, ib, iBas, iCar, iCent, iCnt, iCnttp, iCurCenter, iCurCnttp, iCurMdc, iDCRT(0:7), iGamma, iLoc, &
                     ip, ipA, ipAxyz, ipB, ipBxyz, ipCxyz, ipF1, ipF1a, ipF2, ipF2a, ipIJ, ipK1, ipK2, ipP1, ipP2, ipQ1, iPrim, &
                     iPrint, ipRxyz, ipTmp, ipZ1, ipZ2, ipZI1, ipZI2, iRout, iS, iSbasis, iSEnd, iShll, iSize, iSlocal, iSstart, &
                     iStemp, iStrt, iuvwx(4), iVec, jAng, jBas, jCnt, jCnttp, JndGrd(3,4), jPrim, jS, jSbasis, jShll, jSize, &
                     jSlocal, ld, lDCRT, LmbdT, lOp(4), maxDensSize, mdci, mGrad, mVec, mVecAC, mVecCB, n_Her, nac, ncb, nDAO, &
                     nDCRT, nSkal, nVecAC, nVecCB
logical(kind=iwp) :: JfGrad(3,4), ABeq(3), EQ, EnergyWeight
character(len=80) :: Label
real(kind=wp), external :: DNrm2_
#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(nHer)

iRout = 202
iPrint = nPrint(iRout)

if (iPrint >= 49) then
  call RecPrt(' In FragPGrd: Grad',' ',Grad,1,nGrad)
  call RecPrt(' In FragPGrd: A',' ',A,1,3)
  call RecPrt(' In FragPGrd: RB',' ',RB,1,3)
  call RecPrt(' In FragPGrd: Ccoor',' ',Ccoor,1,3)
  call RecPrt(' In FragPGrd: P',' ',P,nZeta,3)
  call RecPrt(' In FragPGrd: Alpha',' ',Alpha,nAlpha,1)
  call RecPrt(' In FragPGrd: Beta',' ',Beta,nBeta,1)
  write(u6,*) ' In FragPGrd: la,lb=',' ',la,lb
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Setup the fragment shells

call Set_Basis_Mode('Fragments')
call SetUp_iSD()
call Nr_Shells(nSkal)
if (iPrint >= 99) then
  write(u6,*) 'looping over ',nSkal,' shells'
  write(u6,*) 'Shells()%Frag = ',(Shells(i)%Frag,i=1,10)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Reserve space for the largest possible fragment energy weighted
! density matrix
maxDensSize = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%nFragType > 0) maxDensSize = max(maxDensSize,dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over all shells belonging to the fragments
nDAO = nTri_Elem1(la)*nTri_Elem1(lb)
iuvwx(1) = dc(mdc)%nStab
iuvwx(2) = dc(ndc)%nStab
lOp(1) = iOper(kOp(1))
lOp(2) = iOper(kOp(2))

iCurMdc = 0         ! The mdc value of the current fragment placeholder
iCurCnttp = 0       ! The Cnttp of the fragment placeholder
iCurCenter = 999999 ! The index of the fragment in the fragment placeholder list of centers
iSstart = 0         ! The index into the full shells list for the first shell of a fragment
iSbasis = 0         ! The basis function index relative to the start of the fragment
iSEnd = -1          ! Dummy initialize
do iS=1,nSkal
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

  iSize = nTri_Elem1(iAng)
  if (Shells(iShll)%Transf .and. Shells(iShll)%Prjct) iSize = 2*iAng+1
  if (abs(dbsc(iCnttp)%nFragCoor) /= iCurMdc) then
    ! update fragment related quantities
    iCurMdc = abs(dbsc(iCnttp)%nFragCoor)
    iSstart = iS
    iSend = nSkal
    do iStemp=iSstart+1,nSkal
      if (abs(dbsc(iSD(13,iStemp))%nFragCoor) /= iCurMdc) then
        iSend = iStemp-1
        exit
      end if
    end do
    iSbasis = 1
    iCurCenter = iCurCenter+1
    if (iCurCenter > dbsc(iCurCnttp)%nCntr) then
      iCurCenter = 1
      do jCnttp=iCurCnttp+1,nCnttp
        if (dbsc(jCnttp)%nFragType > 0) then
          iCurCnttp = jCnttp
          exit
        end if
      end do
      ! update the energy weighted density matrix of the current fragment
      EnergyWeight = .true.
      call MakeDens(dbsc(iCurCnttp)%nFragDens,dbsc(iCurCnttp)%nFragEner,dbsc(iCurCnttp)%FragCoef,dbsc(iCurCnttp)%FragEner, &
                    EnergyWeight,Array)
      if (iPrint >= 49) call TriPrt('Energy weighted fragment dens',' ',Array,dbsc(iCurCnttp)%nFragDens)
      ! include the minus sign of -2eta_i
      call DScal_(dbsc(iCurCnttp)%nFragDens*(dbsc(iCurCnttp)%nFragDens+1)/2,-One,Array,1)
      if (maxDensSize < dbsc(iCurCnttp)%nFragDens*(dbsc(iCurCnttp)%nFragDens+1)/2) stop 'maxIJSize'
    end if
  end if
  !write(u6,*) '  iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp',iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp
  !write(u6,*) '  iPrim,iBas =',iPrim,iBas

  ! extra derivative stuff
  iuvwx(3) = dc(mdci)%nStab
  iuvwx(4) = dc(mdci)%nStab
  JndGrd(:,1:2) = IndGrd
  JfGrad(:,1:2) = IfGrad

  JfGrad(:,3:4) = .false.
  ! always equivalent of pChrg's
  JndGrd(:,3:4) = 0
  mGrad = 0
  do iCar=1,3
    do i=1,2
      if (JfGrad(iCar,i)) mGrad = mGrad+1
    end do
  end do
  if (iPrint >= 99) write(u6,*) ' mGrad=',mGrad
  if (mGrad == 0) cycle
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Loop over all other shells belonging to the same fragment
  jSbasis = 1
  do jS=iSstart,iSend
    jShll = iSD(0,jS)
    jAng = iSD(1,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jCnttp = iSD(13,jS)
    jCnt = iSD(14,jS)
    jSize = nTri_Elem1(jAng)
    if (Shells(jShll)%Transf .and. Shells(jShll)%Prjct) jSize = 2*jAng+1
    B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
    !write(u6,*) '    jShll,jAng,jCnttp =',jShll,jAng,jCnttp
    !write(u6,*) '    jPrim,jBas =',jPrim,jBas
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Create a rectangular matrix sized (iBas*nTri_Elem1(iAng),jBas*nTri_Elem1(jAng))
    ! from the energy weighted density matrix (desymmetrized)
    ! contains values from iSbasis to iSbasis + iBas*nTri_Elem1(iAng) - 1
    !             and from jSbasis to jSbasis + jBas*nTri_Elem1(jAng) - 1
    ipIJ = 1+maxDensSize
    !write(u6,*) '    extracting values from',iSbasis,' to',iSbasis+iBas*iSize-1,', and from',jSbasis,' to',jSbasis+jBas*jSize-1
    do iSlocal=iSbasis,iSbasis+iBas*iSize-1
      do jSlocal=jSbasis,jSbasis+jBas*jSize-1
        iLoc = ipIJ+(jSlocal-jSbasis)*iBas*iSize+iSlocal-iSbasis
        Array(iLoc) = Array(iTri(iSlocal,jSlocal))
        if (iSlocal /= jSlocal) Array(iLoc) = Array(iLoc)/Two
        !write(u6,*) 'Filling (',iSlocal-iSbasis+1,',',jSlocal-jSbasis+1,') from (',iSlocal,',',jSlocal,')'
      end do
    end do
    if (iPrint >= 99) call RecPrt('W(KC,LD)',' ',Array(ipIJ),iBas*iSize,jBas*jSize)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! DCR stuff (iS and jS have always the same symmetry character)

    call DCR(LmbdT,iStabM,nStabM,dc(mdci)%iStab,dc(mdci)%nStab,iDCRT,nDCRT)
    Fact = real(nStabM,kind=wp)/real(LmbdT,kind=wp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Loop over symmetry operations acting on the basis.

    do lDCRT=0,nDCRT-1
      lOp(3) = iDCRT(lDCRT)
      lOp(4) = lOp(3)
      call OA(iDCRT(lDCRT),C,TC)
      call OA(iDCRT(lDCRT),B,TB)
      if (EQ(A,RB) .and. EQ(A,TC)) cycle
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Calculate the overlap integral < alpha | is > and its derivative
      !
      !*** Storage

      ip = ipIJ+maxDensSize
      ipF1 = ip
      nac = nTri_Elem1(la)*nTri_Elem1(iAng)*4
      ip = ip+nAlpha*nac*iPrim
      ipP1 = ip
      ip = ip+3*nAlpha*iPrim
      ipZ1 = ip
      ip = ip+nAlpha*iPrim
      ipK1 = ip
      ip = ip+nAlpha*iPrim
      ipZI1 = ip
      ip = ip+nAlpha*iPrim
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(1) in FragPGrd'
        call Abend()
      end if

      !*** Effective center and exponent

      call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,iPrim,Alpha,Shells(iShll)%Exp)
      call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

      !*** Overlap and derivative

      n_Her = ((la+1)+iAng+2)/2
      ipAxyz = ip
      ip = ip+nAlpha*iPrim*3*n_Her*(la+2)
      ipCxyz = ip
      ip = ip+nAlpha*iPrim*3*n_Her*(iAng+1)
      ipRxyz = ip
      ip = ip+nAlpha*iPrim*3*n_Her*(nOrdOp+1)
      ipQ1 = ip
      ip = ip+nAlpha*iPrim*3*(la+2)*(iAng+1)*(nOrdOp+1)
      ipA = ip
      ip = ip+nAlpha*iPrim
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(1b) in FragPGrd'
        call Abend()
      end if
      ABeq(1) = A(1) == TC(1)
      ABeq(2) = A(2) == TC(2)
      ABeq(3) = A(3) == TC(3)
      call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,A,Array(ipAxyz),la+1,HerR(iHerR(n_Her)),n_Her,ABeq)
      call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,TC,Array(ipCxyz),iAng,HerR(iHerR(n_Her)),n_Her,ABeq)
      ABeq(1) = .false.
      ABeq(2) = .false.
      ABeq(3) = .false.
      call CrtCmp(Array(ipZ1),Array(ipP1),nAlpha*iPrim,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(n_Her)),n_Her,ABeq)
      if (iPrint >= 49) then
        write(u6,*) ' Array(ipAxyz)=',DNrm2_(nAlpha*iPrim*3*n_Her*(la+2),Array(ipAxyz),1)
        write(u6,*) ' Array(ipCxyz)=',DNrm2_(nAlpha*iPrim*3*n_Her*(iAng+1),Array(ipCxyz),1)
        write(u6,*) ' Array(ipRxyz)=',DNrm2_(nAlpha*iPrim*3*n_Her*(nOrdOp+1),Array(ipRxyz),1)
      end if
      call Assmbl(Array(ipQ1),Array(ipAxyz),la+1,Array(ipRxyz),nOrdOp,Array(ipCxyz),iAng,nAlpha*iPrim,HerW(iHerW(n_Her)),n_Her)
      iStrt = ipA
      do iGamma=1,iPrim
        call dcopy_(nAlpha,Alpha,1,Array(iStrt),1)
        iStrt = iStrt+nAlpha
      end do
      if (iPrint >= 49) write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*iPrim,Array(ipA),1)
      call rKappa_Zeta(Array(ipK1),Array(ipZ1),iPrim*nAlpha)
      ld = 1
      call CmbnAC(Array(ipQ1),nAlpha*iPrim,la,iAng,Array(ipK1),Array(ipF1),Array(ipA),JfGrad(1,1),ld,nVecAC)
      if (iPrint >= 49) then
        write(u6,*) ' Array(ipQ1)=',DNrm2_(nAlpha*iPrim*3*(la+2)*(iAng+1)*(nOrdOp+1),Array(ipQ1),1)
        write(u6,*) ' Array(ipA)=',DNrm2_(nAlpha*iPrim,Array(ipA),1)
      end if
      ip = ip-nAlpha*iPrim*(6+3*n_Her*(la+2)+3*n_Her*(iAng+1)+3*n_Her*(nOrdOp+1)+3*(la+2)*(iAng+1)*(nOrdOp+1)+1)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Calculate the overlap integral < js | beta > and its derivative
      !
      !*** Storage

      ipF2 = ip
      ncb = nTri_Elem1(jAng)*nTri_Elem1(lb)*4
      ip = ip+jPrim*nBeta*ncb
      ipP2 = ip
      ip = ip+3*jPrim*nBeta
      ipZ2 = ip
      ip = ip+jPrim*nBeta
      ipK2 = ip
      ip = ip+jPrim*nBeta
      ipZI2 = ip
      ip = ip+jPrim*nBeta
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(2) in FragPGrd'
        call Abend()
      end if

      !*** Effective center and exponent

      call ZXia(Array(ipZ2),Array(ipZI2),jPrim,nBeta,Shells(jShll)%Exp,Beta)
      call SetUp1(Shells(jShll)%Exp,jPrim,Beta,nBeta,TB,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

      !*** Overlap and derivative

      n_Her = (jAng+(lb+1)+2)/2
      ipCxyz = ip
      ip = ip+nBeta*jPrim*3*n_Her*(jAng+1)
      ipBxyz = ip
      ip = ip+nBeta*jPrim*3*n_Her*(lb+2)
      ipRxyz = ip
      ip = ip+nBeta*jPrim*3*n_Her*(nOrdOp+1)
      ipQ1 = ip
      ip = ip+nBeta*jPrim*3*(jAng+1)*(lb+2)*(nOrdOp+1)
      ipB = ip
      ip = ip+nBeta*jPrim
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(2b) in FragPGrd'
        call Abend()
      end if
      ABeq(1) = TB(1) == RB(1)
      ABeq(2) = TB(2) == RB(2)
      ABeq(3) = TB(3) == RB(3)
      call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,TB,Array(ipCxyz),jAng,HerR(iHerR(n_Her)),n_Her,ABeq)
      call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,RB,Array(ipBxyz),lb+1,HerR(iHerR(n_Her)),n_Her,ABeq)
      ABeq(1) = .false.
      ABeq(2) = .false.
      ABeq(3) = .false.
      call CrtCmp(Array(ipZ2),Array(ipP2),jPrim*nBeta,Ccoor,Array(ipRxyz),nOrdOp,HerR(iHerR(n_Her)),n_Her,ABeq)
      if (iPrint >= 49) then
        write(u6,*) ' Array(ipCxyz)=',DNrm2_(nBeta*jPrim*3*n_Her*(jAng+1),Array(ipCxyz),1)
        write(u6,*) ' Array(ipBxyz)=',DNrm2_(nBeta*jPrim*3*n_Her*(lb+2),Array(ipBxyz),1)
        write(u6,*) ' Array(ipRxyz)=',DNrm2_(nBeta*jPrim*3*n_Her*(nOrdOp+1),Array(ipRxyz),1)
      end if
      call Assmbl(Array(ipQ1),Array(ipCxyz),jAng,Array(ipRxyz),nOrdOp,Array(ipBxyz),lb+1,jPrim*nBeta,HerW(iHerW(n_Her)),n_Her)
      iStrt = ipB
      do iGamma=1,jPrim
        call dcopy_(nBeta,Beta,1,Array(iStrt),jPrim)
        iStrt = iStrt+1
      end do
      if (iPrint >= 49) write(u6,*) ' Array(ipB)=',DNrm2_(jPrim*nBeta,Array(ipB),1)
      call rKappa_Zeta(Array(ipK2),Array(ipZ2),jPrim*nBeta)
      ld = 1
      call CmbnCB(Array(ipQ1),jPrim*nBeta,jAng,lb,Array(ipK2),Array(ipF2),Array(ipB),JfGrad(1,2),ld,nVecCB)
      if (iPrint >= 49) then
        write(u6,*) ' Array(ipQ1)=',DNrm2_(jPrim*nBeta*3*(la+2)*(jAng+1)*(nOrdOp+1),Array(ipQ1),1)
        write(u6,*) ' Array(ipB)=',DNrm2_(JPrim*nBeta,Array(ipB),1)
      end if
      ip = ip-nBeta*jPrim*(6+3*n_Her*(lb+2)+3*n_Her*(jAng+1)+3*n_Her*(nOrdOp+1)+3*(lb+2)*(jAng+1)*(nOrdOp+1)+1)
      nac = nTri_Elem1(la)*nTri_Elem1(iAng)*nVecAC
      ncb = nTri_Elem1(jAng)*nTri_Elem1(lb)*nVecCB
      ipTmp = ip
      ip = ip+max(nAlpha*max(iPrim,jBas)*nac,nBeta*ncb*jBas)
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(3) in FragPGrd'
        call Abend()
      end if
      nac = nTri_Elem1(la)*nTri_Elem1(iAng)
      ncb = nTri_Elem1(jAng)*nTri_Elem1(lb)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Assemble the calculated quantities and contract
      !
      ! Calculate Contraction over components of the fragment
      ! orbitals of type <A|iS>coef<jS|B> where we now have in
      ! Array(ipF1) the cartesian components of <A|iS>, and
      ! similarily, in Array(ipF2), we have stored the cartesian
      ! components of <jS|B>. Observe that the fragment orbitals are
      ! orthonomal atomic orbitals. Hence, the transformation
      ! to the spherical harmonics has to be for normalized
      ! spherical harmonics.
      !
      ! nAlpha = i               nTri_Elem1(la) = a
      ! nBeta  = j               nTri_Elem1(lb) = b
      ! iPrim = k (iBas = K)     nTri_Elem1(iAng) = c (iSize = C)
      ! jPrim = l (jBas = L)     nTri_Elem1(jAng) = d (jSize = D)
      !
      !---From the lefthandside overlap, form iKaC from ikac by
      !   1) i,kac -> k,aci

      call DgeTMo(Array(ipF1),nAlpha,nAlpha,iPrim*nac*nVecAC,Array(ipTmp),iPrim*nac*nVecAC)

      !---2) aciK =  k,aci * k,K (Contract over core orbital)

      call DGEMM_('T','N',nac*nVecAC*nAlpha,iBas,iPrim,One,Array(ipTmp),iPrim,Shells(iShll)%pCff,iPrim,Zero,Array(ipF1), &
                  nac*nVecAC*nAlpha)

      !---3) a,ciK -> ciKa

      call DgeTMo(Array(ipF1),nTri_Elem1(la),nTri_Elem1(la),nTri_Elem1(iAng)*nVecAC*nAlpha*iBas,Array(ipTmp), &
                  nTri_Elem1(iAng)*nVecAC*nAlpha*iBas)

      !---4) iKa,C = c,iKa * c,C

      if (Shells(iShll)%Transf .and. Shells(iShll)%Prjct) then
        call DGEMM_('T','N',nVecAC*nAlpha*iBas*nTri_Elem1(la),iSize,nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng), &
                    RSph(ipSph(iAng)),nTri_Elem1(iAng),Zero,Array(ipF1),nVecAC*nAlpha*iBas*nTri_Elem1(la))
      else
        call DgeTMo(Array(ipTmp),nTri_Elem1(iAng),nTri_Elem1(iAng),nVecAC*iBas*nTri_Elem1(la)*nAlpha,Array(ipF1), &
                    nVecAC*iBas*nTri_Elem1(la)*nAlpha)
      end if
      ! what does this do and is it needed? (from PrjGrd)
      call DgeTMo(Array(ipF1),nVecAC,nVecAC,nAlpha*iBas*nTri_Elem1(la)*iSize,Array(ipTmp),nAlpha*iBas*nTri_Elem1(la)*iSize)
      call dcopy_(nVecAC*nAlpha*iBas*nTri_Elem1(la)*iSize,Array(ipTmp),1,Array(ipF1),1)

      !---And (almost) the same thing for the righthand side, form
      !   LjDb from ljdb
      !   1) jdb,L = l,jdb * l,L

      call DGEMM_('T','N',nBeta*ncb*nVecCB,jBas,jPrim,One,Array(ipF2),jPrim,Shells(jShll)%pCff,jPrim,Zero,Array(ipTmp), &
                  nBeta*ncb*nVecCB)

      !---2)  j,dbL -> dbL,j

      call DgeTMo(Array(ipTmp),nBeta,nBeta,ncb*nVecCB*jBas,Array(ipF2),ncb*nVecCB*jBas)

      !---3) bLj,D = d,bLj * d,D

      if (Shells(jShll)%Transf .and. Shells(jShll)%Prjct) then
        call DGEMM_('T','N',nTri_Elem1(lb)*nVecCB*jBas*nBeta,jSize,nTri_Elem1(jAng),One,Array(ipF2),nTri_Elem1(jAng), &
                    RSph(ipSph(jAng)),nTri_Elem1(jAng),Zero,Array(ipTmp),nTri_Elem1(lb)*nVecCB*jBas*nBeta)
      else
        call DgeTMo(Array(ipF2),nTri_Elem1(jAng),nTri_Elem1(jAng),nVecCB*jBas*nTri_Elem1(lb)*nBeta,Array(ipTmp), &
                    nVecCB*jBas*nTri_Elem1(lb)*nBeta)
      end if

      !---4) b,LjD -> LjD,b

      call DgeTMo(Array(ipTmp),nTri_Elem1(lb)*nVecCB,nTri_Elem1(lb)*nVecCB,jBas*nBeta*jSize,Array(ipF2),jBas*nBeta*jSize)

      !---Next Contract (iKaC)*W(KLCD)*(LjDb) producing ijab

      rFinal(:,:,:,1,:) = Zero

      if (iPrint >= 99) then
        call RecPrt('ipF1 (nVecAC x X)',' ',Array(ipF1),nVecAC,iBas*nAlpha*iSize)
        call RecPrt('ipF2 (nVecCB x Y)',' ',Array(ipF2),nVecCB,jBas*nBeta*jSize)
      end if

      mVec = 0
      mVecAC = 1
      mVecCB = 1
      do iCar=1,3
        do iCent=1,2
          !write(u6,*) 'iCar, iCent = ',iCar,iCent
          if (JfGrad(iCar,iCent)) then
            mVec = mVec+1
            if (iCent == 1) then
              mVecAC = mVecAC+1
              ipF1a = ipF1+(mVecAC-1)*nAlpha*jBas*nTri_Elem1(la)*iSize
              ipF2a = ipF2
            else
              ipF1a = ipF1
              mVecCB = mVecCB+1
              ipF2a = ipF2+(mVecCB-1)*jBas*nBeta*jSize*nTri_Elem1(lb)
            end if
            if (iPrint >= 99) then
              write(u6,*) 'mVecAC, mVecCB = ',mVecAC,mVecCB
              call RecPrt('ipF1a (nAlpha*aAng x iBas*iSize)',' ',Array(ipF1a),nAlpha*nTri_Elem1(la),iBas*iSize)
              call RecPrt('ipF2a (nBeta*bAng x jBas*jSize)',' ',Array(ipF2a),nBeta*nTri_Elem1(lb),jBas*jSize)
            end if

            call FragPCont(Array(ipF1a),nAlpha,iBas,nTri_Elem1(la),iSize,Array(ipF2a),jBas,nBeta,jSize,nTri_Elem1(lb),Array(ipIJ), &
                           rFinal(:,:,:,1,mVec),Fact*Half)
          end if
        end do !iCent
      end do !iCar

      if (iPrint >= 49) then
        do iVec=1,mVec
          write(u6,*) iVec,sqrt(DNrm2_(nZeta*nTri_Elem1(la)*nTri_Elem1(lb),rFinal(:,:,:,1,iVec),1))
        end do
      end if
      if (iPrint >= 99) then
        write(u6,*) ' Result in FragPGrd'
        do ia=1,nTri_Elem1(la)
          do ib=1,nTri_Elem1(lb)
            do iVec=1,mVec
              write(Label,'(A,I2,A,I2,A,I2,A)') ' rFinal(',ia,',',ib,',',iVec,')'
              call RecPrt(Label,' ',rFinal(:,ia,ib,1,iVec),nAlpha,nBeta)
            end do
          end do
        end do
      end if

      !---Distribute contributions to the gradient

      call Distg1X(rFinal,DAO,nZeta,nDAO,mVec,Grad,nGrad,JfGrad,JndGrd,iuvwx,lOp)

    end do !lDCRT
    jSbasis = jSbasis+jBas*jSize
  end do !jS
  iSbasis = iSbasis+iBas*iSize
end do !iS

! Revert to the valence shells

call Free_iSD()

return

end subroutine FragPGrd

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
!               2016, Liviu Ungur                                      *
!***********************************************************************

subroutine FragPInt( &
#                   define _CALLING_
#                   include "int_interface.fh"
                   )
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of Fragment AIEMP         *
!         projection integrals                                         *
!                                                                      *
!     Author: Ben Swerts                                               *
!     based on seward/prjint.f                                         *
!   Modified: Liviu Ungur                                              *
!                                                                      *
! additional info:                                                     *
!       mdci, mdci = show the type number (i.e. from 1 to nCnttp) of   *
!                    the corresponding "fragment" atom. Must never take*
!                    values of the atoms from the "central region"     *
!   iCnttp, jCnttp = index number of the the type of the fragment atom *
!  nFragType(mdci) = number of different basis sets (kinds of atoms)   *
!                    are in the "fragment"                             *
!  nFragCoor(mdci) = number of atoms of the fragment "mdci"            *
!                    must be 1, as we expanded atoms in one atom/group *
!                    check also the "FragExpand" subroutine            *
!  nFragDens(mdci) = nBas of the entire fragment type "mdci"           *
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, nCnttp, Shells
use Center_Info, only: dc
use Symmetry_Info, only: nIrrep, iChTbl
use Index_Functions, only: iTri, nTri_Elem1
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
! Local variables
integer(kind=iwp) :: iAng, iBas, iCnttp, iComp, iCurCenter, iCurCnttp, iCurMdc, iIC, iIrrep, iLoc, iPrim, ip, ipF1, ipF2, ipIJ, &
                     ipK1, ipK2, ipP1, ipP2, ipTmp, ipZ1, ipZ2, ipZI1, ipZI2, iS, iSbasis, iSend, iShll, iSize, iSlocal, iSstart, &
                     iStemp, jAng, jBas, jCnttp, jPrim, jS, jShll, jSize, jSlocal, lDCRT, llOper, LmbdT, mArr, maxDensSize, mdci, &
                     nac, ncb, nDCRT, nOp, nSkal, jSbasis, iCnt, jCnt, iDCRT(0:7)
real(kind=wp) :: C(3), TC(3), B(3), TB(3), Fact, Factor, Xg
logical(kind=iwp) :: EnergyWeight
! external functions:
integer(kind=iwp), external :: NrOpr
!real(kind=wp), external :: DNRM2_
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ia, ib
character(len=24) :: Label
#endif
#include "macros.fh"
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(nHer)
unused_var(iCho)
unused_var(PtChrg)
unused_var(iAddPot)

#ifdef _DEBUGPRINT_
! data for individual fragments:
write(u6,*) ' In FragPInt:    nCnttp          = ',nCnttp
write(u6,*) ' In FragPInt: dbsc(nCnttp)%mdci  = ',(dbsc(i)%mdci,i=1,nCnttp)
write(u6,*) ' In FragPInt:     dbsc(nCnttp)%nCntr  = ',(dbsc(i)%nCntr,i=1,nCnttp)
write(u6,*) ' In FragPInt: nFragType(nCnttp)  = ',(dbsc(i)%nFragType,i=1,nCnttp)
write(u6,*) ' In FragPInt: nFragCoor(nCnttp)  = ',(dbsc(i)%nFragCoor,i=1,nCnttp)
write(u6,*) ' In FragPInt: nFragDens(nCnttp)  = ',(dbsc(i)%nFragDens,i=1,nCnttp)
write(u6,*) ' In FragPInt: nFragDens(nCnttp)  = ',(dbsc(i)%nFragDens,i=1,nCnttp)

write(u6,*) ' In FragPInt: nAlpha,nBeta=',' ',nAlpha,nBeta
write(u6,*) ' In FragPInt: Alpha=',' ',Alpha
write(u6,*) ' In FragPInt: Beta=',' ',Beta

write(u6,*) ' In FragPInt: nZeta=',' ',nZeta
write(u6,*) ' In FragPInt:  nArr=',' ',nArr
write(u6,*) ' In FragPInt:   nIC=',' ',nIC
write(u6,*) ' In FragPInt: la,lb=',' ',la,lb
write(u6,*) ' In FragPInt: nTri_Elem1(la)=',' ',nTri_Elem1(la)
write(u6,*) ' In FragPInt: nTri_Elem1(lb)=',' ',nTri_Elem1(lb)
call RecPrt(' In FragPInt: A     ',' ',A,1,3)
call RecPrt(' In FragPInt: RB    ',' ',RB,1,3)
call RecPrt(' In FragPInt: Ccoor ',' ',Ccoor,1,3)
call RecPrt(' In FragPInt: P     ',' ',P,nZeta,3)
call RecPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
call TrcPrt(' In FragPInt: Array ',' ',Array,nZeta,nArr)
#endif

rFinal(:,:,:,:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Setup the fragment shells

call Free_iSD()
call Set_Basis_Mode('Fragments')
call SetUp_iSD()
call Nr_Shells(nSkal)
#ifdef _DEBUGPRINT_
write(u6,*) 'nSkal_Fragment,nAlpha,nBeta = ',nSkal,nAlpha,nBeta
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Reserve space for the largest possible fragment energy weighted
! density matrix and one combination of shells of it.
maxDensSize = 0
do iCnttp=1,nCnttp
  if (dbsc(iCnttp)%nFragType > 0) then
    maxDensSize = max(maxDensSize,dbsc(iCnttp)%nFragDens*(dbsc(iCnttp)%nFragDens+1)/2)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,i2,A,i6)') 'nFragDens(',iCnttp,')=',dbsc(iCnttp)%nFragDens
#   endif
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over all shells belonging to the fragments

llOper = lOper(1)
iComp = 1
iCurMdc = 0         ! The mdc value of the current fragment placeholder
iCurCnttp = 1       ! The Cnttp of the fragment placeholder
iCurCenter = 999999 ! The index of the fragment in the fragment placeholder list of centers
iSstart = 0         ! The index into the full shells list for the first shell of a fragment
isEnd = nSkal       ! The index into the full shells list for the last shell of a fragment
iSbasis = 0         ! The basis function index relative to the start of the fragment
do iS=1,nSkal
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  mdci = iSD(10,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  iSize = nTri_Elem1(iAng)
  C(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
  ! some printouts:
# ifdef _DEBUGPRINT_
  write(u6,*) 'In FragPInt: iS=',iS,' iShll =',iShll
  write(u6,*) 'In FragPInt: iS=',iS,' iAng  =',iAng
  write(u6,*) 'In FragPInt: iS=',iS,' iBas  =',iBas
  write(u6,*) 'In FragPInt: iS=',iS,' iPrim =',iPrim
  write(u6,*) 'In FragPInt: iS=',iS,' mdci  =',mdci
  write(u6,*) 'In FragPInt: iS=',iS,' iCnttp=',iCnttp
  write(u6,*) 'In FragPInt: iS=',iS,' iSize =',iSize
  write(u6,*) 'In FragPInt: iS=',iS,' iCurMdc =',iCurMdc
  write(u6,*) 'In FragPInt: nFragCoor(',iCnttp,') =',dbsc(iCnttp)%nFragCoor
# endif

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

#   ifdef _DEBUGPRINT_
    write(u6,*) 'start of new fragment encountered'
    write(u6,*) 'iSstart,iSend,iCurCnttp,iCurCenter =',iSstart,iSend,iCurCnttp,iCurCenter
#   endif

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
#     ifdef _DEBUGPRINT_
      call TriPrt('Energy weighted fragment dens',' ',Array,dbsc(iCurCnttp)%nFragDens)
#     endif
      ! include the minus sign of -2eta_i
      call DScal_(dbsc(iCurCnttp)%nFragDens*(dbsc(iCurCnttp)%nFragDens+1)/2,-One,Array,1)
#     ifdef _DEBUGPRINT_
      call TriPrt('-1 Energy weighted fragment dens',' ',Array,dbsc(iCurCnttp)%nFragDens)
#     endif
      if (maxDensSize < dbsc(iCurCnttp)%nFragDens*(dbsc(iCurCnttp)%nFragDens+1)/2) call Abend() !'maxIJSize'
    end if
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) '  iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp',iShll,iAng,mdci,iCnttp,iCurMdc,iCurCnttp
  write(u6,*) '  iPrim,iBas =',iPrim,iBas
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Loop over all other shells belonging to the same fragment
  jSbasis = 1
# ifdef _DEBUGPRINT_
  write(u6,'(3(A,i4))') 'iS = ',iS,' iSstart=',iSstart,' iSEnd=',iSend
# endif
  do jS=iSstart,iSend
    jShll = iSD(0,jS)
    jAng = iSD(1,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jCnttp = iSD(13,jS)
    jCnt = iSD(14,jS)
    jSize = nTri_Elem1(jAng)
    B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)
#   ifdef _DEBUGPRINT_
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jShll =',jShll
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jAng  =',jAng
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jBas  =',jBas
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jPrim =',jPrim
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jCnttp=',jCnttp
    write(u6,'(A,i6,A,i16)') 'In FragPInt: jS=',jS,' jSize =',jSize
#   endif

    if (Shells(jShll)%Transf .and. Shells(jShll)%Prjct) jSize = 2*jAng+1
#   ifdef _DEBUGPRINT_
    write(u6,*) '    jShll,jAng,jCnttp =',jShll,jAng,jCnttp
    write(u6,*) '    jPrim,jBas =',jPrim,jBas
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Create a rectangular matrix sized (iBas*nTri_Elem1(iAng),jBas*nTri_Elem1(jAng))
    ! from the energy weighted density matrix (desymmetrized)
    ! contains values from iSbasis to iSbasis + iBas*nTri_Elem1(iAng) - 1
    !             and from jSbasis to jSbasis + jBas*nTri_Elem1(jAng) - 1
    ipIJ = 1+maxDensSize
#   ifdef _DEBUGPRINT_
    write(u6,*) '    ipIJ=',ipIJ
    write(u6,*) '    extracting values from',iSbasis,' to',iSbasis+iBas*iSize-1,', and from',jSbasis,' to',jSbasis+jBas*jSize-1
#   endif
    do iSlocal=iSbasis,iSbasis+iBas*iSize-1
      do jSlocal=jSbasis,jSbasis+jBas*jSize-1
        iLoc = ipIJ+(jSlocal-jSbasis)*iBas*iSize+iSlocal-iSbasis
#       ifdef _DEBUGPRINT_
        write(u6,'(A,i3,A,i3,A,i4,A,i8)') 'iTri(',iSlocal,',',jSlocal,')=',iTri(iSlocal,jSlocal),' iLoc=',iLoc
#       endif
        Array(iLoc) = Array(iTri(iSlocal,jSlocal))
        if (iSlocal /= jSlocal) Array(iLoc) = Array(iLoc)/Two
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Filling (',iSlocal-iSbasis+1,',',jSlocal-jSbasis+1,') from (',iSlocal,',',jSlocal,')'
#       endif
      end do
    end do
#   ifdef _DEBUGPRINT_
    call RecPrt('W(KC,LD)',' ',Array(ipIJ),iBas*iSize,jBas*jSize)
#   endif
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
      call OA(iDCRT(lDCRT),C,TC)
      call OA(iDCRT(lDCRT),B,TB)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Calculate the overlap integral < alpha | is >

      ip = ipIJ+maxDensSize
      ipF1 = ip
      nac = nTri_Elem1(la)*nTri_Elem1(iAng)
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
        write(u6,*) '  ip-1 > nArr*nZeta(1) in FragPInt'
        call Abend()
      end if
      mArr = (nArr*nZeta-(ip-1))/nZeta

      call ZXia(Array(ipZ1),Array(ipZI1),nAlpha,iPrim,Alpha,Shells(iShll)%Exp)
      call SetUp1(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,A,TC,Array(ipK1),Array(ipP1),Array(ipZI1))

      nHer = (la+iAng+2)/2
      call MltPrm(Alpha,nAlpha,Shells(iShll)%Exp,iPrim,Array(ipZ1),Array(ipZI1),Array(ipK1),Array(ipP1),Array(ipF1),nAlpha*iPrim, &
                  iComp,la,iAng,A,TC,nHer,Array(ip),mArr,CCoor,nOrdOp)
#     ifdef _DEBUGPRINT_
      call RecPrt('<alpha|iS> (aBas x X)',' ',Array(ipF1),nAlpha*iPrim,nac)
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Calculate the overlap integral < jS | beta >

      ip = ip-6*nAlpha*iPrim
      ipF2 = ip
      ncb = nTri_Elem1(jAng)*nTri_Elem1(lb)
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
        write(u6,*) '  ip-1 > nArr*nZeta(2) in FragPInt'
        call Abend()
      end if
      mArr = (nArr*nZeta-(ip-1))/nZeta

      call ZXia(Array(ipZ2),Array(ipZI2),jPrim,nBeta,Shells(jShll)%Exp,Beta)
      call SetUp1(Shells(jShll)%Exp,jPrim,Beta,nBeta,TB,RB,Array(ipK2),Array(ipP2),Array(ipZI2))

      nHer = (jAng+lb+2)/2
      call MltPrm(Shells(jShll)%Exp,jPrim,Beta,nBeta,Array(ipZ2),Array(ipZI2),Array(ipK2),Array(ipP2),Array(ipF2),jPrim*nBeta, &
                  iComp,jAng,lb,TB,RB,nHer,Array(ip),mArr,CCoor,nOrdOp)
      ip = ip-6*jPrim*nBeta
      ipTmp = ip
      ip = ip+max(nAlpha*nac*max(iPrim,jBas),nBeta*ncb*jBas)
      if (ip-1 > nArr*nZeta) then
        write(u6,*) '  ip-1 > nArr*nZeta(3) in FragPInt'
        call Abend()
      end if
#     ifdef _DEBUGPRINT_
      call RecPrt('<jS|beta> (bBas x Y)',' ',Array(ipF2),nBeta*jPrim,ncb)
#     endif
      if ((ipTmp-ipF2) < nBeta*jPrim*ncb) stop 'sizetest 1'
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

      call DgeTMo(Array(ipF1),nAlpha,nAlpha,iPrim*nac,Array(ipTmp),iPrim*nac)
#     ifdef _DEBUGPRINT_
      call RecPrt('<alpha|iS>^T (X x aBas)',' ',Array(ipTmp),iPrim*nac,nAlpha)
#     endif
      if ((ip-ipTmp) < nAlpha*iPrim*nac) stop 'sizetest 2'

      !---2) aciK =  k,aci * k,K

      call DGEMM_('T','N',nac*nAlpha,iBas,iPrim,One,Array(ipTmp),iPrim,Shells(iShll)%pCff,iPrim,Zero,Array(ipF1),nAlpha*nac)
#     ifdef _DEBUGPRINT_
      call RecPrt('<alpha|iS>(regrouped, X x iPrim)',' ',Array(ipTmp),nAlpha*nac,iPrim)
      call RecPrt('Coeffs of iS (iPrim x iBas)',' ',Shells(iShll)%pCff,iPrim,iBas)
      call RecPrt('<alpha|iS>(re) * Coeffs of iS (X x iBas)',' ',Array(ipF1),nAlpha*nac,iBas)
#     endif
      if ((ipF2-ipF1) < nAlpha*iBas*nac) stop 'sizetest 3'

      !---4) a,ciK -> ciKa

      call DgeTMo(Array(ipF1),nTri_Elem1(la),nTri_Elem1(la),nTri_Elem1(iAng)*nAlpha*iBas,Array(ipTmp),nTri_Elem1(iAng)*nAlpha*iBas)
#     ifdef _DEBUGPRINT_
      call RecPrt('result (regrouped, nTri_Elem1(la) x X)',' ',Array(ipF1),nTri_Elem1(la),nTri_Elem1(iAng)*nAlpha*iBas)
      call RecPrt('transpose of result (X x nTri_Elem1(la))',' ',Array(ipTmp),nTri_Elem1(iAng)*nAlpha*iBas,nTri_Elem1(la))
#     endif
      if ((ip-ipTmp) < nAlpha*iBas*nTri_Elem1(iAng)*nTri_Elem1(la)) stop 'sizetest 6'

      !---5) iKa,C = c,iKa * c,C

      if (Shells(iShll)%Transf .and. Shells(iShll)%Prjct) then
        call DGEMM_('T','N',iBas*nTri_Elem1(la)*nAlpha,iSize,nTri_Elem1(iAng),One,Array(ipTmp),nTri_Elem1(iAng),RSph(ipSph(iAng)), &
                    nTri_Elem1(iAng),Zero,Array(ipF1),nAlpha*iBas*nTri_Elem1(la))
#       ifdef _DEBUGPRINT_
        call RecPrt('result (regrouped, X x nTri_Elem1(iAng))',' ',Array(ipTmp),nTri_Elem1(la)*nAlpha*iBas,nTri_Elem1(iAng))
        call RecPrt('Spher of iS (nTri_Elem1(iAng) x (2*iAng+1))',' ',RSph(ipSph(iAng)),nTri_Elem1(iAng),(2*iAng+1))
        call RecPrt('result in spherical gaussians (X x iSize',' ',Array(ipF1),nAlpha*iBas*nTri_Elem1(la),iSize)
#       endif
      else
        ! in this case nTri_Elem1(iAng) = iSize
        call DgeTMo(Array(ipTmp),nTri_Elem1(iAng),nTri_Elem1(iAng),iBas*nTri_Elem1(la)*nAlpha,Array(ipF1), &
                    iBas*nTri_Elem1(la)*nAlpha)
      end if
      if ((ipF2-ipF1) < nAlpha*iBas*nTri_Elem1(la)*iSize) stop 'sizetest 7'

      !---And (almost) the same thing for the righthand side, form
      !   LjDb from ljdb
      !   1) jdb,L = l,jdb * l,L
      call DGEMM_('T','N',nBeta*ncb,jBas,jPrim,One,Array(ipF2),jPrim,Shells(jShll)%pCff,jPrim,Zero,Array(ipTmp),nBeta*ncb)
#     ifdef _DEBUGPRINT_
      call RecPrt('<jS|beta>(regrouped, X x jPrim)',' ',Array(ipF2),nBeta*ncb,jPrim)
      call RecPrt('Coeffs of jS (jPrim x jBas)',' ',Shells(jShll)%pCff,jPrim,jBas)
      call RecPrt('<jS|beta>(re) * Coeffs of jS (Y x jBas)',' ',Array(ipTmp),nBeta*ncb,jBas)
#     endif
      if ((ip-ipTmp) < nBeta*jBas*ncb) stop 'sizetest 8'

      !---2)  j,dbL -> dbL,j

      call DgeTMo(Array(ipTmp),nBeta,nBeta,jBas*ncb,Array(ipF2),jBas*ncb)
#     ifdef _DEBUGPRINT_
      call RecPrt('transposed right 1 (Y x bBas)',' ',Array(ipF2),jBas*ncb,nBeta)
#     endif
      if ((ipTmp-ipF2) < nBeta*jBas*ncb) stop 'sizetest 9'

      !---3) bLj,D = d,bLj * d,D

      if (Shells(jShll)%Transf .and. Shells(jShll)%Prjct) then
        call DGEMM_('T','N',nTri_Elem1(lb)*jBas*nBeta,jSize,nTri_Elem1(jAng),One,Array(ipF2),nTri_Elem1(jAng),RSph(ipSph(jAng)), &
                    nTri_Elem1(jAng),Zero,Array(ipTmp),nTri_Elem1(lb)*jBas*nBeta)
#       ifdef _DEBUGPRINT_
        call RecPrt('multiply right 2 (Y x jSize)',' ',Array(ipTmp),nBeta*jBas*nTri_Elem1(lb),jSize)
#       endif
      else
        ! in this case nTri_Elem1(jAng) = jSize
        call DgeTMo(Array(ipF2),nTri_Elem1(jAng),nTri_Elem1(jAng),jBas*nTri_Elem1(lb)*nBeta,Array(ipTmp),jBas*nTri_Elem1(lb)*nBeta)
      end if
      if ((ip-ipTmp) < nBeta*jBas*nTri_Elem1(lb)*jSize) stop 'sizetest 10'

      !---4) b,LjD -> LjD,b

      call DgeTMo(Array(ipTmp),nTri_Elem1(lb),nTri_Elem1(lb),jBas*nBeta*jSize,Array(ipF2),jBas*nBeta*jSize)
#     ifdef _DEBUGPRINT_
      call RecPrt('transposed right 2 (Y x nTri_Elem1(lb)',' ',Array(ipF2),jBas*nBeta*jSize,nTri_Elem1(lb))
#     endif

      if ((ipTmp-ipF2) < nBeta*jBas*nTri_Elem1(lb)*jSize) stop 'sizetest 11'

      !---Next Contract (iKaC)*W(KLCD)*(LjDb) producing ijab,
      !   by the following procedure:
      !   Loop over a and b
      !     Loop over C
      !       Contract iK(aC)*Kj(Cb), over K producing ij(aCb),
      !         accumulate to ij(ab)
      !     End loop C
      !   End Loop b and a
      !
      ! Total size of ipF1 = nAlpha*nTri_Elem1(la) * iBas*iSize  -> ordered as (nAlpha, iBas,  nTri_Elem1(la), iSize)
      !               ipF2 = nBeta*nTri_Elem1(lb)  * jBas*jSize                (jBas,   nBeta, jSize,          nTri_Elem1(lb))
      !                  W = iBas*iSize * jBas*jSize                           (iBas,   iSize, jBas,           jSize)

#     ifdef _DEBUGPRINT_
      write(u6,*) ' Current contents of rFinal():'
      do ia=1,nTri_Elem1(la)
        do ib=1,nTri_Elem1(lb)
          write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
          call RecPrt(Label,' ',rFinal(:,ia,ib,:),nAlpha,nBeta)
        end do
      end do
#     endif

      iIC = 0
      do iIrrep=0,nIrrep-1
        if (btest(llOper,iIrrep)) then
          iIC = iIC+1
          nOp = NrOpr(iDCRT(lDCRT))
          Xg = real(iChTbl(iIrrep,nOp),kind=wp)
          ! Half is needed because we do a complete loop over iS,jS
          Factor = Xg*Fact*Half
          !write(u6,'(A,i24)') 'FragPInt:  ipIJ=',ipIJ
          !print(u6,*) 'CALL FragPCont'
          call xFlush(u6)

          call FragPCont(Array(ipF1),nAlpha,iBas,nTri_Elem1(la),iSize,Array(ipF2),jBas,nBeta,jSize,nTri_Elem1(lb),Array(ipIJ), &
                         rFinal(:,:,:,iIC),Factor)
        end if
      end do
      if (iIC /= nIC) stop 'iIC /= nIC'

#     ifdef _DEBUGPRINT_
      write(u6,*) ' After contraction:'
      do ia=1,nTri_Elem1(la)
        do ib=1,nTri_Elem1(lb)
          write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
          call RecPrt(Label,' ',rFinal(:,ia,ib,:),nAlpha,nBeta)
        end do
      end do
#     endif

    end do ! end loop over lDCR

    jSbasis = jSbasis+jBas*jSize
  end do ! end loop over jS
  call xFlush(u6)

  iSbasis = iSbasis+iBas*iSize
end do ! end loop over iS
call xFlush(u6)

#ifdef _DEBUGPRINT_
write(u6,*) ' Result in FragPInt'
do ia=1,nTri_Elem1(la)
  do ib=1,nTri_Elem1(lb)
    write(Label,'(A,I2,A,I2,A)') ' rFinal(',ia,',',ib,')'
    call RecPrt(Label,' ',rFinal(:,ia,ib,:),nAlpha,nBeta)
  end do
end do
#endif
! add some verification data
! this check is to ensure that Add_Info is called only on the Master node
! The reason being that the execution of this Kernel routine is split on
! diffeent nodes eariler in the code.
! OneEl -> OneEl_ -> OneEl_IJ -> this routine.
! Normally, the Add_Info must be called after the parallelization is finalized,
! i.e. in the OneEl function.
!if (MyRank == 0) then
!  do ia=1,nTri_Elem1(la)
!    do ib=1,nTri_Elem1(lb)
!      dA = Zero
!      dA = dnrm2_(nAlpha*nBeta,rFinal(1,ia,ib,1),1)
!      if (dA > 1.0e-6_wp) then
!        write(label,'(A,i2,A,i2)') 'Fragpint: ',ia,' ib ',ib
!        call Add_Info(label,dA,1,6)
!      end if
!      call xFlush(u6)
!    end do
!  end do
!end if
call Free_iSD()
call Set_Basis_Mode('Valence')
call SetUp_iSD()
call Nr_Shells(nSkal)

return

end subroutine FragPInt

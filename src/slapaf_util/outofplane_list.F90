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
! Copyright (C) 2004, Roland Lindh                                     *
!***********************************************************************

subroutine OutOfPlane_List(nq,nsAtom,iIter,nIter,Cx,Process,Valu,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv,Proc_dB,iTabBonds, &
                           nBonds,iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)
!***********************************************************************
!     This is a quick and possibly dirty implementation of the out-    *
!     of-plane angle. RL, Tokyo June, 2004.                            *
!***********************************************************************

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: ANr, AtomLbl, Fragments_Bond, jStab, Magic_Bond, nStab, vdW_Bond
use ddvdt, only: aAV, f_Const_Min, rAV, rko
use Constants, only: Zero, Pi, deg2rad
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, iRef, LuIC, iPrv, nBonds, iTabBonds(3,nBonds), mAtoms, &
                                 iTabAI(2,mAtoms), nMax, iTabAtoms(2,0:nMax,mAtoms), nB_Tot, ndB_tot
integer(kind=iwp), intent(inout) :: nq, Indq(3,nB), mB_Tot, mdB_Tot, iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter)
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Valu(nB,nIter), fconst(nB), rMult(nB), BM(nB_Tot), dBM(ndB_Tot)
character(len=14), intent(inout) :: qLbl(nB)
#include "Molcas.fh"
integer(kind=iwp), parameter :: mB = 4*3
integer(kind=iwp) :: iAtom, iAtom_, iCase, iDCR(4), iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iDCRX(0:7), iDCRY(0:7), iDeg, iE1, iE2, &
                     iE3, iE4, iF1, iF2, iF3, iF4, ij, ijDCR, Ind(4), ir, iStabM(0:7), iStabN(0:7), iStabO(0:7), jAtom, jAtom_, &
                     jBond, jBondType, jr, kAtom, kAtom_, kBond, kBondType, kDCRR, kDCRS, kDCRT, kDCRTS, kl, kNeighbor, kr, &
                     Lambda, lAtom, lAtom_, lBond, lBondType, lNeighbor, lr, mCent, nCent, nCoBond_j, nDCRR, nDCRS, nDCRT, nDCRX, &
                     nDCRY, nFgBond_j, nNeighbor_i, nqO, nStabM, nStabN, nStabO
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: A(3,4), Alpha, Deg, delta, delta0, f_Const, f_Const_ij, f_Const_ij_Ref, f_Const_ijk, f_Const_ijk_Ref, &
                 f_Const_Ref, Fi2, Fi3, Fi4, Grad(mB), Grad_ref(9), Hess(mB**2), Prv(3,4), r0, Ref(3,4), rij2, rij2_Ref, rik2, &
                 rik2_Ref, ril2, ril2_Ref, RX4Y(3,3), Val
logical(kind=iwp) :: Help
character(len=LenIn4) :: Lbls(4)
character(len=14) :: Label
integer(kind=iwp), parameter :: iChOp(0:7) = [1,1,1,2,1,2,2,3]
character(len=*), parameter :: ChOp(0:7) = ['E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ']
integer(kind=iwp), external :: iTabRow, nCoBond, nFgBond
logical(kind=iwp), external :: R_Stab_A

!                                                                      *
!***********************************************************************
!                                                                      *

if (nBonds < 3) return
nqO = 0
Hess(:) = Zero

! Loop over out-of-plane angles.

nCent = 4
!
!***********************************************************************
!     Notes relevant to the out-of-plane implementation                *
!                                                                      *
!    Connection is 1-4, 2-4, and 3-4. We renumber them according to    *
!    j-i, k-i, and l-i                                                 *
!                                                                      *
!***********************************************************************
!
! Order will play a role here. That is, the torsion
! A-R(B)-T(C)-TS(D) is NOT identical to
! A-R(B)-TS(D)-T(C). Hence we put no restriction on the
! pairs AB and CD. However, for the pair of pairs we have
! that order is irrelevant, i.e. ABCD is identical to
! DCBA. To garantee this we limit the pairs to the unique
! combinations.

do jBond=1,nBonds
  jBondType = iTabBonds(3,jBond)
  if (jBondType == vdW_Bond) cycle
  if (jBondType > Magic_Bond) cycle

  do iCase=1,2

    if (iCase == 1) then
      iAtom_ = iTabBonds(1,jBond)
      jAtom_ = iTabBonds(2,jBond)
    else
      iAtom_ = iTabBonds(2,jBond)
      jAtom_ = iTabBonds(1,jBond)
    end if
    iAtom = iTabAI(1,iAtom_)
    jAtom = iTabAI(1,jAtom_)
    ir = iTabRow(ANr(iAtom))
    jr = iTabRow(ANr(jAtom))
    Ind(1) = jAtom
    Ind(4) = iAtom

    Help = (ir > 3) .or. (jr > 3)
    iDCR(4) = iTabAI(2,iAtom_)
    iDCR(1) = iTabAI(2,jAtom_)
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'E,R=',AtomLbl(iAtom),ChOp(iDCR(4)),AtomLbl(jAtom),ChOp(iDCR(1))
#   endif
    nCoBond_j = nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nFgBond_j = nFgBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    if ((nCoBond_j > 1) .and. (nFgBond_j == 0)) cycle
    if (iDCR(4) /= iOper(0)) cycle

    ! R

    if (R_Stab_A(iDCR(1),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(1) /= iOper(0))) cycle

    A(:,4) = Cx(:,iAtom,iIter)
    Ref(:,4) = Cx(:,iAtom,iRef)
    Prv(:,4) = Cx(:,iAtom,iPrv)

    ! Form double coset representatives for (iAtom,jAtom)

    call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
    kDCRR = iDCR(1)
#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
    write(u6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
    write(u6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
    write(u6,'(2A)') 'R=',ChOp(kDCRR)
#   endif

    call OA(kDCRR,Cx(:,jAtom,iIter),A(:,1))
    call OA(kDCRR,Cx(:,jAtom,iRef),Ref(:,1))
    call OA(kDCRR,Cx(:,jAtom,iPrv),Prv(:,1))

    ! Form stabilizer for (iAtom,jAtom)

    call Inter(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iStabM,nStabM)
#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#   endif

    if (Help) then
      f_Const_ij_Ref = rko
      f_Const_ij = rko
    else
      r0 = rAV(ir,jr)
      Alpha = aAv(ir,jr)
      rij2_Ref = (Ref(1,4)-Ref(1,1))**2+(Ref(2,4)-Ref(2,1))**2+(Ref(3,4)-Ref(3,1))**2
      f_Const_ij_Ref = rko*exp(Alpha*(r0**2-rij2_Ref))
      rij2 = (A(1,4)-A(1,1))**2+(A(2,4)-A(2,1))**2+(A(3,4)-A(3,1))**2
      f_Const_ij = rko*exp(Alpha*(r0**2-rij2))
    end if

    nNeighbor_i = iTabAtoms(1,0,iAtom_)
    do kNeighbor=1,nNeighbor_i
      kAtom_ = iTabAtoms(1,kNeighbor,iAtom_)
      if (kAtom_ == jAtom_) cycle
      kBond = iTabAtoms(2,kNeighbor,iAtom_)
      kBondType = iTabBonds(3,kBond)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'kBond,kBondType=',kBond,kBondType
#     endif
      if (kBondType == vdW_Bond) cycle
      if (kBondType > Magic_Bond) cycle
      if (kBond == jBond) cycle

      kAtom = iTabAI(1,kAtom_)
      kr = iTabRow(ANr(kAtom))
      Ind(2) = kAtom
      iDCR(2) = iTabAI(2,kAtom_)

      if (iDCR(1) == iOper(0)) then
        if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom)) .and. &
            (iDCR(2) /= iOper(0))) cycle
      else
        if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) cycle
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'T=',AtomLbl(kAtom),ChOp(iDCR(2))
      write(u6,*) 'kAtom=',kAtom
#     endif

      do lNeighbor=1,nNeighbor_i
        lAtom_ = iTabAtoms(1,lNeighbor,iAtom_)
        if (lAtom_ == jAtom_) cycle
        if (lAtom_ <= kAtom_) cycle
        lBond = iTabAtoms(2,lNeighbor,iAtom_)
        lBondType = iTabBonds(3,lBond)
        if (lBondType == vdW_Bond) cycle
        if (lBondType > Magic_Bond) cycle
        if (lBond == jBond) cycle
        if (lBond == kBond) cycle
        lAtom = iTabAI(1,lAtom_)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'lBond,lBondType=',lBond,lBondType
        write(u6,*) 'lAtom=',lAtom
#       endif

        lr = iTabRow(ANr(lAtom))
        Ind(3) = lAtom
        iDCR(3) = iTabAI(2,lAtom_)
        !if (kAtom > lAtom) cycle

        if (iDCR(1) == iOper(0)) then
          if (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)) .and. &
              (iDCR(3) /= iOper(0)) .and. (iDCR(2) /= iOper(0))) cycle
        else
          if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(3) /= iOper(0)) .and. (iDCR(2) /= iOper(0))) cycle
        end if

        if (kAtom == lAtom) then
          if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom)) .and. &
              (iDCR(2) /= iOper(0))) cycle
          if (iDCR(3) == iOper(0)) cycle
          if (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)) .and. &
              (iDCR(3) /= iOper(0)) .and. (iDCR(2) /= iOper(0))) cycle
        end if
#       ifdef _DEBUGPRINT_
        write(u6,*)
        write(u6,*) 'TS=',AtomLbl(lAtom),ChOp(iDCR(3))
#       endif

        Help = (ir > 3) .or. (jr > 3) .or. (kr > 3) .or. (lr > 3)

        write(Label,'(A,I2,A,I2,A,I2,A,I2,A)') 'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'

        ! Form double coset representatives for (kAtom,lAtom)

        call DCR(Lambda,jStab(0,kAtom),nStab(kAtom),jStab(0,lAtom),nStab(lAtom),iDCRS,nDCRS)
        kDCRS = ieor(iDCR(2),iDCR(3))

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'W={',(ChOp(jStab(i,kAtom)),i=0,nStab(kAtom)-1),'}  '
        write(u6,'(10A)') 'X={',(ChOp(jStab(i,lAtom)),i=0,nStab(lAtom)-1),'}  '
        write(u6,'(10A)') 'S={',(ChOp(iDCRS(i)),i=0,nDCRS-1),'}  '
        write(u6,'(2A)') 'S=',ChOp(kDCRS)
#       endif

        Ref(:,2) = Cx(:,kAtom,iRef)
        call OA(kDCRS,Cx(:,lAtom,iRef),Ref(:,3))
        Prv(:,2) = Cx(:,kAtom,iPrv)
        call OA(kDCRS,Cx(:,lAtom,iPrv),Prv(:,3))

        ! Form stabilizer for (kAtom,lAtom)

        call Inter(jStab(0,kAtom),nStab(kAtom),jStab(0,lAtom),nStab(lAtom),iStabN,nStabN)

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
#       endif

        ! Form double coset representatives for
        ! ((iAtom,jAtom),(kAtom,lAtom))

        call DCR(Lambda,iSTabM,nStabM,iStabN,nStabN,iDCRT,nDCRT)

        ! Take care of some special cases which normally
        ! are not included. If A=B we will normally exclude
        ! the pairs R(A)-A and TS(C)-T(C).

        iDCRX(0:nDCRT-1) = iDCRT(0:nDCRT-1)
        iDCRY(0:nDCRT-1) = iDCRT(0:nDCRT-1)
        nDCRX = nDCRT
        nDCRY = nDCRT
        if (iAtom == jAtom) then
          !write(u6,*) ' Special fix'
          call Union(iDCRX,nDCRX,iDCRY,nDCRY,kDCRR,iDCRT,nDCRT)
        else if (kAtom == lAtom) then
          !write(u6,*) ' Special fix'
          call Union(iDCRX,nDCRX,iDCRY,nDCRY,kDCRS,iDCRT,nDCRT)
        end if

        kDCRT = iDCR(2)
        kDCRTS = iDCR(3)

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'T={',(ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
        write(u6,'(2A)') 'T=',ChOp(kDCRT)
#       endif

        call OA(kDCRT,Cx(:,kAtom,iIter),A(:,2))
        call OA(kDCRT,Cx(:,kAtom,iRef),Ref(:,2))
        call OA(kDCRT,Cx(:,kAtom,iPrv),Prv(:,2))
        call OA(kDCRTS,Cx(:,lAtom,iIter),A(:,3))
        call OA(kDCRTS,Cx(:,lAtom,iRef),Ref(:,3))
        call OA(kDCRTS,Cx(:,lAtom,iPrv),Prv(:,3))

        ! Form the stabilizer for the out-of-plane

        if ((iAtom == lAtom) .and. (jAtom == kAtom) .and. (kDCRR == kDCRS)) then
          call Union(iStabM,nStabM,iStabN,nStabN,kDCRTS,iStabO,nStabO)
        else
          call Inter(iStabM,nStabM,iStabN,nStabN,iStabO,nStabO)
        end if

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
        write(u6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
        write(u6,'(10A)') 'O={',(ChOp(iStabO(i)),i=0,nStabO-1),'}  '
        write(u6,*) 'jAtom,iAtom,kAtom,lAtom=',jAtom,iAtom,kAtom,lAtom
#       endif

        ! Compute the degeneracy of the torsion

        iDeg = nIrrep/nStabO
        Deg = sqrt(real(iDeg,kind=wp))

        ! Test if coordinate should be included

        if (Help) then
          f_Const_ijk_Ref = f_Const_ij_Ref
          f_Const_ijk = f_Const_ij
          f_Const_Ref = f_Const_ijk_Ref
          f_Const = f_Const_ijk
        else

          ! Test the ik-pair

          r0 = rAV(ir,kr)
          Alpha = aAv(ir,kr)
          rik2_Ref = (Ref(1,4)-Ref(1,2))**2+(Ref(2,4)-Ref(2,2))**2+(Ref(3,4)-Ref(3,2))**2
          f_Const_ijk_Ref = f_Const_ij_Ref*exp(Alpha*(r0**2-rik2_Ref))
          rik2 = (A(1,4)-A(1,2))**2+(A(2,4)-A(2,2))**2+(A(3,4)-A(3,2))**2
          f_Const_ijk = f_Const_ij*exp(Alpha*(r0**2-rik2))

          ! Test the il-pair

          r0 = rAV(ir,lr)
          Alpha = aAv(ir,lr)
          ril2_Ref = (Ref(1,4)-Ref(1,3))**2+(Ref(2,4)-Ref(2,3))**2+(Ref(3,4)-Ref(3,3))**2
          f_Const_Ref = f_Const_ijk_Ref*exp(Alpha*(r0**2-ril2_Ref))
          ril2 = (A(1,4)-A(1,3))**2+(A(2,4)-A(2,3))**2+(A(3,4)-A(3,3))**2
          f_Const = f_Const_ijk*exp(Alpha*(r0**2-ril2))
        end if
        if ((f_Const_Ref < f_Const_Min) .and. (jBondtype /= Fragments_Bond) .and. (kBondtype /= Fragments_Bond) .and. &
            (lBondtype /= Fragments_Bond)) cycle

        ! Check that valence angles are above threshold

        mCent = 3
        delta0 = 45.0_wp*deg2rad
        !if ((jBondType == Fragments_Bond) .or. (kBondType == Fragments_Bond) .or. (lBondType == Fragments_Bond)) delta = Zero

        ! 1-4-2

        RX4Y(:,1) = Ref(:,1)
        RX4Y(:,2) = Ref(:,4)
        RX4Y(:,3) = Ref(:,2)
        call Bend(RX4Y,mCent,Fi2,Grad_ref,.false.,.false.,'        ',Hess,.false.)
#       ifdef _DEBUGPRINT_
        write(u6,*) '1-4-2: Fi2=',Fi2
#       endif
        delta = delta0
        if ((jBondType == Fragments_Bond) .or. (kBondType == Fragments_Bond)) delta = Zero
        if (Fi2 > Pi-delta) cycle
        if (Fi2 < delta) cycle

        ! 1-4-3

        RX4Y(:,3) = Ref(:,3)
        call Bend(RX4Y,mCent,Fi3,Grad_ref,.false.,.false.,'        ',Hess,.false.)
#       ifdef _DEBUGPRINT_
        write(u6,*) '1-4-3: Fi3=',Fi3
#       endif
        delta = delta0
        if ((jBondType == Fragments_Bond) .or. (lBondType == Fragments_Bond)) delta = Zero
        if (Fi3 > Pi-delta) cycle
        if (Fi3 < delta) cycle

        ! 2-4-3

        RX4Y(:,1) = Ref(:,2)
        call Bend(RX4Y,mCent,Fi4,Grad_ref,.false.,.false.,'        ',Hess,.false.)
#       ifdef _DEBUGPRINT_
        write(u6,*) '2-4-3: Fi4=',Fi4
#       endif
        delta = delta0
        if ((kBondType == Fragments_Bond) .or. (lBondType == Fragments_Bond)) delta = Zero
        if (Fi4 > Pi-delta) cycle
        if (Fi4 < delta) cycle

        call OutofP(Ref,nCent,Val,Grad,.false.,.false.,'        ',Hess,.false.)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Val=',Val/deg2rad
#       endif

        if (abs(Val) > 35.0_wp*deg2rad) cycle

        call OutofP(A,nCent,Val,Grad,.false.,.false.,'        ',Hess,Proc_dB)

        nq = nq+1
        if (.not. Process) mB_Tot = mB_Tot+mB
        if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2
#       ifdef _DEBUGPRINT_
        write(u6,*) 'nq=',nq
#       endif

        nqO = nqO+1
        iF1 = 1
        call NxtWrd(AtomLbl(iAtom),iF1,iE1)
        Lbls(1) = AtomLbl(iAtom)(iF1:iE1)

        iF2 = 1
        call NxtWrd(AtomLbl(jAtom),iF2,iE2)
        Lbls(2) = AtomLbl(jAtom)(iF2:iE2)
        if (kDCRR /= 0) then
          Lbls(2)(iE2+1:iE2+2+iChOp(kDCRR)) = '('//ChOp(kDCRR)(1:iChOp(kDCRR))//')'
          call NxtWrd(Lbls(2),iF2,iE2)
        end if

        iF3 = 1
        call NxtWrd(AtomLbl(kAtom),iF3,iE3)
        Lbls(3) = AtomLbl(kAtom)(iF3:iE3)
        if (kDCRT /= 0) then
          Lbls(3)(iE3+1:iE3+2+iChOp(kDCRT)) = '('//ChOp(kDCRT)(1:iChOp(kDCRT))//')'
          call NxtWrd(Lbls(3),iF3,iE3)
        end if

        iF4 = 1
        call NxtWrd(AtomLbl(lAtom),iF4,iE4)
        Lbls(4) = AtomLbl(lAtom)(iF4:iE4)
        if (kDCRTS /= 0) then
          Lbls(4)(iE4+1:iE4+2+iChOp(kDCRTS)) = '('//ChOp(kDCRTS)(1:iChOp(kDCRTS))//')'
          call NxtWrd(Lbls(4),iF4,iE4)
        end if
        write(LuIC,'(A,I3.3,8A)') 'o',nqO,' = Outofp   ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3),' ',Lbls(4)(iF4:iE4),' ', &
                                  Lbls(1)(iF1:iE1)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,I3.3,8A)') 'o',nqO,' = Outofp   ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3),' ',Lbls(4)(iF4:iE4),' ', &
                                Lbls(1)(iF1:iE1)
        write(u6,*) 'iDeg=',iDeg
#       endif
        Label = ' '
        write(Label,'(A,I3.3)') 'o',nqO

        if (Process) then

          Indq(1,nq) = 6
          ij = (jAtom-1)*nsAtom+iAtom
          kl = (lAtom-1)*nsAtom+kAtom
          Indq(2,nq) = (kl-1)*nsAtom**2+ij
          ijDCR = kDCRT*8+kDCRR+1
          Indq(3,nq) = kDCRS*8**2+ijDCR

          !f_Const = Max(f_Const,f_Const_Min)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'f_const=',f_const
#         endif
          fconst(nq) = sqrt(f_Const)
          rMult(nq) = Deg

          Valu(nq,iIter) = Val
          qLbl(nq) = Label

          ! Project the gradient vector

          call ProjSym(nCent,Ind,A,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq,rMult(nq))

        end if

      end do            ! lNeighbor
    end do              ! kNeighbor
  end do                ! iCase
end do                  ! jBond

return

end subroutine OutOfPlane_List

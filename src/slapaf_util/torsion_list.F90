!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Torsion_List(nq,nsAtom,iIter,nIter,Cx,Process,Valu,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv,Proc_dB,iTabBonds,nBonds, &
                        iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: ANr, AtomLbl, Fragments_Bond, jStab, Magic_Bond, nStab, vdW_Bond
use ddvdt, only: A_Trsn, aAV, f_Const_Min, rAV, rkt
use Constants, only: Zero, One, Two, Ten, Pi, Angstrom, deg2rad
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Slapaf_Info, only: BondType
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, iRef, LuIC, iPrv, nBonds, iTabBonds(3,nBonds), mAtoms, &
                                 iTabAI(2,mAtoms), nMax, iTabAtoms(2,0:nMax,mAtoms), nB_Tot, ndB_Tot
integer(kind=iwp), intent(inout) :: nq, Indq(3,nB), mB_Tot, mdB_Tot, iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter)
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Valu(nB,nIter), fconst(nB), rMult(nB), BM(nB_Tot), dBM(ndB_Tot)
character(len=14), intent(inout) :: qLbl(nB)
#include "Molcas.fh"
integer(kind=iwp), parameter :: mB = 4*3
integer(kind=iwp) :: iAtom, iAtom_, iBond, iBondType, iCase, iDCR(4), iDCRR(0:7), iDCRS(0:7), iDCRT(0:7), iDCRX(0:7), iDCRY(0:7), &
                     iDeg, iE1, iE2, iE3, iE4, iF1, iF2, iF3, iF4, ij, ijDCR, iMagic, Ind(4), iNeighbor, ir, iStabM(0:7), &
                     iStabN(0:7), iStabO(0:7), jAtom, jAtom_, jBond, jBondType, jr, kAtom, kAtom_, kBond, kBondType, kDCRR, kDCRS, &
                     kDCRT, kDCRTS, kl, kr, Lambda, lAtom, lAtom_, lNeighbor, lr, mCent, mE, nCent, nCoBond_j, nCoBond_k, nDCRR, &
                     nDCRS, nDCRT, nDCRX, nDCRY, nE, nFgBond_j, nFgBond_k, nNeighbor_j, nNeighbor_k, nqT, nStabM, nStabN, nStabO
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: A(3,4), Alpha, CosFact, CosFi, CosThr, Deg, delta, Diff, f_Const, f_Const_ij, f_Const_ij_Ref, f_Const_ijk, &
                 f_Const_ijk_Ref, f_Const_Ref, Fact, Fi2, Fi3, Grad(mB), Grad_ref(9), Hess(mB**2), Prv(3,4), r0, Range1, Range2, &
                 Range3, Rbc, RbcCov, Ref(3,4), rij2, rij2_Ref, rjk2, rjk2_Ref, rkl2, rkl2_Ref, Val, Val_Prv
logical(kind=iwp) :: Help, MinBas
character(len=LenIn4) :: Lbls(4)
character(len=14) :: Label
integer(kind=iwp), parameter :: iChOp(0:7) = [1,1,1,2,1,2,2,3]
real(kind=wp), parameter :: f_Const_Min2 = 1.0e-1_wp
character(len=*), parameter :: ChOp(0:7) = ['E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ']
integer(kind=iwp), external :: iTabRow, nCoBond, nFgBond
real(kind=wp), external :: CovRadT
logical(kind=iwp), external :: R_Stab_A, Torsion_Check

!                                                                      *
!***********************************************************************
!                                                                      *
if (nBonds < 3) return

nqT = 0
Hess(:) = Zero

! Loop over dihedrals

MinBas = .false.
if (MinBas) then
  Fact = 1.3_wp
else
  Fact = One
end if
nCent = 4

! Order will play a role here. That is, the torsion
! A-R(B)-T(C)-TS(D) is NOT identical to
! A-R(B)-TS(D)-T(C). Hence we put no restriction on the
! pairs AB and CD. However, for the pair of pairs we have
! that order is irrelevant, i.e. ABCD is identical to
! DCBA. To guarantee this we limit the pairs to the unique
! combinations.

! Start with the center bond: B-C

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'List of all available bonds'
write(u6,*)
do iBond=1,nBonds
  iBondType = iTabBonds(3,iBond)
  write(u6,*)
  jAtom_ = iTabBonds(1,iBond)
  kAtom_ = iTabBonds(2,iBond)
  jAtom = iTabAI(1,jAtom_)
  kAtom = iTabAI(1,kAtom_)
  write(u6,*) 'Atoms-pair:  ',AtomLbl(jAtom),AtomLbl(kAtom)
  write(u6,*) 'iBond,iBondType=',iBond,Bondtype(min(3,iBondType))
end do
#endif
do iBond=1,nBonds

  ! The center bond may be a "magic" bond
  iBondType = iTabBonds(3,iBond)
# ifdef _DEBUGPRINT_
  write(u6,*)
  jAtom_ = iTabBonds(1,iBond)
  kAtom_ = iTabBonds(2,iBond)
  jAtom = iTabAI(1,jAtom_)
  kAtom = iTabAI(1,kAtom_)
  write(u6,*) 'Atoms-pair:',AtomLbl(jAtom),AtomLbl(kAtom)
  write(u6,*) 'iBond,iBondType=',iBond,Bondtype(min(3,iBondType))
# endif

  ! Center bond should not be a van der Waals bond,
  ! anything else goes!

  if (iBondType == vdW_Bond) cycle

  ! Extract index to the center atom in an "Magic" bond if this
  ! is a magic bond.

  if (iBondType > Magic_Bond) then
    iMagic = iBondType-3
  else
    iMagic = 0
  end if

  ! cases: BC or CB

  do iCase=1,2

    if (iCase == 1) then
      jAtom_ = iTabBonds(1,iBond)
      kAtom_ = iTabBonds(2,iBond)
    else
      jAtom_ = iTabBonds(2,iBond)
      kAtom_ = iTabBonds(1,iBond)
    end if

    jAtom = iTabAI(1,jAtom_)
    kAtom = iTabAI(1,kAtom_)
    jr = iTabRow(ANr(jAtom))
    kr = iTabRow(ANr(kAtom))
    Ind(2) = jAtom
    Ind(3) = kAtom
    iDCR(2) = iTabAI(2,jAtom_)
    if (R_Stab_A(iDCR(2),jStab(0,jAtom),nStab(jAtom))) iDCR(2) = iOper(0)
    iDCR(3) = iTabAI(2,kAtom_)
    if (R_Stab_A(iDCR(3),jStab(0,kAtom),nStab(kAtom))) iDCR(3) = iOper(0)
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'R,T=',AtomLbl(jAtom),ChOp(iDCR(2)),AtomLbl(kAtom),ChOp(iDCR(3))
    write(u6,*)
#   endif

    nNeighbor_j = iTabAtoms(1,0,jAtom_)
    nCoBond_j = nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nFgBond_j = nFgBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nNeighbor_k = iTabAtoms(1,0,kAtom_)
    nCoBond_k = nCoBond(kAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nFgBond_k = nFgBond(kAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    if ((nCoBond_j < 2) .and. (nFgBond_j == 0)) cycle
    if ((nCoBond_k < 2) .and. (nFgBond_k == 0)) cycle

#   ifdef _DEBUGPRINT_
    write(u6,*) 'nNeighbor_j,nNeighbor_k=',nNeighbor_j,nNeighbor_k
    write(u6,*)
#   endif

    do iNeighbor=1,nNeighbor_j
      iAtom_ = iTabAtoms(1,iNeighbor,jAtom_)
      !nCoBond_i = nCoBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
      !nFgBond_i = nFgBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
      jBond = iTabAtoms(2,iNeighbor,jAtom_)
      if (jBond == iBond) cycle
      jBondType = iTabBonds(3,jBond)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'jBond,jBondType=',jBond,BondType(min(3,jBondType))
#     endif
      if ((jBondType == vdW_Bond) .or. (jBondType == Magic_Bond)) cycle
      !if ((nCoBond_j > 2) .and. (nCoBond_i >= 4) .and. (nFgBond_i == 0)) cycle
      !if ((nCoBond_i >= 8) .and. (nCoBond_j >= 8) .and. (nCoBond_k >= 8)) cycle
      iAtom = iTabAI(1,iAtom_)
      if ((iBondType > Magic_Bond) .and. (iAtom == iMagic)) cycle
      ir = iTabRow(ANr(iAtom))
      Ind(1) = iAtom
      iDCR(1) = iTabAI(2,iAtom_)
      if (R_Stab_A(iDCR(1),jStab(0,iAtom),nStab(iAtom))) iDCR(1) = iOper(0)

      ! Torsion should be A-..., eliminate P(A)-...

      if (iDCR(1) /= iOper(0)) cycle

      ! Eliminate A-R(B)-T(C)-TS(C) over A-B-RT(C)-RTS(D)
      ! Proceed if A-R-T(C)-TS(C)

      if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) cycle

      ! Eliminate A-R(B)-T(C)-TS(D) over A-R(B)-C-S(D)
      ! Proceed if A-R(B)-C-S(D)

      if (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)) .and. &
          (iDCR(3) /= iOper(0))) cycle
#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'E=',AtomLbl(iAtom),ChOp(iDCR(1))
#     endif
      A(:,1) = Cx(:,iAtom,iIter)
      Ref(:,1) = Cx(:,iAtom,iRef)
      Prv(:,1) = Cx(:,iAtom,iPrv)
      Help = (ir > 3) .or. (jr > 3)

      ! Form double coset representatives for (iAtom,jAtom)

      call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
      kDCRR = iDCR(2)
#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
      write(u6,'(2A)') 'R=',ChOp(kDCRR)
#     endif
      call OA(kDCRR,Cx(:,jAtom,iIter),A(:,2))
      call OA(kDCRR,Cx(:,jAtom,iRef),Ref(:,2))
      call OA(kDCRR,Cx(:,jAtom,iPrv),Prv(:,2))
#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
      write(u6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
#     endif

      ! Form stabilizer for (iAtom,jAtom)

      call Inter(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iStabM,nStabM)
#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#     endif
      if (Help) then
        rij2 = (Ref(1,1)-Ref(1,2))**2+(Ref(2,1)-Ref(2,2))**2+(Ref(3,1)-Ref(3,2))**2
        r0 = Zero
        Alpha = Zero
        f_Const_ij = Zero
        f_Const_ij_Ref = f_Const_ij
      else
        r0 = rAV(ir,jr)
        Alpha = aAv(ir,jr)
        rij2_Ref = (Ref(1,1)-Ref(1,2))**2+(Ref(2,1)-Ref(2,2))**2+(Ref(3,1)-Ref(3,2))**2
        f_Const_ij_Ref = rkt*exp(Alpha*(r0**2-rij2_Ref))
        rij2 = (A(1,1)-A(1,2))**2+(A(2,1)-A(2,2))**2+(A(3,1)-A(3,2))**2
        f_Const_ij = rkt*exp(Alpha*(r0**2-rij2))
      end if

      do lNeighbor=1,nNeighbor_k
        lAtom_ = iTabAtoms(1,lNeighbor,kAtom_)
        !nCoBond_l = nCoBond(lAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
        !nFgBond_l = nFgBond(lAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
        lAtom = iTabAI(1,lAtom_)
        Ind(4) = lAtom
        iDCR(4) = iTabAI(2,lAtom_)
        if (R_Stab_A(iDCR(4),jStab(0,lAtom),nStab(lAtom))) iDCR(4) = iOper(0)
        kBond = iTabAtoms(2,lNeighbor,kAtom_)
        if (kBond == iBond) cycle
        if (lAtom_ == iAtom_) cycle
        kBondType = iTabBonds(3,kBond)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'kBond,kBondType=',kBond,Bondtype(min(3,kBondType))
#       endif
        if ((kBondType == vdW_Bond) .or. (kBondType == Magic_Bond)) cycle
        !if ((nCoBond_k > 2) .and. (nCoBond_l >= 4) .and. (nFgBond_l == 0)) cycle
        !if ((nCoBond_j >= 8) .and. (nCoBond_k >= 8) .and. (nCoBond_l >= 8)) cycle

        if ((iBondType > Magic_Bond) .and. (lAtom == iMagic)) cycle
        lr = iTabRow(ANr(lAtom))
        kDCRT = iDCR(3)
        kDCRTS = iDCR(4)
        kDCRS = ieor(kDCRTS,kDCRT)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'i,j,k,l=',AtomLbl(iAtom),ChOp(iDCR(1)),AtomLbl(jAtom),ChOp(iDCR(2)),AtomLbl(kAtom),ChOp(iDCR(3)), &
                    AtomLbl(lAtom),ChOp(iDCR(4))
#       endif

        ! Eliminate A-R(B)-T(C)-TS(D) over A-TSR(B)-S(C)-D

        if (R_Stab_A(iDCR(4),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(4),jStab(0,jAtom),nStab(jAtom)) .and. &
            R_Stab_A(iDCR(4),jStab(0,kAtom),nStab(kAtom)) .and. (iDCR(4) /= iOper(0))) cycle

        nE = 1
        if (iDCR(2) == iOper(0)) nE = nE+1
        if (iDCR(3) == iOper(0)) nE = nE+1
        if (iDCR(4) == iOper(0)) nE = nE+1
        mE = 1
        if (R_Stab_A(iDCR(4),jStab(0,iAtom),nStab(iAtom))) mE = mE+1
        if (R_Stab_A(ieor(iDCR(4),iDCR(3)),jStab(0,kAtom),nStab(kAtom))) mE = mE+1
        if (R_Stab_A(ieor(iDCR(4),iDCR(2)),jStab(0,jAtom),nStab(jAtom))) mE = mE+1
        if (nE < mE) cycle
        if ((nE == mE) .and. (iAtom > lAtom)) cycle

#       ifdef _DEBUGPRINT_
        write(u6,*)
        write(u6,*) 'TS=',AtomLbl(lAtom),ChOp(iDCR(4))
#       endif

        Help = (ir > 3) .or. (jr > 3) .or. (kr > 3) .or. (lr > 3)

        write(Label,'(A,I2,A,I2,A,I2,A,I2,A)') 'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#       ifdef _DEBUGPRINT_
        write(u6,'(A,I2,A,I2,A,I2,A,I2,A)') 'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#       endif

        ! Form double coset representatives for (kAtom,lAtom)

        call DCR(Lambda,jStab(0,kAtom),nStab(kAtom),jStab(0,lAtom),nStab(lAtom),iDCRS,nDCRS)
#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'S={',(ChOp(iDCRS(i)),i=0,nDCRS-1),'}  '
        write(u6,'(10A)') 'W={',(ChOp(jStab(i,kAtom)),i=0,nStab(kAtom)-1),'}  '
        write(u6,'(10A)') 'X={',(ChOp(jStab(i,lAtom)),i=0,nStab(lAtom)-1),'}  '
        write(u6,'(2A)') 'S=',ChOp(kDCRS)
#       endif

        Ref(1:3,3) = Cx(:,kAtom,iRef)
        call OA(kDCRS,Cx(:,lAtom,iRef),Ref(:,4))
        Prv(1:3,3) = Cx(:,kAtom,iPrv)
        call OA(kDCRS,Cx(:,lAtom,iPrv),Prv(:,4))

        if (Help) then
          rkl2 = (Ref(1,3)-Ref(1,4))**2+(Ref(2,3)-Ref(2,4))**2+(Ref(3,3)-Ref(3,4))**2
        else
          r0 = rAV(kr,lr)
          rkl2 = (Ref(1,3)-Ref(1,4))**2+(Ref(2,3)-Ref(2,4))**2+(Ref(3,3)-Ref(3,4))**2
        end if

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

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'T={',(ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
        write(u6,'(2A)') 'kDCRT=',ChOp(kDCRT)
        write(u6,'(2A)') 'T=',ChOp(kDCRT)
#       endif

        kDCRTS = ieor(kDCRT,kDCRS)

        call OA(kDCRT,Cx(:,kAtom,iIter),A(:,3))
        call OA(kDCRT,Cx(:,kAtom,iRef),Ref(:,3))
        call OA(kDCRT,Cx(:,kAtom,iPrv),Prv(:,3))
        call OA(kDCRTS,Cx(:,lAtom,iIter),A(:,4))
        call OA(kDCRTS,Cx(:,lAtom,iRef),Ref(:,4))
        call OA(kDCRTS,Cx(:,lAtom,iPrv),Prv(:,4))

        ! Form the stabilizer for the torsion

        if ((iAtom == lAtom) .and. (jAtom == kAtom) .and. (kDCRR == kDCRS)) then
          call Union(iStabM,nStabM,iStabN,nStabN,kDCRTS,iStabO,nStabO)
        else
          call Inter(iStabM,nStabM,iStabN,nStabN,iStabO,nStabO)
        end if

#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
        write(u6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
        write(u6,'(10A)') 'O={',(ChOp(iStabO(i)),i=0,nStabO-1),'}  '
#       endif

        ! Compute the degeneracy of the torsion

        iDeg = nIrrep/nStabO
        Deg = sqrt(real(iDeg,kind=wp))

        ! Test if coordinate should be included

        if (Help) then
          rjk2 = (Ref(1,2)-Ref(1,3))**2+(Ref(2,2)-Ref(2,3))**2+(Ref(3,2)-Ref(3,3))**2
          Rbc = sqrt(rjk2)
          RbcCov = (CovRadT(ANr(jAtom))+CovRadT(ANr(kAtom)))/Angstrom
          Diff = RbcCov-Rbc
          if (Diff < Zero) Diff = Zero
          f_Const = A_Trsn(1)+A_Trsn(2)*Diff
          f_Const = f_Const*Fact
          r0 = Zero
          Alpha = Zero
          f_Const_Ref = f_Const
        else
          r0 = rAV(jr,kr)
          Alpha = aAv(jr,kr)
          rjk2_Ref = (Ref(1,2)-Ref(1,3))**2+(Ref(2,2)-Ref(2,3))**2+(Ref(3,2)-Ref(3,3))**2
          f_Const_ijk_Ref = f_Const_ij_Ref*exp(Alpha*(r0**2-rjk2_Ref))
          rjk2 = (A(1,2)-A(1,3))**2+(A(2,2)-A(2,3))**2+(A(3,2)-A(3,3))**2
          f_Const_ijk = f_Const_ij*exp(Alpha*(r0**2-rjk2))

          r0 = rAV(kr,lr)
          Alpha = aAv(kr,lr)
          rkl2_Ref = (Ref(1,3)-Ref(1,4))**2+(Ref(2,3)-Ref(2,4))**2+(Ref(3,3)-Ref(3,4))**2
          f_Const_Ref = f_Const_ijk_Ref*exp(Alpha*(r0**2-rkl2_Ref))
          rkl2 = (A(1,3)-A(1,4))**2+(A(2,3)-A(2,4))**2+(A(3,3)-A(3,4))**2
          f_Const = f_Const_ijk*exp(Alpha*(r0**2-rkl2))
        end if
        if (Torsion_Check(iAtom,jAtom,kAtom,lAtom,Ref,iTabAtoms,nMax,mAtoms)) f_Const_Ref = max(f_Const_Ref,Ten*F_Const_Min)

        if ((f_Const_Ref < f_Const_Min) .and. (iBondType /= Fragments_Bond) .and. (iBondType <= Fragments_Bond) .and. &
            (jBondType /= Fragments_Bond) .and. (kBondType /= Fragments_Bond)) cycle

        ! Check that valence angles are above threshold

        mCent = 3
        delta = 15.0_wp*deg2rad
        if (nsAtom == 4) delta = -Ten
        call Bend(Ref(1,1),mCent,Fi2,Grad_ref,.false.,.false.,'        ',Hess,.false.)
        if (Fi2 > Pi-delta) cycle
        if (Fi2 < delta) cycle
        call Bend(Ref(1,2),mCent,Fi3,Grad_ref,.false.,.false.,'        ',Hess,.false.)
        if (Fi3 > Pi-delta) cycle
        if (Fi3 < delta) cycle
        !write(u6,*) ' T Force Constant:',f_Const

        nq = nq+1
        if (.not. Process) mB_Tot = mB_Tot+mB
        if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2
        !write(u6,*) 'nq=',nq

        nqT = nqT+1
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
        write(LuIC,'(A,I3.3,8A)') 't',nqT,' = Dihedral ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3),' ', &
                                  Lbls(4)(iF4:iE4)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,I3.3,8A)') 't',nqT,' = Dihedral ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3),' ', &
                                Lbls(4)(iF4:iE4)
        write(u6,*) 'iDeg=',iDeg
#       endif
        Label = ' '
        write(Label,'(A,I3.3)') 't',nqT

        call Trsn(A,nCent,Val,Grad,.false.,.false.,'        ',Hess,Proc_dB)
        if (iIter == iPrv) then
          Val_Prv = Val
        else
          Val_Prv = Valu(nq,iPrv)
        end if

        ! correct for 2Pi flip relative to reference.

        Range1 = Two*Pi*0.8_wp
        Range2 = Pi*0.8_wp
        Range3 = Pi*1.2_wp
        if (abs(Val_Prv-Val) > Range1) then
          if (sign(One,Val_Prv) == One) then
            Val = Val+Two*Pi
          else
            Val = Val-Two*Pi
          end if
        else if ((abs(Val_Prv-Val) > Range2) .and. (abs(Val_Prv-Val) < Range3)) then
          if (Val_Prv-Val > Zero) then
            Val = Val+Pi
          else
            Val = Val-Pi
          end if
        end if
#       ifdef _DEBUGPRINT_
        call RecPrt('Trsns:  B',' ',Grad,3,4)
        call RecPrt('Trsns: dB',' ',Hess,12,12)
#       endif

        if (Process) then

          Indq(1,nq) = 5
          ij = (jAtom-1)*nsAtom+iAtom
          kl = (lAtom-1)*nsAtom+kAtom
          Indq(2,nq) = (kl-1)*nsAtom**2+ij
          ijDCR = kDCRT*8+kDCRR+1
          Indq(3,nq) = kDCRS*8**2+ijDCR

          if (iMagic /= 0) then
            f_Const = max(f_Const,f_Const_Min2)
          else if ((.not. Help) .and. &
                   ((iBondType == Fragments_Bond) .or. (jBondType == Fragments_Bond) .or. (kBondType == Fragments_Bond))) then
            f_Const = max(f_Const,f_Const_Min*1.0e3_wp)
          else
            f_Const = max(f_Const,f_Const_Min)
          end if
          fconst(nq) = sqrt(f_Const)
          rMult(nq) = Deg
          ! Scale down fconst if angles are close to linear
          CosFi = max(abs(cos(Fi2)),abs(cos(Fi3)))
          CosThr = 0.97_wp
          if (CosFi > CosThr) then
            CosFact = (CosFi-CosThr)/(One-CosThr)
            CosFact = One-(One-1.0e-3_wp)*CosFact**2
            fconst(nq) = CosFact*fconst(nq)
          end if

          Valu(nq,iIter) = Val
          qLbl(nq) = Label

          ! Project the gradient vector

          call ProjSym(nCent,Ind,A,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq,rMult(nq))

        end if

      end do    ! iNeighbor_k
    end do      ! iNeighbor_j
  end do        ! iCase
end do          ! iBonds
#ifdef _DEBUGPRINT_
write(u6,*) 'nqT=',nqT
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Torsion_List

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

subroutine Torsion_List(nq,nsAtom,iIter,nIter,Cx,Process,value,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv,Proc_dB,iTabBonds,nBonds, &
                        iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)

use Symmetry_Info, only: nIrrep, iOper
use Slapaf_Info, only: jStab, nStab, AtomLbl, ANr

implicit real*8(a-h,o-z)
#include "real.fh"
parameter(mB=4*3)
real*8 Cx(3,nsAtom,nIter), A(3,4), Grad(mB), Hess(mB**2), fconst(nB), value(nB,nIter), Ref(3,4), Prv(3,4), rMult(nB), Grad_ref(9), &
       BM(nB_Tot), dBM(ndB_Tot)
integer iDCRR(0:7), iStabM(0:7), Ind(4), iDCR(4), iDCRT(0:7), iDCRS(0:7), iStabN(0:7), iStabO(0:7), iChOp(0:7), Indq(3,nB), &
        iDCRX(0:7), iDCRY(0:7), nqB(nB), iTabBonds(3,nBonds), iTabAI(2,mAtoms), iTabAtoms(2,0:nMax,mAtoms), iBM(nB_Tot), &
        idBM(2,ndB_Tot)
logical Process, MinBas, Help, Proc_dB, R_Stab_A, Torsion_Check
character*14 Label, qLbl(nB)
character*3 ChOp(0:7)
#include "Molcas.fh"
character*(LenIn4) Lbls(4)
#include "bondtypes.fh"
#define _FMIN_
#include "ddvdt.fh"
#include "ddvdt_trsn.fh"
data ChOp/'E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ'/
data iChOp/1,1,1,2,1,2,2,3/
data f_Const_Min2/1.0D-1/
#include "constants.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (nBonds < 3) return

nqT = 0
call FZero(Hess,144)

! Loop over dihedrals

bohr = CONST_BOHR_RADIUS_IN_SI_*1.0D+10
MinBas = .false.
if (MinBas) then
  Fact = 1.3d0
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
write(6,*)
write(6,*) 'List of all available bonds'
write(6,*)
do iBond=1,nBonds
  iBondType = iTabBonds(3,iBond)
  write(6,*)
  jAtom_ = iTabBonds(1,iBond)
  kAtom_ = iTabBonds(2,iBond)
  jAtom = iTabAI(1,jAtom_)
  kAtom = iTabAI(1,kAtom_)
  write(6,*) 'Atoms-pair:  ',AtomLbl(jAtom),AtomLbl(kAtom)
  write(6,*) 'iBond,iBondType=',iBond,Bondtype(min(3,iBondType))
end do
#endif
do iBond=1,nBonds

  ! The center bond may be a "magic" bond
  iBondType = iTabBonds(3,iBond)
# ifdef _DEBUGPRINT_
  write(6,*)
  jAtom_ = iTabBonds(1,iBond)
  kAtom_ = iTabBonds(2,iBond)
  jAtom = iTabAI(1,jAtom_)
  kAtom = iTabAI(1,kAtom_)
  write(6,*) 'Atoms-pair:',AtomLbl(jAtom),AtomLbl(kAtom)
  write(6,*) 'iBond,iBondType=',iBond,Bondtype(min(3,iBondType))
# endif

  ! Center bond should not be a van der Waals bond,
  ! anything else goes!

  if (iBondType == vdW_Bond) Go To 201

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
    write(6,*)
    write(6,*) 'R,T=',AtomLbl(jAtom),ChOp(iDCR(2)),AtomLbl(kAtom),ChOp(iDCR(3))
    write(6,*)
#   endif

    nNeighbor_j = iTabAtoms(1,0,jAtom_)
    nCoBond_j = nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nFgBond_j = nFgBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nNeighbor_k = iTabAtoms(1,0,kAtom_)
    nCoBond_k = nCoBond(kAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    nFgBond_k = nFgBond(kAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    if ((nCoBond_j < 2) .and. (nFgBond_j == 0)) Go To 250
    if ((nCoBond_k < 2) .and. (nFgBond_k == 0)) Go To 250

#   ifdef _DEBUGPRINT_
    write(6,*) 'nNeighbor_j,nNeighbor_k=',nNeighbor_j,nNeighbor_k
    write(6,*)
#   endif

    do iNeighbor=1,nNeighbor_j
      iAtom_ = iTabAtoms(1,iNeighbor,jAtom_)
      !nCoBond_i = nCoBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
      !nFgBond_i = nFgBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
      jBond = iTabAtoms(2,iNeighbor,jAtom_)
      if (jBond == iBond) Go To 301
      jBondType = iTabBonds(3,jBond)
#     ifdef _DEBUGPRINT_
      write(6,*) 'jBond,jBondType=',jBond,BondType(min(3,jBondType))
#     endif
      if ((jBondType == vdW_Bond) .or. (jBondType == Magic_Bond)) Go To 301
      !if ((nCoBond_j > 2) .and. (nCoBond_i >= 4) .and. (nFgBond_i == 0)) Go To 301
      !if ((nCoBond_i >= 8) .and. (nCoBond_j >= 8) .and. (nCoBond_k >= 8)) Go To 301
      iAtom = iTabAI(1,iAtom_)
      if ((iBondType > Magic_Bond) .and. (iAtom == iMagic)) Go To 301
      ir = iTabRow(ANr(iAtom))
      Ind(1) = iAtom
      iDCR(1) = iTabAI(2,iAtom_)
      if (R_Stab_A(iDCR(1),jStab(0,iAtom),nStab(iAtom))) iDCR(1) = iOper(0)

      ! Torsion should be A-..., eliminate P(A)-...

      if (iDCR(1) /= iOper(0)) Go To 301

      ! Eliminate A-R(B)-T(C)-TS(C) over A-B-RT(C)-RTS(D)
      ! Proceed if A-R-T(C)-TS(C)

      if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) Go To 301

      ! Eliminate A-R(B)-T(C)-TS(D) over A-R(B)-C-S(D)
      ! Proceed if A-R(B)-C-S(D)

      if (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(3),jStab(0,jAtom),nStab(jAtom)) .and. &
          (iDCR(3) /= iOper(0))) Go To 301
#     ifdef _DEBUGPRINT_
      write(6,*)
      write(6,*) 'E=',AtomLbl(iAtom),ChOp(iDCR(1))
#     endif
      call dcopy_(3,Cx(1,iAtom,iIter),1,A,1)
      call dcopy_(3,Cx(1,iAtom,iRef),1,Ref,1)
      call dcopy_(3,Cx(1,iAtom,iPrv),1,Prv,1)
      Help = (ir > 3) .or. (jr > 3)

      ! Form double coset representatives for (iAtom,jAtom)

      call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRR,nDCRR)
      kDCRR = iDCR(2)
#     ifdef _DEBUGPRINT_
      write(6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
      write(6,'(2A)') 'R=',ChOp(kDCRR)
#     endif
      call OA(kDCRR,Cx(1:3,jAtom,iIter),A(1:3,2))
      call OA(kDCRR,Cx(1:3,jAtom,iRef),Ref(1:3,2))
      call OA(kDCRR,Cx(1:3,jAtom,iPrv),Prv(1:3,2))
#     ifdef _DEBUGPRINT_
      write(6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
      write(6,'(10A)') 'V={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
#     endif

      ! Form stabilizer for (iAtom,jAtom)

      call Inter(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iStabM,nStabM)
#     ifdef _DEBUGPRINT_
      write(6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
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
        if (kBond == iBond) Go To 401
        if (lAtom_ == iAtom_) Go To 401
        kBondType = iTabBonds(3,kBond)
#       ifdef _DEBUGPRINT_
        write(6,*) 'kBond,kBondType=',kBond,Bondtype(min(3,kBondType))
#       endif
        if ((kBondType == vdW_Bond) .or. (kBondType == Magic_Bond)) Go To 401
        !if ((nCoBond_k > 2) .and. (nCoBond_l >= 4) .and. (nFgBond_l == 0)) Go To 401
        !if ((nCoBond_j >= 8) .and. (nCoBond_k >= 8) .and. (nCoBond_l >= 8)) Go To 401

        if ((iBondType > Magic_Bond) .and. (lAtom == iMagic)) Go To 401
        lr = iTabRow(ANr(lAtom))
        kDCRT = iDCR(3)
        kDCRTS = iDCR(4)
        kDCRS = ieor(kDCRTS,kDCRT)
#       ifdef _DEBUGPRINT_
        write(6,*) 'i,j,k,l=',AtomLbl(iAtom),ChOp(iDCR(1)),AtomLbl(jAtom),ChOp(iDCR(2)),AtomLbl(kAtom),ChOp(iDCR(3)), &
                   AtomLbl(lAtom),ChOp(iDCR(4))
#       endif

        ! Eliminate A-R(B)-T(C)-TS(D) over A-TSR(B)-S(C)-D

        if (R_Stab_A(iDCR(4),jStab(0,iAtom),nStab(iAtom)) .and. R_Stab_A(iDCR(4),jStab(0,jAtom),nStab(jAtom)) .and. &
            R_Stab_A(iDCR(4),jStab(0,kAtom),nStab(kAtom)) .and. (iDCR(4) /= iOper(0))) Go To 401

        nE = 1
        if (iDCR(2) == iOper(0)) nE = nE+1
        if (iDCR(3) == iOper(0)) nE = nE+1
        if (iDCR(4) == iOper(0)) nE = nE+1
        mE = 1
        if (R_Stab_A(iDCR(4),jStab(0,iAtom),nStab(iAtom))) mE = mE+1
        if (R_Stab_A(ieor(iDCR(4),iDCR(3)),jStab(0,kAtom),nStab(kAtom))) mE = mE+1
        if (R_Stab_A(ieor(iDCR(4),iDCR(2)),jStab(0,jAtom),nStab(jAtom))) mE = mE+1
        if (nE < mE) Go To 401
        if ((nE == mE) .and. (iAtom > lAtom)) Go To 401

#       ifdef _DEBUGPRINT_
        write(6,*)
        write(6,*) 'TS=',AtomLbl(lAtom),ChOp(iDCR(4))
#       endif

        Help = (ir > 3) .or. (jr > 3) .or. (kr > 3) .or. (lr > 3)

        write(Label,'(A,I2,A,I2,A,I2,A,I2,A)') 'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#       ifdef _DEBUGPRINT_
        write(6,'(A,I2,A,I2,A,I2,A,I2,A)') 'D(',iAtom,',',jAtom,',',kAtom,',',lAtom,')'
#       endif

        ! Form double coset representatives for (kAtom,lAtom)

        call DCR(Lambda,jStab(0,kAtom),nStab(kAtom),jStab(0,lAtom),nStab(lAtom),iDCRS,nDCRS)
#       ifdef _DEBUGPRINT_
        write(6,'(10A)') 'S={',(ChOp(iDCRS(i)),i=0,nDCRS-1),'}  '
        write(6,'(10A)') 'W={',(ChOp(jStab(i,kAtom)),i=0,nStab(kAtom)-1),'}  '
        write(6,'(10A)') 'X={',(ChOp(jStab(i,lAtom)),i=0,nStab(lAtom)-1),'}  '
        write(6,'(2A)') 'S=',ChOp(kDCRS)
#       endif

        Ref(1:3,3) = Cx(1:3,kAtom,iRef)
        call OA(kDCRS,Cx(1:3,lAtom,iRef),Ref(1:3,4))
        Prv(1:3,3) = Cx(1:3,kAtom,iPrv)
        call OA(kDCRS,Cx(1:3,lAtom,iPrv),Prv(1:3,4))

        if (Help) then
          rkl2 = (Ref(1,3)-Ref(1,4))**2+(Ref(2,3)-Ref(2,4))**2+(Ref(3,3)-Ref(3,4))**2
        else
          r0 = rAV(kr,lr)
          rkl2 = (Ref(1,3)-Ref(1,4))**2+(Ref(2,3)-Ref(2,4))**2+(Ref(3,3)-Ref(3,4))**2
        end if

        ! Form stabilizer for (kAtom,lAtom)

        call Inter(jStab(0,kAtom),nStab(kAtom),jStab(0,lAtom),nStab(lAtom),iStabN,nStabN)

#       ifdef _DEBUGPRINT_
        write(6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
#       endif

        ! Form double coset representatives for
        ! ((iAtom,jAtom),(kAtom,lAtom))

        call DCR(Lambda,iSTabM,nStabM,iStabN,nStabN,iDCRT,nDCRT)

        ! Take care of some special cases which normally
        ! are not included. If A=B we will normally exclude
        ! the pairs R(A)-A and TS(C)-T(C).

        call iCopy(nDCRT,iDCRT,1,iDCRX,1)
        call iCopy(nDCRT,iDCRT,1,iDCRY,1)
        nDCRX = nDCRT
        nDCRY = nDCRT
        if (iAtom == jAtom) then
          !write(6,*) ' Special fix'
          call Union(iDCRX,nDCRX,iDCRY,nDCRY,kDCRR,iDCRT,nDCRT)
        else if (kAtom == lAtom) then
          !write(6,*) ' Special fix'
          call Union(iDCRX,nDCRX,iDCRY,nDCRY,kDCRS,iDCRT,nDCRT)
        end if

#       ifdef _DEBUGPRINT_
        write(6,'(10A)') 'T={',(ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
        write(6,'(2A)') 'kDCRT=',ChOp(kDCRT)
        write(6,'(2A)') 'T=',ChOp(kDCRT)
#       endif

        kDCRTS = ieor(kDCRT,kDCRS)

        call OA(kDCRT,Cx(1:3,kAtom,iIter),A(1:3,3))
        call OA(kDCRT,Cx(1:3,kAtom,iRef),Ref(1:3,3))
        call OA(kDCRT,Cx(1:3,kAtom,iPrv),Prv(1:3,3))
        call OA(kDCRTS,Cx(1:3,lAtom,iIter),A(1:3,4))
        call OA(kDCRTS,Cx(1:3,lAtom,iRef),Ref(1:3,4))
        call OA(kDCRTS,Cx(1:3,lAtom,iPrv),Prv(1:3,4))

        ! Form the stabilizer for the torsion

        if ((iAtom == lAtom) .and. (jAtom == kAtom) .and. (kDCRR == kDCRS)) then
          call Union(iStabM,nStabM,iStabN,nStabN,kDCRTS,iStabO,nStabO)
        else
          call Inter(iStabM,nStabM,iStabN,nStabN,iStabO,nStabO)
        end if

#       ifdef _DEBUGPRINT_
        write(6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
        write(6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
        write(6,'(10A)') 'O={',(ChOp(iStabO(i)),i=0,nStabO-1),'}  '
#       endif

        ! Compute the degeneracy of the torsion

        iDeg = nIrrep/nStabO
        Deg = sqrt(dble(iDeg))

        ! Test if coordinate should be included

        if (Help) then
          rjk2 = (Ref(1,2)-Ref(1,3))**2+(Ref(2,2)-Ref(2,3))**2+(Ref(3,2)-Ref(3,3))**2
          Rbc = sqrt(rjk2)
          RbcCov = (CovRadT(ANr(jAtom))+CovRadT(ANr(kAtom)))/bohr
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
        if (Torsion_Check(iAtom,jAtom,kAtom,lAtom,Ref,iTabAtoms,nMax,mAtoms)) then
          f_Const_Ref = max(f_Const_Ref,10.d0*F_Const_Min)
        end if

        if ((f_Const_Ref < f_Const_Min) .and. (iBondType /= Fragments_Bond) .and. (iBondType <= Fragments_Bond) .and. &
            (jBondType /= Fragments_Bond) .and. (kBondType /= Fragments_Bond)) Go To 401

        ! Check that valence angles are above threshold

        mCent = 3
        delta = (15.0d0/180.d0)*Pi
        if (nsAtom == 4) delta = -Ten
        call Bend(Ref(1,1),mCent,Fi2,Grad_ref,.false.,.false.,'        ',Hess,.false.)
        if (Fi2 > Pi-delta) Go To 401
        if (Fi2 < delta) Go To 401
        call Bend(Ref(1,2),mCent,Fi3,Grad_ref,.false.,.false.,'        ',Hess,.false.)
        if (Fi3 > Pi-delta) Go To 401
        if (Fi3 < delta) Go To 401
        !write(6,*) ' T Force Constant:',f_Const

        nq = nq+1
        if (.not. Process) mB_Tot = mB_Tot+mB
        if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2
        !write(6,*) 'nq=',nq

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
        write(6,'(A,I3.3,8A)') 't',nqT,' = Dihedral ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3),' ', &
                               Lbls(4)(iF4:iE4)
        write(6,*) 'iDeg=',iDeg
#       endif
        Label = ' '
        write(Label,'(A,I3.3)') 't',nqT

        call Trsn(A,nCent,Val,Grad,.false.,.false.,'        ',Hess,Proc_dB)
        if (iIter == iPrv) then
          Val_Prv = Val
        else
          Val_Prv = value(nq,iPrv)
        end if

        ! correct for 2Pi flip relative to reference.

        Range1 = Two*Pi*0.8d0
        Range2 = Pi*0.80d0
        Range3 = Pi*1.20d0
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
            f_Const = max(f_Const,f_Const_Min*1.d3)
          else
            f_Const = max(f_Const,f_Const_Min)
          end if
          fconst(nq) = sqrt(f_Const)
          rMult(nq) = Deg
          ! Scale down fconst if angles are close to linear
          CosFi = max(abs(cos(Fi2)),abs(cos(Fi3)))
          CosThr = 0.97d0
          if (CosFi > CosThr) then
            CosFact = (CosFi-CosThr)/(One-CosThr)
            CosFact = One-(One-1.0D-3)*CosFact**2
            fconst(nq) = CosFact*fconst(nq)
          end if

          value(nq,iIter) = Val
          qLbl(nq) = Label

          ! Project the gradient vector

          call ProjSym(nCent,Ind,A,iDCR,Grad,Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq,rMult(nq))

        end if

401     continue
      end do              ! iNeighbor_k
301   continue
    end do                ! iNeighbor_j
250 continue
  end do                  ! iCase
201 continue
end do                    ! iBonds
#ifdef _DEBUGPRINT_
write(6,*) 'nqT=',nqT
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Torsion_List

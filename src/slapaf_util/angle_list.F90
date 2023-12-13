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

subroutine Angle_List(nq,nsAtom,iIter,nIter,Cx,Process,Valu,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,Grad_all,iGlow,iGhi,iPrv,Proc_dB, &
                      iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB,Thr_small)

use Symmetry_Info, only: iOper, nIrrep
use Slapaf_Info, only: ANr, AtomLbl, Fragments_Bond, jStab, Magic_Bond, nStab, vdW_Bond
use ddvdt, only: A_Bend, aAV, f_Const_Min, rAV, rkf
use Constants, only: Zero, One, Two, Pi, deg2rad
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, iRef, LuIC, iGlow, iGhi, iPrv, nBonds, iTabBonds(3,nBonds), mAtoms, &
                                 iTabAI(2,mAtoms), nMax, iTabAtoms(2,0:nMax,mAtoms), nB_Tot, ndB_Tot
integer(kind=iwp), intent(inout) :: nq, Indq(3,nB), mB_Tot, mdB_Tot, iBM(nB_Tot), idBM(2,ndB_Tot), nqB(nB)
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter), Thr_small
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Valu(nB,nIter), fconst(nB), rMult(nB), Grad_all(9,iGlow:iGhi,nIter), BM(nB_Tot), dBM(ndB_Tot)
character(len=14), intent(inout) :: qLbl(nB)
#include "Molcas.fh"
integer(kind=iwp), parameter :: mB = 3*3
integer(kind=iwp) :: iAtom, iAtom_, iBond, iBondType, iDCR(3), iDCRR(0:7), iDCRT(0:7), ideg, iE1, iE2, iE3, iF1, iF2, iF3, Ind(3), &
                     iNeighbor, ir, iStabM(0:7), iStabN(0:7), jAtom, jAtom_, jBond, jBondType, jNeighbor, jr, k, kDCR, kDCRR, &
                     kDCRT, Lambda, mAtom, mAtom_, mi, mr, nCent, nCoBond_i, nCoBond_j, nCoBond_m, nDCRR, nDCRT, nk, nNeighbor_m, &
                     nqA, nStabM, nStabN
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: A(3,3), Alpha, Axis(3), BB, Deg, delta, f_Const, f_Const_Ref, Fact, Grad(mB), Grad_Ref(9), Hess(mB**2), &
                 Perp_Axis(3,2), Prv(3,3), r0, Ref(3,3), rim2, rim2_Ref, rmj2, rmj2_Ref, Val, Val_Ref
logical(kind=iwp) :: Help, MinBas
character(len=LenIn4) :: Lbls(3)
character(len=14) :: Label
integer(kind=iwp), parameter :: iChOp(0:7) = [1,1,1,2,1,2,2,3]
character(len=*), parameter :: ChOp(0:7) = ['E  ','X  ','Y  ','XY ','Z  ','XZ ','YZ ','XYZ']
integer(kind=iwp), external :: iTabRow, nCoBond
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: R_Stab_A

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (nBonds < 2) return

nqA = 0
#ifdef _DEBUGPRINT_
write(u6,*) ' Enter Bends.'
#endif
Hess(:) = Zero

! Loop over bends

MinBas = .false.
if (MinBas) then
  Fact = 1.3_wp
else
  Fact = One
end if
nCent = 3
!write(u6,*)

do mAtom_=1,mAtoms
  mAtom = iTabAI(1,mAtom_)
  mr = iTabRow(ANr(mAtom))
  Ind(2) = mAtom
  iDCR(2) = iTabAI(2,mAtom_)

  nNeighbor_m = iTabAtoms(1,0,mAtom_)
  nCoBond_m = nCoBond(mAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
  if (nNeighbor_m < 2) cycle

  do iNeighbor=1,nNeighbor_m
    iAtom_ = iTabAtoms(1,iNeighbor,mAtom_)
    iAtom = iTabAI(1,iAtom_)
    nCoBond_i = nCoBond(iAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
    ir = iTabRow(ANr(iAtom))
    Ind(1) = iAtom
    iDCR(1) = iTabAI(2,iAtom_)

    if (iDCR(1) /= iOper(0)) cycle
    if (R_Stab_A(iDCR(2),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(2) /= iOper(0))) cycle

    iBond = iTabAtoms(2,iNeighbor,mAtom_)
    iBondType = iTabBonds(3,iBond)
    if ((iBondType == vdW_Bond) .or. (iBondType > Magic_Bond)) cycle
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'iAtom,mAtom=',iAtom,mAtom
    write(u6,*) 'iBond,iBondType=',iBond,iBondType
    write(u6,*) 'E,R=',ChOp(iDCR(1)),ChOp(iDCR(2))
#   endif

    A(:,1) = Cx(:,iAtom,iIter)
    Ref(:,1) = Cx(:,iAtom,iRef)
    Prv(:,1) = Cx(:,iAtom,iPrv)

    do jNeighbor=1,nNeighbor_m
      jAtom_ = iTabAtoms(1,jNeighbor,mAtom_)
      jAtom = iTabAI(1,jAtom_)
      nCoBond_j = nCoBond(jAtom_,mAtoms,nMax,iTabBonds,nBonds,iTabAtoms)
      if ((nCoBond_i >= 8) .and. (nCoBond_j >= 8) .and. (nCoBond_m >= 8)) cycle

      jr = iTabRow(ANr(jAtom))
      Ind(3) = jAtom
      iDCR(3) = iTabAI(2,jAtom_)
      if (R_Stab_A(iDCR(3),jStab(0,iAtom),nStab(iAtom)) .and. (iDCR(3) /= iOper(0))) cycle
      if ((iDCR(3) == iOper(0)) .and. (iAtom >= jAtom)) cycle
      kDCR = ieor(iDCR(2),iDCR(3))
      if (R_Stab_A(kDCR,jStab(0,mAtom),nStab(mAtom)) .and. (iDCR(2) /= iOper(0))) cycle

      Help = (mr > 3) .or. (ir > 3) .or. (jr > 3)

      jBond = iTabAtoms(2,jNeighbor,mAtom_)
      jBondType = iTabBonds(3,jBond)
      if ((jBondType == vdW_Bond) .or. (jBondType > Magic_Bond)) cycle
#     ifdef _DEBUGPRINT_
      write(u6,*)
      write(u6,*) 'jAtom,mAtom=',jAtom,mAtom
      write(u6,*) 'jBond,jBondType=',jBond,jBondType
      write(u6,*) 'T=',ChOp(iDCR(3))
#     endif

      write(Label,'(A,I2,A,I2,A,I2,A)') 'A(',iAtom,',',mAtom,',',jAtom,')'

#     ifdef _DEBUGPRINT_
      call RecPrt('A',' ',Cx(:,iAtom,iIter),1,3)
      call RecPrt('B',' ',Cx(:,mAtom,iIter),1,3)
      call RecPrt('C',' ',Cx(:,jAtom,iIter),1,3)
#     endif

      ! Form double coset representatives for (iAtom,jAtom)

      call DCR(Lambda,jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCRT,nDCRT)
#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'T={',(ChOp(iDCRT(i)),i=0,nDCRT-1),'}  '
#     endif
      kDCRT = iDCR(3)

#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'U={',(ChOp(jStab(i,iAtom)),i=0,nStab(iAtom)-1),'}  '
      write(u6,'(10A)') 'V={',(ChOp(jStab(i,mAtom)),i=0,nStab(mAtom)-1),'}  '
      write(u6,'(10A)') 'X={',(ChOp(jStab(i,jAtom)),i=0,nStab(jAtom)-1),'}  '
      write(u6,'(2A)') 'T=',ChOp(kDCRT)
#     endif

      call OA(kDCRT,Cx(:,jAtom,iIter),A(:,3))
      call OA(kDCRT,Cx(:,jAtom,iRef),Ref(:,3))
      call OA(kDCRT,Cx(:,jAtom,iPrv),Prv(:,3))

      ! Form the stabilizer for (iAtom,jAtom)

      if (iAtom == jAtom) then
        call Union(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iDCR(3),iStabN,nStabN)
      else
        call Inter(jStab(0,iAtom),nStab(iAtom),jStab(0,jAtom),nStab(jAtom),iStabN,nStabN)
      end if

#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'N={',(ChOp(iStabN(i)),i=0,nStabN-1),'}  '
#     endif

      ! Form double coset representatives for ((iAtom,mAtom),jAtom)

      call DCR(Lambda,jStab(0,mAtom),nStab(mAtom),iStabN,nStabN,iDCRR,nDCRR)
      kDCRR = iDCR(2)

#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'R={',(ChOp(iDCRR(i)),i=0,nDCRR-1),'}  '
      write(u6,'(2A)') 'R=',ChOp(kDCRR)
#     endif

      call OA(kDCRR,Cx(:,mAtom,iIter),A(:,2))
      call OA(kDCRR,Cx(:,mAtom,iRef),Ref(:,2))
      call OA(kDCRR,Cx(:,mAtom,iPrv),Prv(:,2))

      ! Form the stabilizer for ((iAtom,mAtom),jAtom)

      call Inter(jStab(0,mAtom),nStab(mAtom),iStabN,nStabN,iStabM,nStabM)

#     ifdef _DEBUGPRINT_
      write(u6,'(10A)') 'M={',(ChOp(iStabM(i)),i=0,nStabM-1),'}  '
#     endif

      ! Compute the degeneracy of the angle

      ideg = nIrrep/nStabM
      Deg = sqrt(real(iDeg,kind=wp))
#     ifdef _DEBUGPRINT_
      write(u6,*) ' nIrrep,nStabM=',nIrrep,nStabM
#     endif

      ! Test if coordinate should be included
      ! All angles which are lower than the explicit threshold
      ! are rejected.

      if (Help) then
        rim2 = (Ref(1,1)-Ref(1,2))**2+(Ref(2,1)-Ref(2,2))**2+(Ref(3,1)-Ref(3,2))**2
        rmj2 = (Ref(1,2)-Ref(1,3))**2+(Ref(2,2)-Ref(2,3))**2+(Ref(3,2)-Ref(3,3))**2
        if ((ir == 1) .or. (jr == 1)) then
          f_Const = A_Bend(1)
        else
          f_Const = A_Bend(2)
        end if
        f_Const = f_Const*Fact
        f_Const_Ref = f_Const
      else
        r0 = rAV(ir,mr)
        Alpha = aAv(ir,mr)
        rim2_Ref = (Ref(1,1)-Ref(1,2))**2+(Ref(2,1)-Ref(2,2))**2+(Ref(3,1)-Ref(3,2))**2
        f_Const_Ref = rkf*exp(Alpha*(r0**2-rim2_Ref))
        rim2 = (A(1,1)-A(1,2))**2+(A(2,1)-A(2,2))**2+(A(3,1)-A(3,2))**2
        f_Const = rkf*exp(Alpha*(r0**2-rim2))

        r0 = rAV(mr,jr)
        Alpha = aAv(mr,jr)
        rmj2_Ref = (Ref(1,2)-Ref(1,3))**2+(Ref(2,2)-Ref(2,3))**2+(Ref(3,2)-Ref(3,3))**2
        f_Const_Ref = f_Const_Ref*exp(Alpha*(r0**2-rmj2_Ref))
        rmj2 = (A(1,2)-A(1,3))**2+(A(2,2)-A(2,3))**2+(A(3,2)-A(3,3))**2
        f_Const = f_Const*exp(Alpha*(r0**2-rmj2))

      end if
      if ((f_Const_Ref < f_Const_Min) .and. (iBondType /= Fragments_Bond) .and. (jBondType /= Fragments_Bond)) cycle
#     ifdef _DEBUGPRINT_
      write(u6,*) ' A Force Constant:',f_Const
      write(u6,*) iAtom,mAtom,jAtom,f_Const
#     endif

      call Bend(Ref,nCent,Val_Ref,Grad_Ref,.false.,.false.,'        ',Hess,Proc_dB)
      call Bend(A,nCent,Val,Grad,.false.,.false.,'        ',Hess,Proc_dB)

      ! Skip cases with a too small angle.

      if ((abs(Val_Ref) < Thr_Small) .and. (iBondType /= Fragments_Bond) .and. (jBondType /= Fragments_Bond)) cycle

      !                                                                *
      !*****************************************************************
      !                                                                *
      ! We'd like to avoid the problem with angles that have the
      ! value Pi. We introduce "linear" angles under two
      ! conditions.
      ! (1) the reference angle is within delta of Pi.
      ! (2) the actual angle is within 1.0e-11 of Pi.(?)
      ! Delta is set to 1.0e-11 if there are 3 atoms.(?)

      delta = 45.0_wp*deg2rad
      ! if (mAtoms == 3) delta = 1.0e-11_wp ! I do not understand
      ! although it is probably me who introduced it!

      if (((abs(Val_Ref-Pi) < Delta) .or. (abs(Val-Pi) < 1.0e-11_wp)) .and. &
          (.not. (((iBondType == Fragments_Bond) .or. (jBondType == Fragments_Bond)) .and. (mAtoms <= 4)))) then

        ! Reference is linear(a) or
        ! reference is NOT linear but the new structure is(b).

        if ((abs(Val-Pi) < 1.0e-11_wp) .and. (.not. (abs(Val_Ref-Pi) < Delta))) then
          ! Case b
          nk = 1
        else
          ! Case a
          nk = 2
        end if
        call CoSys(Prv,Axis,Perp_Axis)

        Label = ' '
        Label(1:1) = 'L'
        !do k=1,2
        do k=1,nk
          nq = nq+1
          if (.not. Process) mB_Tot = mB_Tot+mB
          if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2

          nqA = nqA+1
          iF1 = 1
          call NxtWrd(AtomLbl(iAtom),iF1,iE1)
          Lbls(1) = AtomLbl(iAtom)(iF1:iE1)
          iF2 = 1
          call NxtWrd(AtomLbl(mAtom),iF2,iE2)
          Lbls(2) = AtomLbl(mAtom)(iF2:iE2)
          if (kDCRR /= 0) then
            Lbls(2)(iE2+1:iE2+2+iChOp(kDCRR)) = '('//ChOp(kDCRR)(1:iChOp(kDCRR))//')'
            call NxtWrd(Lbls(2),iF2,iE2)
          end if
          iF3 = 1
          call NxtWrd(AtomLbl(jAtom),iF3,iE3)
          Lbls(3) = AtomLbl(jAtom)(iF3:iE3)
          if (kDCRT /= 0) then
            Lbls(3)(iE3+1:iE3+2+iChOp(kDCRT)) = '('//ChOp(kDCRT)(1:iChOp(kDCRT))//')'
            call NxtWrd(Lbls(3),iF3,iE3)
          end if
          write(LuIC,'(A,I3.3,A,I1.1,6A)') 'a',nqA,' = LAngle(',k,') ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,I3.3,A,I1.1,6A)') 'a',nqA,' = LAngle(',k,') ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3)
          write(u6,*) 'iDeg=',iDeg
#         endif
          Label = ' '
          write(Label,'(A,I3.3)') 'a',nqA

          if (Process) then

            call LBend(A,nCent,Val,Grad_all(:,nq,iIter),.false.,'        ',Hess,Proc_dB,Axis,Perp_Axis(1,k),(k == 2))

            ! Flip Angle value and gradient if needed!

            BB = DDot_(9,Grad_all(:,nq,iPrv),1,Grad_all(:,nq,iIter),1)
            if (BB < Zero) then
              !write(u6,*) ' Angle flips, corrected!'
              Val = Two*Pi-Val
              Grad_all(:,nq,iIter) = -Grad_all(:,nq,iIter)
              Hess(:) = -Hess(:)
            end if

            Indq(1,nq) = 2+k
            mi = (iAtom-1)*nsAtom+mAtom
            Indq(2,nq) = (jAtom-1)*nsAtom**2+mi
            Indq(3,nq) = kDCRT*8+kDCRR+1

            f_Const = max(f_Const,f_Const_Min)
            if ((.not. Help) .and. ((iBondType == Fragments_Bond) .or. (jBondType == Fragments_Bond))) f_Const = f_Const*1.0e3_wp
            fconst(nq) = sqrt(f_Const)
            rMult(nq) = Deg

            Valu(nq,iIter) = Val
            qLbl(nq) = Label

            ! Project the gradient vector

            call ProjSym(nCent,Ind,A,iDCR,Grad_all(:,nq,iIter),Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB, &
                         nq,rMult(nq))

          end if

        end do

      else      ! Non-linear case

        nq = nq+1
        if (.not. Process) mB_Tot = mB_Tot+mB
        if (.not. Proc_dB) mdB_Tot = mdB_Tot+mB**2

        nqA = nqA+1
        iF1 = 1
        call NxtWrd(AtomLbl(iAtom),iF1,iE1)
        Lbls(1) = AtomLbl(iAtom)(iF1:iE1)
        iF2 = 1
        call NxtWrd(AtomLbl(mAtom),iF2,iE2)
        Lbls(2) = AtomLbl(mAtom)(iF2:iE2)
        if (kDCRR /= 0) then
          Lbls(2)(iE2+1:iE2+2+iChOp(kDCRR)) = '('//ChOp(kDCRR)(1:iChOp(kDCRR))//')'
          call NxtWrd(Lbls(2),iF2,iE2)
        end if
        iF3 = 1
        call NxtWrd(AtomLbl(jAtom),iF3,iE3)
        Lbls(3) = AtomLbl(jAtom)(iF3:iE3)
        if (kDCRT /= 0) then
          Lbls(3)(iE3+1:iE3+2+iChOp(kDCRT)) = '('//ChOp(kDCRT)(1:iChOp(kDCRT))//')'
          call NxtWrd(Lbls(3),iF3,iE3)
        end if
        write(LuIC,'(A,I3.3,6A)') 'a',nqA,' = Angle ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3)
#       ifdef _DEBUGPRINT_
        write(u6,'(A,I3.3,6A)') 'a',nqA,' = Angle ',Lbls(1)(iF1:iE1),' ',Lbls(2)(iF2:iE2),' ',Lbls(3)(iF3:iE3)
#       endif
        Label = ' '
        write(Label,'(A,I3.3)') 'a',nqA

        if (Process) then

          call Bend(A,nCent,Val,Grad_all(:,nq,iIter),.false.,.false.,'        ',Hess,Proc_dB)

          ! Flip Angle value and gradient if needed!

          BB = DDot_(9,Grad_all(:,nq,iPrv),1,Grad_all(:,nq,iIter),1)
          if (BB < Zero) then
            !write(u6,*) ' Angle flips, corrected!'
            !write(u6,*) ' iRef,iIter=', iRef,iIter
            Val = Two*Pi-Val
            Grad_all(:,nq,iIter) = -Grad_all(:,nq,iIter)
            Hess(:) = -Hess(:)
          end if

          Indq(1,nq) = 2
          mi = (iAtom-1)*nsAtom+mAtom
          Indq(2,nq) = (jAtom-1)*nsAtom**2+mi
          Indq(3,nq) = kDCRT*8+kDCRR+1

          if ((.not. Help) .and. ((iBondType == Fragments_Bond) .or. (jBondType == Fragments_Bond))) f_Const = f_Const*1.0e2_wp
          fconst(nq) = sqrt(f_Const)
          rMult(nq) = Deg

          Valu(nq,iIter) = Val
          qLbl(nq) = Label

          ! Project the gradient vector

          call ProjSym(nCent,Ind,A,iDCR,Grad_all(:,nq,iIter),Hess,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,Proc_dB,nqB,nB,nq, &
                       rMult(nq))

        end if

      end if
      !                                                                *
      !*****************************************************************
      !                                                                *

    end do             ! End loop over jNeighbor
  end do               ! End loop over iCase
end do                 ! End loop over iBond
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Angle_List

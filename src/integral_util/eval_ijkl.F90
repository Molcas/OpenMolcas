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
! Copyright (C) 1991,1993,1999,2023, Roland Lindh                      *
!               1995, Martin Schuetz                                   *
!***********************************************************************

!#define _DEBUGPRINT_
!#define _DEBUGBREIT_
subroutine Eval_ijkl(iiS,jjS,kkS,llS,TInt,nTInt)
!***********************************************************************
!                                                                      *
!  Object: driver for two-electron integrals, parallel region          *
!          contains memory partitioning and loops over uncontracted    *
!          functions...                                                *
!                                                                      *
!  Input:                                                              *
!          iiS,jjS,kkS,llS     : shell indices                         *
!          TInt                : Computed Integrals                    *
!          nTInt               : dimension of TInt                     *
!          Integ_Proc          : subroutine for post processing        *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh / Martin Schuetz,                           *
!             Dept. of Theoretical Chemistry, University of Lund,      *
!             SWEDEN.                                                  *
!             Modified for k2 loop. August '91                         *
!             Modified for direct SCF. January '93                     *
!             Modified to minimize overhead for calculations with      *
!             small basis sets and large molecules. Sept. '93          *
!             parallel region split off in drvtwo.f, April '95         *
!             Total rehack May '99                                     *
!             Total rehack Aug '23                                     *
!***********************************************************************

use Index_Functions, only: iTri
use setup, only: mSkal, nAux, nSOs
use k2_structure, only: IndK2, k2data
use k2_arrays, only: Aux, Create_BraKet, DeDe, Destroy_Braket, FT, ipDijS, iSOSym, nDeDe, nFT, Sew_Scr
use iSD_data, only: iSD
use Basis_Info, only: Shells
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: DoFock, DoIntegrals, Map4
use Integral_interfaces, only: Int_PostProcess, twoel_kernel
#ifdef _DEBUGBREIT_
use Breit, only: nOrdOp
use UnixInfo, only: SuperName
#endif
use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_maxDBLE
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iiS, jjS, kkS, llS, nTInt
real(kind=wp), intent(inout) :: TInt(nTInt)
#include "ibas_ricd.fh"
integer(kind=iwp) :: iAngV(4), iAOst(4), iAOV(4), iBasAO, iBasi, iBasn, iBsInc, iCmpV(4), ijS, ik2, ikS, ilS, ipDDij, ipDDik, &
                     ipDDil, ipDDjk, ipDDjl, ipDDkl, ipDij, ipDik, ipDil, ipDjk, ipDjl, ipDkl, ipDum, ipMem1, ipMem2, iPrimi, &
                     iPrInc, ipTmp, iS, iS_, iShelV(4), iShllV(4), iStabs(4), iTmp, jBasAO, jBasj, jBasn, jBsInc, jk2, jkS, jlS, &
                     jPrimj, jPrInc, jS, jS_, kBasAO, kBask, kBasn, kBsInc, klS, kOp(4), kPrimk, kPrInc, kS, kS_, lBasAO, lBasl, &
                     lBasn, lBsInc, lPriml, lPrInc, lS, lS_, mDCRij, mDCRik, mDCRil, mDCRjk, mDCRjl, mDCRkl, mDij, mDik, mDil, &
                     mDjk, mDjl, mDkl, Mem1, Mem2, MemMax, MemPrm, n, nDCRR, nDCRS, nEta, nIJKL, Nr_of_D, nSO, nZeta
real(kind=wp) :: Coor(3,4), Tmax
logical(kind=iwp) :: IJeqKL, NoInts, Shijij
real(kind=wp), pointer :: SOInt(:), AOInt(:)
integer(kind=iwp), external :: iDAMax_, MemSO2
procedure(twoel_kernel) :: TwoEl_NoSym, TwoEl_Sym
procedure(twoel_kernel), pointer :: Do_TwoEl

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGBREIT_
! use the Breit option computing 1/r^3 integralas but convert to
! conventional 1/r integrals
if ((.not. DoFock) .and. (SuperName /= 'gateway') .and. (nIrrep == 1)) call Set_Breit(1)
#endif
mDCRij = 1
mDCRkl = 1
if (nIrrep == 1) then
  Do_TwoEl => TwoEl_NoSym
else
  Do_TwoEl => TwoEl_Sym
end if
!if (.not. associated(Int_PostProcess)) call Abend()
!                                                                      *
!***********************************************************************
!                                                                      *
if (.not. allocated(iSOSym)) then
  call WarningMessage(2,'Eval_Ints_: Integral environment is not set up!')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
NoInts = .true.
Tmax = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! If memory not allocated already at this point allocate!

if (.not. allocated(Sew_Scr)) then
  !write(u6,*) 'Eval_ints: Allocate memory'
  call mma_MaxDBLE(MemMax)
  if (MemMax > 8000) MemMax = MemMax-8000
  call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
else
  !write(u6,*) 'Eval_ints: Memory already allocated'
  MemMax = size(Sew_Scr)
end if
!write(u6,*) 'Eval_ints: MemMax=',MemMax
ipMem1 = 1

Map4(1) = 1
Map4(2) = 2
Map4(3) = 3
Map4(4) = 4
iS_ = max(iiS,jjS)
jS_ = min(iiS,jjS)
kS_ = max(kkS,llS)
lS_ = min(kkS,llS)
if (iiS /= iS_) then
  iTmp = Map4(1)
  Map4(1) = Map4(2)
  Map4(2) = iTmp
end if
if (kkS /= kS_) then
  iTmp = Map4(3)
  Map4(3) = Map4(4)
  Map4(4) = iTmp
end if
!write(u6,*) ' -->',iS_,jS_,kS_,lS_,'<--'
!                                                                      *
!***********************************************************************
!                                                                      *
call Int_Setup(iSD,mSkal,iS_,jS_,kS_,lS_,Coor,Shijij,iAngV,iCmpV,iShelV,iShllV,iAOV,iStabs)
!                                                                      *
!***********************************************************************
!                                                                      *
iPrimi = Shells(iShllV(1))%nExp
jPrimj = Shells(iShllV(2))%nExp
kPrimk = Shells(iShllV(3))%nExp
lPriml = Shells(iShllV(4))%nExp
iBasi = Shells(iShllV(1))%nBasis
jBasj = Shells(iShllV(2))%nBasis
kBask = Shells(iShllV(3))%nBasis
lBasl = Shells(iShllV(4))%nBasis
nZeta = iPrimi*jPrimj
nEta = kPrimk*lPriml
mDij = nZeta+1 ! Dummy initialize
mDkl = nEta+1  ! Dummy initialize
!                                                                      *
!***********************************************************************
!                                                                      *
! partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q

call Create_BraKet(nZeta,nEta)
!                                                                      *
!***********************************************************************
!                                                                      *
! No SO block in direct construction of the Fock matrix.
nSO = MemSO2(iCmpV(1),iCmpV(2),iCmpV(3),iCmpV(4),iShelV(1),iShelV(2),iShelV(3),iShelV(4),iAOV(1),iAOV(2),iAOV(3),iAOV(4))
if (nSO == 0) return

iS = iShelV(1)
jS = iShelV(2)
kS = iShelV(3)
lS = iShelV(4)
ijS = iTri(iS,jS)
klS = iTri(kS,lS)
ikS = iTri(iS,kS)
ilS = iTri(iS,lS)
jkS = iTri(jS,kS)
jlS = iTri(jS,lS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up pointers to k2 entities.

nDCRR = IndK2(2,ijS)
ik2 = IndK2(3,ijS)
nDCRS = IndK2(2,klS)
jk2 = IndK2(3,klS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up pointers to desymmetrized 1st order density
! matrices. Observe that the desymmetrized 1st order
! density matrices follows the contraction index.

if (DoFock) then
  ipTmp = ipDijs
  Nr_of_D = 1
  call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,Nr_of_D)
  call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,Nr_of_D)
  call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,Nr_of_D)
  call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,Nr_of_D)
  call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,Nr_of_D)
  call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,Nr_of_D)

  !write(u6,*) ' Pointers to D=',ipDij,ipDkl,ipDik,ipDil,ipDjk,ipDjl

end if
!                                                                      *
!***********************************************************************
!                                                                      *
!#ifdef _DEBUGPRINT_
!write(u6,*) ' *** Centers ***'
!write(u6,'(3F7.3,6X,3F7.3)') ((Coor(i,j),i=1,3),j=1,2)
!write(u6,'(3F7.3,6X,3F7.3)') ((Coor(i,j),i=1,3),j=3,4)
!#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute memory request for the primitives, i.e.
! how much memory is needed up to the transfer equation.
call MemRys(iAngV,MemPrm)
!                                                                      *
!***********************************************************************
!                                                                      *
! Decide on the partioning of the shells based on the
! available memory and the requested memory.
call PSOAO0(nSO,MemPrm,MemMax,iAngV,iCmpV,iBasi,iBsInc,jBasj,jBsInc,kBask,kBsInc,lBasl,lBsInc,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk, &
            kPrInc,lPriml,lPrInc,ipMem1,ipMem2,Mem1,Mem2,DoFock)
!#ifdef _DEBUGPRINT_
!write(u6,*) ' ************** Memory partioning **************'
!write(u6,*) ' ipMem1=',ipMem1
!write(u6,*) ' ipMem2=',ipMem2
!write(u6,*) ' Mem1=',Mem1
!write(u6,*) ' Mem2=',Mem2
!write(u6,*) ' iBasi,iBsInc=',iBasi,iBsInc
!write(u6,*) ' jBasj,jBsInc=',jBasj,jBsInc
!write(u6,*) ' kBasi,kBsInc=',kBask,kBsInc
!write(u6,*) ' lBasl,lBsInc=',lBasl,lBsInc
!write(u6,*) ' iPrimi,iPrInc=',iPrimi,iPrInc
!write(u6,*) ' jPrimj,jPrInc=',jPrimj,jPrInc
!write(u6,*) ' kPrimk,kPrInc=',kPrimk,kPrInc
!write(u6,*) ' lPriml,lPrInc=',lPriml,lPrInc
!write(u6,*) ' ***********************************************'
!#endif
SOInt(1:Mem1) => Sew_Scr(ipMem1:ipMem1+Mem1-1)
AOInt(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem1-1)
!                                                                      *
!***********************************************************************
!                                                                      *
jbas_ = jBasj
lbas_ = lBasl
!                                                                      *
!***********************************************************************
!                                                                      *
! These loops will partition the contraction loops if there is not
! enough memory to store the whole SO/AO-block simultaneously. The
! memory partitioning is determined by PSOAO0.

do iBasAO=1,iBasi,iBsInc
  iBasn = min(iBsInc,iBasi-iBasAO+1)
  iAOst(1) = iBasAO-1

  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iAOst(2) = jBasAO-1

    ! Move appropiate portions of the desymmetrized 1st
    ! order density matrix.

    if (DoFock) call Picky(iBasi,iBsInc,iPrimi,iBasAO,iBasn,jBasj,jBsInc,jPrimj,jBasAO,jBasn,iCmpV(1),iCmpV(2),iShelV(1), &
                           iShelV(2),mDCRij,ipDij,ipDDij,mDij,DeDe,nDeDe)

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iAOst(3) = kBasAO-1

      if (DoFock) then
        call Picky(iBasi,iBsInc,iPrimi,iBasAO,iBasn,kBask,kBsInc,kPrimk,kBasAO,kBasn,iCmpV(1),iCmpV(3),iShelV(1),iShelV(3),mDCRik, &
                   ipDik,ipDDik,mDik,DeDe,nDeDe)

        call Picky(jBasj,jBsInc,jPrimj,jBasAO,jBasn,kBask,kBsInc,kPrimk,kBasAO,kBasn,iCmpV(2),iCmpV(3),iShelV(2),iShelV(3),mDCRjk, &
                   ipDjk,ipDDjk,mDjk,DeDe,nDeDe)
      end if

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iAOst(4) = lBasAO-1

        if (DoFock) then
          call Picky(kBask,kBsInc,kPrimk,kBasAO,kBasn,lBasl,lBsInc,lPriml,lBasAO,lBasn,iCmpV(3),iCmpV(4),iShelV(3),iShelV(4), &
                      mDCRkl,ipDkl,ipDDkl,mDkl,DeDe,nDeDe)

          call Picky(iBasi,iBsInc,iPrimi,iBasAO,iBasn,lBasl,lBsInc,lPriml,lBasAO,lBasn,iCmpV(1),iCmpV(4),iShelV(1),iShelV(4), &
                      mDCRil,ipDil,ipDDil,mDil,DeDe,nDeDe)

          call Picky(jBasj,jBsInc,jPrimj,jBasAO,jBasn,lBasl,lBsInc,lPriml,lBasAO,lBasn,iCmpV(2),iCmpV(4),iShelV(2),iShelV(4), &
                      mDCRjl,ipDjl,ipDDjl,mDjl,DeDe,nDeDe)
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        !         Compute SO/AO-integrals

        call Do_TwoEl(iS_,jS_,kS_,lS_,Coor,iAngV,iCmpV,iShelV,iShllV,iAOV,iAOst,NoInts,iStabs,iPrimi,iPrInc,jPrimj,jPrInc,kPrimk, &
                      kPrInc,lPriml,lPrInc,nDCRR,nDCRS,k2Data(:,ik2),k2Data(:,jk2),IJeqKL,kOp,DeDe(ipDDij),mDij,mDCRij, &
                      DeDe(ipDDkl),mDkl,mDCRkl,DeDe(ipDDik),mDik,mDCRik,DeDe(ipDDil),mDil,mDCRil,DeDe(ipDDjk),mDjk,mDCRjk, &
                      DeDe(ipDDjl),mDjl,mDCRjl,Shells(iShllV(1))%pCff(1,iBasAO),iBasn,Shells(iShllV(2))%pCff(1,jBasAO),jBasn, &
                      Shells(iShllV(3))%pCff(1,kBasAO),kBasn,Shells(iShllV(4))%pCff(1,lBasAO),lBasn,FT,nFT,nZeta,nEta,SOInt,nSO, &
                      AOInt,Mem2,Shijij,Aux,nAux)
        !                                                              *
        !***************************************************************
        !                                                              *

        nijkl = iBasn*jBasn*kBasn*lBasn
#       ifdef _DEBUGBREIT_
        if (nOrdOp /= 0) then
          if (nIrrep == 1) then
            n = iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
            call ReSort_Int(AOInt,nijkl,6,n)
          else
            call ReSort_Int(SOInt,nijkl,6,nSO)
          end if
        end if
#       endif
#       ifdef _DEBUGPRINT_
        if (nIrrep == 1) then
          n = iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
          call RecPrt('AOInt',' ',AOInt,nijkl,n)
        else
          call RecPrt('SOInt',' ',SOInt,nijkl,nSO)
        end if
#       endif
        !                                                              *
        !***************************************************************
        !                                                              *
        ! Process SO/AO-integrals

        if (DoIntegrals .and. (.not. NoInts)) then
          ! Get max AO/SO integrals
          if (nIrrep == 1) then
            n = nijkl*iCmpV(1)*iCmpV(2)*iCmpV(3)*iCmpV(4)
            Tmax = max(Tmax,abs(AOInt(iDAMax_(n,AOInt,1))))
          else
            n = nijkl*nSO
            Tmax = max(Tmax,abs(SOInt(iDAMax_(n,SOInt,1))))
          end if
          if (Tmax > CutInt) then
            call Int_PostProcess(iCmpV,iShelV,iBasn,jBasn,kBasn,lBasn,kOp,Shijij,iAOV,iAOst,nijkl,AOInt,SOInt,nSO,iSOSym,nSOs, &
                                 TInt,nTInt,nIrrep)
          else
            Tmax = Zero
          end if
        end if

      end do
    end do
  end do
end do
nullify(SOInt,AOInt)
call Destroy_BraKet()
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGBREIT_
call Set_Breit(0)

contains

subroutine ReSort_Int(IntRaw,nijkl,nComp,nA)

# ifdef _DEBUGPRINT_
  use Definitions, only: u6
# endif

  integer(kind=iwp), intent(in) :: nijkl, nComp, nA
  real(kind=wp), target, intent(inout) :: IntRaw(nijkl*nComp*nA)
  real(kind=wp), pointer :: IntIn(:,:,:), IntOut(:,:)
  integer(kind=iwp) :: i_ijkl, iA

  IntIn(1:nijkl,1:nComp,1:nA) => IntRaw(:)
  IntOut(1:nijkl,1:nA) => IntRaw(1:nijkl*nA)
# ifdef _DEBUGPRINT_
  write(u6,*) 'nijkl,nComp,nA=',nijkl,nComp,nA
  call RecPrt('IntRaw',' ',IntRaw,nijkl,nComp*nA)
# endif

  IntOut(:,:) = IntIn(:,:)+IntIn(:,4,:)+IntIn(:,6,:)

  nullify(IntIn,IntOut)

end subroutine ReSort_Int
#endif

end subroutine Eval_ijkl

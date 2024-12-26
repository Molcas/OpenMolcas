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
use setup, only: mSkal, nSOs
use k2_arrays, only: Create_BraKet, Destroy_Braket, iSOSym, Sew_Scr
use iSD_data, only: iSD, nSD
use Breit, only: nComp
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: DoFock, DoIntegrals, Map4
use Integral_interfaces, only: Int_PostProcess, twoel_kernel
use RI_glob, only: jBas_, lBas_
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
integer(kind=iwp) :: iBasAO, iBasi, iBasn, iBsInc, ijS, ikS, ilS, ipDum, ipMem1, ipMem2, &
                     iS, iS_, iTmp, jBasAO, jBasj, jBasn, jBsInc, jkS, jlS, &
                     jS, jS_, kBasAO, kBask, kBasn, kBsInc, klS, kOp(4), kS, kS_, lBasAO, lBasl, &
                     lBasn, lBsInc, lS, lS_, Mem1, Mem2, MemMax, MemPrm, n, nEta, nIJKL, nSO, nZeta
integer(kind=iwp) :: iSD4(0:nSD,4)
real(kind=wp) :: Coor(3,4), Tmax
logical(kind=iwp) :: NoInts, Shijij, No_batch
real(kind=wp), pointer :: SOInt(:), AOInt(:)
integer(kind=iwp), external :: iDAMax_
procedure(twoel_kernel) :: TwoEl_NoSym, TwoEl_Sym
procedure(twoel_kernel), pointer :: Do_TwoEl
TInt(:)=Zero
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGBREIT_
! use the Breit option computing 1/r^3 integralas but convert to
! conventional 1/r integrals
if ((.not. DoFock) .and. (SuperName /= 'gateway') .and. (nIrrep == 1)) call Set_Breit(1)
#endif
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
call Gen_iSD4(iS_,jS_,kS_,lS_,iSD,nSD,iSD4)
call Int_Setup(iSD,mSkal,iS_,jS_,kS_,lS_,Coor,Shijij)
!                                                                      *
!***********************************************************************
!                                                                      *

nZeta = iSD4(5,1)*iSD4(5,2)
nEta  = iSD4(5,3)*iSD4(5,4)
!                                                                      *
!***********************************************************************
!                                                                      *
! partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q

call Create_BraKet(nZeta,nEta)
!                                                                      *
!***********************************************************************
!                                                                      *
! No SO block in direct construction of the Fock matrix.
Call Size_SOb(iSD4,nSD,nSO,No_batch)
if (No_Batch) return

iS = iSD4(11,1)
jS = iSD4(11,2)
kS = iSD4(11,3)
lS = iSD4(11,4)
ijS = iTri(iS,jS)
klS = iTri(kS,lS)
ikS = iTri(iS,kS)
ilS = iTri(iS,lS)
jkS = iTri(jS,kS)
jlS = iTri(jS,lS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up pointers to desymmetrized 1st order density
! matrices. Observe that the desymmetrized 1st order
! density matrices follows the contraction index.

if (DoFock) Call Dens_Infos()
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
call MemRys(iSD4(1,:),MemPrm)
!                                                                      *
!***********************************************************************
!                                                                      *
! Decide on the partioning of the shells based on the
! available memory and the requested memory.
call PSOAO0(nSO,MemPrm,MemMax,ipMem1,ipMem2,Mem1,Mem2,DoFock,nSD,iSD4)

iBsInc = iSD4(4,1)
jBsInc = iSD4(4,2)
kBsInc = iSD4(4,3)
lBsInc = iSD4(4,4)

SOInt(1:Mem1) => Sew_Scr(ipMem1:ipMem1+Mem1-1)
AOInt(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem2-1)
!                                                                      *
!***********************************************************************
!                                                                      *
iBasi = iSD4(3,1)
jBasj = iSD4(3,2)
kBask = iSD4(3,3)
lBasl = iSD4(3,4)

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
  iSD4( 8,1) = iBasAO-1
  iSD4(19,1) = iBasn

  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iSD4( 8,2) = jBasAO-1
    iSD4(19,2) = jBasn

    ! Move appropiate portions of the desymmetrized 1st
    ! order density matrix.

    if (DoFock) call Picky(nSD,iSD4,1,2)

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iSD4( 8,3) = kBasAO-1
      iSD4(19,3) = kBasn

      if (DoFock) then
        call Picky(nSD,iSD4,1,3)
        call Picky(nSD,iSD4,2,3)
      end if

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iSD4( 8,4) = lBasAO-1
        iSD4(19,4) = lBasn

        if (DoFock) then
          call Picky(nSD,iSD4,3,4)
          call Picky(nSD,iSD4,1,4)
          call Picky(nSD,iSD4,2,4)
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
        !         Compute SO/AO-integrals

        nijkl = iBasn*jBasn*kBasn*lBasn*nComp
        call Do_TwoEl(iS_,jS_,kS_,lS_,Coor,NoInts,kOp, &
                      nZeta,nEta,SOInt,nijkl,nSO, &
                      AOInt,Mem2,Shijij,iSD4)
        !                                                              *
        !***************************************************************
        !                                                              *

#       ifdef _DEBUGBREIT_
        if (nOrdOp /= 0) then
          if (nIrrep == 1) then
            n = iSD4(2,1)*iSD4(2,2)*iSD4(2,3)*iSD4(2,4)
            call ReSort_Int(AOInt,nijkl,6,n)
          else
            call ReSort_Int(SOInt,nijkl,6,nSO)
          end if
        end if
#       endif
#       ifdef _DEBUGPRINT_
        if (nIrrep == 1) then
          n = iSD4(2,1)*iSD4(2,2)*iSD4(2,3)*iSD4(2,4)
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
            n = nijkl*iSD4(2,1)*iSD4(2,2)*iSD4(2,3)*iSD4(2,4)
            Tmax = max(Tmax,abs(AOInt(iDAMax_(n,AOInt,1))))
          else
            n = nijkl*nSO
            Tmax = max(Tmax,abs(SOInt(iDAMax_(n,SOInt,1))))
          end if
          if (Tmax > CutInt) then
            call Int_PostProcess(kOp,Shijij,nijkl,AOInt,SOInt,nSO,iSOSym,nSOs, &
                                 TInt,nTInt,nIrrep,nSD,iSD4)
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
# endif

contains

Subroutine Dens_Infos()
use Dens_stuff, only: mDCRij,mDCRkl,mDCRik,mDCRil,mDCRjk,mDCRjl,&
                       ipDij, ipDkl, ipDik, ipDil, ipDjk, ipDjl,&
                      ipDDij,ipDDkl,ipDDik,ipDDil,ipDDjk,ipDDjl
use k2_arrays, only: ipDijS
Implicit None
integer(kind=iwp), parameter:: Nr_of_D = 1
integer(kind=iwp) ipTmp
ipTmp = ipDijs
call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,Nr_of_D)
call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,Nr_of_D)
call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,Nr_of_D)
call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,Nr_of_D)
call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,Nr_of_D)
call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,Nr_of_D)
End Subroutine Dens_Infos

#ifdef _DEBUGBREIT_
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

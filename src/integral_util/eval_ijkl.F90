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
subroutine Eval_ijkl(iS,jS,kS,lS,TInt,nTInt)
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
use setup, only: nSOs
use k2_arrays, only: Create_BraKet, Destroy_Braket, iSOSym, Sew_Scr
use iSD_data, only: iSD, nSD
use Breit, only: nComp
use Gateway_Info, only: CutInt
use Symmetry_Info, only: nIrrep
use Int_Options, only: DoFock, DoIntegrals
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
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, nTInt
real(kind=wp), intent(inout) :: TInt(nTInt)
integer(kind=iwp) :: iAng(4), iBasAO, iBasi, iBasn, iBsInc, ipDum, ipMem1, ipMem2, iSD4(0:nSD,4), jBasAO, jBasj, jBasn, jBsInc, &
                     kBasAO, kBask, kBasn, kBsInc, lBasAO, lBasl, lBasn, lBsInc, Mem1, Mem2, MemMax, MemPrm, n, nAO, nIJKL, nSO
real(kind=wp) :: Coor(3,4), Tmax
logical(kind=iwp) :: NoInts
real(kind=wp), pointer :: SOInt(:), AOInt(:)
procedure(twoel_kernel) :: TwoEl_NoSym, TwoEl_Sym
procedure(twoel_kernel), pointer :: Do_TwoEl
integer(kind=iwp), parameter :: SCF = 1
integer(kind=iwp), external :: iDAMax_, MemSO2
!                                                                      *
!***********************************************************************
!                                                                      *
TInt(:) = Zero
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

!write(u6,*) ' -->',iS,jS,kS,lS,'<--'
!                                                                      *
!***********************************************************************
!                                                                      *
call Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)
!                                                                      *
!***********************************************************************
!                                                                      *
! No SO block in direct construction of the Fock matrix.
if (nIrrep > 1) then
  nSO = MemSO2(nSD,iSD4)
else
  nSO = 0
end if
if ((nIrrep > 1) .and. (nSO == 0)) return
nAO = product(iSD4(2,:))

call Coor_Setup(iSD4,nSD,Coor)
call Int_Setup(Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! partition memory for K2(ij)/K2(kl) temp spaces zeta,eta,kappa,P,Q

call Create_BraKet(iSD4(5,1)*iSD4(5,2),iSD4(5,3)*iSD4(5,4))
!                                                                      *
!***********************************************************************
!                                                                      *
! Pick up pointers to desymmetrized 1st order density
! matrices. Observe that the desymmetrized 1st order
! density matrices follows the contraction index.

if (DoFock) call Dens_Infos(SCF)
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
iAng(:) = iSD4(1,:)
call MemRys(iAng,MemPrm)
!                                                                      *
!***********************************************************************
!                                                                      *
! Decide on the partitioning of the shells based on the
! available memory and the requested memory.
call PSOAO0(nSO,MemPrm,MemMax,ipMem1,ipMem2,Mem1,Mem2,DoFock,nSD,iSD4)
SOInt(1:Mem1) => Sew_Scr(ipMem1:ipMem1+Mem1-1)
AOInt(1:Mem2) => Sew_Scr(ipMem2:ipMem2+Mem2-1)
!                                                                      *
!***********************************************************************
!                                                                      *

iBsInc = iSD4(4,1)
jBsInc = iSD4(4,2)
kBsInc = iSD4(4,3)
lBsInc = iSD4(4,4)

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
  iSD4(8,1) = iBasAO-1
  iSD4(19,1) = iBasn

  do jBasAO=1,jBasj,jBsInc
    jBasn = min(jBsInc,jBasj-jBasAO+1)
    iSD4(8,2) = jBasAO-1
    iSD4(19,2) = jBasn

    ! Move appropiate portions of the desymmetrized 1st
    ! order density matrix.

    if (DoFock) call Picky(nSD,iSD4,1,2)

    do kBasAO=1,kBask,kBsInc
      kBasn = min(kBsInc,kBask-kBasAO+1)
      iSD4(8,3) = kBasAO-1
      iSD4(19,3) = kBasn

      if (DoFock) then
        call Picky(nSD,iSD4,1,3)
        call Picky(nSD,iSD4,2,3)
      end if

      do lBasAO=1,lBasl,lBsInc
        lBasn = min(lBsInc,lBasl-lBasAO+1)
        iSD4(8,4) = lBasAO-1
        iSD4(19,4) = lBasn

        if (DoFock) then
          call Picky(nSD,iSD4,3,4)
          call Picky(nSD,iSD4,1,4)
          call Picky(nSD,iSD4,2,4)
        end if

        nijkl = iBasn*jBasn*kBasn*lBasn*nComp ! *nComp is a fix for BP integrals

        !                                                              *
        !***************************************************************
        !                                                              *
        ! Compute SO/AO-integrals

        call Do_TwoEl(Coor,NoInts,SOInt,nijkl,nSO,AOInt,Mem2,iSD4)

        nijkl = iBasn*jBasn*kBasn*lBasn

#       ifdef _DEBUGBREIT_
        if (nOrdOp /= 0) then
          if (nIrrep == 1) then
            call ReSort_Int(AOInt,nijkl,6,nAO)
          else
            call ReSort_Int(SOInt,nijkl,6,nSO)
          end if
        end if
#       endif
#       ifdef _DEBUGPRINT_
        if (nIrrep == 1) then
          call RecPrt('AOInt',' ',AOInt,nijkl,nAO)
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
            n = nijkl*nAO
            Tmax = max(Tmax,abs(AOInt(iDAMax_(n,AOInt,1))))
          else
            n = nijkl*nSO
            Tmax = max(Tmax,abs(SOInt(iDAMax_(n,SOInt,1))))
          end if
          if (Tmax > CutInt) then
            call Int_PostProcess(nijkl,AOInt,SOInt,nSO,iSOSym,nSOs,TInt,nTInt,nIrrep,nSD,iSD4)
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

subroutine Dens_Infos(nMethod)

  use Dens_stuff, only: ipDDij, ipDDik, ipDDil, ipDDjk, ipDDjl, ipDDkl, ipDij, ipDik, ipDil, ipDjk, ipDjl, ipDkl, mDCRij, mDCRik, &
                        mDCRil, mDCRjk, mDCRjl, mDCRkl
  use k2_arrays, only: ipDijS

  integer(kind=iwp), intent(in) :: nMethod
  integer(kind=iwp) :: Dum1, Dum2, Dum3, ijS, ikS, ilS, ipTmp, iS, jkS, jlS, jS, klS, kS, lS
  integer(kind=iwp), parameter :: Nr_of_D = 1

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
  ipTmp = ipDijs
  call Dens_Info(ijS,ipDij,ipDum,mDCRij,ipDDij,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)
  call Dens_Info(klS,ipDkl,ipDum,mDCRkl,ipDDkl,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)
  call Dens_Info(ikS,ipDik,ipDum,mDCRik,ipDDik,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)
  call Dens_Info(ilS,ipDil,ipDum,mDCRil,ipDDil,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)
  call Dens_Info(jkS,ipDjk,ipDum,mDCRjk,ipDDjk,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)
  call Dens_Info(jlS,ipDjl,ipDum,mDCRjl,ipDDjl,ipTmp,Nr_of_D,nMethod,Dum1,Dum2,Dum3)

end subroutine Dens_Infos

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

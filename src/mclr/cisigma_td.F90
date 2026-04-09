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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

#include "compiler_features.h"
#ifdef _WARNING_WORKAROUND_
#  ifdef _LIFETIME_BUG_
#    define _SAVE_TARGET_ , save
#  else
#    define _SAVE_TARGET_
#  endif
#else
#  define _SAVE_TARGET_
#endif

subroutine CISigma_td(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,NT,Have_2_el)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use ipPage, only: ipin, ipin1, ipnout, opout, W
use MCLR_Data, only: i12, ipCM, ipMat, iRefSM, ist, KAIN1, KINT2, KINT2A, nConf1, nDens, pInt1, Square, XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: nBas, nCSF, nSym, ntAsh, Page, State_Sym, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real(kind=wp), target, intent(in) :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
character, intent(in) :: NT
logical(kind=iwp), intent(in) :: Have_2_el
integer(kind=iwp) :: i, ij, ijkl, iOp, iS, j, ji, jilk, jS, k, kic(2), kl, l, lk, nDet
real(kind=wp), allocatable :: CIDET(:)
real(kind=wp), allocatable, target _SAVE_TARGET_ :: TI1(:), TI2(:)

! For the timeindep case ipS1 and ipS2 will be half as long
! Avoid sigmavec calls. 95% of the time in mclr is spent in sigmavec

! Interface Anders to Jeppe
! This interface initiates Jeppe's common block
! and will make it easier to use Anders modifications
! of the CI routines

! OK first tell Jeppe where the integrals are.

if (nconf1 == 0) return

! One electron integrals

KAIN1(1:nInt1) => Int1(:)

! Two electron integrals
! symmetric in particle one and two

KINT2(1:nInt2s) => Int2s(:)
KINT2a(1:nInt2a) => Int2a(:)

! Two electron integrals
! anti symmetric in particle one and two

irefsm = iCSym

! Do we have any two-electron integrals?

if (Have_2_el) then
  i12 = 2
else
  i12 = 1
end if

! Symmetry of Sigma vector

iSSM = iSSym
kic(2) = 2
if (issm == State_sym) kic(2) = 1

! Symmetry of CI vector

iCSM = iCSym
kic(1) = 2
if (icsm == State_sym) kic(1) = 1

! Symmetry properties of operator

ndet = nint(max(xispsm(iSSym,1),xispsm(iCSym,1)))
ndet = max(ndet,ncsf(icsym),ncsf(issym))
if (ndet == 0) return
iOP = Mul(iCSM,iSSm)
if (iOp == 1) then
  pInt1(1:nSym) = ipCM(1:nSym)
else
  do iS=1,nSym
    pInt1(is) = ipMat(is,Mul(iS,iOp))
  end do
end if

! Triplet/Singlet operator

ist = iispin+1
square = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
if (TIMEDEP) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (NT == 'T') square = .true.  ! The operator is not sym

  if (page) then
    write(u6,*) 'Page not implemented for Time-dependent perturbations'
    call Abend()
  end if

  ! CIDET is here because sigmavec will destroy the first input vector.
  call mma_allocate(CIDET,nDet,Label='CIDET')
  call ipin(ipCI1)
  CIDET(1:nCSF(iCSM)) = W(ipCI1)%A(1:nCSF(iCSM))

  call ipin(ipci2)
  call SigmaVec(CIDET,W(ipci2)%A,kic)

  if (NT == 'N') then
    call mma_deallocate(CIDET)
    return

  else if (NT == 'S') then

    ! Symmetric operator, no transpose of integrals needed!
    call ipin(ipCI1)
    CIDET(1:nCSF(iCSM)) = W(ipCI1)%A(nConf1+1:nConf1+nCSF(iCSM))

    call ipin(ipci2)
    call SigmaVec(CIDET,W(ipci2)%A(1+nconf1),kic)

  else  ! NT /= 'S'

    ! The operator is not sym --> transpose integrals! NT /= S
    call ipin(ipCI1)
    CIDET(1:nCSF(iCSM)) = W(ipCI1)%A(1:nCSF(iCSM))

    call mma_allocate(TI1,nDens,Label='TI1')
    call mma_allocate(TI2,ntash**4,Label='TI2')

    do i=1,ntash
      do j=1,ntash
        ij = i+ntash*(j-1)
        ji = j+ntash*(i-1)
        do k=1,ntash
          do l=1,ntash
            kl = k+ntash*(l-1)
            lk = l+ntash*(k-1)
            if (ij >= kl) then
              ijkl = iTri(ij,kl)
              jilk = iTri(ji,lk)
              TI2(jilk) = int2s(ijkl)
            end if
          end do
        end do
      end do
    end do

    do is=1,nSym
      js = Mul(Mul(icsym,issym),is)
      if (nbas(js)*nbas(is) /= 0) call DGETMO(Int1(ipmat(is,js)),nbas(is),nbas(is),nbas(js),TI1(ipmat(js,is)),nbas(js))
    end do

    KAIN1 => TI1
    KINT2 => TI2

    call ipin(ipci2)
    call SigmaVec(CIDET,W(ipci2)%A(1+nconf1),kic)

    nullify(KAIN1,KINT2)
    call mma_deallocate(TI1)
    call mma_deallocate(TI2)

  end if  ! End the transpose of integrals.

  call mma_deallocate(CIDET)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else   ! If not timedep
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  if (.not. page) then
    call mma_allocate(CIDET,nDet,Label='CIDET')
    call ipin(ipCI1)
    CIDET(1:nCSF(iCSM)) = W(ipCI1)%A(1:nCSF(iCSM))
    call ipin(ipci2)
    call SigmaVec(CIDET,W(ipci2)%A,kic)
    call mma_deallocate(CIDET)
  else
    call ipnout(ipci2)
    call ipin1(ipCI1,ndet)
    call ipin(ipci2)
    call SigmaVec(W(ipCI1)%A,W(ipci2)%A,kic)
    call opout(ipci1)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine CISigma_td

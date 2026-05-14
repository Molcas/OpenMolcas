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

subroutine CISigma(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)

use Symmetry_Info, only: Mul
use ipPage, only: ipin, ipin1, ipnout, opout, W
use MCLR_Data, only: i12, ipCM, ipMat, iRefSM, iST, KAIN1, KINT2, KINT2A, nConf1, pInt1, Square, XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: nCSF, nSym, Page, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real(kind=wp), target, intent(in) :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
logical(kind=iwp), intent(in) :: Have_2_el
integer(kind=iwp) :: nDet, iOp, iS, kic(2)
real(kind=wp), allocatable :: CIDET(:)

! Interface Anders to Jeppe
! This interface initiates Jeppe's common block
! and will make it easier to use Anders modifications
! of the CI routines

! OK first tell Jeppe where the integrals are.

! yma: notice the nconf1 if DMRG

if (nconf1 == 0) return

! One electron integrals

KAIN1(1:nInt1) => Int1(:)

! Two electron integrals
! symmetric in particle one and two

KINT2(1:nInt2s) => Int2s(:)

! Two electron integrals
! anti symmetric in particle one and two

KINT2A(1:nInt2a) => Int2a(:)

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

if (.not. page) then

  call mma_allocate(CIDET,nDet,Label='CIDET')
# ifdef _MS_
  call ipin(ipCI1)
  call ipin(ipci2)
  do i=1,nroots
    CIDET(1:nCSF(iCSM)) = W(ipCI1)%A((i-1)*ncsf(icsm)+1:i*ncsf(icsm))
    call SigmaVec(CIDET,W(ipci2)%A((i-1)*ncsf(issm)+1),kic)
  end do
# else
  CIDET(1:nCSF(iCSM)) = W(ipCI1)%A(1:nCSF(iCSM))

  call SigmaVec(CIDET,W(ipci2)%A,kic)

# endif
  call mma_deallocate(CIDET)
else
  call ipnout(ipci2)

  call ipin1(ipCI1,ndet)
  call ipin(ipci2)
  call SigmaVec(W(ipCI1)%A,W(ipci2)%A,kic)
  call opout(ipci1)

end if

end subroutine CISigma

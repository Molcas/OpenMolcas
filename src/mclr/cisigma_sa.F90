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
!#error "This file must be compiled inside a module"
!#endif
#else

subroutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)

use Symmetry_Info, only: Mul
use ipPage, only: ipin, ipin1, ipnout, opout, W
use MCLR_Data, only: KAIN1, KINT2, KINT2A, pInt1
use MCLR_Data, only: nConf1, ipCM, ipMat
use MCLR_Data, only: i12, iST, Square
use MCLR_Data, only: iRefSM
use MCLR_Data, only: XISPSM
use CandS, only: ICSM, ISSM
use input_mclr, only: State_Sym, nSym, Page, nCSF, nRoots, Weight
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real*8, target :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
logical Have_2_el
integer kic(2)
real*8, allocatable :: CIDET(:)
integer nDet, iOp, iS, jS, i

! Interface Anders to Jeppe
! This interface initiates Jeppes common block
! and will make it easier to use Anders modifications
! of the CI routines

! OK first tell Jeppe where the integrals are.

! yma: notice the nconf1 if DMRG

if (nconf1 == 0) return

! One electron integrals

KAIN1 => Int1

! Two electron integrals
! symmetric in perticle one and two

KINT2 => Int2s

! Two electron integrals
! anti symmetric in perticle one and two

KINT2A => Int2a

irefsm = iCSym

! Do we have any twoelectron integrals?

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
    jS = Mul(iS,iOp)
    pInt1(is) = ipMat(is,jS)
  end do
end if

! Triplet/Singlet operator

ist = iispin+1
square = .false.

if (.not. page) then
  call mma_allocate(CIDET,nDet,Label='CIDET')
  call ipin(ipCI1)
  call ipin(ipci2)
  do i=1,nroots
    CIDET(1:nCSF(iCSM)) = W(ipCI1)%A((i-1)*ncsf(icsm)+1:i*ncsf(icsm))
    call SigmaVec(CIDET,W(ipci2)%A(1+(i-1)*ncsf(issm)),kic)
    W(ipci2)%A((i-1)*ncsf(issm)+1:i*ncsf(issm)) = weight(i)*W(ipci2)%A((i-1)*ncsf(issm)+1:i*ncsf(issm))
  end do
  call mma_deallocate(CIDET)
else
  call ipnout(ipci2)
  call ipin1(ipCI1,ndet)
  call ipin(ipci2)
  call SigmaVec(W(ipCI1)%A,W(ipci2)%A,kic)
  call opout(ipci1)
end if

end subroutine CISigma_sa
#endif

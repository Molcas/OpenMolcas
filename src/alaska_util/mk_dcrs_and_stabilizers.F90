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
! Copyright (C) 1990,1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine mk_DCRs_and_Stabilizers(Fact,iuvwx,nDCRR,nDCRS,nDCRT,iDCRR,iDCRS,iDCRT,nSD,iSD4)

use Symmetry_Info, only: nIrrep
use Basis_Info, only: MolWgh
use Center_Info, only: dc
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Symmetry_Info, only: ChOper
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(out) :: Fact
integer(kind=iwp), intent(out) :: iuvwx(4), nDCRR, nDCRS, nDCRT, iDCRR(0:7), iDCRS(0:7), iDCRT(0:7)
integer(kind=iwp), intent(in) :: nSD, iSD4(0:nSD,4)
integer(kind=iwp) :: iStabM(0:7), iStabN(0:7), iStb, jStb, kStb, LmbdR, LmbdS, LmbdT, lStabM, lStabN, lStb
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif
real(kind=wp) :: u, v, w, x

iStb = iSD4(10,1)
jStb = iSD4(10,2)
kStb = iSD4(10,3)
lStb = iSD4(10,4)
iuvwx(1) = dc(iStb)%nStab
iuvwx(2) = dc(jStb)%nStab
iuvwx(3) = dc(kStb)%nStab
iuvwx(4) = dc(lStb)%nStab
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center A and B

if (nIrrep == 1) then
  nDCRR = 1
  iDCRR(0) = 0
  LmbdR = 1
else
  call DCR(LmbdR,dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iDCRR,nDCRR)
end if
#ifdef _DEBUGPRINT_
write(u6,'(20A)') ' {R}=(',(ChOper(iDCRR(i)),',',i=0,nDCRR-1),')'
#endif
u = real(dc(iStb)%nStab,kind=wp)
v = real(dc(jStb)%nStab,kind=wp)

! Find stabilizer for center A and B

if (nIrrep == 1) then
  lStabM = 1
  iStabM(0) = 0
else
  call Inter(dc(iStb)%iStab,dc(iStb)%nStab,dc(jStb)%iStab,dc(jStb)%nStab,iStabM,lStabM)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for center C and D.
! Take care of redundancy if {f(aA)f(bB)}={f(cC)f(dD)}. Hence
! we will only use unique combinations of operators from the
! double coset representatives {R} and {S}.

if (nIrrep == 1) then
  nDCRS = 1
  iDCRS(0) = 0
  LmbdS = 1
else
  call DCR(LmbdS,dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iDCRS,nDCRS)
end if
#ifdef _DEBUGPRINT_
write(u6,'(20A)') ' {S}=(',(ChOper(iDCRS(i)),',',i=0,nDCRS-1),')'
#endif
w = real(dc(kStb)%nStab,kind=wp)
x = real(dc(lStb)%nStab,kind=wp)

! Find stabilizer for center C and D

if (nIrrep == 1) then
  lStabN = 1
  iStabN(0) = 0
else
  call Inter(dc(kStb)%iStab,dc(kStb)%nStab,dc(lStb)%iStab,dc(lStb)%nStab,iStabN,lStabN)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Double Coset Representatives for the two charge distributions.

if (nIrrep == 1) then
  nDCRT = 1
  iDCRT(0) = 0
  LmbdT = 1
else
  call DCR(LmbdT,iStabM,lStabM,iStabN,lStabN,iDCRT,nDCRT)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Factor due to summation over DCR

if (MolWgh == 1) then
  Fact = real(nIrrep,kind=wp)/real(LmbdT,kind=wp)
else if (MolWgh == 0) then
  Fact = u*v*w*x/real(nIrrep**3*LmbdT,kind=wp)
else
  Fact = sqrt(u*v*w*x)/real(nIrrep*LmbdT,kind=wp)
end if

end subroutine mk_DCRs_and_Stabilizers

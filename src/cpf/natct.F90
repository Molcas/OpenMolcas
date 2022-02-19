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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine NATCT(H,LIC0)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: LIC0
real(kind=wp) :: H(LIC0)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: ICMO, IDISK, iDum, iDummy(7,8), iiRC, IOCC, iOpt, iSYM, LW91A, LW91B, M, N2SUM, n2Tri, nbMax, NSUM
real(kind=wp) :: dum, Dummy(1), EREL, ErelDC, ErelMV
character(len=72) :: Header

NSUM = 0
N2SUM = 0
n2Tri = 0
nbMax = 0
do ISYM=1,NSYM
  nbMax = max(nbMax,nBas(iSym))
  NSUM = NSUM+NBAS(ISYM)
  N2SUM = N2SUM+NBAS(ISYM)**2
  n2Tri = n2Tri+nBas(iSym)*(nBas(iSym)+1)/2
end do

! Read MO coefficients
IDISK = ITOC17(1)
call dDAFILE(Lu_TraOne,2,H(LW(87)),N2SUM,IDISK)
if (LW(87)+N2SUM-1 >= LW(88)) then
  write(u6,*)
  write(u6,'(6X,A)') '*** ERROR IN SUBROUTINE NATCT ***'
  write(u6,'(6X,A)') 'NO SPACE LEFT TO GENERATE FINAL ORBITALS'
  write(u6,*)
  call XFLUSH(u6)
  call Abend()
end if

! Loop over irreps and compute natural orbitals

IOCC = LW(90)
ICMO = LW(87)
do M=1,NSYM
  ! set occupation number of orbitals prefrozen in MOTRA
  call DCOPY_(NBAS(M),[Zero],0,H(IOCC),1)
  ! skip orbitals prefrozen in MOTRA
  call DCOPY_(NPFRO(M),[Two],0,H(IOCC),1)
  call NATORB_CPF(H(LW(62)),H(ICMO+NBAS(M)*NPFRO(M)),H(LW(88)),H(LW(89)),H(LW(89)),H(IOCC+NPFRO(M)),M)
  call DCOPY_(NORB(M)*NBAS(M),H(LW(89)),1,H(ICMO+NBAS(M)*NPFRO(M)),1)
  ICMO = ICMO+NBAS(M)**2
  IOCC = IOCC+NBAS(M)
end do

LW91A = LW(91)
LW91B = LW91A+n2Sum
if (LW91B+n2Tri-1 > Lic) then
  write(u6,*) ' Not enough core in NATCT'
  call XFLUSH(u6)
  call ErrTra()
  call Abend()
end if
call RelEne(ErelMV,ErelDC,nSym,nBas,H(LW(87)),H(LW(90)),H(LW91A),H(LW91B))

EREL = ERELMV+ERELDC
write(u6,'(/,5X,A)') 'FIRST ORDER RELATIVISTIC CORRECTIONS'
write(u6,'(5X,A,F17.8)') 'MASS-VELOCITY        ',ErelMV
write(u6,'(5X,A,F17.8)') '1-EL DARWIN CONTACT  ',ErelDC
write(u6,'(5X,A,F17.8)') 'TOTAL REL. CORRECTION',Erel
if (ISDCI == 1) then
  write(u6,'(5X,A,F17.8)') 'REL. CI ENERGY       ',ETOT+Erel
  write(u6,'(5X,A,F17.8)') 'REL. CI+Q ENERGY     ',DETOT+Erel
else
  write(u6,'(5X,A,F17.8)') 'TOTAL REL. ENERGY    ',ETOT+Erel
end if
call XFLUSH(u6)
call dPRWF(H)
if (iCPF == 1) then
  Header = ' CPF natural orbitals'
else if (iSDCI == 1) then
  Header = ' SDCI natural orbitals'
else if (iNCPF == 1) then
  Header = ' ACPF natural orbitals'
else
  Header = ' MCPF natural orbitals'
end if
call Primo(Header,.true.,.false.,1.0e-4_wp,dum,nSym,nBas,nBas,Name,Dummy,H(LW(90)),H(LW(87)),-1)

! Read the overlap matrix in ao basis
iiRC = -1
iOpt = 6
call RdOne(iiRC,iOpt,'MLTPL  0',1,H(LW(91)),iDum)
if (iiRC /= 0) then
  write(u6,*) 'Natct: Error reading overlap matrix!'
  call Abend()
end if
call Charge(nSym,nBas,Name,H(LW(87)),H(LW(90)),H(LW(91)),2,.true.,.true.)
call Prpt_old(nSym,nBas,nSum,n2Sum,H(LW(87)),H(LW(90)))

if (iCPF == 1) then
  Header = '* CPF NO COEFS'
else if (iSDCI == 1) then
  Header = '* SDCI NO COEFS'
else if (iNCPF == 1) then
  Header = '* ACPF NO COEFS'
else
  Header = '* MCPF NO COEFS'
end if
call WrVec('CPFORB',Lu_CPFORB,'CO',nSym,nBas,nBas,H(LW(87)),H(LW(90)),Dummy,iDummy,Header)

return

! This is to allow type punning without an explicit interface
contains

subroutine dPRWF(H)

  real(kind=wp), target :: H(*)
  integer(kind=iwp), pointer :: iH1(:), iH2(:), iH3(:)

  call c_f_pointer(c_loc(H(LW(1))),iH1,[1])
  call c_f_pointer(c_loc(H(LW(2))),iH2,[1])
  call c_f_pointer(c_loc(H(LW(3))),iH3,[1])
  call PRWF_CPF(iH1,iH2,iH3,H(LW(26)))
  nullify(iH1,iH2,iH3)

end subroutine dPRWF

end subroutine NATCT

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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_Setup_Add_Info()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Data for RPA tests, checking the setup.
! Testing these data is as much a test of the preceding SCF run as
! a test of the RPA setup.

use RPA_globals, only: nDel, nFro, nOcc, nSym, NuclearRepulsionEnergy, nVir, OccEn, Reference, VirEn
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: Tol, iUHF, l_orbitals, iSpin, iSym, ipO, ipV, i
real(kind=wp) :: Tst(8)
character(len=13) :: orbitals
character(len=*), parameter :: SecNam = 'RPA_Setup_Add_Info'
integer(kind=iwp), external :: RPA_iUHF, Cho_X_GetTol
real(kind=wp), external :: Cho_dSumElm, dDot_

! Check that molecular geometry is the expected one:
! nuclear repulsion energy
Tol = 12
call Add_Info('PotNuc',NuclearRepulsionEnergy,1,Tol)

! Check orbital spaces: sums and norms of orbital energies.
Tol = min(Cho_X_GetTol(2),2)
iUHF = RPA_iUHF()
if (iUHF == 1) then
  orbitals = ' orbital'
  l_orbitals = 8
else if (iUHF == 2) then
  orbitals = ' spin-orbital'
  l_orbitals = 13
else
  write(u6,'(A,I6)') 'iUHF=',iUHF
  call RPA_Warn(3,SecNam//': iUHF error')
  orbitals = ' '
  l_orbitals = 1
end if
Tst(:) = Zero
do iSpin=1,iUHF
  ipO = 1
  ipV = 1
  do iSym=1,nSym
    Tst(1) = Tst(1)+Cho_dSumElm(OccEn(ipO,iSpin),nFro(iSym,iSpin))
    Tst(2) = Tst(2)+dDot_(nFro(iSym,iSpin),OccEn(ipO,iSpin),1,OccEn(ipO,iSpin),1)
    ipO = ipO+nFro(iSym,iSpin)
    Tst(3) = Tst(3)+Cho_dSumElm(OccEn(ipO,iSpin),nOcc(iSym,iSpin))
    Tst(4) = Tst(4)+dDot_(nOcc(iSym,iSpin),OccEn(ipO,iSpin),1,OccEn(ipO,iSpin),1)
    ipO = ipO+nOcc(iSym,iSpin)
    Tst(5) = Tst(5)+Cho_dSumElm(VirEn(ipV,iSpin),nVir(iSym,iSpin))
    Tst(6) = Tst(6)+dDot_(nVir(iSym,iSpin),VirEn(ipV,iSpin),1,VirEN(ipV,iSpin),1)
    ipV = ipV+nVir(iSym,iSpin)
    Tst(7) = Tst(7)+Cho_dSumElm(VirEn(ipV,iSpin),nDel(iSym,iSpin))
    Tst(8) = Tst(8)+dDot_(nDel(iSym,iSpin),VirEn(ipV,iSpin),1,VirEn(ipV,iSpin),1)
    ipV = ipV+nDel(iSym,iSpin)
  end do
end do
do i=2,8,2
  Tst(i) = sqrt(Tst(i))
end do
call Add_Info(Reference//orbitals(1:l_orbitals)//' energy',Tst,8,Tol)

end subroutine RPA_Setup_Add_Info

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

subroutine espf_mltp(natom,MltOrd,nMult,nGrdPt,TTT,Mltp,Grid,IsMM,Ext,iPL)
! Compute the expected values of the ESPF operators
! i.e. charges or charges+dipoles

use espf_global, only: MxExtPotComp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom, MltOrd, nMult, nGrdPt, IsMM(natom), iPL
real(kind=wp), intent(in) :: TTT(nGrdPt,nMult), Grid(3,nGrdPt), Ext(MxExtPotComp,natom)
real(kind=wp), intent(out) :: Mltp(nMult)
#include "LenIn.fh"
integer(kind=iwp) :: iAddPot, iAtom, iMult, jMlt, kOrd, kPnt, ncmp
real(kind=wp) :: opnuc(1), SumOfChg, TotElecInt
real(kind=wp), allocatable :: Charge(:), D2(:), EI(:)
character(len=LenIn), allocatable :: CName(:)
character(len=3), parameter :: Axis(3) = [' x ',' y ',' z ']

if (iPL >= 5) then
  write(u6,*) ' In espf_mltp:',MltOrd,nMult,nGrdPt
  call RecPrt('TTT',' ',TTT,nGrdPt,nMult)
end if
call mma_allocate(Charge,natom,label='Charge')
call Get_Nuc_Charge_All(Charge,natom)
iMult = 1
do iAtom=1,natom
  if (IsMM(iAtom) == 0) then
    Mltp(iMult) = Charge(iAtom)
    if (MltOrd /= 1) Mltp(iMult+1:iMult+MltOrd-1) = Zero
    iMult = iMult+MltOrd
  end if
end do
call mma_deallocate(Charge)
opnuc = Zero
ncmp = 1
iAddPot = -2
call mma_allocate(D2,nGrdPt,label='dESPF2')
call DrvPot(Grid,opnuc,ncmp,D2,nGrdPt,iAddPot)
if (iPL >= 5) call RecPrt('PV',' ',D2,nGrdPt,1)
do jMlt=1,nMult
  do kPnt=1,nGrdPt
    Mltp(jMlt) = Mltp(jMlt)+D2(kPnt)*TTT(kPnt,jMlt)
  end do
end do
call mma_deallocate(D2)

if (iPL >= 3) then
  write(u6,'(/,A,/)') '      Expectation values of the ESPF operators:'
  call mma_allocate(EI,natom,label='ElecInt')
  call mma_allocate(CName,natom,label='CName')
  call Get_CArray('Unique Atom Names',CName,LenIn*natom)
  SumOfChg = Zero
  jMlt = 1
  TotElecInt = Zero
  do iAtom=1,natom
    EI(iAtom) = Zero
    if (IsMM(iAtom) == 1) cycle
    do kOrd=0,MltOrd-1
      if (kOrd == 0) then
        write(u6,1000) CName(iAtom),Mltp(jMlt)
        SumOfChg = SumOfChg+Mltp(jMlt)
      else
        write(u6,1001) Axis(kOrd),Mltp(jMlt+kOrd)
      end if
      EI(iAtom) = EI(iAtom)+Mltp(jMlt+kOrd)*Ext(kOrd+1,iAtom)
    end do
    jMlt = jMlt+MltOrd
    TotElecInt = TotElecInt+EI(iAtom)
  end do
  write(u6,1002) SumOfChg
  write(u6,1003) TotElecInt
  do iAtom=1,natom
    if (IsMM(iAtom) == 0) write(u6,1004) CName(iAtom),EI(iAtom)
  end do
  write(u6,'(A)')
  call mma_deallocate(EI)
  call mma_deallocate(CName)
end if

return

1000 format('        Charge on ',A,'      = ',F10.4)
1001 format('        + Dipole component ',A3,'= ',F10.4)
1002 format(/,'      Total ESPF charge     = ',F10.4,/)
1003 format(/,'      Total ESPF QM/MM interaction energy = ',F10.6,/)
1004 format('        ',A,' individual contribution =',F10.6)

end subroutine espf_mltp

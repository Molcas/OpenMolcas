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

subroutine LoProp_Print(rMP,nij,nElem,nAtoms,Q_Nuc,LblCnt,lSave)

use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
#include "WrkSpc.fh"
integer(kind=iwp), intent(in) :: nij, nElem, nAtoms
real(kind=wp), intent(in) :: rMP(nij,nElem), Q_Nuc(nAtoms)
character(len=LenIn4), intent(in) :: LblCnt(nAtoms)
logical(kind=iwp), intent(in) :: lSave
integer(kind=iwp) :: i, iAtom, iEnd, ij, Inc, iPL, ipLPChg, iSt, mAtoms
real(kind=wp) :: E_Charge(MxAtom), Q_Charge(MxAtom)
character(len=LenIn) :: Lbl(MxAtom)
character(len=120) :: Banner_Line(3)
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

!                                                                      *
!***********************************************************************
!                                                                      *
! Get the print level

iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
if (iPL < 2) return
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Binom()
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over all domains

write(u6,*)
Banner_Line(1) = 'LoProp Charges per center'
#ifdef _DEBUGPRINT_
Banner_Line(2) = ' '
Banner_Line(3) = 'Note that this charge analysis only makes sense if the orbital basis is of true AO type!'
call Banner(Banner_Line,3,120)
#else
write(u6,'(6X,A)') trim(Banner_Line(1))
#endif

! Collect data

mAtoms = 0
ij = 0
do iAtom=1,nAtoms
  ij = iAtom*(iAtom+1)/2
  if ((LblCnt(iAtom)(LENIN1:LENIN4) == ':E  ') .or. (LblCnt(iAtom)(LENIN1:LENIN4) == '    ')) then
    mAtoms = mAtoms+1
    Q_Charge(mAtoms) = Q_nuc(iAtom)
    E_Charge(mAtoms) = rMP(ij,1)
    Lbl(mAtoms) = LblCnt(iAtom)(1:LENIN)
  end if
end do
if (lSave) then
  call GetMem('LoProp Chg','Allo','Real',ipLPChg,mAtoms)
  call dCopy_(mAtoms,Q_Charge,1,Work(ipLPChg),1)
  call daxpy_(mAtoms,One,E_Charge,1,Work(ipLPChg),1)
  call Put_dArray('LoProp Charge',Work(ipLPChg),mAtoms)
  call GetMem('LoProp Chg','Free','Real',ipLPChg,mAtoms)
end if

! Print out the stuff!

Inc = 10
do iSt=1,mAtoms,Inc
  iEnd = min(mAtoms,iSt+Inc-1)
  write(u6,*)
  write(u6,'(/16X,10(3X,A))') (Lbl(i),i=iSt,iEnd)
  write(u6,'(6X,A,10F9.4)') 'Nuclear   ',(Q_Charge(i),i=iSt,iEnd)
  write(u6,'(6X,A,10F9.4)') 'Electronic',(E_Charge(i),i=iSt,iEnd)
  write(u6,*)
  write(u6,'(6X,A,10F9.4)') 'Total     ',(Q_Charge(i)+E_Charge(i),i=iSt,iEnd)
end do

#ifdef _DEBUGPRINT_
write(u6,*)
call Banner(Banner_Line,0,120)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine LoProp_Print

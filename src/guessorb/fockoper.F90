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
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine creates a model fock operator based on atomic orbital   *
! energies.                                                            *
!                                                                      *
! This is a very preliminary routine for testing purposes. It relies   *
! on the basis set being of ANO type.                                  *
!                                                                      *
! Absolutely NOT to be used for production!!!!                         *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************

subroutine FockOper(RC,Fock)

use GuessOrb_Global, only: Label, MxAtom, AtName, nBas, nNuc, nSym, xCharge
use Constants, only: Zero, One, Three, Six, Seven, Eight
use Definitions, only: wp, iwp, u6

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp), intent(out) :: RC
real(kind=wp), intent(out) :: Fock(*)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp), parameter :: MxComp = 4
logical(kind=iwp) :: Debug, Trace, Found
integer(kind=iwp) :: iSym, iBas, iOff, iNuc, nBasTot, iUse(MxAtom,MxComp), nData, i, k, lenName
real(kind=wp) :: energy
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
Debug = .false.
Trace = .false.
lenName = len(AtName)
if (Trace) write(u6,*) '>>> Entering fockoper'
RC = 0
!----------------------------------------------------------------------*
! Setup various counters.                                              *
!----------------------------------------------------------------------*
nBasTot = 0
do iSym=1,nSym
  nBasTot = nBasTot+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Is Fock operator on disk?                                            *
!----------------------------------------------------------------------*
call Qpg_dArray('Eorb',Found,nData)
if (Found) then
  call Get_dArray('Eorb',Fock,nData)
  if (Debug) then
    write(u6,*)
    write(u6,*) 'Found Eorb'
    write(u6,'(10f12.6)') (Fock(i),i=1,nData)
    write(u6,*)
  end if
  if (.true.) return
end if
if (.true.) then
  RC = 1
  return
end if
write(u6,*) '***'
write(u6,*) '*** Warning: using built-in fock operator'
write(u6,*) '***'
!----------------------------------------------------------------------*
! Create model Fock operator.                                          *
!----------------------------------------------------------------------*
iOff = 0
do iSym=1,nSym
  if (Debug) then
    write(u6,*) '***'
    write(u6,*) '*** Symmetry',iSym
    write(u6,*) '***'
  end if
  do i=1,MxAtom
    do k=1,MxComp
      iUse(i,k) = 0
    end do
  end do
  do iBas=1,nBas(iSym)
    iNuc = 0
    do i=1,nNuc
      if (AtName(i) == Label(iBas+iOff)(1:lenName)) iNuc = i
    end do
    if (Debug) then
      write(u6,'(2(a,i3),3a,i3,f6.2)') 'iSym:',iSym,' iBas:',iBas,' = ',Label(iBas+iOff)(1:lenName),Label(iBas+iOff)(lenName+1:), &
                                       iNuc,xCharge(iNuc)
    end if
    if (iNuc == 0) then
      call SysAbendMsg('fockoper','Fatal','001')
    end if
    energy = Zero
    if (abs(xCharge(iNuc)-One) < 1.0e-3_wp) then
      if (Label(iBas+iOff)(lenName+1:) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -0.50000_wp
      end if
    else if (abs(xCharge(iNuc)-Three) < 1.0e-3_wp) then
      if (Label(iBas+iOff)(lenName+1:) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -2.47773_wp
        if (iUse(iNuc,1) == 2) energy = -0.19632_wp
      end if
    else if (abs(xCharge(iNuc)-Six) < 1.0e-3_wp) then
      if (Label(iBas+iOff)(lenName+1:) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -11.32554_wp
        if (iUse(iNuc,1) == 2) energy = -0.70563_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.43335_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.43335_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.43335_wp
      end if
    else if (abs(xCharge(iNuc)-Seven) < 1.0e-3_wp) then
      if (Label(iBas+iOff)(lenName+1:) == '01s ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -15.62909_wp
        if (iUse(iNuc,1) == 2) energy = -0.94531_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.56758_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.56758_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.56758_wp
      end if
    else if (abs(xCharge(iNuc)-Eight) < 1.0e-3_wp) then
      if (Label(iBas+iOff)(lenName+1:) == '01s ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -20.66866_wp
        if (iUse(iNuc,1) == 2) energy = -1.24433_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.63192_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.63192_wp
      else if (Label(iBas+iOff)(lenName+1:) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.63192_wp
      end if
    else
      call SysAbendMsg('fockoper','Fatal','002')
    end if
    Fock(iOff+iBas) = energy
  end do
  iOff = iOff+nBas(iSym)
end do
!----------------------------------------------------------------------*
! Done, deallocate the rest.                                           *
!----------------------------------------------------------------------*
if (trace) write(u6,*) '<<< Exiting fockoper'

return

end subroutine FockOper

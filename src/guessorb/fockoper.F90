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

use GuessOrb_Global, only: Label, LenIn, LenIn1, LenIn8, MxAtom, Name, nBas, nNuc, nSym, xCharge

implicit none
!----------------------------------------------------------------------*
! Parameters                                                           *
!----------------------------------------------------------------------*
integer MxComp
parameter(MxComp=4)
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 Fock(*)
integer RC
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical Debug
logical Trace
logical Found
integer iSym
integer iBas
integer iOff
integer iNuc
integer nBasTot
integer iUse(MxAtom,MxComp)
integer nData
integer i, k
real*8 energy
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
Debug = .false.
Trace = .false.
if (Trace) write(6,*) '>>> Entering fockoper'
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
    write(6,*)
    write(6,*) 'Found Eorb'
    write(6,'(10F12.6)')(Fock(i),i=1,nData)
    write(6,*)
  end if
  if (.true.) return
end if
if (.true.) then
  RC = 1
  return
end if
write(6,*) '***'
write(6,*) '*** Warning: using built in fock operator'
write(6,*) '***'
!----------------------------------------------------------------------*
! Create model Fock operator.                                          *
!----------------------------------------------------------------------*
iOff = 0
do iSym=1,nSym
  if (Debug) then
    write(6,*) '***'
    write(6,*) '*** Symmetry',iSym
    write(6,*) '***'
  end if
  do i=1,MxAtom
    do k=1,MxComp
      iUse(i,k) = 0
    end do
  end do
  do iBas=1,nBas(iSym)
    iNuc = 0
    do i=1,nNuc
      if (Name(i) == Label(iBas+iOff)(1:LENIN)) iNuc = i
    end do
    if (Debug) then
      write(6,'(2(a,i3),3a,i3,f6.2)') 'iSym:',iSym,' iBas:',iBas,' = ',Label(iBas+iOff)(1:LENIN),Label(iBas+iOff)(LENIN1:LENIN8), &
                                      iNuc,xCharge(iNuc)
    end if
    if (iNuc == 0) then
      call SysAbendMsg('fockoper','Fatal','001')
    end if
    energy = 0.0d0
    if (abs(xCharge(iNuc)-1.0d0) < 1.0d-3) then
      if (Label(iBas+iOff)(LENIN1:LENIN8) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -0.50000d0
      end if
    else if (abs(xCharge(iNuc)-3.0d0) < 1.0d-3) then
      if (Label(iBas+iOff)(LENIN1:LENIN8) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -2.47773d0
        if (iUse(iNuc,1) == 2) energy = -0.19632d0
      end if
    else if (abs(xCharge(iNuc)-6.0d0) < 1.0d-3) then
      if (Label(iBas+iOff)(LENIN1:LENIN8) == '01s     ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -11.32554d0
        if (iUse(iNuc,1) == 2) energy = -0.70563d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.43335d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.43335d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.43335d0
      end if
    else if (abs(xCharge(iNuc)-7.0d0) < 1.0d-3) then
      if (Label(iBas+iOff)(LENIN1:LENIN8) == '01s ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -15.62909d0
        if (iUse(iNuc,1) == 2) energy = -0.94531d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.56758d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.56758d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.56758d0
      end if
    else if (abs(xCharge(iNuc)-8.0d0) < 1.0d-3) then
      if (Label(iBas+iOff)(LENIN1:LENIN8) == '01s ') then
        iUse(iNuc,1) = iUse(iNuc,1)+1
        if (iUse(iNuc,1) == 1) energy = -20.66866d0
        if (iUse(iNuc,1) == 2) energy = -1.24433d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02px    ') then
        iUse(iNuc,2) = iUse(iNuc,2)+1
        if (iUse(iNuc,2) == 1) energy = -0.63192d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02py    ') then
        iUse(iNuc,3) = iUse(iNuc,3)+1
        if (iUse(iNuc,3) == 1) energy = -0.63192d0
      else if (Label(iBas+iOff)(LENIN1:LENIN8) == '02pz    ') then
        iUse(iNuc,4) = iUse(iNuc,4)+1
        if (iUse(iNuc,4) == 1) energy = -0.63192d0
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
if (trace) write(6,*) '<<< Exiting fockoper'

return

end subroutine FockOper

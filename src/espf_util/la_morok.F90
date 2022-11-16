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

subroutine LA_Morok(nAtom,CorG,iMode)
! Morokuma's scaling scheme:
!   k = (q_LA - q_QM)/(q_MM - q_QM), k = constant
!
! iMode = 1 => the LA gradient is distributed
! on the frontier QM and MM atoms according to:
! dE/dq_QM = dE/dq_QM + dE/dq_LA * dq_LA/dq_QM
! dE/dq_MM = dE/dq_MM + dE/dq_LA * dq_LA/dq_MM
!
! iMode = 2 => the LA position is updated
! q_LA = q_QM + k * (q_MM - q_QM)

implicit real*8(a-h,o-z)
real*8 CorG(3,nAtom)
#include "espf.fh"
#include "stdalloc.fh"
logical Exist, Exist2, isOkLA, isOkMM, isOkQM, lMorok
logical DoTinker, DoGromacs
character*10 ESPFKey
character*180 Line
character*180 Get_Ln
character*256 Message
external Get_Ln
integer, dimension(:), allocatable :: AT, GroToMol
integer, dimension(:,:), allocatable :: DefLA
real*8, dimension(:), allocatable :: FactLA

iPL = iPL_espf()
lMorok = .false.
DoTinker = .false.
DoGromacs = .false.
call F_Inquire('ESPF.DATA',Exist)
if (Exist) then
  IPotFl = IsFreeUnit(1)
  call Molcas_Open(IPotFl,'ESPF.DATA')
10 Line = Get_Ln(IPotFl)
  ESPFKey = Line(1:10)
  if (ESPFKey == 'LA_MOROK  ') then
    lMorok = .true.
  else if (ESPFKey == 'TINKER    ') then
    DoTinker = .true.
  else if (ESPFKey == 'GROMACS   ') then
    DoGromacs = .true.
  else if (ESPFKey == 'ENDOFESPF ') then
    goto 11
  end if
  goto 10
11 close(IPotFl)
end if
if (.not. lMorok) return
#ifdef _DEBUGPRINT_
iPL = 4
call RecPrt('LA_Morok: coord or grad:',' ',CorG,3,nAtom)
#endif

! Tinker part

Exist = .false.
if (DoTinker) then
  call F_Inquire('QMMM',Exist)
end if
if (Exist) then
  ITkQMMM = IsFreeUnit(25)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'LAH') /= 0) then
      call Get_I1(2,iLA)
      call Get_I1(3,iMM)
      call Get_I1(4,iQM)
      call Get_F1(5,Fact)
      if ((iMM < 1) .or. (iQM < 1)) then
        write(6,*) 'LA_Morok: link atoms badly defined'
        write(6,*) '          check each LA connectivity'
        call Quit_OnUserError()
      end if
#     ifdef _DEBUGPRINT_
      write(6,*)
      write(6,*) 'LA_Morok: LAH ',iLA,' between ',iQM,' and ',iMM
      write(6,*) '          scaling factor : ',Fact
#     endif
      if (iMode == 1) then
        if (iPL >= 2) write(6,*) 'LA_Morok: scaling gradients'
        do iXYZ=1,3
          CorG(iXYZ,iQM) = CorG(iXYZ,iQM)+CorG(iXYZ,iLA)*(One-Fact)
          CorG(iXYZ,iMM) = CorG(iXYZ,iMM)+CorG(iXYZ,iLA)*Fact
          CorG(iXYZ,iLA) = Zero
        end do
      else if (iMode == 2) then
        if (iPL >= 2) write(6,*) 'LA_Morok: updating positions'
        do iXYZ=1,3
          CorG(iXYZ,iLA) = CorG(iXYZ,iQM)+(CorG(iXYZ,iMM)-Corg(iXYZ,iQM))*Fact
        end do
      else
        write(6,*) 'LA_Morok: wrong iMode'
        call Quit_OnUserError()
      end if
    end if
  end do
  close(ITkQMMM)
end if

! Gromacs part

Exist = .false.
if (DoGromacs) then
  call Qpg_iArray('LA Def',Exist,nLink)
end if
if (Exist) then
  nLink = nLink/3
  call mma_allocate(DefLA,3,nLink)
  call mma_allocate(FactLA,nLink)
  call Get_iArray('LA Def',DefLA,3*nLink)
  call Get_dArray('LA Fact',FactLA,nLink)
  ! Check for consistency
  call Qpg_iArray('Atom Types',Exist2,nTot)
  if (.not. Exist2) then
    Message = 'LA_Morok: no atom type info on runfile'
    call WarningMessage(2,Message)
    call Abend()
  end if
  call mma_allocate(AT,nTot)
  call Get_iArray('Atom Types',AT,nTot)
  do iLink=1,nLink
    isOkLA = AT(DefLA(1,iLink)) == QM
    isOkQM = AT(DefLA(2,iLink)) == QM
    isOkMM = AT(DefLA(3,iLink)) == MMI
    if (.not. (isOkLa .and. isOkQM .and. isOkMM)) then
      Message = 'Link atoms badly defined. Check input!'
      call WarningMessage(2,Message)
      call Abend()
    end if
  end do
  ! Generate vector for translating from Gromacs to Molcas numbering
  call mma_allocate(GroToMol,nTot)
  iAtIn = 1
  iAtOut = 1
  do iAt=1,nTot
    if ((AT(iAt) == QM) .or. (AT(iAt) == MMI)) then
      GroToMol(iAt) = iAtIn
      iAtIn = iAtIn+1
    else if (AT(iAt) == MMO) then
      GroToMol(iAt) = iAtOut
      iAtOut = iAtOut+1
    else
      Message = 'LA_Morok: unknown atom type'
      call WarningMessage(2,Message)
      call Abend()
    end if
  end do
  if (iMode == 1) then
    ! Apply Morokuma scheme to gradient...
    if (iPL >= 2) then
      write(6,*) 'Applying Morokuma scheme to gradient'
    end if
    do iLink=1,nLink
      iLA = GroToMol(DefLA(1,iLink))
      iQM = GroToMol(DefLA(2,iLink))
      iMM = GroToMol(DefLA(3,iLink))
      Fact = FactLA(iLink)
      do ixyz=1,3
        CorG(ixyz,iQM) = CorG(ixyz,iQM)+CorG(ixyz,iLA)*(1-Fact)
        CorG(ixyz,iMM) = CorG(ixyz,iMM)+CorG(ixyz,iLA)*Fact
        CorG(ixyz,iLA) = Zero
      end do
    end do
  else if (iMode == 2) then
    ! ...or to position
    if (iPL >= 2) then
      write(6,*) 'Applying Morokuma scheme to positions'
    end if
    do iLink=1,nLink
      iLA = GroToMol(DefLA(1,iLink))
      iQM = GroToMol(DefLA(2,iLink))
      iMM = GroToMol(DefLA(3,iLink))
      Fact = FactLA(iLink)
      do ixyz=1,3
        CorG(ixyz,iLA) = CorG(ixyz,iQM)+(CorG(ixyz,iMM)-CorG(ixyz,iQM))*Fact
      end do
    end do
  else
    Message = 'LA_Morok: wrong iMode'
    call WarningMessage(2,Message)
    call Abend()
  end if
  call mma_deallocate(DefLA)
  call mma_deallocate(FactLA)
  call mma_deallocate(AT)
  call mma_deallocate(GroToMol)
end if

#ifdef _DEBUGPRINT_
call RecPrt('LA_Morok: coord or grad:',' ',CorG,3,nAtom)
#endif

end subroutine LA_Morok

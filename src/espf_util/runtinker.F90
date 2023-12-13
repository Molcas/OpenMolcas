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

subroutine RunTinker(nAtom,Cord,Mltp,First,IsMM,MltOrd,DynExtPot,iQMChg,nAtMM,StandAlone,DoDirect)

use espf_global, only: MxExtPotComp
use Para_Info, only: MyRank
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtom, IsMM(nAtom), MltOrd, iQMChg
real(kind=wp), intent(in) :: Cord(3,nAtom), Mltp(*)
logical(kind=iwp), intent(in) :: First, StandAlone, DoDirect
logical(kind=iwp), intent(inout) :: DynExtPot
integer(kind=iwp), intent(inout) :: nAtMM
integer(kind=iwp) :: iAtom, iLast, iMlt, iPL, iq, iRelax, iSomething, istatus, ITkPot, ITkQMMM, jLast, Lu, mLine, nLine, nMMq, RC
character(len=256) :: TkLine
character(len=180) :: Line
character(len=12) :: ExtPotFormat
real(kind=wp), allocatable :: ESPF(:,:), LPC(:), MMqx(:,:), Mull(:)
integer(kind=iwp), external :: iCLast, iPL_espf, IsFreeUnit
character(len=180), external :: Get_Ln

iPL = iPL_espf()
write(ExtPotFormat,'(a4,i2,a6)') '(I4,',MxExtPotComp,'F13.8)'

! Always update the coordinates of the tinker xyz file
! WARNING: coordinates are converted to Angstroms
! This is done through a communication file: Project.qmmm

! Only call Tinker on the master node

ITkQMMM = 1
if (MyRank == 0) then
  ITkQMMM = IsFreeUnit(ITkQMMM)
  call Molcas_Open(ITkQMMM,'QMMM')

  ! The MM subsystem can relax (microiterations, MD, ...) unless:
  ! 1) there are no QM multipoles
  ! 2) this is a call to retrieve MM energy/gradient/electrostatic potential only

  iRelax = 1
  if (First .or. (.not. StandAlone)) iRelax = 0
  if (DoDirect) then
    write(ITkQMMM,1000) iRelax,-1
  else
    write(ITkQMMM,1000) iRelax,MltOrd/4
  end if
  do iAtom=1,nAtom
    write(ITkQMMM,1010) Cord(:,iAtom)*Angstrom
  end do
  if (.not. First) then
    write(ITkQMMM,'(A)') 'Multipoles'
    if (iQMChg == 0) then
      if (iPL >= 3) write(u6,'(A)') ' Multipoles passed to Tinker'
      iMlt = 1
      do iAtom=1,nAtom
        if (IsMM(iAtom) == 0) then
          if (MltOrd == 1) then
            write(ITkQMMM,'(I6,4F15.8)') iAtom,Mltp(iMlt),Zero,Zero,Zero
          else
            write(ITkQMMM,'(I6,4F15.8)') iAtom,Mltp(iMlt:iMlt+3)
          end if
          iMlt = iMlt+MltOrd
        else
          write(ITkQMMM,'(I6,4F15.8)') iAtom,Zero,Zero,Zero,Zero
        end if
      end do
    else if (iQMChg == 1) then
      if ((StandAlone .and. (iPL >= 2)) .or. ((.not. StandAlone) .and. (iPL >= 3))) &
        write(u6,'(A)') ' ESPF multipoles passed to Tinker'
      iMlt = 1
      do iAtom=1,nAtom
        if (IsMM(iAtom) == 0) then
          if (MltOrd == 1) then
            write(ITkQMMM,'(I6,4F15.8)') iAtom,Mltp(iMlt),Zero,Zero,Zero
          else
            write(ITkQMMM,'(I6,4F15.8)') iAtom,Mltp(iMlt:iMlt+3)
          end if
          iMlt = iMlt+MltOrd
        else
          write(ITkQMMM,'(I6,4F15.8)') iAtom,Zero,Zero,Zero,Zero
        end if
      end do
    else if (iQMChg == 2) then
      if ((StandAlone .and. (iPL >= 2)) .or. ((.not. StandAlone) .and. (iPL >= 3))) &
        write(u6,'(A)') ' Mulliken charges passed to Tinker'
      call mma_allocate(Mull,nAtom,Label='Mull')
      call Get_dArray('Mulliken Charge',Mull,nAtom)
      do iAtom=1,nAtom
        if (IsMM(iAtom) == 0) write(ITkQMMM,'(I6,4F15.8)') iAtom,Mull(iAtom),Zero,Zero,Zero
      end do
      call mma_deallocate(Mull)
    else if (iQMChg == 3) then
      if ((StandAlone .and. (iPL >= 2)) .or. ((.not. StandAlone) .and. (iPL >= 3))) &
        write(u6,'(A)') ' LoProp charges passed to Tinker'
      call mma_allocate(LPC,nAtom,Label='LPC')
      call Get_dArray('LoProp Charge',LPC,nAtom)
      do iAtom=1,nAtom
        if (IsMM(iAtom) == 0) write(ITkQMMM,'(I6,4F15.8)') iAtom,LPC(iAtom),Zero,Zero,Zero
      end do
      call mma_deallocate(LPC)
    end if
  end if
  close(ITkQMMM)

  ! Tinker is running

  call Getenvf('Project ',Line)
  mLine = len(Line)
  iLast = iCLast(Line,mLine)
  Line = Line(1:iLast)//'.xyz'
  Line = Line(1:iLast)//'.key'
  Line = '/tkr2qm_s ${Project}.xyz>${Project}.Tinker.log'
  call Getenvf('TINKER ',TkLine)
  mLine = len(TkLine)
  iLast = iCLast(TkLine,mLine)
  if (iLast == 0) then
    call Getenvf('MOLCAS',TkLine)
    mLine = len(TkLine)
    iLast = iCLast(TkLine,mLine)
    TkLine = TkLine(1:iLast)//'/tinker/bin'
  end if
  iLast = iCLast(TkLine,mLine)
  nLine = len(Line)
  jLast = iCLast(Line,nLine)
  Line = TkLine(1:iLast)//Line(1:jLast)
  call StatusLine(' espf:',' Calling Tinker')
  RC = 0
  call Systemf(Line(1:iLast+jLast),RC)
end if
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  call GA_Sync()
  call PFGET_ASCII('TINKER.LOG')
  call PFGET_ASCII('QMMM')
  call PFGET_ASCII('ESPF.EXTPOT')
  call GA_Sync()
end if
#endif
if ((StandAlone .and. (iPL >= 2)) .or. ((.not. StandAlone) .and. (iPL >= 3))) then
  iSomething = 0
  Lu = 55
  Lu = IsFreeUnit(Lu)
  call Molcas_Open(Lu,'TINKER.LOG')
  do
    read(Lu,'(A)',iostat=istatus) Line
    if (istatus < 0) exit
    iSomething = iSomething+1
    nLine = len(Line)
    iLast = iCLast(Line,nLine)
    write(u6,*) Line(1:iLast)
  end do
  close(Lu)
  if (iSomething == 0) then
    write(u6,*) ' Something bad with Tinker: no output !'
    call Quit_OnUserError()
  end if
end if

! Tinker post-processing
! WARNING: all Tinker results must be converted to atomic units !!!
! Convert the ESPF external potential and derivatives to something
! understandable by molcas, stored in the ESPF.EXTPOT file

if (iPL >= 2) then
  write(u6,*) 'Back from Tinker'
  write(u6,*)
end if
call mma_allocate(ESPF,MxExtPotComp,nAtom,Label='ESPF')
ESPF(:,:) = Zero
ITkQMMM = IsFreeUnit(ITkQMMM)
call Molcas_Open(ITkQMMM,'QMMM')
Line = Get_Ln(ITkQMMM)
if (index(Line,'MMisOK') == 0) then
  write(u6,*) 'Something wrong happend with Tinker'
  call Abend()
end if
do while (index(Line,'TheEnd ') == 0)
  Line = Get_Ln(ITkQMMM)
  if (index(Line,'NMM ') /= 0) then
    call Get_I1(2,nAtMM)
  else if (index(Line,'ESPF1 ') /= 0) then
    call Get_I1(2,iAtom)
    call Get_F(3,ESPF(1:4,iAtom),4)
  else if (index(Line,'ESPF21 ') /= 0) then
    call Get_I1(2,iAtom)
    call Get_F(3,ESPF(5:7,iAtom),3)
  else if (index(Line,'ESPF22 ') /= 0) then
    call Get_I1(2,iAtom)
    call Get_F(3,ESPF(8:10,iAtom),3)
  else if (index(Line,'FullCoupling ') /= 0) then
    DynExtPot = .true.
  else if (index(Line,'MMq ') /= 0) then
    call Get_I1(2,nMMq)
    call mma_allocate(MMqx,4,nMMq,Label='MMqx')
    do iq=1,nMMq
      Line = Get_Ln(ITkQMMM)
      call Get_F(1,MMqx(1:4,iq),4)
    end do
  end if
end do
close(ITkQMMM)
ITkPot = IsFreeUnit(ITkQMMM)
call Molcas_Open(ITkPot,'ESPF.EXTPOT')
if (DoDirect) then
  write(ITkPot,'(I10,1X,I2)') nMMq,0
  do iq=1,nMMq
    MMqx(1:3,iq) = MMqx(1:3,iq)*Angstrom
    write(ITkPot,'(4F15.8)') MMqx(1:4,iq)
  end do
  call mma_deallocate(MMqx)
else
  write(ITkPot,'(I1)') 0
  do iAtom=1,nAtom
    ESPF(1,iAtom) = ESPF(1,iAtom)*Angstrom
    ESPF(2:4,iAtom) = ESPF(2:4,iAtom)*Angstrom**2
    ESPF(5:10,iAtom) = ESPF(5:10,iAtom)*Angstrom**3
    write(ITkPot,ExtPotFormat) iAtom,ESPF(:,iAtom)
  end do
end if
close(ITkPot)
call mma_deallocate(ESPF)

return

1000 format('Molcas  ',i2,2x,i2)
1010 format(3F15.8)

end subroutine RunTinker

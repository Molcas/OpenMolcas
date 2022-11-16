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

subroutine ReadIn_ESPF(natom,ipCord,ipExt,MltOrd,iRMax,DeltaR,Forces,Show_espf,ipIsMM,StandAlone,iGrdTyp,DoTinker,DoGromacs, &
                       DynExtPot,ipMltp,natMM,lMorok,DoDirect,ipGradCl,EnergyCl)

use external_centers
implicit real*8(a-h,o-z)
#include "espf.fh"
#include "opt_mmo.fh"
#include "stdalloc.fh"
#include "print.fh"
character*180 Key, Line, PotFile, UpKey
character*10 ESPFKey
character*12 ExtPotFormat
logical Convert, DoTinker_old, DoTinker, DoGromacs_old, DoGromacs, Exist, Forces, Show_espf, StandAlone, DynExtPot, lMorok_old, &
        lMorok, NoExt, DoDirect_old, DoDirect
save fift
data fift/1.5d1/
character*180 Get_Ln
external Get_Ln
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine RunTinker(nAtom,Cord,ipMltp,IsMM,MltOrd,DynExtPot,iQMChg,nAtMM,StandAlone,DoDirect)
    integer, intent(In) :: nAtom
    real*8, intent(In) :: Cord(3,nAtom)
    integer, intent(In) :: ipMltp
    integer, intent(In) :: IsMM(nAtom)
    integer, intent(In) :: MltOrd
    logical, intent(InOut) :: DynExtPot
    integer, intent(In) :: iQMChg
    integer, intent(InOut) :: nAtMM
    logical, intent(In) :: StandAlone
    logical, intent(In) :: DoDirect
  end subroutine RunTinker
end interface

#ifndef _GROMACS_
#include "macros.fh"
unused_var(ipGradCl)
unused_var(EnergyCl)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! If some keywords are not given, what are the defauts ?
! 3 cases:
!   1) ESPF.DATA does not exist:
!      MULT: MltOrd = 0 (monopole)
!      GRID: Type = PNT ; iRMax = 4 shells ; DeltaR = 1 angstrom
!      EXTE: MANDATORY
!   2) ESPF.DATA exists:
!      Get back all values from ESPF.DATA
!      a) If the keyword "Forces" is read, only the old values
!      are retained, so it is clever not to include any other
!      ESPF keyword
!      b) If the keyword "Forces" is not read, all the new values
!      are compared to the old ones.

! Initialize values

write(ExtPotFormat,'(a4,i2,a6)') '(I4,',MxExtPotComp,'F10.5)'
MltOrd = 1
MltOrd_old = MltOrd
nChg = -1
iGrdTyp = 1
iRMax = 4
!iGrdTyp = 2
!iRMax = 1
DeltaR = One/Angstrom
Convert = .false.
DoTinker = .false.
DoTinker_old = DoTinker
DoGromacs = .false.
DoGromacs_old = DoGromacs
Forces = .not. StandAlone
Show_espf = .false.
nMult = 0
DynExtPot = .false.
iQMChg = 1
natMM = 0
lMorok = .false.
lMorok_old = lMorok
NoExt = .false.
DoDirect = .false.
DoDirect_old = DoDirect
nOrd_ext = 0
MMIterMax = 0
ConvF = 2.0D-4*AuToKjPerMolNm

! Print level

iPL = iPL_espf()

! If the ESPF.DATA file exists, retrieve all data
! from it in "*_old" variables and arrays.

IPotFl = 15
PotFile = '***'
call F_Inquire('ESPF.DATA',Exist)
if (Exist) then
  IPotFl = IsFreeUnit(IPotFl)
  call Molcas_Open(IPotFl,'ESPF.DATA')
  do
    Line = Get_Ln(IPotFl)
    ESPFKey = Line(1:10)
    if (ESPFKey == 'MLTORD    ') then
      call Get_I1(2,MltOrd_old)
      ibla = 0
      do ii=0,MltOrd_old
        ibla = ibla+(ii+2)*(ii+1)/2
      end do
      MltOrd_old = ibla
    else if (ESPFKey == 'IRMAX     ') then
      call Get_I1(2,iRMax_old)
    else if (ESPFKey == 'DELTAR    ') then
      call Get_F1(2,DeltaR_old)
    else if (ESPFKey == 'GRIDTYPE  ') then
      call Get_I1(2,iGrdTyp_old)
    else if (ESPFKey == 'MULTIPOLE ') then
      call Get_I1(2,nMult)
      call GetMem('ESPFMltp','ALLO','REAL',ipMltp,nMult)
      do iMlt=1,nMult,MltOrd_old
        Line = Get_Ln(IPotFl)
        call Get_I1(1,iAt)
        call Get_F(2,Work(ipMltp+iMlt-1),MltOrd_old)
      end do
    else if (ESPFKey == 'TINKER    ') then
      DoTinker_old = .true.
    else if (ESPFKey == 'GROMACS   ') then
      DoGromacs_old = .true.
    else if (ESPFKey == 'DIRECT    ') then
      DoDirect_old = .true.
    else if (ESPFKey == 'LA_MOROK  ') then
      lMorok_old = .true.
    else if (ESPFKey == 'ENDOFESPF ') then
      exit
    end if
  end do
  close(IPotFl)
  iRMax = iRMax_old
  DeltaR = DeltaR_old
  MltOrd = MltOrd_old
  iGrdTyp = iGrdTyp_old
  DoTinker = DoTinker_old
  DoGromacs = DoGromacs_old
  lMorok = lMorok_old
end if

if (StandAlone) then

  ! Copy input from standard input to a local scratch file

  LuSpool = isFreeUnit(IPotFl)
  call SpoolInp(LuSpool)

  ! Locate "start of input"
  rewind(LuSpool)
  call RdNLst(LuSpool,'espf')

  do
    Key = Get_Ln(LuSpool)
    Line = Key
    call UpCase(Line)
    select case (Line(1:4))
      case ('MULT')
        !>>>>>>>>>>>>> MULT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        Key = Get_Ln(LuSpool)
        call Get_I1(1,MltOrd)
        if (MltOrd < 0) then
          write(6,'(A)') ' Error in espf/readin: MltOrd < 0!'
          call Quit_OnUserError()
        end if
        if (DoGromacs .and. (MltOrd > 0)) then
          write(6,'(A)') ' Error in espf/readin: Gromacs calculation requested with MltOrd > 0'
          write(6,'(A)') ' Only MltOrd = 0 is currently allowed'
          call Quit_OnUserError()
        end if
        if (MltOrd > 1) then
          write(6,'(A)') ' Error in espf/readin: MltOrd > 1 NYI!'
          call Quit_OnUserError()
        end if
        ibla = 0
        do ii=0,MltOrd
          ibla = ibla+(ii+2)*(ii+1)/2
        end do
        MltOrd = ibla

      case ('EXTE')
        !>>>>>>>>>>>>> EXTE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (Forces) cycle
        Key = Get_Ln(LuSpool)
        UpKey = Key
        call Upcase(UpKey)
        call Get_iNumber(Key(1:(index(Key,' ')-1)),ibla,iErr)
        if (iErr == 0) then
          PotFile = '* *'
          nChg = ibla

          ! nChg < 0: error
          ! nChg > 0: external potential is given as point charges and dipoles,
          ! like for the seward xfield keyword
          ! nchg = 0: external potential directly given on atom centers as:
          ! pot field_x field_y field_z dfield_xx dfield_yy dfield_zz
          ! dfield_xy dfield_xz dfield_yz (ONE LINE per CENTER)

          if (nChg < 0) then
            write(6,*) 'Error in readin_espf: nChg < 0!'
            call Quit_OnUserError()
          else if (nChg > 0) then
            Convert = (index(UpKey,'ANGSTROM') /= 0)
            nXF = nChg
            nData_XF = 7
            call mma_allocate(XF,nData_XF,nXF,Label='XF')
            do iChg=1,nChg
              Key = Get_Ln(LuSpool)
              call Get_F(1,XF(1,iChg),7)
              if (Convert) then
                XF(1:3,iChg) = XF(1:3,iChg)/Angstrom
                XF(5:7,iChg) = XF(5:7,iChg)/Angstrom
              end if
            end do
            Convert = .false.
          else
            do iAt=1,natom
              Key = Get_Ln(LuSpool)
              call Get_I1(1,jAt)
              if ((jAt < 1) .or. (jAt > natom)) then
                write(6,'(A)') ' Error in espf/readin: atom out of range.'
                call Quit_OnUserError()
              end if
              call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
            end do
          end if
        else
          iAt = index(UpKey,'NONE ')
          NoExt = (iAt /= 0)

          ! Is it a QM/MM computation ?

          iAt = index(UpKey,'TINKER ')
          DoTinker = (iAt /= 0)

          iAt = index(UpKey,'GROMACS ')
          if (iAt /= 0) then
#           ifdef _GROMACS_
            DoGromacs = .true.
#           else
            write(6,*) 'Interface to Gromacs not installed'
            call Quit_OnUserError()
#           endif
          end if

          DoDirect = (index(UpKey(iAt+7:120),'DIRECT') /= 0)
          ! tmp
          if (DoDirect) then
            write(6,*) 'Direct not yet implemented, abort.'
            call Quit_OnUserError()
          end if
          ! tmp

          ! What kind of charges Tinker or Gromacs will use in the microiterations

          if (NoExt) then
            PotFile = '* *'
          else if (DoTinker .or. DoGromacs) then
            PotFile = 'ESPF.EXTPOT'
            if (index(UpKey(iAt+7:120),'MULL') /= 0) then
              iQMChg = 2
            else if (index(UpKey(iAt+7:120),'LOPR') /= 0) then
              iQMChg = 3
            end if
            if (DoDirect .and. (iQMChg == 1)) iQMChg = 2
          else
            PotFile = Key(1:len(Key))
          end if
        end if

      case ('GRID')
        !>>>>>>>>>>>>> GRID <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        Key = Get_Ln(LuSpool)
        call Upcase(Key)
        if (index(Key,'GEPOL') /= 0) then
          iGrdTyp = 2
          call Get_I1(2,iRMax)
          if ((iRMax <= 0) .or. (iRMax > 4)) then
            write(6,'(A)') 'Error in readin_espf: 1 <= iRMax <= 4 !!!'
            call Quit_OnUserError()
          end if
        else if (index(Key,'PNT') /= 0) then
          iGrdTyp = 1
          call Get_I1(2,iRMax)
          if (iRMax <= 0) then
            write(6,'(A)') 'Error in espf/readin: iRMax < 1 !!!'
            call Quit_OnUserError()
          end if
          call Get_F1(3,DeltaR)
          if (DeltaR <= Zero) then
            write(6,'(A)') 'Error in espf/readin: DeltaR < 0.0 !!!'
            call Quit_OnUserError()
          end if
          if (index(Key,'ANGSTROM') /= 0) Convert = .true.
          if (Convert) DeltaR = DeltaR/Angstrom
          Convert = .false.
        else
          write(6,'(A)') 'Unrecognized GRID: GEPOL or PNT(default)'
          call Quit_OnUserError()
        end if

      case ('FORC')
        !>>>>>>>>>>>>> FORC <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (.not. Exist) then
          write(6,*) 'Error! Forces: the ESPF data are missing'
          call Quit_OnUserError()
        end if
        if (DoTinker) then
          write(6,*) 'Please erase the @Tinker call together with Forces'
          call Quit_OnUserError()
        end if
        Forces = .true.
        write(6,'(/,A,/,A)') ' This ESPF run will compute energy gradient',' Any other keyword is ignored !'

        ! Here I assume all I need can be retrieved from the $Project.espf file

      case ('SHOW')
        !>>>>>>>>>>>>> SHOW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        Show_espf = .true.

      case ('LAMO')
        !>>>>>>>>>>>>> LAMO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        lMorok = .true.
        if (iPL >= 2) write(6,'(A)') ' Morokuma scheme on'

      case ('MMIT')
        !>>>>>>>>>>>>> MMIT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (.not. DoGromacs) then
          write(6,'(A)') ' MM microiterations only available with Gromacs'
          call Quit_OnUserError()
        end if
        Line = Get_Ln(LuSpool)
        call Get_I1(1,MMIterMax)

      case ('MMCO')
        !>>>>>>>>>>>>> MMCO <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        Line = Get_Ln(LuSpool)
        call Get_F1(1,ConvF)
        ConvF = ConvF*AuToKjPerMolNm

      case ('END ')
        !>>>>>>>>>>>>> END  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        exit

      case default
        if (.not. Exist) then
          write(6,*) ' Unidentified keyword:',Key
          call FindErrorLine()
          call Quit_OnUserError()
        end if

    end select
  end do
end if

! Remove local copy of standard input

if (StandAlone) close(LuSpool)

! "Forces" case: retrieve all the data and update the MM gradient

if (Forces) then
  MltOrd = MltOrd_old
  iRMax = iRMax_old
  DeltaR = DeltaR_old
  iGrdTyp = iGrdTyp_old
  DoTinker = DoTinker_old
  DoGromacs = DoGromacs_old
  iQMChg = 0
  if (DoTinker) call RunTinker(natom,Work(ipCord),ipMltp,iWork(ipIsMM),MltOrd,DynExtPot,iQMchg,natMM,StandAlone,DoDirect)
# ifdef _GROMACS_
  if (DoGromacs) call RunGromacs(natom,Work(ipCord),ipMltp,MltOrd,Forces,ipGradCl,EnergyCl)
# endif
  if (nAtMM /= 0) write(6,*) 'MM gradients have been updated'
  lMorok = lMorok_old
  IPotFl = IsFreeUnit(IPotFl)
  call Molcas_Open(IPotFl,'ESPF.EXTPOT')
  Line = Get_Ln(IPotFl)
  call Get_I1(1,nChg)
  do iAt=1,natom
    Line = Get_Ln(IPotFl)
    call Get_I1(1,jAt)
    call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
  end do
  close(IPotFl)

else if (NoExt) then

  ! No external potential

  nChg = -1
  write(6,'(/,A)') ' No external electrostatic potential'

else if (nChg == -1) then

  ! External potential read from a file

  if (DoTinker) call RunTinker(natom,Work(ipCord),ipMltp,iWork(ipIsMM),MltOrd,DynExtPot,iQMChg,natMM,StandAlone,DoDirect)
# ifdef _GROMACS_
  if (DoGromacs) call RunGromacs(natom,Work(ipCord),ipMltp,MltOrd,Forces,ipGradCl,EnergyCl)
# endif
  LuSpool = IsFreeUnit(1)
  call Molcas_Open(LuSpool,PotFile(1:(index(PotFile,' ')-1)))
  if (iPL >= 3) write(6,'(/,A,A)') ' External potential read in ',PotFile(1:(index(PotFile,' ')-1))
  Key = Get_Ln(LuSpool)
  UpKey = Key
  call Upcase(UpKey)
  call Get_I1(1,nChg)
  if (nChg < 0) then
    write(6,*) 'Error in readin_espf: nChg < 0!'
    call Quit_OnUserError()
  else if (nChg > 0) then
    call Get_I1(2,nOrd_ext)
    Convert = (index(UpKey,'ANGSTROM') /= 0)
    nData_XF = 4+3*nOrd_ext
    nXF = nChg
    call mma_allocate(XF,nData_XF,nXF,Label='XF')
    do iChg=1,nChg
      Key = Get_Ln(LuSpool)
      call Get_F(1,XF(1,iChg),iShift)
      if (Convert) then
        XF(1:3,iChg) = XF(1:3,iChg)/Angstrom
        if (nOrd_ext /= 0) then
          XF(5:7,iChg) = XF(5:7,iChg)*Angstrom
        end if
      end if
    end do
    Convert = .false.
  else
    do iAt=1,natom
      Key = Get_Ln(LuSpool)
      call Get_I1(1,jAt)
      if ((jAt < 1) .or. (jAt > natom)) then
        write(6,'(A)') ' Error in espf/readin: atom out of range.'
        call Quit_OnUserError()
      end if
      call Get_F(2,Work(ipExt+(jAt-1)*MxExtPotComp),MxExtPotComp)
    end do
  end if
  close(LuSpool)
end if

! If nChg > 0, 2 possibilities:
! a) read external point charges (only, no dipoles) for a direct
!    QM/MM electrostatic coupling
! b) external potential calculated from point charges and dipoles

if ((nChg > 0) .and. DoDirect) then
  nXF = nChg
  nOrd_XF = nOrd_ext
  iXPolType = 0
  nXMolnr = 0
  call mma_deallocate(XF)
else if (nChg > 0) then
  do iAt=1,natom
    do iChg=1,nChg
      dx = Work(ipCord+(iAt-1)*3)-XF(1,iChg)
      dy = Work(ipCord+(iAt-1)*3+1)-XF(2,iChg)
      dz = Work(ipCord+(iAt-1)*3+2)-XF(3,iChg)
      qChg = XF(4,iChg)
      dpxChg = XF(5,iChg)
      dpyChg = XF(6,iChg)
      dpzChg = XF(7,iChg)
      rAtChg = sqrt(dx*dx+dy*dy+dz*dz)

      rAC2 = rAtChg*rAtChg
      rAC3 = rAtChg*rAC2
      rAC5 = rAC2*rAC3
      rAC7 = rAC2*rAC5
      iStart = ipExt+(iAt-1)*MxExtPotComp
      ! Potential E
      Work(iStart) = Work(iStart)+qChg/rAtChg-(dpxChg*dx+dpyChg*dy+dpzChg*dz)/rAC3
      ! Field F / x
      Work(iStart+1) = Work(iStart+1)-qChg*dx/rAC3+(dpxChg*(three*dx*dx-rAC2)+dpyChg*(three*dx*dy)+dpzChg*(three*dx*dz))/rAC5
      ! Field F / y
      Work(iStart+2) = Work(iStart+2)-qChg*dy/rAC3+(dpxChg*(three*dy*dx)+dpyChg*(three*dy*dy-rAC2)+dpzChg*(three*dy*dz))/rAC5
      ! Field F / z
      Work(iStart+3) = Work(iStart+3)-qChg*dz/rAC3+(dpxChg*(three*dz*dx)+dpyChg*(three*dz*dy)+dpzChg*(three*dz*dz-rAC2))/rAC5
      ! Gradient G / xx
      Work(iStart+4) = Work(iStart+4)+qChg*(three*dx*dx-rAC2)/rAC5-(dpxChg*(fift*dx*dx*dx-nine*dx*rAC2)+ &
                       dpyChg*(fift*dx*dx*dy-three*dy*rAC2)+dpzChg*(fift*dx*dx*dz-three*dz*rAC2))/rAC7
      ! Gradient G / yy
      Work(iStart+5) = Work(iStart+5)+qChg*(three*dy*dy-rAC2)/rAC5-(dpxChg*(fift*dy*dy*dx-three*dx*rAC2)+ &
                       dpyChg*(fift*dy*dy*dy-nine*dy*rAC2)+dpzChg*(fift*dy*dy*dz-three*dz*rAC2))/rAC7
      ! Gradient G / zz
      Work(iStart+6) = Work(iStart+6)+qChg*(three*dz*dz-rAC2)/rAC5-(dpxChg*(fift*dz*dz*dx-three*dx*rAC2)+ &
                       dpyChg*(fift*dz*dz*dy-three*dy*rAC2)+dpzChg*(fift*dz*dz*dz-nine*dz*rAC2))/rAC7
      ! Gradient G / xy
      Work(iStart+7) = Work(iStart+7)+qChg*(three*dx*dy)/rAC5-(dpxChg*(fift*dx*dy*dx-three*dx*rAC2)+ &
                       dpyChg*(fift*dx*dy*dy-three*dy*rAC2)+dpzChg*(fift*dx*dy*dz))/rAC7
      ! Gradient G / xz
      Work(iStart+8) = Work(iStart+8)+qChg*(three*dx*dz)/rAC5-(dpxChg*(fift*dx*dz*dx-three*dx*rAC2)+ &
                       dpyChg*(fift*dx*dz*dy)+dpzChg*(fift*dx*dz*dz-three*dz*rAC2))/rAC7
      ! Gradient G / yz
      Work(iStart+9) = Work(iStart+9)+qChg*(three*dy*dz)/rAC5-(dpxChg*(fift*dy*dz*dx)+ &
                       dpyChg*(fift*dy*dz*dy-three*dy*rAC2)+dpzChg*(fift*dy*dz*dz-three*dz*rAC2))/rAC7
    end do
  end do
  call mma_deallocate(XF)
end if

! Check the compatibility between old and new keywords

if (Exist) then
  if ((.not. Forces) .and. ((MltOrd /= MltOrd_old) .or. (iRMax /= iRMax_old) .or. (iGrdTyp /= iGrdTyp_old) .or. &
                            (lMorok .neqv. lMorok_old) .or. (DoDirect .neqv. DoDirect_old) .or. &
                            (abs(DeltaR-DeltaR_old) > 1.0d-6))) then
    write(6,*) 'Conficts between some old and new ESPF keywords'
    write(6,*) '      ','MltOrd     iRMax    DeltaR    iGrdTyp   lMorok   DoDirect'
    write(6,*) ' OLD: ',MltOrd_old,iRMax_old,DeltaR_old,iGrdTyp_old,lMorok_old,DoDirect_Old
    write(6,*) ' NEW: ',MltOrd,iRMax,DeltaR,iGrdTyp,lMorok,DoDirect
    write(6,'(A)') ' Check these values or erase ESPF.DATA'
  end if
else
  if (PotFile(1:(index(PotFile,' ')-1)) == '***') then
    write(6,*) 'Error! The EXTE data are missing'
    call Quit_OnUserError()
  end if
end if

! Some output

if (iPL >= 2) then
  if (DoDirect) then
    write(6,'(A)') ' DIRECT keyword found',' The ESPF scheme is switched off'
    write(6,'(A,I5,A)') ' External potential due to',nChg,' point charges'
  else
    if (nChg == 0) then
      write(6,'(A)') ' External potential:'
    else if (nChg > 0) then
      write(6,'(A,I5,A)') ' External potential due to',nChg,' point charges:'
    end if
    if (nChg >= 0) write(6,'(A)') ' Atom     E         Fx        Fy        Fz        Gxx       Gyy       Gzz       Gxy       '// &
                                  'Gxz       Gyz'
  end if
end if

! Write the external potential in ESPF.EXTPOT for later use

IPotFl = IsFreeUnit(IPotFl)
call Molcas_Open(IPotFl,'ESPF.EXTPOT')
write(IPotFl,'(I1)') 0
do iAt=1,natom
  if ((.not. DoDirect) .and. (nChg >= 0) .and. (iPL >= 2)) &
    write(6,ExtPotFormat) iAt,(Work(ipExt+(iAt-1)*MxExtPotComp+j),j=0,MxExtPotComp-1)
  write(IPotFl,ExtPotFormat) iAt,(Work(ipExt+(iAt-1)*MxExtPotComp+j),j=0,MxExtPotComp-1)
end do
close(IPotFl)
write(6,*)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine ReadIn_ESPF

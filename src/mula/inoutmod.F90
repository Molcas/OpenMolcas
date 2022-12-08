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
! Copyright (C) 1995, Niclas Forsberg                                  *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

!module InOutMod

!  Contains:
!    ReadInp
!    WriteLog
!    WriteDip
!    IntCalcHeader
!    ExpPointHeader
!    ISCHeader
!    WriteHeader
!    WrMold
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!contains

subroutine ReadInp(Title,AtomLbl,Mass,InterVec,Bond,nBond,NumInt,NumOfAt,trfName1,trfName2,m_max,n_max,max_dip,max_term,MatEl, &
                   ForceField,Cartesian,lExpan,lISC,iCode,dMinWind)
!  Purpose:
!    Read the input file.
!
!  Output:
!    Title    : String - title of the project.
!    AtomLbl  : Array of character - contains the labels for the atoms.
!    Mass     : Real array - contains the mass of the atoms.
!    InterVec : Integer array - containis the atoms that are used in the calculations of each internal coordinate.
!    Bond     : Integer array - contains atom pairs that are to be bonded together in a plot.
!    nBond    : Integer - dim of Bond, i.e. 2*(number of bonds).
!    NumInt   : Integer - the total number of internal coordinates.
!    NumOfAt  : Integer - the number of atoms.
!    trfName  : Character array - type of transformation of variables.
!    m_max    : Integer - maximum level for the first state.
!    n_max    : Integer - maximum level for the second state.
!    m_plot   : Integer array - level(s) to plot for the first state.
!    n_plot   : Integer array - level(s) to plot for the second state.
!    max_dip  : Integer - highest order of term in transition dipole.
!    max_term : Integer - highest power of a term in polynomial fitted to energy values.
!    MatEl    : Logical
!    lISC     : Logical to calculate InterSystem Crossing
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.
!
!  Modified to use isotopes module
!    Ignacio Fdez. Galvan, 2017

use mula_global, only: AtCoord1, AtCoord2, broadplot, cmstart, cmend, energy1, energy2, hbarcm, Hess1, Hess2, Huge_Print, inpUnit, &
                       ipow, LifeTime, m_plot, MaxNumAt, n_plot, ndata, NormModes, nPolyTerm, nvar, OscStr, plotwindow, t_dipin1, &
                       t_dipin2, TranDip, TranDipGrad, Use_cm, Use_nm, var, VibModPlot, WriteVibLevels, yin1, yin2
use Isotopes, only: Isotope
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight, Quart, Angstrom, auTocm, auToeV, deg2rad, uToau
use Definitions, only: wp, iwp, u6

implicit none
character(len=80), intent(out) :: Title, trfName1(MaxNumAt), trfName2(MaxNumAt)
character(len=4), intent(out) :: AtomLbl(MaxNumAt)
real(kind=wp), intent(out) :: Mass(MaxNumAt), dMinWind
integer(kind=iwp), intent(out) :: InterVec(MaxNumAt*15), Bond(MaxNumAt*2), nBond, NumInt, NumOfAt, m_max, n_max, max_dip, max_term
logical(kind=iwp), intent(out) :: MatEl, ForceField, Cartesian, lExpan, lISC
integer(kind=iwp), intent(inout) :: iCode
integer(kind=iwp), parameter :: max_len_xvec = 150, max_plot_temp = 100
integer(kind=iwp) :: i, iAtom, icoord, idata, iend, ierror, ii, iInt, inpUnit1, inpUnit2, IntVal, iPrint, iStart, istatus, iStop, &
                     iterm, j, jvar, k, l, l_a, l_n, Length, m, modeTemp, n, nAtom, nFcart, nsum, Num, plot_temp(max_plot_temp)
                     !, AtData(120,4)
real(kind=wp) :: const, coord(1000,10), error, FWHM, GrdVal(1000,10), m1, m12, m2, m23, m3, r1, r2, rfact, sign1_1, sign1_2, &
                 sign2_1, sign2_2, theta, TmpCoord(3), xvec(max_len_xvec) !, MassData(400)
logical(kind=iwp) :: exists, XnotFound, YnotFound, ZnotFound
character(len=80) :: InLine, OutLine
character(len=9) :: CoordType
character(len=7) :: AtInp
character(len=6) :: PlUnit
character(len=4) :: Atom, Atom1, Atom2, Atom3, Atom4
character(len=3) :: DegOrRad
character(len=2) :: AAOrAu, AtName, Eunit !, AtList(0:120)
character :: c
integer(kind=iwp), allocatable :: GeoVec(:), tempmodes(:)
real(kind=wp), allocatable :: Fcart(:,:), ScaleParam(:), Sinv(:,:), SS(:,:), Temp(:,:)
character(len=*), parameter :: frmt = '(A80)'
! User-defined functions called:
integer(kind=iwp), external :: iPrintLevel, isfreeunit, iStrToInt
real(kind=wp), external :: StrToDble

! Read from MassFile:
! - name of atoms into array of string - AtList.
! - atomic number, default mass number, smallest mass number and
!   table offset into two dimensional array - AtData.
! - masses of all isotopes into array - MassData.

!i = 0
!call Molcas_Open(massUnit,'MASSUNIT')
!read(massUnit,frmt) InLine
! Read sequential lines from atomic data file.
!call Normalize(InLine,OutLine)
!do while (OutLine(1:4) /= 'END ')
!  l = Index(OutLine,'*')
!  if (l == 0) then
!    i = i+1
!    AtList(i) = OutLine(1:2)
!    read(OutLine(3:80),*) (AtData(i,j),j=1,4)
!  end if
!  read(massUnit,frmt) InLine
!  call Normalize(InLine,OutLine)
!end do
!iAtList = i

!i = 0
!read(massUnit,frmt) InLine
!call Normalize(InLine,OutLine)
!do while (OutLine(1:4) /= 'END ')
!  l = Index(OutLine,'*')
!  if (l == 0) then
!    i = i+1
!    read(OutLine,*) MassData(i)
!  end if
!  read(massUnit,frmt) InLine
!  call Normalize(InLine,OutLine)
!end do
!close(massUnit)
iPrint = iPrintLevel(-1)

! ----------------------------------------------------------------------
! TITLe: Read the title into a string - Title

Title = ' '
call KeyWord(inpUnit,'TITL',.true.,exists)
if (exists) read(inpUnit,frmt) Title

! TITLe: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! ATOMs:  Process ATOMS

call KeyWord(inpUnit,'ATOM',.true.,exists)
NumOfAt = 0
read(inpUnit,frmt) InLine
call Normalize(InLine,OutLine)
do while (OutLine(1:4) /= 'END ')
  NumOfAt = NumOfAt+1
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
end do

! Read the labels of atoms into an array - AtomLbl.

call KeyWord(inpUnit,'ATOM',.true.,exists)
do nAtom=1,NumOfAt
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
  k = 1
  call WordPos(k,OutLine,iStart,iStop)
  k = iStop+1
  AtInp = ' '
  AtInp = OutLine(iStart:iStop)
  !- Check if a massnumber is given.
  Num = 0
  do i=1,7
    c = AtInp(i:i)
    IntVal = index('0123456789',c)-1
    if (IntVal >= 0) then
      Num = 10*Num+IntVal
    else
      exit
    end if
  end do
  if ((iStart+i-1) > iStop) then
    call WordPos(k,OutLine,iStart,iStop)
    AtInp = ' '
    AtInp = OutLine(iStart:iStop)
    i = 1
    c = AtInp(i:i)
  end if
  ! Extract atomic name and possible label.
  AtName = ' '
  AtomLbl(nAtom) = ' '
  j = 1
  do ii=i,7
    IntVal = index('0123456789',c)-1
    if ((IntVal < 0) .and. (c /= ' ')) then
      AtName(j:j) = c
      AtomLbl(nAtom)(j:j) = c
      j = j+1
      i = i+1
      c = AtInp(i:i)
    end if
  end do
  do ii=i,7
    IntVal = index('0123456789',c)-1
    if (((iStart+i-1) <= iStop) .and. (IntVal >= 0)) then
      AtomLbl(nAtom)(j:j) = c
      j = j+1
      i = i+1
      c = AtInp(i:i)
    end if
  end do
  ! Remove massnumber, atom and label from OutLine.
  Length = len(OutLine)
  OutLine = OutLine(iStop+1:Length)

  ! If mass is given in input, then read mass, else use MassData.
  k = 1
  call WordPos(k,OutLine,iStart,iStop)
  k = iStop
  if (k /= Length) then
    read(OutLine,*) Mass(nAtom)
  else
    !do i=1,iAtList
    !  if (AtList(i) == AtName ) exit
    !end do
    !if (Num == 0) then
    !  Num = AtData(i,2)
    !end if
    !Mass(nAtom) = MassData(AtData(i,4)+Num+1-AtData(i,3))
    call Isotope(Num,AtName,Mass(nAtom))
    Mass(nAtom) = Mass(nAtom)/uToau
    ! Print error message if the isotope is not listed in MassFile.
    !nOffset = AtData(i+1,4)
    !if ((Num > nOffset) .or. (Mass(nAtom) == -One)) then
    !  write(u6,*)
    !  write(u6,*) ' *************** ERROR *****************'
    !  write(u6,'(A,A2,A,I3,A)') ' The isotope of ',AtName,' with mass number ',Num
    !  write(u6,*) ' is not listed in MassFile.'
    !  write(u6,*) '****************************************'
    !  call Quit_OnUserError()
    !end if
  end if
end do

! ATOMs: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! INTErnal: Resolve the different internal coordinates specified in the input.

call KeyWord(inpUnit,'INTE',.true.,exists)
j = 1
NumInt = 0
nBond = 0
read(inpUnit,frmt) InLine
call Normalize(InLine,OutLine)

do while (OutLine(1:4) /= 'END ')

  ! Bond Stretching.
  l = index(OutLine,'BOND')
  if (l /= 0) then
    k = l+len('BOND')
    call WordPos(k,OutLine,iStart,iStop)
    Atom1 = OutLine(iStart:iStop)
    k = iStop+1
    call WordPos(k,OutLine,iStart,iStop)
    Atom2 = OutLine(iStart:iStop)
    InterVec(j) = 1
    do i=1,NumOfAt
      if (Atom1 == AtomLbl(i)) then
        InterVec(j+1) = i
      else if (Atom2 == AtomLbl(i)) then
        InterVec(j+2) = i
      end if
    end do
    Bond(nBond+1) = InterVec(j+1)
    Bond(nBond+2) = InterVec(j+2)
    nBond = nBond+2
    j = j+3
    NumInt = NumInt+1
  end if
  ! Valence Angle Bending.
  l = index(OutLine,'ANGLE')
  if (l /= 0) then
    k = l+len('ANGLE')
    call WordPos(k,OutLine,iStart,iStop)
    Atom1 = OutLine(iStart:iStop)
    k = iStop+1
    call WordPos(k,OutLine,iStart,iStop)
    Atom2 = OutLine(iStart:iStop)
    k = iStop+1
    call WordPos(k,OutLine,iStart,iStop)
    Atom3 = OutLine(iStart:iStop)
    InterVec(j) = 2
    do i=1,NumOfAt
      if (Atom1 == AtomLbl(i)) then
        InterVec(j+1) = i
      else if (Atom2 == AtomLbl(i)) then
        InterVec(j+2) = i
      else if (Atom3 == AtomLbl(i)) then
        InterVec(j+3) = i
      end if
    end do
    j = j+4
    NumInt = NumInt+1
  end if
  ! Linear Valence Angle.
  l = index(OutLine,'LINANG')
  if (l /= 0) then
    k = l+len('LINANG')
    call WordPos(k,OutLine,iStart,iStop)
    Atom1 = OutLine(iStart:iStop)
    k = iStop+1
    call WordPos(k,OutLine,iStart,iStop)
    Atom2 = OutLine(iStart:iStop)
    k = iStop+1
    call WordPos(k,OutLine,iStart,iStop)
    Atom3 = OutLine(iStart:iStop)
    InterVec(j) = 3
    do i=1,NumOfAt
      if (Atom1 == AtomLbl(i)) then
        InterVec(j+1) = i
      else if (Atom2 == AtomLbl(i)) then
        InterVec(j+2) = i
      else if (Atom3 == AtomLbl(i)) then
        InterVec(j+3) = i
      end if
    end do
    j = j+4
    NumInt = NumInt+2
  end if
  ! Torsion.
  l = index(OutLine,'TORSION')
  if (l /= 0) then
    k = l+len('TORSION')
    call WordPos(k,OutLine,iStart,iStop)
    Atom1 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom2 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom3 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom4 = OutLine(iStart:iStop)
    InterVec(j) = 4
    do i=1,NumOfAt
      if (Atom1 == AtomLbl(i)) then
        InterVec(j+1) = i
      else if (Atom2 == AtomLbl(i)) then
        InterVec(j+2) = i
      else if (Atom3 == AtomLbl(i)) then
        InterVec(j+3) = i
      else if (Atom4 == AtomLbl(i)) then
        InterVec(j+4) = i
      end if
    end do
    j = j+5
    NumInt = NumInt+1
  end if
  ! Out of Plane Angle.
  l = index(OutLine,'OUTOFPL')
  if (l /= 0) then
    k = l+len('OUTOFPL')
    call WordPos(k,OutLine,iStart,iStop)
    Atom1 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom2 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom3 = OutLine(iStart:iStop)
    call WordPos(k,OutLine,iStart,iStop)
    Atom4 = OutLine(iStart:iStop)
    InterVec(j) = 5
    do i=1,NumOfAt
      if (Atom1 == AtomLbl(i)) then
        InterVec(j+1) = i
      else if (Atom2 == AtomLbl(i)) then
        InterVec(j+2) = i
      else if (Atom3 == AtomLbl(i)) then
        InterVec(j+3) = i
      else if (Atom4 == AtomLbl(i)) then
        InterVec(j+4) = i
      end if
    end do
    j = j+5
    NumInt = NumInt+1
  end if

  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
end do

call mma_allocate(AtCoord1,3,NumOfAt,label='AtCoord1')
call mma_allocate(AtCoord2,3,NumOfAt,label='AtCoord2')
call mma_allocate(Hess1,NumInt,NumInt,label='Hess1')
call mma_allocate(Hess2,NumInt,NumInt,label='Hess2')
n = 3*NumOfAt
call mma_allocate(TranDipGrad,3,n,label='TranDipGrad')

! INTErnal: End --------------------------------------------------------

! ----------------------------------------------------------------------
! ISC: InterSystem Crossing. GG-30-Dec-08

lISC = .false.
dMinWind = Zero
call KeyWord(inpUnit,'ISC ',.true.,lISC)
if (lISC) then
  call KeyWord(inpUnit,'EXPF',.true.,exists)
  if (exists) then
    read(InpUnit,frmt) InLine
    call Normalize(InLine,Outline)
    read(OutLine,*) dMinWind
  end if
  call KeyWord(inpUnit,'OLDC',.true.,exists)
  if (exists) iCode = iCode+1
  call KeyWord(inpUnit,'DISK',.true.,exists)
  if (exists) iCode = iCode+10
end if

! ISC: End -------------------------------------------------------------

! ----------------------------------------------------------------------
! MODEs: Chose which modes to use in intensity calculations.

call KeyWord(inpUnit,'MODE',.true.,exists)
if (exists) then
  call mma_allocate(tempmodes,NumInt,label='tempmodes')
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
  i = 1
  do while (OutLine(1:4) /= 'END ')
    k = 1
    call WordPos(k,OutLine,iStart,iStop)
    do while (iStop < len(OutLine))
      tempModes(i) = iStrToInt(OutLine(iStart:iStop))
      if (tempModes(i) > NumInt) then
        write(u6,*)
        write(u6,*) ' *********** ERROR *************'
        write(u6,*) ' Too many normal modes specified'
        write(u6,*) ' *******************************'
        call Quit_OnUserError()
      end if
      i = i+1
      call WordPos(k,OutLine,iStart,iStop)
    end do
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
  end do
  n = i-1
  call mma_allocate(NormModes,n,label='NormModes')
  NormModes(:) = tempModes(1:n)
  call mma_deallocate(tempmodes)
  do i=1,n-1
    do j=i+1,n
      if (NormModes(j) < NormModes(i)) then
        modeTemp = NormModes(i)
        NormModes(i) = NormModes(j)
        NormModes(j) = modeTemp
      end if
    end do
  end do
else
  ! No MODEs specified
  call mma_allocate(NormModes,NumInt,label='NormModes')
  do i=1,NumInt
    NormModes(i) = i
  end do
  if (iPrint >= 1) then
    write(u6,*) ' *** Warning: MODEs not specified !'
    write(u6,*) '     All ',NumInt,' will be calculated.'
    write(u6,*)
  end if
end if

! MODEs: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! MXLEvels: Maximun number of vibr. quanta in first and second state

call KeyWord(inpUnit,'MXLE',.true.,exists)
m_max = 0 ! Default
n_max = 1 ! Default
if (.not. exists) then
  if (iPrint >= 1) then
    write(u6,*) ' *** Warning: MXLEvels not specified !'
    write(u6,*) '     Default values are 0 and 1.'
    write(u6,*)
  end if
else
  read(inpUnit,*) m_max,n_max
end if

! MXLEvels: End --------------------------------------------------------

! ----------------------------------------------------------------------
! VARIational: Check if a full matrix element calculation or a
!              simpler harmonic approximation is wanted.

MatEl = .false.
call KeyWord(inpUnit,'VARI',.true.,MatEl)

! VARIational: End -----------------------------------------------------

! ----------------------------------------------------------------------
! TRANsitions: Check which transitions that are wanted in the output.

call KeyWord(inpUnit,'TRAN',.true.,exists)
if (.not. exists) then
  if (.not. lISC) then
    write(u6,*) ' *** Warning: TRANsitions not specified!'
    write(u6,*) '     All Levels will be printed.'
    write(u6,*)
  end if
  n = m_max+1
  do i=1,n
    plot_temp(i) = i-1
  end do
else
  call KeyWord(inpUnit,'FIRS',.false.,exists)
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
  k = 1
  i = 1
  call WordPos(k,OutLine,iStart,iStop)
  do while (iStop < len(OutLine))
    if (i > max_plot_temp) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' TRANsitions: i  >  max_plot_temp'
      write(u6,*) ' ********************************'
      call Quit_OnUserError()
    end if
    plot_temp(i) = iStrToInt(OutLine(iStart:iStop))
    if (plot_temp(i) > m_max) then
      write(u6,*)
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' A quantum specified for the output'
      write(u6,*) ' is larger than max quantum.'
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    i = i+1
    call WordPos(k,OutLine,iStart,iStop)
  end do
  n = i-1
end if
call mma_allocate(m_plot,n,label='m_plot')
m_plot(:) = plot_temp(1:n)

call KeyWord(inpUnit,'TRAN',.true.,exists)
if (.not. exists) then
  n = n_max+1
  do i=1,n
    plot_temp(i) = i-1
  end do
else
  call KeyWord(inpUnit,'SECO',.false.,exists)
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
  k = 1
  i = 1
  call WordPos(k,OutLine,iStart,iStop)
  do while (iStop < len(OutLine))
    if (i > max_plot_temp) then
      write(u6,*)
      write(u6,*) ' *********** ERROR ************'
      write(u6,*) ' TRANsitions: i > max_plot_temp'
      write(u6,*) ' ******************************'
      call Quit_OnUserError()
    end if
    plot_temp(i) = iStrToInt(OutLine(iStart:iStop))
    if (plot_temp(i) > n_max) then
      write(u6,*)
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' A quantum specified for the output'
      write(u6,*) ' is larger than max quantum.'
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    i = i+1
    call WordPos(k,OutLine,iStart,iStop)
  end do
  n = i-1
end if
call mma_allocate(n_plot,n,label='n_plot')
n_plot(:) = plot_temp(1:n)

! TRANsitions: End -----------------------------------------------------

call KeyWord(inpUnit,'FORC',.true.,Forcefield)
if (Forcefield) then

  !--------------------------------------------------------------------!
  !
  !                              Forcefield part
  !
  !--------------------------------------------------------------------!

  ! --------------------------------------------------------------------
  ! ENERgy: Read minimum energies for the two surfaces.

  call KeyWord(inpUnit,'ENER',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ************** ERROR *************'
    write(u6,*) ' Electronic T_e energies not given.'
    write(u6,*) ' **********************************'
    call Quit_OnUserError()
  end if
  call KeyWord(inpUnit,'FIRS',.false.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ******************** ERROR *********************'
    write(u6,*) ' Electronic T_e energy for first state not given.'
    write(u6,*) ' ************************************************'
    call Quit_OnUserError()
  end if
  !read(inpUnit,*) energy1,Eunit
  !call Upcase(Eunit)
  ! Replace above two lines with:
  read(InpUnit,frmt) InLine
  call Normalize(InLine,Outline)
  read(OutLine,*) Energy1
  EUnit = 'AU'
  if (index(OutLine,'EV') > 0) EUnit = 'EV'
  if (index(OutLine,'CM') > 0) EUnit = 'CM'
  ! End of replacement.

  rfact = One
  if (Eunit == 'AU') rfact = One
  if (Eunit == 'EV') rfact = One/auToeV
  if (Eunit == 'CM') rfact = One/auTocm
  energy1 = energy1*rfact
  call KeyWord(inpUnit,'SECO',.false.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ******************** ERROR **********************'
    write(u6,*) ' Electronic T_e energy for second state not given.'
    write(u6,*) ' *************************************************'
    call Quit_OnUserError()
  end if
  !read(inpUnit,*) energy2,Eunit
  !call Upcase(Eunit)
  ! Replace above two lines with:
  read(InpUnit,frmt) InLine
  call Normalize(InLine,Outline)
  read(OutLine,*) Energy2
  EUnit = 'AU'
  if (index(OutLine,'EV') > 0) EUnit = 'EV'
  if (index(OutLine,'CM') > 0) EUnit = 'CM'
  ! End of replacement.
  if (Eunit == 'AU') rfact = One
  if (Eunit == 'EV') rfact = One/auToeV
  if (Eunit == 'CM') rfact = One/auTocm
  energy2 = energy2*rfact

  ! ENERgy: : End ------------------------------------------------------

  ! --------------------------------------------------------------------
  ! GEOMetries: Read geometries for the two surfaces.

  call KeyWord(inpUnit,'GEOM',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ***** ERROR *******'
    write(u6,*) ' GEOMetry not found.'
    write(u6,*) ' *******************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read geometry ("GEOMETRY").'
  !read(inpUnit,*) CoordType
  !call UpCase(CoordType)
  ! Replace above two lines with:
  read(InpUnit,frmt) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  ! End of replacement.

  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit1 = isfreeunit(31)
    inpUnit2 = isfreeunit(32)
    call molcas_open(inpUnit1,'UNSYM1')
    !open(31,'UNSYM1')
    call KeyWord(inpUnit1,'*LABEL COORDINATES CHARGE',.false.,exists)
    call molcas_open(inpUnit2,'UNSYM2')
    !open(32,'UNSYM2')
    call KeyWord(inpUnit2,'*LABEL COORDINATES CHARGE',.false.,exists)
  else
    inpUnit1 = inpUnit
    inpUnit2 = inpUnit
  end if

  if (CoordType == 'CARTESIAN') then
    Cartesian = .true.
    !D write(u6,*) ' READINP: Read cartesian coordinates for first state.'
    call KeyWord(inpUnit,'FIRS',.false.,exists)
    if (.not. exists) then
      write(u6,*)
      write(u6,*) ' ****** ERROR ********'
      write(u6,*) ' Wrong geometry input.'
      write(u6,*) ' *********************'
      call Quit_OnUserError()
    end if
    ierror = 0
    do iAtom=1,NumOfAt
      read(inpUnit1,frmt) InLine
      call UpCase(InLine)
      Atom = InLine(1:4)
      read(InLine(5:80),*) (TmpCoord(i),i=1,3)
      j = 10000000
      do i=1,NumOfAt
        if (Atom == AtomLbl(i)) j = i
      end do
      if (j > NumOfAt) then
        write(u6,*)
        write(u6,*) '*********** ERROR ***********'
        write(u6,*) ' Unrecognized atom label "'//Atom//'"'
        write(u6,*) '*****************************'
        ierror = ierror+1
      else
        AtCoord1(:,j) = TmpCoord(:)
      end if
    end do
    !D write(u6,*) ' READINP: Read cartesian coordinates for second state.'
    call KeyWord(inpUnit,'SECO',.false.,exists)
    if (.not. exists) then
      write(u6,*)
      write(u6,*) ' ****** ERROR ********'
      write(u6,*) ' Wrong geometry input.'
      write(u6,*) ' *********************'
      call Quit_OnUserError()
    end if
    do iAtom=1,NumOfAt
      read(inpUnit2,frmt) InLine
      call UpCase(InLine)
      Atom = InLine(1:4)
      read(InLine(5:80),*) (TmpCoord(i),i=1,3)
      j = 1
      do while (Atom /= AtomLbl(j))
        j = j+1
      end do
      if (j > NumOfAt) then
        write(u6,*)
        write(u6,*) '*********** ERROR ***********'
        write(u6,*) ' Unrecognized atom label "'//Atom//'"'
        write(u6,*) '*****************************'
        ierror = ierror+1
      else
        AtCoord2(:,j) = TmpCoord(:)
      end if
    end do

  else if (CoordType == 'INTERNAL') then

    Cartesian = .false.
    call KeyWord(inpUnit,'FIRS',.false.,exists)
    !D write(u6,*) ' READINP: Read internal coordinates for first state.'
    n = 5*(3*NumOfAt-5)
    ! Largest possible necessary GeoVec
    call mma_allocate(GeoVec,n,label='GeoVec')
    iInt = 1
    j = 1
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
    do while (OutLine(1:4) /= 'END ')
      ! Bond distance.
      l = index(OutLine,'BOND')
      if (l /= 0) then
        k = l+len('BOND')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        GeoVec(j) = 1
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        AAOrAu = OutLine(iStart:iStop)
        if (AAOrAu == 'AA') then
          xvec(iInt) = xvec(iInt)/Angstrom
        end if
        j = j+3
      end if
      ! Valence Angle.
      l = index(OutLine,'ANGLE')
      if (l /= 0) then
        k = l+len('ANGLE')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        GeoVec(j) = 2
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+4
      end if
      ! Linear Valence Angle.
      l = index(OutLine,'LINANG')
      if (l /= 0) then
        k = l+len('LINANG')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        GeoVec(j) = 3
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+4
      end if
      ! Dihedral angle.
      l = index(OutLine,'TORSION')
      if (l /= 0) then
        k = l+len('TORSION')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom4 = OutLine(iStart:iStop)
        GeoVec(j) = 4
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          else if (Atom4 == AtomLbl(i)) then
            GeoVec(j+4) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+5
      end if
      ! Out of Plane Angle.
      l = index(OutLine,'OUTOFPL')
      if (l /= 0) then
        k = l+len('OUTOFPL')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom4 = OutLine(iStart:iStop)
        GeoVec(j) = 5
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          else if (Atom4 == AtomLbl(i)) then
            GeoVec(j+4) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+5
      end if
      iInt = iInt+1
      read(inpUnit,frmt) InLine
      call Normalize(InLine,OutLine)
    end do
    iInt = iInt-1
    !call Int_to_Cart(GeoVec,xvec,AtCoord1,NumOfAt,iInt,Mass)
    l_a = NumOfAt
    l_n = max_len_xvec
    call Int_to_Cart1(GeoVec,xvec,AtCoord1,l_a,l_n)

    call mma_deallocate(GeoVec)

    call KeyWord(inpUnit,'SECO',.false.,exists)
    !D write(u6,*) ' READINP: Read internal coordinates for second state.'
    n = 5*(3*NumOfAt-5)
    call mma_allocate(GeoVec,n,label='GeoVec')
    iInt = 1
    j = 1
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
    do while (OutLine(1:4) /= 'END ')
      ! Bond distance.
      l = index(OutLine,'BOND')
      if (l /= 0) then
        k = l+len('BOND')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        GeoVec(j) = 1
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        AAOrAu = OutLine(iStart:iStop)
        if (AAOrAu == 'AA') then
          xvec(iInt) = xvec(iInt)/Angstrom
        end if
        j = j+3
      end if
      ! Valence Angle.
      l = index(OutLine,'ANGLE')
      if (l /= 0) then
        k = l+len('ANGLE')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        GeoVec(j) = 2
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+4
      end if
      ! Linear Valence Angle.
      l = index(OutLine,'LINANG')
      if (l /= 0) then
        k = l+len('LINANG')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        GeoVec(j) = 3
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+4
      end if
      ! Dihedral angle.
      l = index(OutLine,'TORSION')
      if (l /= 0) then
        k = l+len('TORSION')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom4 = OutLine(iStart:iStop)
        GeoVec(j) = 4
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          else if (Atom4 == AtomLbl(i)) then
            GeoVec(j+4) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+5
      end if
      ! Out of Plane Angle.
      l = index(OutLine,'OUTOFPL')
      if (l /= 0) then
        k = l+len('OUTOFPL')
        call WordPos(k,OutLine,iStart,iStop)
        Atom1 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom2 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom3 = OutLine(iStart:iStop)
        call WordPos(k,OutLine,iStart,iStop)
        Atom4 = OutLine(iStart:iStop)
        GeoVec(j) = 5
        do i=1,NumOfAt
          if (Atom1 == AtomLbl(i)) then
            GeoVec(j+1) = i
          else if (Atom2 == AtomLbl(i)) then
            GeoVec(j+2) = i
          else if (Atom3 == AtomLbl(i)) then
            GeoVec(j+3) = i
          else if (Atom4 == AtomLbl(i)) then
            GeoVec(j+4) = i
          end if
        end do
        call WordPos(k,OutLine,iStart,iStop)
        xvec(iInt) = StrToDble(OutLine(iStart:iStop))
        k = iStop+1
        call WordPos(k,OutLine,iStart,iStop)
        DegOrRad = OutLine(iStart:iStop)
        if (DegOrRad == 'DEG') then
          xvec(iInt) = xvec(iInt)*deg2rad
        end if
        j = j+5
      end if
      iInt = iInt+1
      read(inpUnit,frmt) InLine
      call Normalize(InLine,OutLine)
    end do
    iInt = iInt-1
    !call Int_to_Cart(GeoVec,xvec,AtCoord2,NumOfAt,iInt,Mass)
    l_a = NumOfAt
    l_n = max_len_xvec
    call Int_to_Cart1(GeoVec,xvec,AtCoord2,l_a,l_n)
    call mma_deallocate(GeoVec)
  end if

  if (inpUnit1 /= inpUnit) close(inpUnit1)
  if (inpUnit2 /= inpUnit) close(inpUnit2)

  ! GEOMetries: End ----------------------------------------------------

  ! --------------------------------------------------------------------
  ! MXORder: Read maximum order for the transition dipole moment.

  call KeyWord(inpUnit,'MXOR',.true.,exists)
  if (.not. exists) then
    if (.not. lISC) then
      write(u6,*) ' *** Warning: MXORder not specified !'
      write(u6,*) '     Transition Dipole is assumed as constant.'
      write(u6,*)
    end if
    max_dip = 0
  else
    read(inpUnit,*) max_dip
  end if

  ! MXORder: End -------------------------------------------------------

  ! --------------------------------------------------------------------
  ! OSCStr: Determines if you get the oscillator strength
  ! or the transition intensity in the output.

  OscStr = .false.
  if (.not. lISC) then
    call KeyWord(inpUnit,'OSCS',.true.,OscStr)
    if (OscStr) then
      write(u6,*) '     The Oscillator Strength will be calculated.'
      write(u6,*)
    else
      write(u6,*) '     The Intensity will be calculated.'
      write(u6,*)
    end if
  end if

  ! OSCStr: End --------------------------------------------------------

  ! --------------------------------------------------------------------
  ! NANOmeters: Generate plot file in nanometers.

  Use_nm = .false.
  if (.not. lISC) call KeyWord(inpUnit,'NANO',.true.,Use_nm)

  ! NANOmeters: End ----------------------------------------------------

  ! --------------------------------------------------------------------
  ! CM-1: Generate plot file in cm-1.

  Use_cm = .false.
  if (.not. lISC) call KeyWord(inpUnit,'CM-1',.true.,Use_cm)

  ! CM-1: End ----------------------------------------------------------

  ! --------------------------------------------------------------------
  ! PLOT: Enter the limits (in eV, cm-1 or in nm) for the plot file.

  cmstart = Zero
  cmend = Zero
  plotwindow = .false.
  if (.not. lISC) then
    call KeyWord(inpUnit,'PLOT',.true.,plotwindow)
    if (plotwindow) then
      read(inpUnit,*) cmstart,cmend
    end if
  end if

  ! PLOT : End ---------------------------------------------------------

  PlUnit = ' '
  broadplot = .false.
  LifeTime = 130.0e-15_wp
  if (.not. lISC) then
    if (Use_nm .and. Use_cm) then
      write(u6,*)
      write(u6,*) ' ***************** ERROR *******************'
      write(u6,*) ' NANOmeters and CM-1 are mutually exclusive.'
      write(u6,*) ' *******************************************'
      call Quit_OnUserError()
    end if
    if (Use_nm) then
      PlUnit = ' nm.'
      write(u6,*) '     The plot file will be in nanometers and the'
    else if (Use_cm) then
      PlUnit = ' cm-1.'
      write(u6,*) '     The plot file will be in cm-1 and the'
    else
      PlUnit = ' eV.'
      write(u6,*) '     The plot file will be in eV and the'
    end if
    if (plotwindow) then
      write(u6,*) '     interval is from ',cmstart,' to ',cmend,PlUnit
    else
      write(u6,*) '     interval automatically defined.'
    end if
    write(u6,*)

    ! --------------------------------------------------------------------
    ! BROAplot: Gives the peaks in the spectrum an artificial halfwidth.

    call KeyWord(inpUnit,'BROA',.true.,broadplot)
    call KeyWord(inpUnit,'LIFE',.true.,exists)
    if (exists) read(inpUnit,*) LifeTime
    write(u6,*) '     Life time : ',LifeTime,' sec'
    FWHM = hbarcm/(Two*LifeTime)
    write(u6,'(1X,A,F8.1,A,F9.6,A)') '     FWHM : ',FWHM,' cm-1 /',FWHM*auToeV/auTocm,' eV'

    ! BROAplot: End ------------------------------------------------------
  end if

  ! --------------------------------------------------------------------
  ! VIBWrite: Writes vibrational levels to logfile.

  WriteVibLevels = .false.
  call KeyWord(inpUnit,'VIBW',.true.,WriteVibLevels)

  ! VIBWrite: End ------------------------------------------------------

  ! --------------------------------------------------------------------
  ! VibPlot: Check for keyword VibPlot.

  VibModPlot = .false.
  call KeyWord(inpUnit,'VIBP',.true.,VibModPlot)
  if (VibModPlot) then
    call KeyWord(inpUnit,'CYCL',.false.,exists)
    if (exists) then
      call KeyWord(inpUnit,'VIBP',.true.,exists)
      read(inpUnit,frmt) InLine
      call Normalize(InLine,OutLine)
      k = 1
      call WordPos(k,OutLine,iStart,iStop)
      call WordPos(k,OutLine,iStart,iStop)
      Bond(nBond+1) = iStrToInt(OutLine(iStart:iStop))
      call WordPos(k,OutLine,iStart,iStop)
      Bond(nBond+2) = iStrToInt(OutLine(iStart:iStop))
      nBond = nBond+2
    end if
  end if

  ! VibPlot: End -------------------------------------------------------

  ! --------------------------------------------------------------------
  ! HUGElog: If exists, then the logfile will be more detailed.

  Huge_Print = .false.
  call KeyWord(inpUnit,'HUGE',.true.,Huge_Print)

  ! HUGElog: End -------------------------------------------------------

  ! --------------------------------------------------------------------
  ! EXPAnsion: Program will be aborted after expansion point geometry
  !            is calculated.

  lExpan = .false.
  call KeyWord(inpUnit,'EXPA',.true.,lExpan)

  ! HUGElog: End -------------------------------------------------------

  ! --------------------------------------------------------------------
  ! FORCe: Here really read force constants.

  call KeyWord(inpUnit,'FORC',.true.,exists)
  max_term = 2
  call KeyWord(inpUnit,'FIRS',.false.,exists)
  !D write(u6,*) ' READINP: Read force constants for first state.'
  read(InpUnit,frmt) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  !D write(u6,*) ' READINP: CoordType=',CoordType

  ! Read Hessian. If Hessian is given in cartesian coordinates,
  ! first remove total translation and total rotation and then
  ! transform it to internal coordinates.
  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit1 = isfreeunit(31)
    call molcas_open(inpUnit1,'UNSYM1')
    !open(31,'UNSYM1')
    call KeyWord(inpUnit1,'UNSYMMETRIZED HESSIAN',.false.,exists)
    read(inpUnit1,'(a17)') Inline
    read(inpUnit1,'(a17)') Inline
  else
    inpUnit1 = inpUnit
  end if

  if (CoordType == 'INTERNAL') then
    do i=1,NumInt
      read(inpUnit,*) Hess1(i,:)
    end do
  else if (CoordType == 'CARTESIAN') then
    call mma_allocate(SS,3*NumOfAt,NumInt,label='SS')
    SS(:,:) = Zero
    call CalcS(AtCoord1,InterVec,SS,NumInt,NumOfAt)

    ! Invert S matrix and remove total translation and total rotation.
    call mma_allocate(Sinv,3*NumOfAt,NumInt,label='Sinv')

    call RotTranRem(Sinv,SS,Mass,AtCoord1,NumOfAt,NumInt)
    call mma_deallocate(SS)

    ! Read cartesian force constant matrix.
    n = 3*NumOfAt
    nFcart = n
    call mma_allocate(Fcart,nFcart,nFcart,label='Fcart')

    do m=1,3*NumOfAt
      read(inpUnit1,'(a17)',iostat=istatus) Inline
      if (istatus <= 0) read(inpUnit1,*,iostat=istatus) (Fcart(m,n),n=1,3*NumOfAt)
      if (istatus > 0) then
        write(u6,*)
        write(u6,*) ' ***************** ERROR *******************'
        write(u6,*) ' Error trying to read cartesian force cnsts.'
        write(u6,*) ' for the first state. Probably wrong input'
        write(u6,*) ' structure, perhaps too few values.'
        write(u6,*) ' *******************************************'
        call Quit_OnUserError()
      end if
    end do

    ! Check if Hessian is symmetric.
    error = Zero
    do m=1,(3*NumOfAt)
      do n=1,(3*NumOfAt)
        error = error+(Fcart(m,n)-Fcart(n,m))**2
      end do
    end do
    if (error > 1.0e-10_wp) then
      write(u6,*)
      write(u6,*) ' ******************** ERROR **********************'
      write(u6,*) ' Non-symmetric error in cartesian Hessian:',error
      write(u6,*) ' *************************************************'
      call Quit_OnUserError()
    end if

    ! Transform cartesian F matrix to internal coordinates.
    !PAM01: Here follows an n**4-scaling 2-index transformation
    !       It has been replaced by two n**3 matrix multiplies
    !do j=1,NumInt
    !  do i=1,NumInt
    !    Fsum = Zero
    !    do n=1,(3*NumOfAt)
    !      do m=1,(3*NumOfAt)
    !        Fsum = Fsum+Fcart(m,n)*Sinv(1+Mod((m+2),3),Int((m+2)/3),i)*Sinv(1+Mod((n+2),3),Int((n+2)/3),j)
    !      end do
    !    end do
    !    Hess1(i,j) = Fsum
    !  end do
    !end do
    !PAM01 Here follows replacement code:
    call mma_allocate(Temp,3*NumOfAt,NumInt,label='Temp')
    call DGEMM_('N','N',3*NumOfAt,NumInt,3*NumOfAt,One,FCart,3*NumOfAt,SInv,3*NumOfAt,Zero,Temp,3*NumOfAt)
    call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,SInv,3*NumOfAt,Temp,3*NumOfAt,Zero,Hess1,NumInt)
    call mma_deallocate(Temp)
    call mma_deallocate(Sinv)
    call mma_deallocate(Fcart)

  end if

  ! Scale Hessian if scaling factors were given.
  call KeyWord(inpUnit,'SCAL',.true.,exists)
  if (exists) then
    call mma_allocate(ScaleParam,NumInt,label='Scale1')

    call KeyWord(inpUnit,'FIRS',.false.,exists)
    do i=1,NumInt
      read(inpUnit,*) ScaleParam(i)
      ScaleParam(i) = sqrt(ScaleParam(i))
    end do
    do j=1,Numint
      Hess1(:,j) = Hess1(:,j)*ScaleParam*ScaleParam(j)
    end do
    call mma_deallocate(ScaleParam)
  end if
  !D write(u6,*) ' Scaled.'

  ! Hessian for second surface.
  call KeyWord(inpUnit,'FORC',.true.,exists)
  call KeyWord(inpUnit,'SECO',.false.,exists)
  !D write(u6,*) ' READINP: Read force constants for second state.'
  !read(inpUnit,*) CoordType
  !call UpCase(CoordType)
  ! Replace above two lines with:
  read(InpUnit,frmt) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  ! End of replacement
  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit2 = isfreeunit(32)
    call molcas_open(inpUnit2,'UNSYM2')
    !open (32,'UNSYM2')
    call KeyWord(inpUnit2,'UNSYMMETRIZED HESSIAN',.false.,exists)
    read(inpUnit2,'(a17)') Inline
    read(inpUnit2,'(a17)') Inline
  else
    inpUnit2 = inpUnit
  end if

  if (CoordType == 'INTERNAL') then
    do i=1,NumInt
      read(inpUnit,*) Hess2(i,:)
    end do
  else if (CoordType == 'CARTESIAN') then
    call mma_allocate(SS,3*NumOfAt,NumInt,label='SS')
    SS(:,:) = Zero
    call CalcS(AtCoord2,InterVec,SS,NumInt,NumofAt)

    ! Invert S matrix and remove total translation and total rotation.
    call mma_allocate(Sinv,3*NumOfAt,NumInt,label='Sinv')

    call RotTranRem(Sinv,SS,Mass,AtCoord2,NumOfAt,NumInt)
    call mma_deallocate(SS)

    ! Read cartesian force constant matrix.
    n = 3*NumOfAt
    nFcart = n
    call mma_allocate(Fcart,nFcart,nFcart,label='Fcart')

    do m=1,3*NumOfAt
      read(inpUnit2,'(a17)',iostat=istatus) Inline
      if (istatus <= 0) read(inpUnit2,*,iostat=istatus) (Fcart(m,n),n=1,3*NumOfAt)
      if (istatus > 0) then
        write(u6,*)
        write(u6,*) ' ***************** ERROR *******************'
        write(u6,*) ' Error trying to read cartesian force cnsts.'
        write(u6,*) ' for the second state. Probably wrong input'
        write(u6,*) ' structure, perhaps too few values.'
        write(u6,*) ' *******************************************'
        call Quit_OnUserError()
      end if
    end do

    ! Check if Hessian is symmetric.
    error = Zero
    do n=1,(3*NumOfAt)
      do m=1,(3*NumOfAt)
        error = error+(Fcart(m,n)-Fcart(n,m))**2
      end do
    end do
    if (error > 1.0e-2_wp) then
      write(u6,*)
      write(u6,*) ' ******************** ERROR **********************'
      write(u6,*) ' Non-symmetric error in cartesian Hessian:',error
      write(u6,*) ' *************************************************'
      call Quit_OnUserError()
    end if

    ! Transform cartesian F matrix to internal coordinates.
    !PAM01: Here follows an n**4-scaling 2-index transformation
    !       It has been replaced by two n**3 matrix multiplies
    !do j=1,NumInt
    !  do i=1,NumInt
    !    Fsum = Zero
    !    do n=1,(3*NumOfAt)
    !      do m=1,(3*NumOfAt)
    !        Fsum = Fsum+Fcart(m,n)*Sinv(1+Mod((m+2),3),Int((m+2)/3),i)*Sinv(1+Mod((n+2),3),Int((n+2)/3),j)
    !      end do
    !    end do
    !    Hess2(i,j) = Fsum
    !  end do
    !end do
    !PAM01 Here follows replacement code:
    call mma_allocate(Temp,3*NumOfAt,NumInt,label='Temp')
    call DGEMM_('N','N',3*NumOfAt,NumInt,3*NumOfAt,One,FCart,3*NumOfAt,SInv,3*NumOfAt,Zero,Temp,3*NumOfAt)
    call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,SInv,3*NumOfAt,Temp,3*NumOfAt,Zero,Hess2,NumInt)
    call mma_deallocate(Temp)
    call mma_deallocate(Sinv)
    call mma_deallocate(Fcart)

  end if

  if (inpUnit1 /= inpUnit) close(inpUnit1)
  if (inpUnit2 /= inpUnit) close(inpUnit2)

  ! FORCe: End ---------------------------------------------------------

  ! --------------------------------------------------------------------
  ! SCALe: Scale Hessian if scaling factors were given.

  call KeyWord(inpUnit,'SCAL',.true.,exists)
  if (exists) then
    call mma_allocate(ScaleParam,NumInt,label='ScaleParam')

    call KeyWord(inpUnit,'SECO',.false.,exists)
    do i=1,NumInt
      read(inpUnit,*) ScaleParam(i)
      ScaleParam(i) = sqrt(ScaleParam(i))
    end do
    do j=1,NumInt
      Hess2(:,j) = Hess2(:,j)*ScaleParam*ScaleParam(j)
    end do
    call mma_deallocate(ScaleParam)
  end if

  ! SCALe: End ---------------------------------------------------------

  ! --------------------------------------------------------------------
  ! DIPOles: Read Transition Dipoles or Spin-Orbit Coupling.

  call KeyWord(inpUnit,'DIPO',.true.,exists)
  if (.not. exists) call KeyWord(inpUnit,'SOC ',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ************ ERROR ****************'
    write(u6,*) ' The DIPOLES/SOC keyword is missing.'
    write(u6,*) ' ***********************************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read dipoles ("DIPOLES").'
  read(inpUnit,frmt) Inline
  call Normalize(InLine,OutLine)
  if (OutLine(1:4) == 'FILE') then
    !D write(u6,*) ' Read trans dips from file UNSYM21'
    call molcas_open(31,'UNSYM21')
    !open(31,'UNSYM21')
    if (max_dip == 1) then
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT X',.false.,exists)
      read(31,*) TranDipGrad(1,:)
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Y',.false.,exists)
      read(31,*) TranDipGrad(2,:)
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Z',.false.,exists)
      read(31,*) TranDipGrad(3,:)
    end if
    XnotFound = .true.
    YnotFound = .true.
    ZnotFound = .true.
    call KeyWord(31,'*BEGIN TRANSITION PROPERTIES',.false.,exists)
    do while (XnotFound .or. YnotFound .or. ZnotFound)
      read(31,'(A37)') InLine
      if (InLine(1:20) == 'MLTPL  1 COMPONENT1 ') then
        XnotFound = .false.
        read(InLine(22:37),*) TranDip(1)
      end if
      if (InLine(1:20) == 'MLTPL  1 COMPONENT2 ') then
        YnotFound = .false.
        read(InLine(22:37),*) TranDip(2)
      end if
      if (InLine(1:20) == 'MLTPL  1 COMPONENT3 ') then
        ZnotFound = .false.
        read(InLine(22:37),*) TranDip(3)
      end if
    end do
    close(31)
  else
    call KeyWord(inpUnit,'DIPO',.true.,exists)
    !D write(u6,*) ' READINP: Read dipoles ("DIPOLES").'
    read(inpUnit,*) TranDip(:)
    if (max_dip == 1) then
      read(inpUnit,*) TranDipGrad(1,:)
      read(inpUnit,*) TranDipGrad(2,:)
      read(inpUnit,*) TranDipGrad(3,:)
    end if
  end if

  !                         End of Forcefield part
  ! --------------------------------------------------------------------

else
  write(u6,*) ' Forcefield not found. Will use Energy Surface.'

  !--------------------------------------------------------------------!
  !
  !                       Energy surface input part
  !
  !--------------------------------------------------------------------!

  ! Read transformations.
  call KeyWord(inpUnit,'NONL',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' *********** ERROR **************'
    write(u6,*) ' Cannot find keyword  NONLINEAR !'
    write(u6,*) ' ********************************'
    call Quit_OnUserError()
  end if
  do icoord=1,NumInt
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
    trfName1(icoord) = OutLine
  end do
  do icoord=1,NumInt
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
    trfName2(icoord) = OutLine
  end do

  ! Read terms in polynomial.
  call KeyWord(inpUnit,'POLY',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ************ ERROR **************'
    write(u6,*) ' Cannot find keyword  POLYNOMIAL !'
    write(u6,*) ' *********************************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read polynomial.'
  k = 1
  nvar = 0
  read(inpUnit,frmt) InLine
  call WordPos(k,InLine,iStart,iStop)
  do while (iStop < len(Inline))
    k = iStop+1
    nvar = nvar+1
    call WordPos(k,InLine,iStart,iStop)
  end do
  nPolyTerm = 1
  read(inpUnit,frmt) InLine
  call Normalize(InLine,OutLine)
  do while (OutLine(1:4) /= 'END ')
    nPolyTerm = nPolyTerm+1
    read(inpUnit,frmt) InLine
    call Normalize(InLine,OutLine)
  end do
  call KeyWord(inpUnit,'POLY',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ************ ERROR **************'
    write(u6,*) ' Cannot find keyword  POLYNOMIAL !'
    write(u6,*) ' *********************************'
    call Quit_OnUserError()
  end if

  call mma_allocate(ipow,nPolyTerm,nvar,label='ipow')
  do iterm=1,nPolyTerm
    read(inpUnit,*) ipow(iterm,:)
  end do

  ! Find highest power of term in polynomial.
  max_term = 0
  do iterm=1,nPolyTerm
    nsum = 0
    do jvar=1,nvar
      nsum = nsum+ipow(iterm,jvar)
      if (nsum > max_term) max_term = nsum
    end do
  end do

  ! Read grid points and energies.
  CoordType = 'INTERNAL'
  call KeyWord(inpUnit,'DATA',.true.,exists)
  if (.not. exists) then
    write(u6,*)
    write(u6,*) ' ********* ERROR ***********'
    write(u6,*) ' Cannot find keyword  DATA !'
    write(u6,*) ' ***************************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read grid points (Key word "DATA").'
  idata = 0
  do while (.true.)
    idata = idata+1
    read(inpUnit,*,iostat=istatus) (coord(idata,jvar),jvar=1,nvar),(GrdVal(idata,j),j=1,10)
    if (istatus < 0) exit
  end do

  ndata = idata-1
  call mma_allocate(var,ndata,nvar,label='var')
  var(:,:) = coord(1:idata,1:nvar)

  ! CASPT2 energies for ground state.
  call mma_allocate(yin1,ndata,label='yin1')
  yin1(:) = GrdVal(1:ndata,2)

  ! CASPT2 energies for excited state.
  call mma_allocate(yin2,ndata,label='yin2')
  yin2(:) = GrdVal(1:ndata,6)

  ! If three atomic molecule, add Watson correction to energy values.
  if (NumOfAt == 3) then
    m1 = Mass(1)*uToAu
    m2 = Mass(2)*uToAu
    m3 = Mass(3)*uToAu
    m12 = (m1*m2)/(m1+m2)
    m23 = (m2*m3)/(m2+m3)
    do i=1,ndata
      r1 = var(i,1)
      r2 = var(i,2)
      theta = var(i,3)*deg2rad
      const = -(One/Eight)*(One/(m12*r1**2)+One/(m23*r2**2))*(One+(One/(sin(theta)**2)))
      const = const+Quart*(cos(theta)/(m2*r1*r2))*(cos(theta)/sin(theta))**2
      yin1(i) = yin1(i)+const
      yin2(i) = yin2(i)+const
    end do
  end if

  call mma_allocate(t_dipin1,ndata,label='t_dipin1')
  t_dipin1(:) = GrdVal(1:ndata,9)

  call mma_allocate(t_dipin2,ndata,label='t_dipin2')
  t_dipin2(:) = GrdVal(:,10)

  ! Make sure that all transition dipole values have the same sign.
  sign1_1 = sign(One,t_dipin1(1))
  sign2_1 = sign(One,t_dipin2(1))
  do idata=2,ndata
    sign1_2 = sign(One,t_dipin1(idata))
    sign2_2 = sign(One,t_dipin2(idata))
    if (int(sign1_1*sign2_1) /= int(sign1_2*sign2_2)) then
      write(u6,*) 'Sign shift in transition dipole data:',idata
    end if
    t_dipin1(idata) = sign1_1*sign1_2*t_dipin1(idata)
    t_dipin2(idata) = sign2_1*sign2_2*t_dipin2(idata)
  end do
  !D write(u6,*) ' READINP: Transition dipoles finished.'

  !                   End of Energy surface input part
  ! --------------------------------------------------------------------

end if

!D write(u6,*) ' Ending READINP.'

end subroutine ReadInp

subroutine WriteLog(PotCoef,AtomLbl,AtCoord,Mass,InterVec,stand_dev,max_err,energy,Hess,G,V,harmfreq,qMat,Bond,nBond,xvec,D3,D4, &
                    PED,x_anharm,anharmfreq,max_term,nState,ForceField,NumOfAt,nOsc)
!  Purpose:
!    Write results to log file.
!
!  Input:
!    AtomLbl    : Array of character - contains the labels for the atoms.
!    AtCoord    : Two dimensional Real array - contains the cartesian coordinates of the atoms.
!    Mass       : Real array - contains the mass of the atoms.
!    InterVec   : Integer array - contains the atoms that are used in the calculations of each internal coordinate.
!    ipow       : Two dimesional integer array - terms of polynomial.
!    PotCoef    : Real two dimensional array - coefficients of terms given in ipow.
!    stand_dev  : Real variable - standard deviation of fitted values of polynomial.
!    max_err    : Real variable - maximum error of fitted values of polynomial.
!    energy     : Real variable - energy in minimum.
!    Hess       : Real array -  contains the force constants expressed in internal coordinates.
!    G          : Real array.
!    V          : Real array - contains the eigenvectors of G*F as columns.
!    harmfreq   : Real array - harmonic frequencies.
!    qMat       : Real array - contains the cartesian displacements of the atoms.
!    Bond       : Integer array - contains atom pairs that are to be bonded together in a plot.
!    nBond      : Integer - dim of Bond, i.e. 2*(number of bonds).
!    xvec       : Real array - contains the internal coordinates at equilibrium.
!    D3         : Real array - contains the third derivatives.
!    D4         : Real array - contains the fourth derivatives.
!    PED        : Real three dimensional array - potential energy distribution.
!    x_anharm   : Real two dimensional array - anharmonicity constants.
!    anharmfreq : Real array - anharmonic frequencies.
!    max_term   : Integer - highest power of term in fitted polynomial.
!    nState     : Integer - 1 or 2, depending upon which state it is.
!
!  Output:
!    Log file
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use mula_global, only: Huge_Print, ipow, MaxNumAt, nPolyTerm, VibModPlot
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: auTocm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: InterVec(MaxNumAt*15), Bond(MaxNumAt*2), nBond, max_term, nState, NumOfAt, nOsc
real(kind=wp), intent(in) :: PotCoef(nPolyTerm,1), AtCoord(3,NumOfAt), Mass(MaxNumAt), stand_dev, max_err, energy, &
                             Hess(nOsc,nOsc), G(nOsc,nOsc), V(nOsc,nOsc), harmfreq(nOsc), qMat(3*NumOfAt,nOsc), xvec(nOsc), &
                             D3(nOsc,nOsc,nOsc), D4(nOsc,nOsc,nOsc,nOsc), PED(nOsc,nOsc,nOsc), x_anharm(nOsc,nOsc), anharmfreq(nOsc)
character(len=4), intent(in) :: AtomLbl(MaxNumAt)
logical(kind=iwp), intent(in) :: ForceField
integer(kind=iwp) :: i, iBond, iInt, imode, j, k, l, m1, m2, maxCol, mInt, nAtom, nCol, n_Int, nRow, NumInt, VibPlotUnit
real(kind=wp) :: vLength, vMax, x1, x2, x3
logical(kind=iwp) :: cont
character(len=11) :: VibPlotFile
character(len=8) :: BondString
integer(kind=iwp), allocatable :: aNormModes(:)
integer(kind=iwp), external :: isfreeunit
real(kind=wp), external :: Dnrm2_

NumInt = nOsc
call mma_allocate(aNormModes,NumInt,label='aNormModes')
do i=1,NumInt
  aNormModes(i) = i
end do

if (nState == 1) then
  write(u6,*)
  write(u6,*)
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ','|           Results for the first state           |'
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,*)
end if
if (nState == 2) then
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ','|           Results for the second state          |'
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,*)
end if

! Write coefficients of polynomial fit, standard deviation and
! maximum error to log.
if (.not. ForceField) then
  write(u6,*)
  write(u6,'(A)') '  Coefficients of different terms in polynomial:'
  write(u6,'(A)') '  ----------------------------------------------'
  do i=1,nPolyTerm
    write(u6,'(1x,f18.12,10(3x,i2))') PotCoef(i,1),ipow(i,:)
  end do
  write(u6,*)
  write(u6,*)
  write(u6,'(A,es14.6)') '  Standard deviation of fitted values:',stand_dev
  write(u6,'(A)') '  ------------------------------------'
  write(u6,*)
  write(u6,'(A,es14.6)') '  Maximum error of fitted values:     ',max_err
  write(u6,'(A)') '  -------------------------------'

  ! Write energy in minimum.
  write(u6,*)
  write(u6,*)
  write(u6,'(A,F15.8,A)') '  Energy in minimum:',energy
  write(u6,'(A)') '  ------------------'
end if

! Write force constant matrix to log file.
if (Huge_Print) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Force Constant Matrix:'
  write(u6,*) ' ','======================'
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (n_Int,n_Int=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (Hess(mInt,n_Int),n_Int=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do

  ! Write matrix G to log file.
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Inverse Mass Tensor :'
  write(u6,*) ' ','====================='
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (n_Int,n_Int=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (G(mInt,n_Int),n_Int=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do
  write(u6,*)

  ! Eigenvectors.
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Eigenvectors :'
  write(u6,*) ' ','=============='
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (n_Int,n_Int=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (V(mInt,n_Int),n_Int=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do
  write(u6,*)
end if

! Third derivatives.
if (max_term > 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Cubic force constants:'
  write(u6,'(A)') '  ----------------------'
  do i=1,NumInt
    do j=i,NumInt
      do k=j,NumInt
        write(u6,'(a,i2,i2,i2,f15.8)') '  ',i,j,k,D3(i,j,k)
      end do
    end do
  end do
end if

! Fourth derivatives.
if (max_term > 3) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Quartic force constants:'
  write(u6,'(A)') '  ------------------------'
  do i=1,NumInt
    do j=i,NumInt
      do k=j,NumInt
        do l=k,NumInt
          write(u6,'(a,i2,i2,i2,i2,f15.8)') '  ',i,j,k,l,D4(i,j,k,l)
        end do
      end do
    end do
  end do
end if

! Anharmonicity constants.
if (max_term > 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Anharmonicity constants:'
  write(u6,'(A)') '  ------------------------'
  do i=1,NumInt
    write(u6,'(a,20f15.8)') ' ',(auTocm*x_anharm(i,j),j=1,i)
  end do
end if

! Write labels, coordinates and masses to log file.
call WriteCartCoord(AtomLbl,AtCoord,Mass,NumOfAt)

! Write internal coordinates to log file.
call WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)

! Harmonic frequencies in reciprocal cm, GHz and hartrees.
call WriteFreq(harmfreq,aNormModes,NumInt,'Harmonic frequencies')

! Fundamental frequencies in reciprocal cm, GHz and hartrees.
if (max_term > 2) call WriteFreq(anharmfreq,aNormModes,NumInt,'Anharmonic frequencies')

! Write Potential Energy Distribution for each mode to log
! if the molecule is not too large.
if (Huge_Print .and. (NumInt <= 10)) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Potential Energy Distribution :'
  write(u6,*) ' ','==============================='
  do imode=1,NumInt
    write(u6,*)
    write(u6,*) ' ','Mode',imode
    write(u6,*) ' ','-------'
    do k=1,NumInt
      if (NumInt < 21) then
        write(u6,'(60F6.2)') (PED(k,l,imode),l=1,NumInt)
      else
        write(u6,'(60F5.2)') (PED(k,l,imode),l=1,NumInt)
      end if
    end do
  end do
end if

! Find length of longest displacement vector.
vMax = Dnrm2_(3,qMat(1,1),1)
do nAtom=1,NumOfAt
  do j=1,NumInt
    vLength = Dnrm2_(3,qMat((nAtom-1)*3+1,j),1)
    if (vLength > vMax) then
      vMax = vLength
    end if
  end do
end do

! Print displacement vectors scaled with length of longest vector
! if the molecule is not too large.
if (Huge_Print .and. (NumInt <= 10)) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Scaled displacement vectors :'
  write(u6,*) ' ','============================='
  write(u6,*) ' ','(length of longest vector set equal to one)'
  write(u6,*)
  maxCol = 10
  do nAtom=1,NumOfAt
    write(u6,*) ' ',' Atom: ',AtomLbl(nAtom)
    write(u6,*) ' ',' ---------'
    nRow = NumInt/maxCol
    nCol = mod(NumInt,maxCol)
    m1 = 1
    m2 = 0
    do i=0,nRow
      if (i == nRow) then
        m2 = m2+nCol
      else
        m2 = m2+maxCol
      end if
      if (m2 /= 0) then
        write(u6,'(10I12)') (j,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  x ',(qMat((nAtom-1)*3+1,j)/vMax,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  y ',(qMat((nAtom-1)*3+2,j)/vMax,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  z ',(qMat((nAtom-1)*3+3,j)/vMax,j=m1,m2)
        write(u6,*)
        m1 = m2+1
      end if
    end do
  end do
  write(u6,*)
  write(u6,*)
end if

! Generate MOLDEN input:

if (nState == 1) then
  call WrMold('MOLDEN-1',NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
else
  call WrMold('MOLDEN-2',NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
end if

! Print geometry of molecule, together with cartesian displacements,
! in a file.
if (VibModPlot) then
  VibPlotUnit = isfreeunit(21)
  if (nState == 1) then
    VibPlotFile = 'plot.modes1'
  end if
  if (nState == 2) then
    VibPlotFile = 'plot.modes2'
  end if
  call molcas_Open(vibplotunit,vibplotfile)
  !open(VibPlotUnit,VibPlotFile)
  write(VibPlotUnit,*) 2*NumOfAt
  write(VibPlotUnit,*) nBond/2
  write(VibPlotUnit,*) '1'
  write(VibPlotUnit,*) NumInt
  iBond = 1
  do while (iBond+1 <= nBond)
    write(BondString,'(i4,i4)') Bond(iBond),Bond(iBond+1)
    write(VibPlotUnit,'(a8)') BondString
    iBond = iBond+2
  end do
  do nAtom=1,NumOfAt
    x1 = AtCoord(1,nAtom)
    x2 = AtCoord(2,nAtom)
    x3 = AtCoord(3,nAtom)
    write(VibPlotUnit,'(a,3f12.3)') '>',x2,x3,x1
  end do
  do iInt=1,NumInt
    do nAtom=1,NumOfAt
      x1 = AtCoord(1,nAtom)
      x2 = AtCoord(2,nAtom)
      x3 = AtCoord(3,nAtom)
      write(VibPlotUnit,'(a,3f12.3)') AtomLbl(nAtom),x2+qMat((nAtom-1)*3+2,iInt)/vMax,x3+qMat((nAtom-1)*3+3,iInt)/vMax, &
                                      x1+qMat((nAtom-1)*3+1,iInt)/vMax
    end do
  end do
  close(VibPlotUnit)
end if

call mma_deallocate(aNormModes)

end subroutine WriteLog

subroutine WriteDip(DipGrad,Modes,Title,nOsc)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOsc, Modes(nOsc)
real(kind=wp), intent(in) :: DipGrad(3,nOsc)
character(len=*), intent(in) :: Title
integer(kind=iwp) :: i

write(u6,*)
write(u6,*)
write(u6,'(a2,a)') ' ',Title
write(u6,'(a2,a)') ' ','=============================================='
write(u6,'(a2,a)') ' ',' mode          X           Y           Z'
write(u6,'(a2,a)') ' ','----------------------------------------------'
do i=1,nOsc
  write(u6,'(a3,i2,a1,a3,3f12.5)') ' ',Modes(i),'.',' ',DipGrad(:,i)
end do
write(u6,'(a2,a)') ' ','=============================================='
write(u6,*)

end subroutine WriteDip

subroutine IntCalcHeader()

use Definitions, only: u6

implicit none

write(u6,*)
write(u6,*)
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ','|             Intensity calculation               |'
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,*)

end subroutine IntCalcHeader

subroutine ExpPointHeader()

use Definitions, only: u6

implicit none

write(u6,*)
write(u6,*)
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ','|            Expansion point geometry             |'
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,*)

end subroutine ExpPointHeader

subroutine ISCHeader()

use Definitions, only: u6

implicit none

write(u6,*)
write(u6,*)
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ','|      InterSystem Crossing rate calculation      |'
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,*)

end subroutine ISCHeader

subroutine WriteHeader(Title)
!  Purpose:
!    Write header and title to logfile.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Definitions, only: u6

implicit none
character(len=80), intent(in) :: Title

write(u6,*)

! Write title of project.
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,'(A,A)') '  Title : ',Title
write(u6,'(A)') '  -------'
write(u6,*)

end subroutine WriteHeader

subroutine WrMold(FName,NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
! Open file named FName, and create a MOLDEN input file.

use Constants, only: Zero, auTocm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NumOfAt, NumInt
character(len=*), intent(in) :: FName, AtomLbl(NumOfAt)
real(kind=wp), intent(in) :: AtCoord(3,NumOfAt), HarmFreq(NumInt), QMat(3,NumOfAt,NumInt)
integer(kind=iwp) :: i, iAtom, iInt
real(kind=wp) :: D2, DispMx, DMax, Factor
character(len=2) :: AtName

call molcas_open(9,Fname)
!open(9,FName,status='UNKNOWN')
write(9,*) '[MOLDEN FORMAT]'

! Harmonic frequencies:
write(9,*) '[N_FREQ]'
write(9,*) NumInt
write(9,*) '[FREQ]'
do iInt=1,NumInt
  write(9,'(1X,F10.3)') auTocm*HarmFreq(iInt)
end do
! VV: Not implemented???
write(9,*) '[INT]'
do iInt=1,NumInt
  write(9,'(1X,F10.3)') Zero
end do
! Atom coordinates:
write(9,*) '[NATOM]'
write(9,*) NumOfAt
write(9,*) '[FR-COORD]'
do iAtom=1,NumOfAt
  AtName = AtomLbl(iAtom)(1:2)
  if ((ichar('0') <= ichar(AtName(2:2))) .and. (ichar(AtName(2:2)) <= ichar('9'))) then
    AtName(2:2) = ' '
  end if
  write(9,'(1x,A2,3F16.8)') AtName,(AtCoord(i,iAtom),i=1,3)
end do

! Cartesian displacement coordinates:
! Scale such that max displacement is 0.2 A.U.
DMax = 0.2_wp
write(9,*) '[FR-NORM-COORD]'
do iInt=1,NumInt
  write(9,*) 'vibration ',iInt
  DispMx = Zero
  do iAtom=1,NumOfAt
    D2 = Zero
    do i=1,3
      D2 = D2+QMat(i,iAtom,iInt)**2
    end do
    DispMx = max(DispMx,sqrt(D2))
  end do
  Factor = DMax/DispMx
  do iAtom=1,NumOfAt
    write(9,'(1X,3F16.8)') (Factor*QMat(i,iAtom,iInt),i=1,3)
  end do
end do

close(9)

end subroutine WrMold

!end module InOutMod

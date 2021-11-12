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

subroutine ReadInp(Title,AtomLbl,Mass,InterVec,Bond,nBond,NumInt,NumOfAt,trfName1,trfName2,m_max,n_max,max_dip,max_term,MatEl, &
                   ForceField,Cartesian,lExpan,lISC,iCode,dMinWind)
!  Purpose:
!    Read the input file.
!
!  Output:
!    Title    : String - title of the project.
!    AtomLbl  : Array of character - contains the labels for the
!               atoms.
!    Mass     : Real*8 array - contains the mass of the
!               atoms.
!    InterVec : Integer array - containis the atoms that are used
!               in the calculations of each internal coordinate.
!    Bond     : Integer array - contains atom pairs that are to be
!               bonded together in a plot.
!    nBond    : Integer - dim of Bond, i.e. 2*(number of bonds).
!    NumInt   : Integer - the total number of internal coordinates.
!    NumOfAt  : Integer - the number of atoms.
!    trfName  : Character array - type of transformation of variables.
!    m_max    : Integer - maximum level for the first state.
!    n_max    : Integer - maximum level for the second state.
!    m_plot   : Integer array - level(s) to plot for the first state.
!    n_plot   : Integer array - level(s) to plot for the second state.
!    max_dip  : Integer - highest order of term in transition dipole.
!    max_term : Integer - highest power of a term in polynomial fitted
!               to energy values.
!    MatEl    : Logical
!    lISC     : Logical to calculate InterSystem Crossing
!
!  Uses:
!    IOTools
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.
!
!  Modified to use isotopes module
!    Ignacio Fdez. Galvan, 2017

use Isotopes, only: Isotope
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Eight, Quart
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
#include "inout.fh"
#include "Constants_mula.fh"
#include "inputdata.fh"
#include "indims.fh"
integer InterVec(MaxNumAt*15)
character*80 trfName1(MaxNumAt), trfName2(MaxNumAt)
real*8 Mass(MaxNumAt)
character*4 AtomLbl(MaxNumAt)
integer Bond(MaxNumAt*2)
character*80 Title
!logical Plot
character*9 CoordType
character*1 c
character*7 AtInp
character*2 AtName
character*2 Eunit
character*3 DegOrRad
character*2 AAOrAu
character*4 Atom
character*6 PlUnit
character*4 Atom1, Atom2, Atom3, Atom4
!character*2 AtList(0:120)
!integer AtData(120,4)
!real*8 MassData(400)
real*8 coord(1000,10)
real*8 GrdVal(1000,10)
logical MatEl, ForceField, lISC
logical Cartesian
logical exist
logical lExpan
real*8 m1, m2, m3, m12, m23
real*8 TmpCoord(3)
parameter(max_len_xvec=150)
real*8 xvec(max_len_xvec)
parameter(max_plot_temp=100)
integer plot_temp(max_plot_temp), iCode
logical XnotFound, YnotFound, ZnotFound
! Format declarations.
character*8 format
! User-defined functions called:
external StrToDble, iStrToInt
#include "WrkSpc.fh"
integer, allocatable :: GeoVec(:), tempmodes(:)
real*8, allocatable :: Fcart(:,:), ScaleParam(:), Sinv(:,:), SS(:,:), Temp(:,:)

format = '(A80)'

! Read from MassFile:
! - name of atoms into array of string - AtList.
! - atomic number, default mass number, smallest mass number and
!   table offset into two dimensional array - AtData.
! - masses of all isotopes into array - MassData.

!i = 0
!call Molcas_Open(massUnit,'MASSUNIT')
!read(massUnit,Format) InLine
! Read sequential lines from atomic data file.
!call Normalize(InLine,OutLine)
!do while (OutLine(1:4) /= 'END ')
!  l = Index(OutLine,'*')
!  if (l == 0) then
!    i = i+1
!    AtList(i) = OutLine(1:2)
!    read(OutLine(3:80),*) (AtData(i,j),j=1,4)
!  end if
!  read(massUnit,Format) InLine
!  call Normalize(InLine,OutLine)
!end do
!iAtList = i

!i = 0
!read(massUnit,Format) InLine
!call Normalize(InLine,OutLine)
!do while (OutLine(1:4) /= 'END ')
!  l = Index(OutLine,'*')
!  if (l == 0) then
!    i = i+1
!    read(OutLine,*) MassData(i)
!  end if
!  read(massUnit,Format) InLine
!  call Normalize(InLine,OutLine)
!end do
!close(massUnit)
iPrint = iPrintLevel(-1)

! ----------------------------------------------------------------------
! TITLe: Read the title into a string - Title

Title = ' '
call KeyWord(inpUnit,'TITL',.true.,exist)
if (exist) read(inpUnit,format) Title

! TITLe: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! ATOMs:  Process ATOMS

call KeyWord(inpUnit,'ATOM',.true.,exist)
NumOfAt = 0
read(inpUnit,format) InLine
call Normalize(InLine,OutLine)
do while (OutLine(1:4) /= 'END ')
  NumOfAt = NumOfAt+1
  read(inpUnit,format) InLine
  call Normalize(InLine,OutLine)
end do

! Read the labels of atoms into an array - AtomLbl.

call KeyWord(inpUnit,'ATOM',.true.,exist)
do nAtom=1,NumOfAt
  read(inpUnit,format) InLine
  call Normalize(InLine,OutLine)
  k = 1
  call WordPos(k,OutLine,iStart,iStop)
  k = iStop+1
  AtInp = '       '
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
    AtInp = '      '
    AtInp = OutLine(iStart:iStop)
    i = 1
    c = AtInp(i:i)
  end if
  ! Extract atomic name and possible label.
  AtName = '  '
  AtomLbl(nAtom) = '    '
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
    !  write(u6,*) ' is not listed in MassFile.             '
    !  write(u6,*) '****************************************'
    !  call Quit_OnUserError()
    !end if
  end if
end do

! ATOMs: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! INTErnal: Resolve the different internal coordinates specified in the input.

call KeyWord(inpUnit,'INTE',.true.,exist)
j = 1
NumInt = 0
nBond = 0
read(inpUnit,format) InLine
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

  read(inpUnit,format) InLine
  call Normalize(InLine,OutLine)
end do

call GetMem('AtCoord1','Allo','Real',ipAtCoord1,3*NumOfAt)
call GetMem('AtCoord2','Allo','Real',ipAtCoord2,3*NumOfAt)
l_Hess1 = NumInt
call GetMem('Hess1','Allo','Real',ipHess1,l_Hess1*l_Hess1)
call GetMem('Hess2','Allo','Real',ipHess2,l_Hess1*l_Hess1)
n = 3*NumOfAt
call GetMem('TranDipGrad','Allo','Real',ipTranDipGrad,3*n)

! INTErnal: End --------------------------------------------------------

! ----------------------------------------------------------------------
! ISC: InterSystem Crossing. GG-30-Dec-08

lISC = .false.
dMinWind = Zero
call KeyWord(inpUnit,'ISC ',.true.,lISC)
if (lISC) then
  call KeyWord(inpUnit,'EXPF',.true.,exist)
  if (exist) then
    read(InpUnit,format) InLine
    call Normalize(InLine,Outline)
    read(OutLine,*) dMinWind
  end if
  call KeyWord(inpUnit,'OLDC',.true.,exist)
  if (exist) iCode = iCode+1
  call KeyWord(inpUnit,'DISK',.true.,exist)
  if (exist) iCode = iCode+10
end if

! ISC: End -------------------------------------------------------------

! ----------------------------------------------------------------------
! MODEs: Chose which modes to use in intensity calculations.

call KeyWord(inpUnit,'MODE',.true.,exist)
if (exist) then
  call mma_allocate(tempmodes,NumInt,label='tempmodes')
  read(inpUnit,format) InLine
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
    read(inpUnit,format) InLine
    call Normalize(InLine,OutLine)
  end do
  n = i-1
  l_NormModes = n
  call GetMem('NormModes','Allo','Inte',ipNormModes,l_NormModes)
  do iv=1,n
    iWork(ipNormModes+iv-1) = tempModes(iv)
  end do
  call mma_deallocate(tempmodes)
  do i=1,n-1
    do j=i+1,n
      if (iWork(ipNormModes+j-1) < iWork(ipNormModes+i-1)) then
        modeTemp = iWork(ipNormModes+i-1)
        iWork(ipNormModes+i-1) = iWork(ipNormModes+j-1)
        iWork(ipNormModes+j-1) = modeTemp
      end if
    end do
  end do
else
  ! No MODEs specified
  l_NormModes = NumInt
  call GetMem('NormModes','Allo','Inte',ipNormModes,l_NormModes)
  do i=1,NumInt
    iWork(ipNormModes+i-1) = i
  end do
  if (iPrint >= 1) then
    write(u6,*) ' *** Warning: MODEs not specified !'
    write(u6,*) '     All ',l_NormModes,' will be calculated.'
    write(u6,*)
  end if
end if

! MODEs: End -----------------------------------------------------------

! ----------------------------------------------------------------------
! MXLEvels: Maximun number of vibr. quanta in first and second state

call KeyWord(inpUnit,'MXLE',.true.,exist)
m_max = 0 ! Defaul
n_max = 1 ! Defaul
if (.not. exist) then
  if (iPrint >= 1) then
    write(u6,*) ' *** Warning: MXLEvels not specified !'
    write(u6,*) '     Default values are 0 and 1.      '
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

call KeyWord(inpUnit,'TRAN',.true.,exist)
if (.not. exist) then
  if (.not. lISC) then
    write(u6,*) ' *** Warning: TRANsitions not specified!'
    write(u6,*) '     All Levels will be printed.        '
    write(u6,*)
  end if
  n = m_max+1
  do i=1,n
    plot_temp(i) = i-1
  end do
else
  call KeyWord(inpUnit,'FIRS',.false.,exist)
  read(inpUnit,format) InLine
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
      write(u6,*) ' is larger than max quantum.       '
      write(u6,*) ' **********************************'
      call Quit_OnUserError()
    end if
    i = i+1
    call WordPos(k,OutLine,iStart,iStop)
  end do
  n = i-1
end if
l_m_plot = n
call GetMem('m_plot','Allo','Inte',ipm_plot,l_m_plot)
do iv=1,n
  iWork(ipm_plot+iv-1) = plot_temp(iv)
end do

call KeyWord(inpUnit,'TRAN',.true.,exist)
if (.not. exist) then
  n = n_max+1
  do i=1,n
    plot_temp(i) = i-1
  end do
else
  call KeyWord(inpUnit,'SECO',.false.,exist)
  read(inpUnit,format) InLine
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
l_n_plot = n
call GetMem('n_plot','Allo','Inte',ipn_plot,l_n_plot)
do iv=1,n
  iWork(ipn_plot+iv-1) = plot_temp(iv)
end do

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

  call KeyWord(inpUnit,'ENER',.true.,exist)
  if (.not. exist) then
    write(u6,*)
    write(u6,*) ' ************** ERROR *************'
    write(u6,*) ' Electronic T_e energies not given.'
    write(u6,*) ' **********************************'
    call Quit_OnUserError()
  end if
  call KeyWord(inpUnit,'FIRS',.false.,exist)
  if (.not. exist) then
    write(u6,*)
    write(u6,*) ' ******************** ERROR *********************'
    write(u6,*) ' Electronic T_e energy for first state not given.'
    write(u6,*) ' ************************************************'
    call Quit_OnUserError()
  end if
  !read(inpUnit,*) energy1,Eunit
  !call Upcase(Eunit)
  ! Replace above two lines with:
  read(InpUnit,format) InLine
  call Normalize(InLine,Outline)
  read(OutLine,*) Energy1
  EUnit = 'AU'
  if (index(OutLine,'EV') > 0) EUnit = 'EV'
  if (index(OutLine,'CM') > 0) EUnit = 'CM'
  ! End of replacement.

  rfact = One
  if (Eunit == 'AU') rfact = One
  if (Eunit == 'EV') rfact = One/auToeV
  if (Eunit == 'CM') rfact = One/HarToRcm
  energy1 = energy1*rfact
  call KeyWord(inpUnit,'SECO',.false.,exist)
  if (.not. exist) then
    write(u6,*)
    write(u6,*) ' ******************** ERROR **********************'
    write(u6,*) ' Electronic T_e energy for second state not given.'
    write(u6,*) ' *************************************************'
    call Quit_OnUserError()
  end if
  !read(inpUnit,*) energy2,Eunit
  !call Upcase(Eunit)
  ! Replace above two lines with:
  read(InpUnit,format) InLine
  call Normalize(InLine,Outline)
  read(OutLine,*) Energy2
  EUnit = 'AU'
  if (index(OutLine,'EV') > 0) EUnit = 'EV'
  if (index(OutLine,'CM') > 0) EUnit = 'CM'
  ! End of replacement.
  if (Eunit == 'AU') rfact = One
  if (Eunit == 'EV') rfact = One/auToeV
  if (Eunit == 'CM') rfact = One/HarToRcm
  energy2 = energy2*rfact

  ! ENERgy: : End ------------------------------------------------------

  ! --------------------------------------------------------------------
  ! GEOMetries: Read geometries for the two surfaces.

  call KeyWord(inpUnit,'GEOM',.true.,exist)
  if (.not. exist) then
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
  read(InpUnit,format) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  ! End of replacement.

  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit1 = 31
    inpUnit2 = 32
    call molcas_open(31,'UNSYM1')
    !open(unit=31,file='UNSYM1')
    call KeyWord(inpUnit1,'*LABEL COORDINATES CHARGE',.false.,exist)
    call molcas_open(32,'UNSYM2')
    !open(unit=32,file='UNSYM2')
    call KeyWord(inpUnit2,'*LABEL COORDINATES CHARGE',.false.,exist)
  else
    inpUnit1 = inpunit
    inpUnit2 = inpunit
  end if

  if (CoordType == 'CARTESIAN') then
    Cartesian = .true.
    !D write(u6,*) ' READINP: Read cartesian coordinates for first state.'
    call KeyWord(inpUnit,'FIRS',.false.,exist)
    if (.not. exist) then
      write(u6,*)
      write(u6,*) ' ****** ERROR ********'
      write(u6,*) ' Wrong geometry input.'
      write(u6,*) ' *********************'
      call Quit_OnUserError()
    end if
    ierror = 0
    do iAtom=1,NumOfAt
      read(inpUnit1,format) InLine
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
        do iv=1,3
          Work(ipAtCoord1+iv+3*(j-1)-1) = TmpCoord(iv)
        end do
      end if
    end do
    !D write(u6,*) ' READINP: Read cartesian coordinates for second state.'
    call KeyWord(inpUnit,'SECO',.false.,exist)
    if (.not. exist) then
      write(u6,*)
      write(u6,*) ' ****** ERROR ********'
      write(u6,*) ' Wrong geometry input.'
      write(u6,*) ' *********************'
      call Quit_OnUserError()
    end if
    do iAtom=1,NumOfAt
      read(inpUnit2,format) InLine
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
        do iv=1,3
          Work(ipAtCoord2+iv+3*(j-1)-1) = TmpCoord(iv)
        end do
      end if
    end do

  else if (CoordType == 'INTERNAL') then

    Cartesian = .false.
    call KeyWord(inpUnit,'FIRS',.false.,exist)
    !D write(u6,*) ' READINP: Read internal coordinates for first state.'
    n = 5*(3*NumOfAt-5)
    ! Largest possible necessary GeoVec
    call mma_allocate(GeoVec,n,label='GeoVec')
    iInt = 1
    j = 1
    read(inpUnit,format) InLine
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
      read(inpUnit,format) InLine
      call Normalize(InLine,OutLine)
    end do
    iInt = iInt-1
    !call Int_to_Cart(GeoVec,xvec,AtCoord1,NumOfAt,iInt,Mass)
    l_a = NumOfAt
    l_n = max_len_xvec
    call Int_to_Cart1(GeoVec,xvec,Work(ipAtCoord1),l_a,l_n)

    call mma_deallocate(GeoVec)

    call KeyWord(inpUnit,'SECO',.false.,exist)
    !D write(u6,*) ' READINP: Read internal coordinates for second state.'
    n = 5*(3*NumOfAt-5)
    call mma_allocate(GeoVec,n,label='GeoVec')
    iInt = 1
    j = 1
    read(inpUnit,format) InLine
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
      read(inpUnit,format) InLine
      call Normalize(InLine,OutLine)
    end do
    iInt = iInt-1
    !call Int_to_Cart(GeoVec,xvec,AtCoord2,NumOfAt,iInt,Mass)
    l_a = NumOfAt
    l_n = max_len_xvec
    call Int_to_Cart1(GeoVec,xvec,Work(ipAtCoord2),l_a,l_n)
    call mma_deallocate(GeoVec)
    call GetMem('GeoVec','Free','Inte',ipGeoVec,5*(3*NumOfAt-5))
  end if

  if (inpunit1 == 31) then
    close(inpUnit1)
    close(inpUnit2)
  end if

  ! GEOMetries: End ----------------------------------------------------

  ! --------------------------------------------------------------------
  ! MXORder: Read maximum order for the transition dipole moment.

  call KeyWord(inpUnit,'MXOR',.true.,exist)
  if (.not. exist) then
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
      read(inpunit,*) cmstart,cmend
    end if
  end if

  ! PLOT : End ---------------------------------------------------------

  PlUnit = '      '
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
      PlUnit = ' nm.  '
      write(u6,*) '     The plot file will be in nanometers and the '
    else if (Use_cm) then
      PlUnit = ' cm-1.'
      write(u6,*) '     The plot file will be in cm-1 and the'
    else
      PlUnit = ' eV.  '
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
    call KeyWord(inpUnit,'LIFE',.true.,exist)
    if (exist) read(inpunit,*) LifeTime
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
    call KeyWord(inpUnit,'CYCL',.false.,exist)
    if (exist) then
      call KeyWord(inpUnit,'VIBP',.true.,exist)
      read(inpUnit,format) InLine
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

  call KeyWord(inpUnit,'FORC',.true.,exist)
  max_term = 2
  call KeyWord(inpUnit,'FIRS',.false.,exist)
  !D write(u6,*) ' READINP: Read force constants for first state.'
  read(InpUnit,format) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  !D write(u6,*) ' READINP: CoordType=',CoordType

  ! Read Hessian. If Hessian is given in cartesian coordinates,
  ! first remove total translation and total rotation and then
  ! transform it to internal coordinates.
  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit1 = 31
    call molcas_open(31,'UNSYM1')
    !open(unit=31,file='UNSYM1')
    call KeyWord(inpUnit1,'UNSYMMETRIZED HESSIAN',.false.,exist)
    read(inpUnit1,'(a17)') Inline
    read(inpUnit1,'(a17)') Inline
  else
    inpUnit1 = inpUnit
  end if

  if (CoordType == 'INTERNAL') then
    do i=1,NumInt
      read(inpUnit,*) (Work(ipHess1+i+l_Hess1*(j-1)-1),j=1,NumInt)
    end do
  else if (CoordType == 'CARTESIAN') then
    call mma_allocate(SS,3*NumOfAt,NumInt,label='SS')
    SS(:,:) = Zero
    call CalcS(Work(ipAtCoord1),InterVec,SS,NumInt,NumOfAt)

    ! Invert S matrix and remove total translation and total rotation.
    call mma_allocate(Sinv,3*NumOfAt,NumInt,label='Sinv')

    call RotTranRem(Sinv,SS,Mass,Work(ipAtCoord1),NumOfAt,NumInt)
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
        write(u6,*) ' for the first state. Probably wrong input  '
        write(u6,*) ' structure, perhaps too few values.         '
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
    call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,SInv,3*NumOfAt,Temp,3*NumOfAt,Zero,Work(ipHess1),NumInt)
    call mma_deallocate(Temp)
    call mma_deallocate(Sinv)
    call mma_deallocate(Fcart)

  end if

  ! Scale Hessian if scaling factors were given.
  call KeyWord(inpUnit,'SCAL',.true.,exist)
  if (exist) then
    call mma_allocate(ScaleParam,NumInt,label='Scale1')

    call KeyWord(inpUnit,'FIRS',.false.,exist)
    do i=1,NumInt
      read(inpUnit,*) ScaleParam(i)
      ScaleParam(i) = sqrt(ScaleParam(i))
    end do
    do j=1,Numint
      do i=1,NumInt
        Work(ipHess1+i+l_Hess1*(j-1)-1) = Work(ipHess1+i+l_Hess1*(j-1)-1)*ScaleParam(i)*ScaleParam(j)
      end do
    end do
    call mma_deallocate(ScaleParam)
  end if
  !D write(u6,*) ' Scaled.'

  ! Hessian for second surface.
  call KeyWord(inpUnit,'FORC',.true.,exist)
  call KeyWord(inpUnit,'SECO',.false.,exist)
  !D write(u6,*) ' READINP: Read force constants for second state.'
  !read(inpUnit,*) CoordType
  !call UpCase(CoordType)
  ! Replace above two lines with:
  read(InpUnit,format) InLine
  call Normalize(InLine,Outline)
  iend = index(OutLine,' ')
  CoordType = OutLine(1:iend-1)
  ! End of replacement
  if (CoordType == 'FILE') then
    Coordtype = 'CARTESIAN'
    inpUnit2 = 32
    call molcas_open(32,'UNSYM2')
    !open (unit=32,file='UNSYM2')
    call KeyWord(inpUnit2,'UNSYMMETRIZED HESSIAN',.false.,exist)
    read(inpUnit2,'(a17)') Inline
    read(inpUnit2,'(a17)') Inline
  else
    inpUnit2 = inpUnit
  end if

  if (CoordType == 'INTERNAL') then
    do i=1,NumInt
      read(inpUnit,*) (Work(ipHess2+i+l_Hess1*(j-1)-1),j=1,NumInt)
    end do
  else if (CoordType == 'CARTESIAN') then
    call mma_allocate(SS,3*NumOfAt,NumInt,label='SS')
    SS(:,:) = Zero
    call CalcS(Work(ipAtCoord2),InterVec,SS,NumInt,NumofAt)

    ! Invert S matrix and remove total translation and total rotation.
    call mma_allocate(Sinv,3*NumOfAt,NumInt,label='Sinv')

    call RotTranRem(Sinv,SS,Mass,Work(ipAtCoord2),NumOfAt,NumInt)
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
        write(u6,*) ' for the second state. Probably wrong input '
        write(u6,*) ' structure, perhaps too few values.         '
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
    call DGEMM_('T','N',NumInt,NumInt,3*NumOfAt,One,SInv,3*NumOfAt,Temp,3*NumOfAt,Zero,Work(ipHess2),NumInt)
    call mma_deallocate(Temp)
    call mma_deallocate(Sinv)
    call mma_deallocate(Fcart)

  end if

  if (inpunit1 == 31) then
    close(inpunit1)
    close(inpunit2)
  end if

  ! FORCe: End ---------------------------------------------------------

  ! --------------------------------------------------------------------
  ! SCALe: Scale Hessian if scaling factors were given.

  call KeyWord(inpUnit,'SCAL',.true.,exist)
  if (exist) then
    call mma_allocate(ScaleParam,NumInt,label='ScaleParam')

    call KeyWord(inpUnit,'SECO',.false.,exist)
    do i=1,NumInt
      read(inpUnit,*) ScaleParam(i)
      ScaleParam(i) = sqrt(ScaleParam(i))
    end do
    do j=1,NumInt
      do i=1,NumInt
        Work(ipHess2+i+l_Hess1*(j-1)-1) = Work(ipHess2+i+l_Hess1*(j-1)-1)*ScaleParam(i)*ScaleParam(j)
      end do
    end do
    call mma_deallocate(ScaleParam)
  end if

  ! SCALe: End ---------------------------------------------------------

  ! --------------------------------------------------------------------
  ! DIPOles: Read Transition Dipoles or Spin-Orbit Coupling.

  call KeyWord(inpUnit,'DIPO',.true.,exist)
  if (.not. exist) call KeyWord(inpUnit,'SOC ',.true.,exist)
  if (.not. exist) then
    write(u6,*)
    write(u6,*) ' ************ ERROR ****************'
    write(u6,*) ' The DIPOLES/SOC keyword is missing.'
    write(u6,*) ' ***********************************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read dipoles ("DIPOLES").'
  read(inpUnit,format) Inline
  call Normalize(InLine,OutLine)
  if (OutLine(1:4) == 'FILE') then
    !D write(u6,*) ' Read trans dips from file UNSYM21'
    call molcas_open(31,'UNSYM21')
    !open(unit=31,file='UNSYM21')
    if (max_dip == 1) then
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT X',.false.,exist)
      read(31,*) (Work(ipTranDipGrad+3*(i-1)),i=1,3*NumOfAt)
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Y',.false.,exist)
      read(31,*) (Work(ipTranDipGrad+1+3*(i-1)),i=1,3*NumOfAt)
      call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Z',.false.,exist)
      read(31,*) (Work(ipTranDipGrad+2+3*(i-1)),i=1,3*NumOfAt)
    end if
    XnotFound = .true.
    YnotFound = .true.
    ZnotFound = .true.
    call KeyWord(31,'*BEGIN TRANSITION PROPERTIES',.false.,exist)
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
    call KeyWord(inpUnit,'DIPO',.true.,exist)
    !D write(u6,*) ' READINP: Read dipoles ("DIPOLES").'
    read(inpUnit,*) (TranDip(i),i=1,3)
    if (max_dip == 1) then
      read(inpUnit,*) (Work(ipTranDipGrad+3*(i-1)),i=1,3*NumOfAt)
      read(inpUnit,*) (Work(ipTranDipGrad+1+3*(i-1)),i=1,3*NumOfAt)
      read(inpUnit,*) (Work(ipTranDipGrad+2+3*(i-1)),i=1,3*NumOfAt)
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
  call KeyWord(inpUnit,'NONL',.true.,exist)
  if (.not. Exist) then
    write(u6,*)
    write(u6,*) ' *********** ERROR **************'
    write(u6,*) ' Cannot find keyword  NONLINEAR !'
    write(u6,*) ' ********************************'
    call Quit_OnUserError()
  end if
  do icoord=1,NumInt
    read(inpUnit,format) InLine
    call Normalize(InLine,OutLine)
    trfName1(icoord) = OutLine
  end do
  do icoord=1,NumInt
    read(inpUnit,format) InLine
    call Normalize(InLine,OutLine)
    trfName2(icoord) = OutLine
  end do

  ! Read terms in polynomial.
  call KeyWord(inpUnit,'POLY',.true.,exist)
  if (.not. Exist) then
    write(u6,*)
    write(u6,*) ' ************ ERROR **************'
    write(u6,*) ' Cannot find keyword  POLYNOMIAL !'
    write(u6,*) ' *********************************'
    call Quit_OnUserError()
  end if
  !D write(u6,*) ' READINP: Read polynomial.'
  k = 1
  nvar = 0
  read(inpUnit,format) InLine
  call WordPos(k,InLine,iStart,iStop)
  do while (iStop < len(Inline))
    k = iStop+1
    nvar = nvar+1
    call WordPos(k,InLine,iStart,iStop)
  end do
  nPolyTerm = 1
  read(inpUnit,format) InLine
  call Normalize(InLine,OutLine)
  do while (OutLine(1:4) /= 'END ')
    nPolyTerm = nPolyTerm+1
    read(inpUnit,format) InLine
    call Normalize(InLine,OutLine)
  end do
  call KeyWord(inpUnit,'POLY',.true.,exist)
  if (.not. Exist) then
    write(u6,*)
    write(u6,*) ' ************ ERROR **************'
    write(u6,*) ' Cannot find keyword  POLYNOMIAL !'
    write(u6,*) ' *********************************'
    call Quit_OnUserError()
  end if

  call GetMem('ipow','Allo','Inte',ipipow,nPolyTerm*nvar)
  do iterm=1,nPolyTerm
    read(inpUnit,*) (iWork(ipipow+iterm+nPolyTerm*(ivar-1)-1),ivar=1,nvar)
  end do

  ! Find highest power of term in polynomial.
  max_term = 0
  do iterm=1,nPolyTerm
    nsum = 0
    do jvar=1,nvar
      nsum = nsum+iWork(ipipow+iterm+nPolyTerm*(jvar-1)-1)
      if (nsum > max_term) max_term = nsum
    end do
  end do

  ! Read grid points and energies.
  CoordType = 'INTERNAL'
  call KeyWord(inpUnit,'DATA',.true.,exist)
  if (.not. Exist) then
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
  call GetMem('var','Allo','Real',ipvar,ndata*nvar)
  do idata=1,ndata
    do jvar=1,nvar
      Work(ipvar+idata+ndata*(jvar-1)-1) = coord(idata,jvar)
    end do
  end do

  ! CASPT2 energies for ground state.
  call GetMem('yin1','Allo','Real',ipyin1,ndata)

  do idata=1,ndata
    Work(ipyin1+idata-1) = GrdVal(idata,2)
  end do

  ! CASPT2 energies for excited state.
  call GetMem('yin2','Allo','Real',ipyin2,ndata)
  do idata=1,ndata
    Work(ipyin2+idata-1) = GrdVal(idata,6)
  end do

  ! If three atomic molecule, add Watson correction to energy values.
  if (NumOfAt == 3) then
    m1 = Mass(1)*uToAu
    m2 = Mass(2)*uToAu
    m3 = Mass(3)*uToAu
    m12 = (m1*m2)/(m1+m2)
    m23 = (m2*m3)/(m2+m3)
    do i=1,ndata
      r1 = Work(ipvar+i-1)
      r2 = Work(ipvar+i+ndata-1)
      theta = Work(ipvar+i+ndata*2-1)*deg2rad
      const = -(One/Eight)*(One/(m12*r1**2)+One/(m23*r2**2))*(One+(One/(sin(theta)**2)))
      const = const+Quart*(cos(theta)/(m2*r1*r2))*(cos(theta)/sin(theta))**2
      Work(ipyin1+i-1) = Work(ipyin1+i-1)+const
      Work(ipyin2+i-1) = Work(ipyin2+i-1)+const
    end do
  end if

  call GetMem('t_dipin1','Allo','Real',ipt_dipin1,ndata)
  do idata=1,ndata
    Work(ipt_dipin1+idata-1) = GrdVal(idata,9)
  end do

  call GetMem('t_dipin2','Allo','Real',ipt_dipin2,ndata)
  do idata=1,ndata
    Work(ipt_dipin2+idata-1) = GrdVal(idata,10)
  end do
  call GetMem('t_dipin3','Allo','Real',ipt_dipin3,ndata)
  call dcopy_(ndata,[Zero],0,Work(ipt_dipin3),1)
  !t_dipin3 = Zero

  ! Make sure that all transition dipole values have the same sign.
  sign1_1 = Work(ipt_dipin1)/abs(Work(ipt_dipin1))
  sign2_1 = Work(ipt_dipin2)/abs(Work(ipt_dipin2))
  do idata=2,ndata
    sign1_2 = Work(ipt_dipin1+idata-1)/abs(Work(ipt_dipin1+idata-1))
    sign2_2 = Work(ipt_dipin2+idata-1)/abs(Work(ipt_dipin2+idata-1))
    if (int(sign1_1*sign2_1) /= int(sign1_2*sign2_2)) then
      write(u6,*) 'Sign shift in transition dipole data:',idata
    end if
    Work(ipt_dipin1+idata-1) = sign1_1*sign1_2*Work(ipt_dipin1+idata-1)
    Work(ipt_dipin2+idata-1) = sign2_1*sign2_2*Work(ipt_dipin2+idata-1)
  end do
  !D write(u6,*) ' READINP: Transition dipoles finished.'

  !                   End of Energy surface input part
  ! --------------------------------------------------------------------

end if

!D write(u6,*) ' Ending READINP.'

end subroutine ReadInp

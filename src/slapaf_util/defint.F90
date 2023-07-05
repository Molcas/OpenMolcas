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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine DefInt(nBVct,BMtrx,nQQ,nAtom,rInt,Lbl,Coor,nDim)
!***********************************************************************
!                                                                      *
! Object: to generate the B matrix which is the transformation matrix  *
!         between an infinitesimal displacement in the symmetry adapted*
!         internal coordinates and the symmetry unique cartesian       *
!         coordinates.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May 1991                                                 *
!***********************************************************************

use Slapaf_Info, only: AtomLbl
use Slapaf_Parameters, only: iRow, Redundant

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
real*8 BMtrx(3*nAtom,nQQ), rInt(nQQ), Coor(3,nAtom)
character type*6, Temp*120, Lbl(nQQ)*8, cNum*4, Line*120, format*8, filnam*16
logical Flip, lPIC(6*nAtom), lAtom(nAtom), Found
logical, save :: First = .true.
logical :: lWrite = .false., lErr = .false.
integer, allocatable :: Ind(:)
real*8, allocatable :: xyz(:), Tmp2(:), Mass(:), TM(:)
real*8, allocatable :: BVct(:,:), value(:), rMult(:)
character(len=8), allocatable :: Labels(:)

call mma_allocate(BVct,3*nAtom,nBVct,Label='BVct')
call mma_allocate(value,nBVct,Label='Value')
call mma_allocate(rMult,nBVct,Label='rMult')
call mma_allocate(Labels,nBVct,Label='Labels')
BVct(:,:) = Zero
value(:) = Zero

iRout = 30
iPrint = nPrint(iRout)

if (iPrint >= 6) lWrite = .true.
do i=1,6*nAtom
  lPIC(i) = .true.
end do
do i=1,nAtom
  lAtom(i) = .true.
end do
do jBVct=1,nBVct
  Labels(jBVct) = ' '
end do

!Lu = 6
nTemp = len(Temp)
write(format,'(A,I3.3,A)') '(F',nTemp,'.0)'

Lu_UDIC = 91
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
!open(Lu_UDIC,File=filnam,Form='Formatted',Status='OLD')
rewind(Lu_UDIC)

call dcopy_(nBVct,[Zero],0,rMult,1)
if (iPrint == 99) First = .true.
if (lWrite) then
  write(6,*)
  write(6,'(80A)') ('*',i=1,60)
  write(6,*) ' User-Defined Internal coordinates'
  write(6,'(80A)') ('-',i=1,60)
  do iLines=1,iRow
    read(Lu_UDIC,'(A)') Temp
    write(6,'(A)') Temp
  end do
  write(6,'(80A)') ('*',i=1,60)
  write(6,*)
  write(6,*)
  write(6,*) '*************************************************************'
  write(6,*) '* Values of primitive internal coordinates                  *'
  write(6,*) '-------------------------------------------------------------'
  rewind(Lu_UDIC)
end if

! Step 1. Set up the b vectors from which we will define the
! internal coordinates.

iBVct = 0
Found = .false.
do iLines=1,iRow
  Flip = .false.
  read(Lu_UDIC,'(A)') Line
  Temp = Line
  call UpCase(Temp)

  if (Temp(1:4) == 'VARY') then
    Found = .true.
    exit
  end if
  iBVct = iBVct+1

  ! Move the label of the internal coordinate

  neq = index(Line,'=')
  if (neq == 0) then
    call WarningMessage(2,'Error in DefInt')
    write(6,*)
    write(6,'(A)') '***********************************'
    write(6,'(A)') ' Syntax error in line :            '
    write(6,'(A)') Line(1:33),'...'
    write(6,'(A)') '***********************************'
    call Quit_OnUserError()
  else
    iFrst = 1
    call NxtWrd(Line,iFrst,iEnd)
    jEnd = iEnd
    if (Line(iEnd:iEnd) == '=') jEnd = jEnd-1
    if (jEnd-iFrst+1 > 8) then
      call WarningMessage(2,'Error in DefInt')
      write(6,*)
      write(6,'(A)') '***********************************'
      write(6,'(A)') ' Syntax error in line :            '
      write(6,'(A)') Line(1:33),'...'
      write(6,'(A,A)') Line(iFrst:jEnd),' has more than 8 character'
      write(6,'(A)') '***********************************'
      call Quit_OnUserError()
    end if
    call ChkLbl(Line(iFrst:jEnd),Labels,iBVct-1)
    Labels(iBVct) = Line(iFrst:jEnd)
  end if

  ! Construct the corresponding transformation vector

  mCntr = 0
  if (index(Temp,'CART') /= 0) then
    nCntr = 1
    nGo = index(Temp,'CART')
    call NxtWrd(Temp,nGo,nGo2)
    nGo = nGo2+1
    call NxtWrd(Temp,nGo,nGo2)
    if (Temp(nGo:nGo2) == 'X') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      type = 'X     '
    else if (Temp(nGo:nGo2) == 'Y') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      type = 'Y     '
    else if (Temp(nGo:nGo2) == 'Z') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      type = 'Z     '
    else
      nGo = -1
      call WarningMessage(2,'Error in DefInt')
      write(6,*)
      write(6,*) '*********** ERROR ************'
      write(6,*) ' DefInt: wrong cartesian type '
      write(6,'(A,A)') ' Temp=',Temp
      write(6,*) '******************************'
      call Quit_OnUserError()
    end if
  else if (index(Temp,'BOND') /= 0) then
    nCntr = 2
    nGo = index(Temp,'BOND')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'STRTCH'
  else if (index(Temp,'LANGLE(2)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(2)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'LBEND2'
  else if (index(Temp,'LANGLE(1)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(1)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'LBEND1'
  else if (index(Temp,'ANGL') /= 0) then
    nCntr = 3
    nGo = index(Temp,'ANGL')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'BEND  '
  else if (index(Temp,'DIHE') /= 0) then
    nCntr = 4
    nGo = index(Temp,'DIHE')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'TRSN  '
    Flip = .not. First
  else if (index(Temp,'OUTO') /= 0) then
    nCntr = 4
    nGo = index(Temp,'OUTO')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'OUTOFP'
  else if (index(Temp,'DISS') /= 0) then
    i1 = index(Line,'(')
    i2 = index(Line,'+')
    i3 = index(Line,')')
    if ((i1 >= i2) .or. (i2 >= i3)) then
      call WarningMessage(2,'Error in DefInt')
      write(6,*)
      write(6,*) '********** ERROR ************'
      write(6,*) ' Line contains syntax error !'
      write(6,'(A)') Line
      write(6,*) i1,i2,i3
      write(6,*) '*****************************'
      call Quit_OnUserError()
    end if
    nGo = i3+1
    Temp = Line(i1+1:i2-1)
    read(Temp,format) Tmp
    nCntr = nint(Tmp)
    Temp = Line(i2+1:i3-1)
    read(Temp,format) Tmp
    mCntr = nint(Tmp)
    type = 'DISSOC'
  else
    nGo = -1
    call WarningMessage(2,'Error in DefInt')
    write(6,*)
    write(6,*) '*********** ERROR ***********'
    write(6,*) ' Line contains syntax error !'
    write(6,'(A)') Line
    write(6,*) '*****************************'
    call Quit_OnUserError()
  end if

  msAtom = nCntr+mCntr
  call mma_allocate(xyz,3*msAtom,Label='xyz')
  call mma_allocate(Tmp2,3*msAtom,Label='Tmp2')
  call mma_allocate(Ind,2*msAtom,Label='Ind')
  call mma_allocate(Mass,2*msAtom,Label='Mass')
  call mma_allocate(TM,9*nAtom*(nCntr+mCntr),Label='TM')

  call Cllct(Line(nGo:nTemp),BVct(1,iBVct),Value_Temp,nAtom,Coor,nCntr,mCntr,xyz,Tmp2,Ind,type,Mass,TM,Labels(iBVct),lWrite, &
             rMult(iBVct),lAtom)

  if ((.not. First) .and. (type == 'TRSN') .and. (abs(Value_Temp) < Pi*Half)) Flip = .false.
  if (Flip .and. (value(iBVct)*Value_Temp < Zero)) then
    !write(Lu,*) 'Flip Sign for ',Labels(iBVct)
    if (value(iBVct) < Zero) then
      value(iBVct) = -Pi-(Pi-Value_Temp)
    else
      value(iBVct) = Pi+(Pi+Value_Temp)
    end if

  else
    value(iBVct) = Value_Temp
  end if

  call mma_deallocate(TM)
  call mma_deallocate(Mass)
  call mma_deallocate(Ind)
  call mma_deallocate(Tmp2)
  call mma_deallocate(xyz)
end do

if (.not. Found) then
  call WarningMessage(2,'Error in DefInt')
  write(6,*) '**********************************************'
  write(6,*) ' ERROR: No internal coordinates are defined ! '
  write(6,*) '**********************************************'
  call Quit_OnUserError()
end if

nDefPICs = iBVct
if (iPrint >= 59) call RecPrt(' The B-vectors',' ',BVct,3*nAtom,nBVct)
if (iPrint >= 19) call RecPrt(' Values of primitive internal coordinates / au or rad',' ',value,nBVct,1)

! Step 2. Define internal coordinates as linear combinations of
! the previously defined primitive internal coordinates.

iBMtrx = 0
jLines = iLines
do
  jLines = jLines+1
  if (jLines > iRow) exit

  ! Jump over the FIX keyword if present

  read(Lu_UDIC,'(A)') Line
  Temp = Line
  call UpCase(Temp)
  if (Temp(1:3) == 'FIX') cycle
  if (Temp(1:4) == 'ROWH') exit

  iBMtrx = iBMtrx+1
  rInt(iBMtrx) = Zero
  RR = Zero

  iFrst = 1
  call NxtWrd(Line,iFrst,iEnd)
  jEnd = iEnd
  if (Line(iEnd:iEnd) == '=') jEnd = iEnd-1
  Lbl(iBMtrx) = Line(iFrst:jEnd)
  neq = index(Line,'=')
  if (neq == 0) then

    ! a single vector (this will only extend over one line)

    iBVct = 0
    do jBVct=1,nBVct
      if (Line(iFrst:jEnd) == Labels(jBVct)) iBVct = jBVct
    end do
    if (iBVct == 0) then
      call WarningMessage(2,'Error in DefInt')
      write(6,*)
      write(6,*) '*******************************'
      write(6,*) ' ERROR: A single vector        '
      write(6,*) ' Undefined internal coordinate '
      write(6,'(A,A)') Line
      write(6,'(A,A)') Line(iFrst:jEnd)
      write(6,*) '*******************************'
      call Quit_OnUserError()
    end if

    lPIC(iBVct) = .false.
    call dcopy_(3*nAtom,BVct(1,iBVct),1,BMtrx(1,iBMtrx),1)
    call DScal_(3*nAtom,rMult(iBVct)**2,BMtrx(1,iBMtrx),1)
    rInt(iBMtrx) = rMult(iBVct)**2*value(iBVct)
    RR = RR+rMult(iBVct)**2

  else

    ! A linear combination of vectors

    call dcopy_(3*nAtom,[Zero],0,BMtrx(1,iBMtrx),1)
    iFrst = neq+1
    Sgn = One

    ! Process the rest of the line and possible extension lines

    do
      ! Get the factor
      call NxtWrd(Line,iFrst,iEnd)
      Temp = Line(iFrst:iEnd)
      read(Temp,format) Fact
      Fact = Fact*Sgn
      iFrst = iEnd+1
      ! Get the label
      call NxtWrd(Line,iFrst,iEnd)
      iBVct = 0
      do jBVct=1,nBVct
        if (Line(iFrst:iEnd) == Labels(jBVct)) iBVct = jBVct
      end do
      if (iBVct == 0) then
        call WarningMessage(2,'Error in DefInt')
        write(6,*)
        write(6,*) '************ ERROR *************'
        write(6,*) ' Linear combinations of vectors '
        write(6,*) ' Undefined internal coordinate  '
        write(6,'(A,A)') Line
        write(6,'(A,A)') Line(iFrst:iEnd)
        write(6,*) '********************************'
        call Quit_OnUserError()
      end if

      lPIC(iBVct) = .false.
      call DaXpY_(3*nAtom,Fact*rMult(iBVct)**2,BVct(1,iBvct),1,BMtrx(1,iBMtrx),1)
      rInt(iBMtrx) = rInt(iBMtrx)+Fact*rMult(iBVct)**2*value(iBVct)
      RR = RR+rMult(iBVct)**2*Fact**2

      iFrst = iEnd+1
      Temp = Line(iFrst:nTemp)
      nPlus = index(Temp,'+')
      nMinus = index(Temp,'-')
      if ((nPlus /= 0) .and. ((nPlus < nMinus) .eqv. (nMinus > 0))) then
        Sgn = One
        iFrst = iFrst+nPlus
        cycle
      end if
      if ((nMinus /= 0) .and. ((nMinus < nPlus) .eqv. (nPlus > 0))) then
        Sgn = -One
        iFrst = iFrst+nMinus
        cycle
      end if

      ! Here if all statements processed of this line

      if (index(Line,'&') == 0) exit
      jLines = jLines+1
      if (jLines > iRow) then
        call WarningMessage(2,'Error in DefInt')
        write(6,*)
        write(6,*) '********** ERROR *********'
        write(6,*) ' DefInt: jLines > iRow '
        write(6,*) '**************************'
        call Quit_OnUserError()
      end if
      read(Lu_UDIC,'(A)') Line
      iFrst = 1
      call NxtWrd(Line,iFrst,iEnd)
      if (Line(iFrst:iEnd) == '+') then
        iFrst = iEnd+1
        Sgn = One
      else if (Line(iFrst:iEnd) == '-') then
        iFrst = iEnd+1
        Sgn = -One
      else
        call WarningMessage(2,'Error in DefInt')
        write(6,*)
        write(6,*) '************** ERROR *************'
        write(6,*) ' Syntax Error: first character in  extension line is not + or -'
        write(6,'(A)') Line
        write(6,'(3A)') '-->',Line(iFrst:iEnd),'<--'
        write(6,*) '**********************************'
        call Quit_OnUserError()
      end if
    end do

  end if

  rInt(iBMtrx) = rInt(iBMtrx)/sqrt(RR)
  call DScal_(3*nAtom,One/sqrt(RR),BMtrx(1,iBMtrx),1)

end do

! Skip the  RowH  section of input ---

do
  jLines = jLines+1
  if (jLines > iRow) exit
  read(Lu_UDIC,'(A)') Line
  Temp = Line
  call UpCase(Temp)
end do

if (iPrint >= 99) call RecPrt(' The B-matrix',' ',BMtrx,3*nAtom,nQQ)

if (lWrite) then
  write(6,*)
  write(6,*)
  write(6,*) '*********************************************'
  write(6,*) '* Value of internal coordinates / au or rad *'
  write(6,*) '---------------------------------------------'
  write(6,'(1X,A,2X,F10.4)') (Lbl(iInt),rInt(iInt),iInt=1,nQQ)
  write(6,*)
  call XFlush(6)
end if

! Some checks: Warnings & Errors ---

iDiff = 3*nAtom-nDim
if (iDiff == 0) then
  cNum = '3N'
else
  write(cNum,'(A,I1)') '3N-',iDiff
end if

if (iBMtrx < nDim) then
  write(6,*) '**************** ERROR **********************'
  write(6,*) ' N.r of Internal Coordinates lower than ',cNum
  write(6,*) '*********************************************'
  write(6,*)
  lErr = .true.
end if

if ((iBMtrx > nDim) .and. (.not. Redundant)) then
  write(6,*) '***************** ERROR ***********************'
  write(6,*) ' N.r of Internal Coordinates greater than ',cNum
  write(6,*) '***********************************************'
  write(6,*)
  lErr = .true.
end if

if (nDefPICs < iBMtrx) then
  write(6,*)
  write(6,*) '****************** ERROR ******************'
  write(6,*) ' Too many Internal Coordinates !           '
  write(6,'(A,I4)') ' N.r of Primitive Internal Coordinates:',nDefPICs
  write(6,'(A,I4)') ' N.r of Internal Coordinates:          ',iBMtrx
  write(6,*) '*******************************************'
  write(6,*)
  lErr = .true.
end if

do i=1,nDefPICs
  if (lPIC(i)) then
    write(6,*)
    write(6,*) '*************** ERROR *******************'
    write(6,*) ' Primitive Internal Coordinate not used: '
    write(6,*) ' ',Labels(i)
    write(6,*) '*****************************************'
    write(6,*)
    lErr = .true.
  end if
end do

do i=1,nAtom
  if (lAtom(i)) then
    call WarningMessage(2,'Error in DefInt')
    write(6,*)
    write(6,*) '*********** ERROR ****************'
    write(6,*) ' No Coordinate defined for atom:  '
    write(6,*) ' ',AtomLbl(i)
    write(6,*) '**********************************'
    write(6,*)
    lErr = .true.
  end if
end do

if (lErr) call Quit_OnUserError()

if (iBMtrx /= nQQ) then
  call WarningMessage(2,'Error in DefInt')
  write(6,*)
  write(6,*) '******************* ERROR *****************'
  if (iBMtrx > nQQ) write(6,*) ' Too many internal coordinates are defined'
  if (iBMtrx < nQQ) write(6,*) ' Too few internal coordinates are defined'
  write(6,*) ' You have defined',iBMtrx
  write(6,*) ' There should be ',nQQ
  write(6,*) '*******************************************'
  call Quit_OnUserError()
end if

close(Lu_UDIC)
First = .false.

call mma_deallocate(Labels)
call mma_deallocate(rMult)
call mma_deallocate(value)
call mma_deallocate(BVct)

return

end subroutine DefInt

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

use Slapaf_Info, only: AtomLbl, iRow, Redundant
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBVct, nQQ, nAtom, nDim
real(kind=wp), intent(out) :: BMtrx(3*nAtom,nQQ), rInt(nQQ)
character(len=8), intent(out) :: Lbl(nQQ)
real(kind=wp), intent(in) :: Coor(3,nAtom)
#include "print.fh"
integer(kind=iwp) :: i, i1, i2, i3, iBMtrx, iBVct, iDiff, iEnd, iFrst, iInt, iLines, iPrint, iRout, jBVct, jEnd, jLines, Lu_UDIC, &
                     mCntr, msAtom, nCntr, nDefPICs, neq, nGo, nGo2, nMinus, nPlus, nTemp
real(kind=wp) :: Fact, RR, Sgn, Tmp, Value_Temp
logical(kind=iwp) :: First = .true., Flip, Found, lErr, lWrite
character(len=120) :: Line, Temp
character(len=16) :: filnam
character(len=8) :: Frmt
character(len=6) :: Typ
character(len=4) :: cNum
integer(kind=iwp), allocatable :: Ind(:)
real(kind=wp), allocatable :: BVct(:,:), Mass(:), rMult(:), TM(:), Tmp2(:), Val(:), xyz(:)
logical(kind=iwp), allocatable :: lAtom(:), lPIC(:)
character(len=8), allocatable :: Labels(:)

call mma_allocate(BVct,3*nAtom,nBVct,Label='BVct')
call mma_allocate(Val,nBVct,Label='Val')
call mma_allocate(rMult,nBVct,Label='rMult')
call mma_allocate(Labels,nBVct,Label='Labels')
BVct(:,:) = Zero
Val(:) = Zero
lWrite = .false.
lErr = .false.

iRout = 30
iPrint = nPrint(iRout)

if (iPrint >= 6) lWrite = .true.
call mma_allocate(lAtom,nAtom,Label='lAtom')
call mma_allocate(lPIC,6*nAtom,Label='lPIC')
lPIC(:) = .true.
lAtom(:) = .true.
Labels(:) = ''

!Lu = u6
nTemp = len(Temp)
write(Frmt,'(A,I3.3,A)') '(F',nTemp,'.0)'

Lu_UDIC = 91
filnam = 'UDIC'
call molcas_open(Lu_UDIC,filnam)
!open(Lu_UDIC,File=filnam,Form='Formatted',Status='OLD')
rewind(Lu_UDIC)

rMult(:) = Zero
if (iPrint == 99) First = .true.
if (lWrite) then
  write(u6,*)
  write(u6,'(80A)') ('*',i=1,60)
  write(u6,*) ' User-Defined Internal coordinates'
  write(u6,'(80A)') ('-',i=1,60)
  do iLines=1,iRow
    read(Lu_UDIC,'(A)') Temp
    write(u6,'(A)') Temp
  end do
  write(u6,'(80A)') ('*',i=1,60)
  write(u6,*)
  write(u6,*)
  write(u6,*) '*************************************************************'
  write(u6,*) '* Values of primitive internal coordinates                  *'
  write(u6,*) '-------------------------------------------------------------'
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
    write(u6,*)
    write(u6,'(A)') '***********************************'
    write(u6,'(A)') ' Syntax error in line :            '
    write(u6,'(A)') Line(1:33),'...'
    write(u6,'(A)') '***********************************'
    call Quit_OnUserError()
  else
    iFrst = 1
    call NxtWrd(Line,iFrst,iEnd)
    jEnd = iEnd
    if (Line(iEnd:iEnd) == '=') jEnd = jEnd-1
    if (jEnd-iFrst+1 > 8) then
      call WarningMessage(2,'Error in DefInt')
      write(u6,*)
      write(u6,'(A)') '***********************************'
      write(u6,'(A)') ' Syntax error in line :            '
      write(u6,'(A)') Line(1:33),'...'
      write(u6,'(A,A)') Line(iFrst:jEnd),' has more than 8 character'
      write(u6,'(A)') '***********************************'
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
      Typ = 'X     '
    else if (Temp(nGo:nGo2) == 'Y') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      Typ = 'Y     '
    else if (Temp(nGo:nGo2) == 'Z') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      Typ = 'Z     '
    else
      nGo = -1
      call WarningMessage(2,'Error in DefInt')
      write(u6,*)
      write(u6,*) '*********** ERROR ************'
      write(u6,*) ' DefInt: wrong cartesian type '
      write(u6,'(A,A)') ' Temp=',Temp
      write(u6,*) '******************************'
      call Quit_OnUserError()
    end if
  else if (index(Temp,'BOND') /= 0) then
    nCntr = 2
    nGo = index(Temp,'BOND')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'STRTCH'
  else if (index(Temp,'LANGLE(2)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(2)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'LBEND2'
  else if (index(Temp,'LANGLE(1)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(1)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'LBEND1'
  else if (index(Temp,'ANGL') /= 0) then
    nCntr = 3
    nGo = index(Temp,'ANGL')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'BEND  '
  else if (index(Temp,'DIHE') /= 0) then
    nCntr = 4
    nGo = index(Temp,'DIHE')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'TRSN  '
    Flip = .not. First
  else if (index(Temp,'OUTO') /= 0) then
    nCntr = 4
    nGo = index(Temp,'OUTO')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    Typ = 'OUTOFP'
  else if (index(Temp,'DISS') /= 0) then
    i1 = index(Line,'(')
    i2 = index(Line,'+')
    i3 = index(Line,')')
    if ((i1 >= i2) .or. (i2 >= i3)) then
      call WarningMessage(2,'Error in DefInt')
      write(u6,*)
      write(u6,*) '********** ERROR ************'
      write(u6,*) ' Line contains syntax error !'
      write(u6,'(A)') Line
      write(u6,*) i1,i2,i3
      write(u6,*) '*****************************'
      call Quit_OnUserError()
    end if
    nGo = i3+1
    Temp = Line(i1+1:i2-1)
    read(Temp,Frmt) Tmp
    nCntr = nint(Tmp)
    Temp = Line(i2+1:i3-1)
    read(Temp,Frmt) Tmp
    mCntr = nint(Tmp)
    Typ = 'DISSOC'
  else
    nGo = -1
    call WarningMessage(2,'Error in DefInt')
    write(u6,*)
    write(u6,*) '*********** ERROR ***********'
    write(u6,*) ' Line contains syntax error !'
    write(u6,'(A)') Line
    write(u6,*) '*****************************'
    call Quit_OnUserError()
  end if

  msAtom = nCntr+mCntr
  call mma_allocate(xyz,3*msAtom,Label='xyz')
  call mma_allocate(Tmp2,3*msAtom,Label='Tmp2')
  call mma_allocate(Ind,2*msAtom,Label='Ind')
  call mma_allocate(Mass,2*msAtom,Label='Mass')
  call mma_allocate(TM,9*nAtom*(nCntr+mCntr),Label='TM')

  call Cllct(Line(nGo:nTemp),BVct(1,iBVct),Value_Temp,nAtom,Coor,nCntr,mCntr,xyz,Tmp2,Ind,Typ,Mass,TM,Labels(iBVct),lWrite, &
             rMult(iBVct),lAtom)

  if ((.not. First) .and. (Typ == 'TRSN') .and. (abs(Value_Temp) < Pi*Half)) Flip = .false.
  if (Flip .and. (Val(iBVct)*Value_Temp < Zero)) then
    !write(Lu,*) 'Flip Sign for ',Labels(iBVct)
    if (Val(iBVct) < Zero) then
      Val(iBVct) = -Pi-(Pi-Value_Temp)
    else
      Val(iBVct) = Pi+(Pi+Value_Temp)
    end if

  else
    Val(iBVct) = Value_Temp
  end if

  call mma_deallocate(TM)
  call mma_deallocate(Mass)
  call mma_deallocate(Ind)
  call mma_deallocate(Tmp2)
  call mma_deallocate(xyz)
end do

if (.not. Found) then
  call WarningMessage(2,'Error in DefInt')
  write(u6,*) '**********************************************'
  write(u6,*) ' ERROR: No internal coordinates are defined ! '
  write(u6,*) '**********************************************'
  call Quit_OnUserError()
end if

nDefPICs = iBVct
if (iPrint >= 59) call RecPrt(' The B-vectors',' ',BVct,3*nAtom,nBVct)
if (iPrint >= 19) call RecPrt(' Values of primitive internal coordinates / au or rad',' ',Val,nBVct,1)

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
      write(u6,*)
      write(u6,*) '*******************************'
      write(u6,*) ' ERROR: A single vector        '
      write(u6,*) ' Undefined internal coordinate '
      write(u6,'(A,A)') Line
      write(u6,'(A,A)') Line(iFrst:jEnd)
      write(u6,*) '*******************************'
      call Quit_OnUserError()
    end if

    lPIC(iBVct) = .false.
    BMtrx(:,iBMtrx) = rMult(iBVct)**2*BVct(:,iBVct)
    rInt(iBMtrx) = rMult(iBVct)**2*Val(iBVct)
    RR = RR+rMult(iBVct)**2

  else

    ! A linear combination of vectors

    BMtrx(:,iBMtrx) = Zero
    iFrst = neq+1
    Sgn = One

    ! Process the rest of the line and possible extension lines

    do
      ! Get the factor
      call NxtWrd(Line,iFrst,iEnd)
      Temp = Line(iFrst:iEnd)
      read(Temp,Frmt) Fact
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
        write(u6,*)
        write(u6,*) '************ ERROR *************'
        write(u6,*) ' Linear combinations of vectors '
        write(u6,*) ' Undefined internal coordinate  '
        write(u6,'(A,A)') Line
        write(u6,'(A,A)') Line(iFrst:iEnd)
        write(u6,*) '********************************'
        call Quit_OnUserError()
      end if

      lPIC(iBVct) = .false.
      BMtrx(:,iBMtrx) = BMtrx(:,iBMtrx)+Fact*rMult(iBVct)**2*BVct(:,iBVct)
      rInt(iBMtrx) = rInt(iBMtrx)+Fact*rMult(iBVct)**2*Val(iBVct)
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
        write(u6,*)
        write(u6,*) '********** ERROR *********'
        write(u6,*) ' DefInt: jLines > iRow '
        write(u6,*) '**************************'
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
        write(u6,*)
        write(u6,*) '************** ERROR *************'
        write(u6,*) ' Syntax Error: first character in  extension line is not + or -'
        write(u6,'(A)') Line
        write(u6,'(3A)') '-->',Line(iFrst:iEnd),'<--'
        write(u6,*) '**********************************'
        call Quit_OnUserError()
      end if
    end do

  end if

  rInt(iBMtrx) = rInt(iBMtrx)/sqrt(RR)
  BMtrx(:,iBMtrx) = BMtrx(:,iBMtrx)/sqrt(RR)

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
  write(u6,*)
  write(u6,*)
  write(u6,*) '*********************************************'
  write(u6,*) '* Value of internal coordinates / au or rad *'
  write(u6,*) '---------------------------------------------'
  write(u6,'(1X,A,2X,F10.4)') (Lbl(iInt),rInt(iInt),iInt=1,nQQ)
  write(u6,*)
  call XFlush(u6)
end if

! Some checks: Warnings & Errors ---

iDiff = 3*nAtom-nDim
if (iDiff == 0) then
  cNum = '3N'
else
  write(cNum,'(A,I1)') '3N-',iDiff
end if

if (iBMtrx < nDim) then
  write(u6,*) '**************** ERROR **********************'
  write(u6,*) ' N.r of Internal Coordinates lower than ',cNum
  write(u6,*) '*********************************************'
  write(u6,*)
  lErr = .true.
end if

if ((iBMtrx > nDim) .and. (.not. Redundant)) then
  write(u6,*) '***************** ERROR ***********************'
  write(u6,*) ' N.r of Internal Coordinates greater than ',cNum
  write(u6,*) '***********************************************'
  write(u6,*)
  lErr = .true.
end if

if (nDefPICs < iBMtrx) then
  write(u6,*)
  write(u6,*) '****************** ERROR ******************'
  write(u6,*) ' Too many Internal Coordinates !           '
  write(u6,'(A,I4)') ' N.r of Primitive Internal Coordinates:',nDefPICs
  write(u6,'(A,I4)') ' N.r of Internal Coordinates:          ',iBMtrx
  write(u6,*) '*******************************************'
  write(u6,*)
  lErr = .true.
end if

do i=1,nDefPICs
  if (lPIC(i)) then
    write(u6,*)
    write(u6,*) '*************** ERROR *******************'
    write(u6,*) ' Primitive Internal Coordinate not used: '
    write(u6,*) ' ',Labels(i)
    write(u6,*) '*****************************************'
    write(u6,*)
    lErr = .true.
  end if
end do

do i=1,nAtom
  if (lAtom(i)) then
    call WarningMessage(2,'Error in DefInt')
    write(u6,*)
    write(u6,*) '*********** ERROR ****************'
    write(u6,*) ' No Coordinate defined for atom:  '
    write(u6,*) ' ',AtomLbl(i)
    write(u6,*) '**********************************'
    write(u6,*)
    lErr = .true.
  end if
end do

call mma_deallocate(lAtom)
call mma_deallocate(lPIC)

if (lErr) call Quit_OnUserError()

if (iBMtrx /= nQQ) then
  call WarningMessage(2,'Error in DefInt')
  write(u6,*)
  write(u6,*) '******************* ERROR *****************'
  if (iBMtrx > nQQ) write(u6,*) ' Too many internal coordinates are defined'
  if (iBMtrx < nQQ) write(u6,*) ' Too few internal coordinates are defined'
  write(u6,*) ' You have defined',iBMtrx
  write(u6,*) ' There should be ',nQQ
  write(u6,*) '*******************************************'
  call Quit_OnUserError()
end if

close(Lu_UDIC)
First = .false.

call mma_deallocate(Labels)
call mma_deallocate(rMult)
call mma_deallocate(Val)
call mma_deallocate(BVct)

return

end subroutine DefInt

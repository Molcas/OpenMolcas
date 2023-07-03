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
! Copyright (C) 1991,1997, Roland Lindh                                *
!***********************************************************************

subroutine DefInt2(BVct,dBVct,nBvct,BMtrx,mInt,nAtom,nLines,value,rInt,rInt0,Lbl,lWrite,rMult,dBMtrx,Value0,lIter,iFlip)
!***********************************************************************
!                                                                      *
! Object: to generate the B matrix for the constraints                 *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             May '91                                                  *
!                                                                      *
!             Modified to constraints, June '97 by R. Lindh            *
!***********************************************************************

use UnixInfo, only: SuperName

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
real*8 BVct(3*nAtom,nBVct), dBVct(3*nAtom,3*nAtom,nBVct), value(nBVct), BMtrx(3*nAtom,mInt), rInt(mInt), rInt0(mInt), &
       rMult(nBVct,nBVct), dBMtrx(3*nAtom,3*nAtom,mInt), Value0(nBVct), MaxErr
character Line*120, type*6, format*8, Temp*120, Lbl(mInt)*8, filnam*16
logical lWrite, Start, rInt0_on_file, rInt0_in_memory, InSlapaf
integer, parameter :: Flip = 1, NoFlip = 0
integer, external :: StrnLn
integer iFlip(nBVct)
integer, allocatable :: Ind(:), tpc(:)
real*8, allocatable :: Hess(:), Mass(:), Grad(:), xyz(:), r0(:)
character(len=8), allocatable :: Labels(:)
#include "angstr.fh"

call mma_allocate(Labels,nBVct,Label='Labels')
Lu = 6

! Initiate some stuff for automatic setting of
! the constraints to be those that the structure
! actually has initially.

rInt0_on_file = .false.
InSlapaf = SuperName == 'slapaf'
if (InSlapaf) call qpg_dArray('rInt0',rInt0_on_file,nrInt0)
if (.not. rInt0_on_File) nrInt0 = mInt
rInt0_in_memory = .false.

iRout = 30
iPrint = nPrint(iRout)
Start = lIter == 1
call ICopy(nBVct,[Flip],0,iFlip,1)

nTemp = len(Temp)
write(format,'(A,I3.3,A)') '(F',nTemp,'.0)'

Lu_UDC = 91
filnam = 'UDC'
call molcas_open(Lu_UDC,filnam)
!open(Lu_UDC,File=filnam,Form='FORMATTED',Status='OLD')
rewind(Lu_UDC)

call dcopy_(nBVct**2,[Zero],0,rMult,1)
if (iPrint >= 99) lWrite = .true.
if ((iPrint >= 99) .or. lWrite) then
  write(Lu,*)
  call CollapseOutput(1,'Constraints section')
  write(Lu,'(34X,A)') 'CONSTRAINTS'
  write(Lu,*)
  write(Lu,'(80A)') ('*',i=1,80)
  do iLines=1,nLines
    read(Lu_UDC,'(A)') Line
    write(Lu,'(A)') Line(:StrnLn(Line))
  end do
  write(Lu,'(80A)') ('*',i=1,80)
  write(Lu,*)
  write(Lu,*)
  write(Lu,*) '*************************************************************'
  write(Lu,*) '* Values of the primitive constraints                       *'
  write(Lu,*) '*************************************************************'
  rewind(Lu_UDC)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Step 1. Set up the b vectors from which we will define the constraints.
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
iBVct = 0
call mma_allocate(tpc,nBVct,Label='tpc')
do iLines=1,nLines
  read(Lu_UDC,'(A)') Line
  Temp = Line
  call UpCase(Temp)
  if (Temp(1:4) == 'VALU') Go To 100
  iBVct = iBVct+1

  ! Move the label of the internal coordinate

  neq = index(Line,'=')
  if (neq == 0) then
    call WarningMessage(2,'Error in DefInt2')
    write(Lu,'(A)') ' Syntax error in line:'
    write(Lu,'(A)') Line
    call Quit_OnUserError()
  else
    iFrst = 1
    call NxtWrd(Line,iFrst,iEnd)
    jEnd = iEnd
    if (Line(iEnd:iEnd) == '=') jEnd = jEnd-1
    if (jEnd-iFrst+1 > 8) then
      call WarningMessage(2,'Error in DefInt2')
      write(Lu,'(A,A)') Line(iFrst:jEnd),' has more than 8 character, syntax error!'
      call Quit_OnUserError()
    end if
    call ChkLbl(Line(iFrst:jEnd),Labels,iBVct-1)
    Labels(iBVct) = Line(iFrst:jEnd)
  end if

  ! Construct the corresponding transformation vector

  mCntr = 0
  iType = 0
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
      iType = 1
    else if (Temp(nGo:nGo2) == 'Y') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      type = 'Y     '
      iType = 2
    else if (Temp(nGo:nGo2) == 'Z') then
      nGo = nGo2+1
      call NxtWrd(Temp,nGo,nGo2)
      type = 'Z     '
      iType = 3
    else
      nGo = -1
      call WarningMessage(2,'Error in DefInt2')
      write(Lu,*) 'DefInt2: wrong cartesian type'
      write(Lu,'(A,A)') 'Temp=',Temp
      call Quit_OnUserError()
    end if
  else if (index(Temp,'BOND') /= 0) then
    nCntr = 2
    nGo = index(Temp,'BOND')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'STRTCH'
    iType = 4
  else if (index(Temp,'LANGLE(2)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(2)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'LBEND2'
    iType = 5
  else if (index(Temp,'LANGLE(1)') /= 0) then
    nCntr = 3
    nGo = index(Temp,'LANGLE(1)')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'LBEND1'
    iType = 6
  else if (index(Temp,'ANGL') /= 0) then
    nCntr = 3
    nGo = index(Temp,'ANGL')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'BEND  '
    iType = 7
  else if (index(Temp,'DIHE') /= 0) then
    nCntr = 4
    nGo = index(Temp,'DIHE')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'TRSN  '
    ! AOM!! Remove flip for Torsions!!!
    iFlip(iBVct) = NoFlip
    iType = 8
  else if (index(Temp,'OUTO') /= 0) then
    nCntr = 4
    nGo = index(Temp,'OUTO')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'OUTOFP'
    iFlip(iBVct) = NoFlip
    iType = 9
  else if (index(Temp,'EDIF') /= 0) then
    nCntr = nAtom
    nGo = index(Temp,'EDIF')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'EDIFF '
    iFlip(iBVct) = NoFlip
    iType = 10
  else if (index(Temp,'SPHE') /= 0) then
    nCntr = nAtom
    nGo = index(Temp,'SPHE')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'SPHERE'
    iType = 11
  else if (index(Temp,'NAC ') /= 0) then
    nCntr = nAtom
    nGo = index(Temp,'NAC ')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'NAC   '
    iType = 12
  else if (index(Temp,'TRAN') /= 0) then
    nCntr = nAtom
    nGo = index(Temp,'TRAN')
    nGo = nGo-1+index(Temp(nGo:nTemp),' ')
    type = 'TRANSV'
    iFlip(iBVct) = NoFlip
    iType = 13
  else if (index(Temp,'DISS') /= 0) then
    i1 = index(Line,'(')
    i2 = index(Line,'+')
    i3 = index(Line,')')
    if ((i1 >= i2) .or. (i2 >= i3)) then
      call WarningMessage(2,'Error in DefInt2')
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
    iType = 14
  else
    nGo = -1
    call WarningMessage(2,'Error in DefInt2')
    write(Lu,*) ' Line contains syntax error!'
    write(Lu,'(A)') Line
    call Quit_OnUserError()
  end if
  tpc(iBVct) = iType

  msAtom = nCntr+mCntr
  call mma_allocate(xyz,3*msAtom,Label='xyz')
  call mma_allocate(Grad,3*msAtom,Label='Grad')
  call mma_allocate(Ind,2*msAtom,Label='Ind')
  call mma_allocate(Mass,msAtom,Label='Mass')
  call mma_allocate(Hess,(3*msAtom)**2,Label='Hess')

  call Cllct2(Line(nGo:nTemp),BVct(1,iBVct),dBVct(1,1,iBVct),value(iBVct),nAtom,nCntr,mCntr,xyz,Grad,Ind,type,Mass,Labels(iBVct), &
              lWrite,rMult(iBVct,iBVct),Hess,lIter)

  if ((type == 'TRSN  ') .and. (abs(value(iBVct)) < Pi*Half)) iFlip(iBVct) = NoFlip

  call mma_deallocate(Hess)
  call mma_deallocate(Mass)
  call mma_deallocate(Ind)
  call mma_deallocate(Grad)
  call mma_deallocate(xyz)

end do
call WarningMessage(2,'Error in DefInt2')
write(Lu,*) 'DefInt2: No internal coordinates are defined!'
call Quit_OnUserError()

100 continue

! Process primitive value to correct for flips

if (.not. Start) then
  do iBVct=1,nBVct

    ! Test if we have a flip in the sign of the value

    if ((value(iBVct)*Value0(iBVct) < Zero) .and. (iFlip(iBVct) == Flip)) then
      !write(Lu,*) 'Flip Sign for ',Labels(iBVct)
      value(iBVct) = -value(iBVct)
    end if
  end do
end if
call dcopy_(nBVct,value,1,Value0,1)

if (iPrint >= 59) call RecPrt(' The B-vectors',' ',BVct,3*nAtom,nBVct)
if (iPrint >= 19) then
  call RecPrt(' Value of primitive internal coordinates / au or rad',' ',value,nBVct,1)
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Step 2. Define constraints as linear combinations of
!         the previously defined primitive internal coordinates.
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
iBMtrx = 0
jLines = iLines
201 continue
jLines = jLines+1
if (jLines > nLines) Go To 200

read(Lu_UDC,'(A)') Line
Temp = Line
call UpCase(Temp)

iBMtrx = iBMtrx+1
rInt(iBMtrx) = Zero
write(Lbl(iBMtrx),'(A,I3.3)') 'Cns',iBMtrx

! Find the label of the first primitive

iFrst = 1
call NxtWrd(Line,iFrst,iEnd)
jEnd = iEnd
if (Line(iEnd:iEnd) == '=') jEnd = iEnd-1

nPlus = index(Line,' + ')
nMinus = index(Line,' - ')
!                                                                      *
!***********************************************************************
!                                                                      *
if ((nPlus == 0) .and. (nMinus == 0) .and. (index(Line,'&') == 0)) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! a single vector (this will only extend over one line)

  if (index(Line,'&') /= 0) then
    call WarningMessage(2,'Error in DefInt2')
    write(Lu,*) 'Single vector lines should not extend'
    write(Lu,*) 'over more than one line!'
    write(Lu,'(A)') Line
    call Quit_OnUserError()
  end if
  iBVct = 0
  do jBVct=1,nBVct
    if (Line(iFrst:jEnd) == Labels(jBVct)) iBVct = jBVct
  end do
  if (iBVct == 0) then
    call WarningMessage(2,'Error in DefInt2')
    write(Lu,*) ' A single vector'
    write(Lu,*) ' Undefined internal coordinate'
    write(Lu,'(A,A)') Line
    write(Lu,'(A,A)') Line(iFrst:jEnd)
    call Quit_OnUserError()
  end if

  call dcopy_(3*nAtom,BVct(1,iBVct),1,BMtrx(1,iBMtrx),1)
  call dcopy_((3*nAtom)**2,dBVct(1,1,iBVct),1,dBMtrx(1,1,iBMtrx),1)
  rInt(iBMtrx) = value(iBVct)

  iFrst = iEnd+1
  call NxtWrd(Line,iFrst,iEnd)
  if (Line(iEnd:iEnd) == '=') then
    iFrst = iEnd+1
    call NxtWrd(Line,iFrst,iEnd)
  end if
  Temp = Line(iFrst:iEnd)

  if (index(Temp,'FIX') /= 0) then

    ! Pick up values from the runfile. Written there on the first iteration.

    if (.not. rInt0_in_memory) then
      rInt0_in_memory = .true.
      call mma_allocate(r0,nrInt0,Label='r0')
      if (rInt0_on_file) then
        call Get_dArray('rInt0',r0,nrInt0)
      else
        r0(:) = Zero
      end if
    end if

    if (rInt0_on_file) then
      rInt0(iBMtrx) = r0(iBMtrx)
    else
      r0(iBmtrx) = rInt(iBMtrx)
      rInt0(iBMtrx) = rInt(iBMtrx)
    end if

  else

    ! Read value from input file.

    read(Temp,format) rInt0(iBMtrx)
    Temp = Line
    call UpCase(Temp)
    if (index(Temp,'ANGSTROM') /= 0) rInt0(iBMtrx) = rInt0(iBMtrx)/angstr
    if (index(Temp,'DEGREE') /= 0) rInt0(iBMtrx) = rInt0(iBMtrx)*Pi/1.800D+02
  end if

  ! AOM: trying to correct torsion...
  if (tpc(iBVct) == 8) then
    n0 = int((rInt(iBMtrx)-rInt0(iBMtrx))/(Two*Pi))
    if (abs(rInt(IBMtrx)-rInt0(iBMtrx)-dble(n0)*Two*Pi) > abs(rInt(IBMtrx)-rInt0(iBMtrx)-dble(n0+1)*Two*Pi)) then
      n0 = n0+1
    else if (abs(rInt(IBMtrx)-rInt0(iBMtrx)-dble(n0)*Two*Pi) > abs(rInt(IBMtrx)-rInt0(iBMtrx)-dble(n0-1)*Two*Pi)) then
      n0 = n0-1
    end if
    rInt(iBMtrx) = rInt(iBMtrx)-dble(n0)*Two*Pi
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! A linear combination of vectors

  call dcopy_(3*nAtom,[Zero],0,BMtrx(1,iBMtrx),1)
  call dcopy_((3*nAtom)**2,[Zero],0,dBMtrx(1,1,iBMtrx),1)
  iFrst = 1
  Sgn = One

  ! Process the rest of the line and possible extension lines

  22 continue

  !*                                                                   *
  !*********************************************************************
  !*                                                                   *
  !> Get the factor
  call NxtWrd(Line,iFrst,iEnd)
  Temp = Line(iFrst:iEnd)
  read(Temp,format) Fact
  Fact = Fact*Sgn
  iFrst = iEnd+1
  !> Get the label
  call NxtWrd(Line,iFrst,iEnd)
  if (Line(iEnd:iEnd) == '=') iEnd = iEnd-1
  iBVct = 0
  do jBVct=1,nBVct
    if (Line(iFrst:iEnd) == Labels(jBVct)) iBVct = jBVct
  end do
  if (iBVct == 0) then
    call WarningMessage(2,'Error in DefInt2')
    write(Lu,*) ' Linear combinations of vectors'
    write(Lu,*) ' Undefined internal coordinate'
    write(Lu,'(A,A)') Line
    write(Lu,'(A,A)') Line(iFrst:iEnd)
    call Quit_OnUserError()
  end if

  call DaXpY_(3*nAtom,Fact,BVct(1,iBvct),1,BMtrx(1,iBMtrx),1)
  call DaXpY_((3*nAtom)**2,Fact,dBVct(1,1,iBvct),1,dBMtrx(1,1,iBMtrx),1)
  rInt(iBMtrx) = rInt(iBMtrx)+Fact*value(iBVct)
  !*                                                                   *
  !*********************************************************************
  !*                                                                   *

  iFrst = iEnd+1
  25 continue
  Temp = Line(iFrst:nTemp)
  nEq = index(Temp,'=')
  nPlus = index(Temp,'+')
  nMinus = index(Temp,'-')

  if ((nEq /= 0) .and. ((nEq < nMinus) .eqv. (nMinus > 0)) .and. ((nEq < nPlus) .eqv. (nPlus > 0))) then
    if (index(Line,'&') /= 0) then
      call WarningMessage(2,'Error in DefInt2')
      write(Lu,*) 'This line should not be extended'
      write(Lu,'(A)') Line
      call Quit_OnUserError()
    end if
    iFrst = iFrst+nEq+1
    call NxtWrd(Line,iFrst,iEnd)
    Temp = Line(iFrst:iEnd)

    if (index(Temp,'FIX') /= 0) then

      ! Pick up values from the runfile. Written there on the first iteration.

      if (.not. rInt0_in_memory) then
        rInt0_in_memory = .true.
        call mma_allocate(r0,nrInt0,Label='r0')
        if (rInt0_on_file) then
          call Get_dArray('rInt0',r0,nrInt0)
        else
          r0(:) = Zero
        end if
      end if

      if (rInt0_on_file) then
        rInt0(iBMtrx) = r0(iBMtrx)
      else
        r0(iBmtrx) = rInt(iBMtrx)
      end if

    else

      ! Read value from input file.

      read(Temp,format) rInt0(iBMtrx)
      Temp = Line
      call UpCase(Temp)
      if (index(Temp,'ANGSTROM') /= 0) rInt0(iBMtrx) = rInt0(iBMtrx)/angstr
      if (index(Temp,'DEGREE') /= 0) rInt0(iBMtrx) = rInt0(iBMtrx)*Pi/1.800D+02
    end if
    Go To 24
  end if

  if ((nPlus /= 0) .and. ((nPlus < nMinus) .eqv. (nMinus > 0))) then
    Sgn = One
    iFrst = iFrst+nPlus
    Go To 22
  end if

  if ((nMinus /= 0) .and. ((nMinus < nPlus) .eqv. (nPlus > 0))) then
    Sgn = -One
    iFrst = iFrst+nMinus
    Go To 22
  end if

  ! Here if all statements processed of this line

  if (index(Line,'&') /= 0) then
    jLines = jLines+1
    if (jLines > nLines) then
      call WarningMessage(2,'Error in DefInt2')
      write(Lu,*) 'DefInt2: jLines > nLines'
      call Quit_OnUserError()
    end if
    read(Lu_UDC,'(A)') Line
    iFrst = 1
    call NxtWrd(Line,iFrst,iEnd)
    if (Line(iFrst:iEnd) == '+') then
      iFrst = iEnd+1
      Sgn = One
      Go To 22
    else if (Line(iFrst:iEnd) == '-') then
      iFrst = iEnd+1
      Sgn = -One
      Go To 22
    else if (Line(iFrst:iEnd) == '=') then
      iFrst = 1
      Go To 25
    else
      call WarningMessage(2,'Error in DefInt2')
      write(Lu,*) ' Syntax Error: first character in  extension line is not + or -'
      write(Lu,'(A)') Line
      write(Lu,'(3A)') '-->',Line(iFrst:iEnd),'<--'
      call Quit_OnUserError()
    end if
  end if

  ! At the end of this line

  24 continue
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
Go To 201

200 continue
if (iPrint >= 99) then
  call RecPrt(' The B-matrix',' ',BMtrx,3*nAtom,mInt)
  do iInt=1,mInt
    call RecPrt(' The dB-matrix',' ',dBMtrx(1,1,iInt),3*nAtom,3*nAtom)
  end do
end if
close(Lu_UDC)
call mma_deallocate(tpc)
if (rint0_in_memory) then
  if (InSlapaf) call Put_dArray('rInt0',r0,mInt)
  call mma_deallocate(r0)
end if

! Compute the maximum error in the constraints

MaxErr = Zero
do iInt=1,mInt
  MaxErr = max(MaxErr,abs(rInt(iInt)-rInt0(iInt)))
end do
call Put_dScalar('Max error',MaxErr)

if ((iPrint >= 99) .or. lWrite) then
  write(Lu,*)
  write(Lu,*)
  write(Lu,*) '*******************************************'
  write(Lu,*) '* Values of the constraints   / au or rad *'
  write(Lu,*) '*******************************************'
  write(Lu,*) '  Label        C         C0'
  write(Lu,'(1X,A,2X,F10.6,F10.6)') (Lbl(iInt),rInt(iInt),rInt0(iInt),iInt=1,mInt)
  write(Lu,*)
  call CollapseOutput(0,'Constraints section')
  write(Lu,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Labels)

return

end subroutine DefInt2

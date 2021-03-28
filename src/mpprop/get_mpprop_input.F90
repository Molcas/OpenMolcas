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

subroutine Get_Mpprop_input(nAtoms,iPol,LNearestAtom,LAllCenters,AveOrb,LLumOrb,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,ThrsMul,iPrint)

implicit real*8(a-h,o-z)

#include "MolProp.fh"
#include "warnings.fh"

character*3 EndKey
! Jose Character*4 TestLabe(0:nAtoms), KWord
character*4 KWord
character*6 TestLabe(0:nAtoms)
character*180 Key, BLine
character*180 Get_Ln
logical Debug, LNearestAtom
logical LAllCenters, AveOrb, Diffuse(3)
logical LLumorb
dimension dLimmo(2)

external Get_Ln

data Debug/.false./

iStdOut = 6

do i=1,180
  BLine(i:i) = ' '
end do
Title = BLine

LuRd = 21
call SpoolInp(LuRd)

rewind(LuRd)
call RdNLst(LuRd,'MPPROP')

! KeyWord directed input

998 Key = Get_Ln(LuRd)
if (Debug) write(iStdOut,*) ' Processing:',Key
KWord = trim(Key)
call UpCase(KWord)
if (KWord(1:1) == '*') Go To 998
if (KWord == BLine) Go To 998
if (KWord(1:4) == 'BOND') Go To 981
!if (KWord(1:4) == 'METH') Go To 982
if (KWord(1:4) == 'TITL') Go To 983
if (KWord(1:4) == 'TYPE') Go To 984
if (KWord(1:4) == 'POLA') Go To 985
if (KWord(1:4) == 'NONE') Go To 986
if (KWord(1:4) == 'ALLC') Go To 987
if (KWord(1:4) == 'LUMO') Go To 988
if (KWord(1:4) == 'DIFF') Go To 989
if (KWord(1:4) == 'PRIN') Go To 991
! Keyword added to handle average orbitals
if (KWord(1:4) == 'AVER') Go To 7912

if (KWord(1:4) == 'END ') Go To 997
iChrct = len(KWord)
Last = iCLast(KWord,iChrct)
write(6,*)
write(6,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
write(6,*) ' Error in keyword.'
call ErrTra
call Quit(_RC_INPUT_ERROR_)
write(6,*) ' Premature end of input file.'
call Quit(_RC_INPUT_ERROR_)
write(6,*) ' Error while reading input file.'
990 call ErrTra
call Quit(_RC_INPUT_ERROR_)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the input for bonds
981 LAllCenters = .true.
do i=1,MxAtomMP
  NUB(i) = 0
  do j=1,MxAtomMP
    NBI(i,j) = 0
    BondMat(i,j) = .false.
  end do
end do
do i=1,nAtoms
  nBonds = 0
  Key = Get_Ln(LuRd)
  m = 1
  do j=1,nAtoms
    TestLabe(j) = ' '
  end do
  do j=1,180
    EndKey = Key(j:j+2)
    call UpCase(EndKey)
    if ((Key(j:j) == ' ') .or. (Key(j:j) == ',') .or. (Key(j:j) == ';')) then
      if ((j == 1) .and. ((Key(j:j) == ';') .or. (Key(j:j) == ','))) then
        write(iStdOut,*) 'Error in input, breaker in first position'
        goto 990
      elseif (m < j) then
        TestLabe(nBonds) = Key(m:j-1)
        nBonds = nBonds+1
        m = j+1
      else
        m = j+1
      end if
    elseif (EndKey == 'END') then
      goto 9811
    end if
  end do
  do j=1,nAtoms
    if (TestLabe(0) == Labe(j)) then
      do k=1,nAtoms
        do l=1,nAtoms
          if (TestLabe(k) == Labe(l)) then
            BondMat(j,l) = .true.
            BondMat(l,j) = .true.
          else

          end if
        end do
      end do
    end if
  end do
end do
9811 continue

do i=1,nAtoms
  do j=1,nAtoms
    if (BondMat(i,j)) then
      NuB(i) = NuB(i)+1
      NBI(i,NuB(i)) = j
    end if
  end do
end do

write(iStdOut,*)
write(iStdOut,'(10X,A)') '************************'
write(iStdOut,'(10X,A)') '**** Bonding matrix ****'
write(iStdOut,'(10X,A)') '************************'
write(iStdOut,*)
write(iStdOut,'(A8,A,A)') 'Atom','  No bonds','   Bonding with'
do i=1,nAtoms
  write(iStdOut,'(A8,I6,A11,1000A8)') LABE(I),NUB(I),(LABE(NBI(I,J)),J=1,NUB(I))
end do
write(iStdOut,*)
write(iStdOut,*)
goto 998

! Set Method level
!982 continue
!Key = Get_Ln(LuRd)
!if (Debug) write (6,*) ' Processing:',Key
!Method = Key
!call UpCase(Method)
!if ((Method == 'SCF') .or. (Method == 'RASSCF')) then
!
!else
!  write(6,*) 'The Method Label is not correct'
!  write(6,*) Method
!  call Quit(20)
!end if
!Goto 998
! Set the Title
983 Key = Get_Ln(LuRd)
if (Debug) write(iStdOut,*) ' Processing:',Key
call UpCase(Key)
Title = Key
goto 998
! Set the specific atom types
984 Key = Get_Ln(LuRd)
call UpCase(Key)
if (Key(1:3) == 'END') then
  Go To 998
else
  read(Key,*) j,iAtomPar(j)
end if
goto 984
! Set that local polarizabilities should be calculated
985 Key = Get_Ln(LuRd)
read(Key,*) iPol
goto 998
! Set the nearest atom calculation to .false.
986 LNearestAtom = .false.
goto 998
! Set the average orbital option to .true.
7912 AveOrb = .true.
iPol = 0
goto 998
! Set the logical variable that defines if all centers or not
987 LAllCenters = .true.
LNearestAtom = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
! Set the bonds to all sites
! Get information from input
!do i=1,mxAtomMP
!  NuB(i) = 0
!  do j=1,mxAtomMP
!    NBI(i,j) = 0
!  end do
!end do
do i=1,nAtoms
  NuB(i) = nAtoms-1
  do j=1,nAtoms-1
    if (j >= i) then
      NBI(i,j) = j+1
      BondMat(i,j+1) = .true.
    else
      NBI(i,j) = j
      BondMat(i,j) = .true.
    end if
  end do
end do
goto 998
! Get the vectors from the INPORB file
988 LLumOrb = .true.
goto 998

! Obtain diffuse stuff to the MpProp decomposition.
989 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
if (Key(1:4) == 'NUME') then
  Diffuse(1) = .true.
  Diffuse(2) = .true.
9891 continue
  Key = Get_Ln(LuRd)
  call UpCase(Key)
  if (Key(1:4) == 'LIMI') then
    Key = Get_Ln(LuRd)
    call Get_F(1,dLimmo,2)
  elseif (Key(1:4) == 'THRE') then
    Key = Get_Ln(LuRd)
    call Get_F1(1,Thrs1)
    call Get_F1(2,Thrs2)
    call Get_I1(3,nThrs)
    call Get_F1(4,ThrsMul)
  elseif (Key(1:4) == 'END ') then
    goto 987 !Not an error, Diffuse implies AllCenters
  else
    write(6,*) 'Undefined option for ''DIFFuse'':',Key
    call FindErrorLine
    call Quit_OnUserError()
  end if
  goto 9891
elseif (Key(1:4) == 'REXT') then
  Diffuse(1) = .true.
  Diffuse(3) = .true.
else
  write(6,*) 'Undefined option for ''DIFFuse'':',Key
  call FindErrorLine
  call Quit_OnUserError()
end if
goto 987  !Not an error, Diffuse implies AllCenters

!PrintLevel.
991 continue
Key = Get_Ln(LuRd)
call UpCase(Key)
call Get_I1(1,iPrint)
goto 998

997 continue
if (Title == Bline) then
  write(iStdOut,*)
  write(iStdOut,*) ' !!WARNING!! The molecule do not have a name'
  write(iStdOut,*)
end if

return

end subroutine Get_Mpprop_input

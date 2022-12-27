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

use MPProp_globals, only: BondMat, iAtomPar, Labe, Title
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms
integer(kind=iwp), intent(inout) :: iPol, nThrs, iPrint
logical(kind=iwp), intent(inout) :: LNearestAtom, LAllCenters, AveOrb, LLumOrb, Diffuse(3)
real(kind=wp), intent(inout) :: dLimmo(2), Thrs1, Thrs2, ThrsMul
#include "warnings.h"
integer(kind=iwp) :: i, iChrct, iStdOut, j, k, l, Last, LuRd, m, nBonds
character(len=3) :: EndKey
character(len=4) :: KWord
character(len=180) :: Key
logical(kind=iwp) :: Skip
integer(kind=iwp), allocatable :: NuB(:), NBI(:,:)
character(len=6), allocatable :: TestLabe(:)
logical(kind=iwp), parameter :: Debug = .false.
integer(kind=iwp), external :: iCLast
character(len=180), external :: Get_Ln

iStdOut = u6

Title = ' '

LuRd = 21
call SpoolInp(LuRd)

rewind(LuRd)
call RdNLst(LuRd,'MPPROP')

call mma_allocate(NuB,nAtoms,label='NuB')
call mma_allocate(NBI,nAtoms,nAtoms,label='NBI')
NuB(:) = 0
NBI(:,:) = 0

! KeyWord directed input

Skip = .false.
do
  if (.not. Skip) then
    Key = Get_Ln(LuRd)
    if (Debug) write(iStdOut,*) ' Processing:',Key
    KWord = trim(Key)
    call UpCase(KWord)
  end if
  Skip = .false.
  if (KWord(1:1) == '*') cycle
  if (KWord == ' ') cycle
  select case (KWord(1:4))
    case ('BOND')
      ! Get the input for bonds
      LAllCenters = .true.
      atoms: do i=1,nAtoms
        nBonds = 0
        Key = Get_Ln(LuRd)
        m = 1
        call mma_allocate(TestLabe,nAtoms,label='TestLabe')
        TestLabe(:) = ' '
        do j=1,180
          EndKey = Key(j:j+2)
          call UpCase(EndKey)
          if ((Key(j:j) == ' ') .or. (Key(j:j) == ',') .or. (Key(j:j) == ';')) then
            if ((j == 1) .and. ((Key(j:j) == ';') .or. (Key(j:j) == ','))) then
              write(iStdOut,*) 'Error in input, breaker in first position'
              call Quit(_RC_INPUT_ERROR_)
            else if (m < j) then
              TestLabe(nBonds) = Key(m:j-1)
              nBonds = nBonds+1
            end if
            m = j+1
          else if (EndKey == 'END') then
            exit atoms
          end if
        end do
        do j=1,nAtoms
          if (TestLabe(0) == Labe(j)) then
            do k=1,nAtoms
              do l=1,nAtoms
                if (TestLabe(k) == Labe(l)) then
                  BondMat(j,l) = .true.
                  BondMat(l,j) = .true.
                end if
              end do
            end do
          end if
        end do
        call mma_deallocate(TestLabe)
      end do atoms

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
        write(iStdOut,'(A8,I6,A11,1000A8)') Labe(i),NuB(i),(Labe(NBI(i,j)),j=1,NuB(i))
      end do
      write(iStdOut,*)
      write(iStdOut,*)

    !case ('METH')
    !  ! Set Method level
    !  Key = Get_Ln(LuRd)
    !  if (Debug) write(iStdOut,*) ' Processing:',Key
    !  Method = Key
    !  call UpCase(Method)
    !  if ((Method == 'SCF') .or. (Method == 'RASSCF')) then
    !
    !  else
    !    write(iStdOut,*) 'The Method Label is not correct'
    !    write(iStdOut,*) Method
    !    call Quit(20)
    !  end if

    case ('TITL')
      ! Set the Title
      Key = Get_Ln(LuRd)
      if (Debug) write(iStdOut,*) ' Processing:',Key
      call UpCase(Key)
      Title = Key

    case ('TYPE')
      ! Set the specific atom types
      do
        Key = Get_Ln(LuRd)
        call UpCase(Key)
        if (Key(1:3) == 'END') exit
        read(Key,*) j,iAtomPar(j)
      end do

    case ('POLA')
      ! Set that local polarizabilities should be calculated
      Key = Get_Ln(LuRd)
      read(Key,*) iPol

    case ('NONE')
      ! Set the nearest atom calculation to .false.
      LNearestAtom = .false.

    case ('ALLC')
      ! Set the logical variable that defines if all centers or not
      LAllCenters = .true.
      LNearestAtom = .false.
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Set the bonds to all sites
      ! Get information from input
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

    case ('LUMO')
      ! Get the vectors from the INPORB file
      LLumOrb = .true.

    case ('DIFF')
      ! Obtain diffuse stuff to the MpProp decomposition.
      Key = Get_Ln(LuRd)
      call UpCase(Key)
      if (Key(1:4) == 'NUME') then
        Diffuse(1) = .true.
        Diffuse(2) = .true.
        do
          Key = Get_Ln(LuRd)
          call UpCase(Key)
          if (Key(1:4) == 'LIMI') then
            Key = Get_Ln(LuRd)
            call Get_F(1,dLimmo,2)
          else if (Key(1:4) == 'THRE') then
            Key = Get_Ln(LuRd)
            call Get_F1(1,Thrs1)
            call Get_F1(2,Thrs2)
            call Get_I1(3,nThrs)
            call Get_F1(4,ThrsMul)
          else if (Key(1:4) == 'END ') then
            exit
          else
            write(iStdOut,*) 'Undefined option for ''DIFFuse'':',Key
            call FindErrorLine()
            call Quit_OnUserError()
          end if
        end do
      else if (Key(1:4) == 'REXT') then
        Diffuse(1) = .true.
        Diffuse(3) = .true.
      else
        write(iStdOut,*) 'Undefined option for ''DIFFuse'':',Key
        call FindErrorLine()
        call Quit_OnUserError()
      end if
      !Not an error, Diffuse implies AllCenters
      KWord = 'ALLC'
      Skip = .true.

    case ('PRIN')
      !PrintLevel.
      Key = Get_Ln(LuRd)
      call UpCase(Key)
      call Get_I1(1,iPrint)

    ! Keyword added to handle average orbitals
    case ('AVER')
      ! Set the average orbital option to .true.
      AveOrb = .true.
      iPol = 0

    case ('END ')
      exit

    case default
      iChrct = len(KWord)
      Last = iCLast(KWord,iChrct)
      write(iStdOut,*)
      write(iStdOut,'(1X,A,A)') KWord(1:Last),' is not a keyword!'
      write(iStdOut,*) ' Error in keyword.'
      call Quit(_RC_INPUT_ERROR_)
      !write(iStdOut,*) ' Premature end of input file.'
      !call Quit(_RC_INPUT_ERROR_)
      !write(iStdOut,*) ' Error while reading input file.'
      !call Quit(_RC_INPUT_ERROR_)
  end select
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (Title == ' ') then
  write(iStdOut,*)
  write(iStdOut,*) ' !!WARNING!! The molecule does not have a name'
  write(iStdOut,*)
end if

call mma_deallocate(NuB)
call mma_deallocate(NBI)

return

end subroutine Get_Mpprop_input

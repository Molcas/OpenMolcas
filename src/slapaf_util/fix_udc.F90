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

subroutine Fix_UDC(iRow_c,nLambda,nsAtom,nStab,Remove)
!***********************************************************************
!                                                                      *
!     Object: Replace the fragment specifications with actual          *
!             constraint specifications.                               *
!             Remove constraints to be ignored.                        *
!                                                                      *
!***********************************************************************

use Slapaf_Info, only: AtomLbl
use UnixInfo, only: SuperName
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: iRow_c, nLambda
integer(kind=iwp), intent(in) :: nsAtom, nStab(nsAtom)
logical(kind=iwp), intent(in) :: Remove
integer(kind=iwp) :: i, iAtom, iEnd, iEq, iFrag, iFrst, iLine, iSoft, iZMat, jAtom, Lu_TMP, Lu_UDC, nDeg, nFrag, nLabel1, nLabel2, &
                     nLabel3, nLabel4, nLines, nSoft, Num, nZMat, ZMatOffset
logical(kind=iwp) :: Ignore, Values
character(len=180) :: Label1, Label2, Label3, Label4, Line, Line2
character(len=16) :: FilNam
integer(kind=iwp), allocatable :: FragZMat(:,:)
character(len=180), allocatable :: FragLabels(:), Lines(:), SoftLabels(:)
character(len=180), external :: Get_Ln
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
if (iRow_c <= 1) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Open files for user defined internal constraints.

Lu_UDC = 91
FilNam = 'UDC'
Lu_UDC = IsFreeUnit(Lu_UDC)
call Molcas_Open(Lu_UDC,FilNam)
rewind(Lu_UDC)

Lu_TMP = Lu_UDC+1
FilNam = 'UDCTMP'
Lu_TMP = IsFreeUnit(Lu_TMP)
call Molcas_Open(Lu_TMP,FilNam)
rewind(Lu_TMP)
!                                                                      *
!***********************************************************************
!                                                                      *
! Process the content of LU_UDC and put it onto LU_TMP

call mma_allocate(FragZMat,iRow_c,2,Label='FragZMat')
call mma_allocate(FragLabels,iRow_c,Label='FragLabels')
call mma_allocate(Lines,iRow_c,Label='Lines')
call mma_allocate(SoftLabels,iRow_c,Label='SoftLabels')

iRow_c = 0 ! Number of lines in the section
nZMat = 0  ! Number of additional lines with constraints
nFrag = 0  ! Number of fragment constraints
nSoft = 0  ! Number of default soft constraints
Values = .false.

iZMat = 0
ZMatOffset = 0
do
  Line = Get_Ln(LU_UDC)   ! Read the line
  call UpCase(Line)

  if (Line(1:4) == 'VALU') Values = .true.
  if (Line(1:4) == 'END') exit

  ! Read all possible continuation lines

  nLines = 1
  Lines(nLines) = Line
  do while (index(Lines(nLines),'&') > 0)
    nLines = nLines+1
    Lines(nLines) = Get_Ln(LU_UDC)
    call UpCase(Lines(nLines))
  end do
  iLine = 1
  Line = Lines(iLine)
  iEq = index(Line,'=')

  ! If there are already some ZMAT constraints in the file,
  ! start numbering from the next one

  if (Line(1:4) == 'ZMAT') then
    read(Line(5:7),*) Num
    ZMatOffset = max(ZMatOffset,Num)
  end if

  if (.not. Values) then

    ! Find out if the primitive defaults to soft

    iFrst = iEq+1
    call NxtWrd(Line,iFrst,iEnd)
    if ((Line(iFrst:iFrst+3) == 'SPHE') .or. (Line(iFrst:iFrst+3) == 'TRAN') .or. (Line(iFrst:iFrst+3) == 'EDIF') .or. &
        (Line(iFrst:iFrst+3) == 'NAC ')) then
      nSoft = nSoft+1
      iFrst = 1
      call NxtWrd(Line,iFrst,iEnd)
      SoftLabels(nSoft) = Line(iFrst:iEnd)
    end if

    ! Modify the lines with Fragments to Z-matrix format

    iFrst = iEq+1
    call NxtWrd(Line,iFrst,iEnd)
    if (Line(iFrst:iFrst+7) == 'FRAGMENT') then

      ! Here we do the expansion

      nFrag = nFrag+1
      FragZMat(nFrag,1) = 0
      FragZMat(nFrag,2) = 0
      iFrst = 1
      call NxtWrd(Line,iFrst,iEnd)
      FragLabels(nFrag) = Line(iFrst:iEnd)
      iFrst = index(Line,'FRAGMENT')
      call NxtWrd(Line,iFrst,iEnd)

      ! Pick up the first atom label

      iFrst = iEnd+1
      call NxtWrd(Line,iFrst,iEnd)
      if (Line(iFrst:iEnd) == '&') then
        iLine = iLine+1
        Line = Lines(iLine)     ! Read the continuation line
        iFrst = 1
        call NxtWrd(Line,iFrst,iEnd)
      end if
      Label1 = Line(iFrst:iEnd)
      nLabel1 = iEnd-iFrst+1

      ! Pick up the second atom label

      iFrst = iEnd+1
      call NxtWrd(Line,iFrst,iEnd)
      if (Line(iFrst:iEnd) == '&') then
        iLine = iLine+1
        Line = Lines(iLine)     ! Read the continuation line
        iFrst = 1
        call NxtWrd(Line,iFrst,iEnd)
      end if
      Label2 = Line(iFrst:iEnd)
      nLabel2 = iEnd-iFrst+1

      ! Case: Bond

      nZMat = nZMat+1
      write(Line2,'(A,I3.3,A,2(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Bond',Label2(1:nLabel2),Label1(1:nLabel1)

      write(Lu_TMP,'(A)') Line2
      iRow_c = iRow_c+1
      iZMat = iZMat+1
      FragZMat(nFrag,1) = iZMat
      FragZMat(nFrag,2) = iZMat

      ! Pick up the third atom label

      iFrst = iEnd+1
      call NxtWrd(Line,iFrst,iEnd)
      if (iEnd /= -1) then
        if (Line(iFrst:iEnd) == '&') then
          iLine = iLine+1
          Line = Lines(iLine)     ! Read the continuation line
          iFrst = 1
          call NxtWrd(Line,iFrst,iEnd)
        end if
        Label3 = Line(iFrst:iEnd)
        nLabel3 = iEnd-iFrst+1
        jAtom = 0
        do iAtom=1,nsAtom
          if (Label3 == AtomLbl(iAtom)) jAtom = iAtom
        end do
        if (nStab(jAtom) == 1) then
          nDeg = 3
        else if (nStab(jAtom) == 2) then
          nDeg = 2
        else if (nStab(jAtom) == 1) then
          nDeg = 1
        else
          nDeg = 0
        end if

        ! Case: Bond-Angle
        if (nDeg > 1) then
          nZMat = nZMat+1
          write(Line2,'(A,I3.3,A,2(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Bond',Label3(1:nLabel3),Label2(1:nLabel2)
          write(Lu_TMP,'(A)') Line2
          iRow_c = iRow_c+1
        end if

        if (nDeg > 2) then
          nZMat = nZMat+1
          write(Line2,'(A,I3.3,A,3(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Angle',Label3(1:nLabel3),Label2(1:nLabel2),Label1(1:nLabel1)
          write(Lu_TMP,'(A)') Line2
          iRow_c = iRow_c+1
          iZMat = iZMat+2
          FragZMat(nFrag,2) = iZMat
        end if

        do

          ! Pick up next atom (fourth)

          iFrst = iEnd+1
          call NxtWrd(Line,iFrst,iEnd)

          ! Branch out if no more labels. Look out for line continuations

          if (iEnd == -1) exit
          if (Line(iFrst:iEnd) == '&') then
            iLine = iLine+1
            Line = Lines(iLine)     ! Read the continuation line
            iFrst = 1
            call NxtWrd(Line,iFrst,iEnd)
          end if

          Label4 = Line(iFrst:iEnd)
          nLabel4 = iEnd-iFrst+1
          jAtom = 0
          do iAtom=1,nsAtom
            if (Label4 == AtomLbl(iAtom)) jAtom = iAtom
          end do
          if (nStab(jAtom) == 1) then
            nDeg = 3
          else if (nStab(jAtom) == 2) then
            nDeg = 2
          else if (nStab(jAtom) == 1) then
            nDeg = 1
          else
            nDeg = 0
          end if

          ! Case: Bond-Angle-Dihedral
          if (nDeg /= 0) then
            nZMat = nZMat+1
            write(Line2,'(A,I3.3,A,2(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Bond',Label4(1:nLabel4),Label3(1:nLabel3)
            write(Lu_TMP,'(A)') Line2
            iRow_c = iRow_c+1

            if (nDeg /= 1) then
              nZMat = nZMat+1
              write(Line2,'(A,I3.3,A,3(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Angle',Label4(1:nLabel4),Label3(1:nLabel3), &
                                                Label2(1:nLabel2)
              write(Lu_TMP,'(A)') Line2
              iRow_c = iRow_c+1

              if (nDeg /= 2) then
                nZMat = nZMat+1
                write(Line2,'(A,I3.3,A,4(1X,A))') 'ZMAT',nZMat+ZMatOffset,' = Dihedral',Label4(1:nLabel4),Label3(1:nLabel3), &
                                                  Label2(1:nLabel2),Label1(1:nLabel1)
                write(Lu_TMP,'(A)') Line2
                iRow_c = iRow_c+1
                iZMat = iZMat+3
                FragZMat(nFrag,2) = iZMat
              end if
            end if
          end if

          ! Push the labels

          Label1 = Label2
          nLabel1 = nLabel2
          Label2 = Label3
          nLabel2 = nLabel3
          Label3 = Label4
          nLabel3 = nLabel4

          ! Get the next label

        end do

      end if

    else

      ! No processing just write to the file.

      write(Lu_TMP,'(A)') Line
      iRow_c = iRow_c+1
    end if

  else ! Values

    ! Remove constraints to be ignored

    Ignore = .false.
    iFrst = 1
    call NxtWrd(Lines(1),iFrst,iEnd)
    iFrag = 0
    do i=1,nFrag
      if (Line(iFrst:iEnd) == trim(FragLabels(i))) iFrag = i
    end do
    ! Some constraints default to soft unless marked as hard
    iSoft = 0
    if (index(Lines(nLines),'HARD') <= iEq) then
      do i=1,nSoft
        if (Line(iFrst:iEnd) == trim(SoftLabels(i))) iSoft = i
      end do
    end if
    iEq = index(Lines(nLines),'=')
    ! Ignore soft in num. grad., phantom in slapaf
    if ((SuperName == 'numerical_gradient') .and. ((index(Lines(nLines),'SOFT') > iEq) .or. (iSoft > 0))) Ignore = Remove
    if ((SuperName == 'slapaf') .and. (index(Lines(nLines),'PHANTOM') > iEq)) Ignore = Remove
    ! Ignore if this is a fragment constraint
    if (iFrag > 0) then
      ! if it is already ignored, ignore all replacement constraints
      if (Ignore) FragZMat(iFrag,2) = 0
      Ignore = .true.
    end if
    if (.not. Ignore) then
      do i=1,nLines
        write(Lu_TMP,'(A)') Lines(i)
      end do
      iRow_c = iRow_c+nLines
    else
      nLambda = nLambda-1
    end if

  end if ! Values

end do ! Get the next line

! Put in the additional constraints here.

do iFrag=1,nFrag
  ! Note that FragZMat(iFrag,2)=0 if the fragment is ignored
  do iZMat=FragZMat(iFrag,1),FragZMat(iFrag,2)
    write(Lu_TMP,'(A,I3.3,A)') 'ZMAT',iZMat+ZMatOffset,' = FIX'
    iRow_c = iRow_c+1
    nLambda = nLambda+1
  end do
end do

call mma_deallocate(FragZMat)
call mma_deallocate(FragLabels)
call mma_deallocate(Lines)
call mma_deallocate(SoftLabels)

!                                                                      *
!***********************************************************************
!                                                                      *
! Write out the last line again

Line = 'END OF CONSTRAINTS'
write(Lu_TMP,'(A)') Line
!                                                                      *
!***********************************************************************
!                                                                      *
! Copy the stuff from LU_TMP back to LU_UDC

rewind(Lu_UDC)
rewind(Lu_TMP)
do
  Line = Get_Ln(Lu_TMP)
  !write(u6,*) Line

  write(Lu_UDC,'(A)') trim(Line)
  if (Line(1:4) == 'END ') exit
end do
!write(u6,*) 'iRow_C,nLambda=',iRow_C,nLambda

!                                                                      *
!***********************************************************************
!                                                                      *
! Close the files.

close(Lu_TMP)
close(Lu_UDC)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Fix_UDC

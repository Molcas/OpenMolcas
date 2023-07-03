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

subroutine Merge_Constraints(FileIn1,FileIn2,FileOut,nLambda,iRow_c)

implicit none
character(Len=*), intent(in) :: FileIn1, FileIn2, FileOut
integer, intent(inout) :: nLambda, iRow_c
integer :: Lu1, Lu2, Lu3, isFreeUnit, AixRm
character(Len=180) :: Line, Get_Ln
character(Len=6) :: Tag
logical :: Found, AuxFile
external :: Get_Ln, isFreeUnit, AixRm
integer :: i, j, Lu, iErr

! Open the input files if they exist
Lu1 = 0
if (FileIn1 /= '') then
  call f_Inquire(FileIn1,Found)
  if (Found) then
    Lu1 = isFreeUnit(30)
    call Molcas_Open(Lu1,FileIn1)
  end if
end if

Lu2 = 0
if ((FileIn2 /= '') .and. (FileIn2 /= FileIn1)) then
  call f_Inquire(FileIn2,Found)
  if (Found) then
    Lu2 = isFreeUnit(30)
    call Molcas_Open(Lu2,FileIn2)
  end if
end if

! Set the output file as a third file, or a temporary one
if (((FileOut == FileIn1) .and. (Lu1 /= 0)) .or. ((FileOut == FileIn2) .and. (Lu2 /= 0))) then
  AuxFile = .true.
else
  AuxFile = .false.
end if

nLambda = 0
iRow_c = 0

! If there are no input files, only delete the output file if it exists
if ((Lu1 == 0) .and. (Lu2 == 0)) then
  if (.not. AuxFile) then
    call f_Inquire(FileOut,Found)
    if (Found) iErr = AixRm(FileOut)
  end if
  return
! Otherwise, open it
else
  Lu3 = isFreeUnit(30)
  if (AuxFile) then
    call Molcas_Open(Lu3,'purge')
  else
    call Molcas_Open(Lu3,FileOut)
  end if
end if

! Copy the constraints to the output file, counting the lines
! Copy only from existing files
do i=1,2
  if (i == 1) Tag = 'VALUES'
  if (i == 2) Tag = 'END'

  do j=1,2
    if (j == 1) Lu = Lu1
    if (j == 2) Lu = Lu2
    if (Lu == 0) cycle
    Line = adjustl(Get_Ln(Lu))
    call UpCase(Line)
    do while (Line(1:4) /= Tag(1:4))
      if (index(Line,'&') == 0) then
        iRow_c = iRow_c+1
        if (i == 2) nLambda = nLambda+1
        if (index(Line,'=') == 0) call FixEqualSign(Line,Lu)
      end if
      write(Lu3,100) trim(Line)
      Line = adjustl(Get_Ln(Lu))
      call UpCase(Line)
    end do
  end do

  write(Lu3,100) trim(Tag)
end do
iRow_c = iRow_c+1

! Close files, and copy temporary if used
if (Lu1 /= 0) close(Lu1)
if (Lu2 /= 0) close(Lu2)
if (Lu3 /= 0) close(Lu3)
if (AuxFile) call fCopy('purge',FileOut,iErr)

100 format(A)

end subroutine Merge_Constraints

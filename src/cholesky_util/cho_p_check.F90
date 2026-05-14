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

subroutine Cho_P_Check(irc)
!
! Purpose: check settings for parallel Cholesky.
!
!-TODO/FIXME: The features not allowed in parallel execution of the
!             Cholesky decomposition should be implemented later. Thus,
!             this subroutine is, in effect, a TODO-list.

use Para_Info, only: Is_Real_Par, nProcs
use Cholesky, only: Cho_AdrVec, Cho_DecAlg, Cho_Fake_Par, Cho_IntChk, Cho_Real_Par, Cho_ReOrd, Cho_SimRI, Cho_SScreen, &
                    Cho_TstScreen, IFCSew, LuPri, MxShPr, RstCho, RstDia
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp) :: WriteBlank

irc = 0
WriteBlank = .true.

if (Cho_Real_Par) then ! TRUE PARALLEL
  if ((Cho_DecAlg /= 4) .and. (Cho_DecAlg /= 5) .and. (Cho_DecAlg /= 6)) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Only possible parallel Cholesky decomposition algorithm is "PARAllel".'
    write(Lupri,'(A,I3,A)') 'Resetting Cho_DecAlg from ',Cho_DecAlg,' to 5 (parallel two-step algorithm),'
    Cho_DecAlg = 5
  end if
  if (MxShPr /= 1) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Max. number of shell pair distributions calculated in each pass is 1 for parallel Cholesky.'
    write(Lupri,'(A,I6,A)') 'Resetting MxShPr from ',MxShPr,' to 1'
    MxShPr = 1
  end if
  if (Cho_IntChk) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'You have requested integral checking.'
    write(Lupri,'(A)') 'Integral checking is not possible for parallel Cholesky.'
    irc = irc+1
  end if
  if (RstDia .or. RstCho) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    if (RstDia) then
      write(Lupri,'(A)') 'You have requested diagonal restart.'
      irc = irc+1
    end if
    if (RstCho) then
      write(Lupri,'(A)') 'You have requested decomposition restart.'
      irc = irc+1
    end if
    write(Lupri,'(A)') 'Restart is not possible for parallel Cholesky.'
  end if
  if (Cho_ReOrd) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Vector reordering is not possible for parallel Cholesky.'
    irc = irc+1
  end if
  if (Cho_AdrVec /= 1) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Address mode for vector I/O must be word-addressable for parallel Cholesky.'
    write(Lupri,'(A,I4,A)') 'Resetting Cho_AdrVec from ',Cho_AdrVec,' to 1'
    Cho_AdrVec = 1
  end if
  if (IfcSew /= 2) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Seward interface must be directly in reduced sets for parallel Cholesky.'
    write(Lupri,'(A,I4,A)') 'Resetting IfcSew from ',IfcSew,' to 2'
    IfcSew = 2
  end if
  if (Cho_TstScreen) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Test of subtraction screening is not possible for parallel Cholesky.'
    write(Lupri,'(A)') 'Turning Cho_TstScreen off.'
    Cho_TstScreen = .false.
  end if
  if (Cho_SScreen) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Subtraction screening is not possible for parallel Cholesky.'
    irc = irc+1
  end if
  if (Cho_SimRI) then
    if (WriteBlank) then
      write(Lupri,*)
      WriteBlank = .false.
    end if
    write(Lupri,'(A)') 'Simulation of RI is not possible for parallel Cholesky.'
    irc = irc+1
  end if
else
  if (CHO_FAKE_PAR .and. (nProcs > 1) .and. Is_Real_Par()) then ! FAKE PARALLEL
    if (Cho_ReOrd) then
      if (WriteBlank) then
        write(Lupri,*)
        WriteBlank = .false.
      end if
      write(Lupri,'(A)') 'Vector reordering is not possible for parallel Cholesky.'
      irc = irc+1
    end if
  end if
end if

end subroutine Cho_P_Check

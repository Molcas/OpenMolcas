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

subroutine CHO_MCA_DRV()
!
! Purpose: MOLCAS interface to Cholesky decomposition driver.

#ifdef _DEBUGPRINT_
use Gateway_Info, only: CutInt, ThrInt
#endif
use Cholesky, only: HaltIt, Lupri, MySP, nShell
use stdalloc, only: mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ICODE, irc
real(kind=wp) :: THRAO
#ifdef _DEBUGPRINT_
real(kind=wp) :: CUTINT1, CUTINT2, THRINT1, THRINT2
#endif
logical(kind=iwp) :: DOFOCK, DOGRAD, INDEXATION
character(len=*), parameter :: SECNAM = 'CHO_MCA_DRV'

call STATUSLINE('Seward: ','Cholesky decomposition of ERIs')

! Initialize integral program (this does some memory
! allocations; thus, DO NOT move this.
! --------------------------------------------------

#ifdef _DEBUGPRINT_
CUTINT1 = CutInt
THRINT1 = ThrInt
#endif

call Set_Basis_Mode('Valence')
call Setup_iSD()
NSHELL = -1
INDEXATION = .true.
THRAO = Zero
DOFOCK = .false.
DOGRAD = .false.
call SETUP_INTS(NSHELL,INDEXATION,THRAO,DOFOCK,DOGRAD)

#ifdef _DEBUGPRINT_
CUTINT2 = CutInt
THRINT2 = ThrInt
write(LUPRI,*) SECNAM,': CutInt before Setup_Ints: ',CUTINT1
write(LUPRI,*) SECNAM,': CutInt after  Setup_Ints: ',CUTINT2
write(LUPRI,*) SECNAM,': ThrInt before Setup_Ints: ',THRINT1
write(LUPRI,*) SECNAM,': ThrInt after  Setup_Ints: ',THRINT2
if ((CUTINT2 /= CUTINT1) .or. (THRINT2 /= THRINT1)) call CHO_QUIT('Initialization error in '//SECNAM,102)
#endif

! Start the Cholesky decomposition program.
! -----------------------------------------

ICODE = 0
call CHO_DRV(ICODE)
if (ICODE /= 0) then
  write(LUPRI,*) SECNAM,': decomposition driver returned code ',ICODE
  call CHO_QUIT('Decomposition failed!',104)
end if

! Finalize integral program.
! --------------------------

call TERM_INTS()

! Halt execution if requested.
! ----------------------------

if (HALTIT) then
  write(LUPRI,*) SECNAM,': halting execution after decomposition as requested...'
  call GASYNC()
  call CHO_QUIT('End of Test (in '//SECNAM//')',100)
end if

call GASYNC()
call Free_iSD()

call mma_deallocate(MySP,safe='*')
call Cho_X_dealloc(irc)

end subroutine CHO_MCA_DRV

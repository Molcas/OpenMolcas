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

subroutine SINGLE_ANISO_OPEN(IReturn)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: input_to_read, nDir, nDirZee, nH, nMult, NSS, NSTATE, nT, nTempMagn
logical(kind=iwp) :: GRAD, ifrestart
character(len=180) :: input_file_name

!-----------------------------------------------------------------------
!* initializations
ireturn = 0

! check for the "restart" option:
#ifdef _DEBUGPRINT_
write(u6,*) 'Enter restart_check'
#endif
call restart_check(ifrestart,input_to_read,input_file_name,nT,nH,nTempMagn,nDir,nDirZee,nMult,GRAD)
#ifdef _DEBUGPRINT_
write(u6,*) 'Exit restart_check'
#endif

if (ifrestart) then
  call restart_sa(input_to_read,input_file_name,nss,nstate)
else ! not a restart job -- take all data form RUNFILE
  call fetch_data_RunFile_init(nss,nstate)
end if !Ifrestart

write(u6,'(A)') repeat('@',95)
write(u6,'(A)') '   SINGLE_ANISO (OPEN)'
write(u6,'(A)') '(last updated on 12-March-2018)'
write(u6,'(A)') '   New features: '
write(u6,*)
write(u6,'(A)') '1.  Calculation of the SIGN of the product gX * gY * gZ for any moment;'
write(u6,'(A)') '2.  Calculation of the parameters of the Crystal-Field for lanthanides (CRYS).'
write(u6,'(A)') '3.  Automatic generation of various plot: (PLOT)'
write(u6,'(A)') '     -- Powder magnetic susceptibilty:  XT=f(T)'
write(u6,'(A)') '     -- Powder molar magnetization:  M=f(H,T)'
write(u6,'(A)') '4.  Support for various restart options: (REST)'
write(u6,'(A)') '     -- from  $Project.rassi.h5 file.'
write(u6,'(A)') '     -- from  $Project.aniso (binary) file.'
write(u6,'(A)') '     -- from  ANISOINPUT (ascii) file.'
write(u6,'(A)') '     -- from  $Project.RunFile file.'
write(u6,'(A)') '5.  RASSI was adjusted to provide more accurate tunnelling gaps between near-degenerate states'
write(u6,'(A)') 'Check the MOLCAS manual for details and input examples.'
write(u6,*)
write(u6,'(A)') repeat('@',95)
call xFlush(u6)

call SINGLE_ANISO2(nH,nT,nTempMagn,nDir,nDirZee,nss,nstate,nMult,input_file_name,Ifrestart,IReturn,GRAD)

return

end subroutine SINGLE_ANISO_OPEN

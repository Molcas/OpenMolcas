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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine ClsFls_SCF()
!***********************************************************************
!                                                                      *
!     purpose: Close files after SCF calculations                      *
!                                                                      *
!***********************************************************************

#ifdef _HDF5_
use mh5, only: mh5_close_file
use SCFWfn, only: wfn_fileid
#endif
use InfSCF, only: DoCholesky, DSCF
use SCFFiles, only: LuDel, LuDGd, LuDSt, LuGrd, LuOSt, LuTSt, Lux, Luy
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iRC

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

!---  close two-electron integral file --------------------------------*
if ((.not. DSCF) .and. (.not. DoCholesky)) then
  iRc = -1
  call ClsOrd(iRc)
  if (iRc /= 0) then
    write(u6,*) 'ClsFls: Error closing ORDINT'
    call Abend()
  end if
end if

!---  close DNSMAT, dVxcdR, TWOHAM and GRADIENT -----------------------*
call DaClos(LuDSt)
call DaClos(LuOSt)
call DaClos(LuTSt)
call DaClos(LuGrd)

!---  close 2nd order updatefiles -------------------------------------*
call DaClos(LuDGd)
call DaClos(Lux)
call DaClos(LuDel)
call DaClos(Luy)

#ifdef _HDF5_
call mh5_close_file(wfn_fileid)
#endif

end subroutine ClsFls_SCF

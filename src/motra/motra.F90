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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1999, Roland Lindh                                     *
!***********************************************************************

subroutine Motra(ireturn)
!***********************************************************************
!                                                                      *
! Objective: AO to MO integral transformation                          *
!                                                                      *
! Modify one-electron integrals to use dynamic memory allocation.      *
! R. Lindh, March 1999.                                                *
!                                                                      *
!**** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

#ifdef _HDF5_QCM_
use hdf5_utils, only: ijklname, file_id, hdf5_close, hdf5_create, hdf5_exit, hdf5_init
use Cholesky, only: tv2disk
#endif
use motra_global, only: CMO, HOne, iCTonly, iDoInt, ihdf5, iOneOnly, iPrint, Kine, nTot2, Ovlp
use stdalloc, only: mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: irc
real(kind=wp) :: tcpu_reo, TCR1, TCR2, TWR1, TWR2
logical(kind=iwp) :: DoCholesky, Do_int
integer(kind=iwp), external :: iPrintLevel

!----------------------------------------------------------------------*
! ( Dynamic work area has been allocated in Start() )                  *
!----------------------------------------------------------------------*
!call IniMem()
!----------------------------------------------------------------------*
! Run through the input section                                        *
!----------------------------------------------------------------------*
call init_motra()
if (iPrintLevel(-1) <= 0) iPrint = -1
call InpCtl_Motra()

!> initialize HDF5 interface
#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  !> enable HDF5 support and open the file ijklname
  call hdf5_init()
  call hdf5_create(ijklname,file_id(1))
end if
#endif

!----------------------------------------------------------------------*
! Cholesky check
call DecideOnCholesky(DoCholesky)
!----------------------------------------------------------------------*
! Use of MOTRA only for AO-->MO transf of the Cholesky vectors
if (iCTonly == 1) then
  if (.not. DoCholesky) then
    write(u6,*) '      Warning! This is not RI/CD calculation: '
    write(u6,*) '                      keyword CTonly ignored! '
  else
#   ifdef _HDF5_QCM_
    if ((ihdf5 == 1) .and. (tv2disk /= 'KPQ')) then
      write(u6,*) ' Transformed Cholesky vectors cannot be  written as (pq,K) in HDF5 file as of now. Activate the KPQ option to ' &
                  //'store them or disable the HDF5 option.'
      call Abend()
    end if
#   endif
    write(u6,*)
    write(u6,*) '      ... Skipping MoTRA of ERIs ...'
    write(u6,*)
    write(u6,*) '      ... but Cholesky vectors will be MoTRA.'
    write(u6,*)
    ! Cholesky vectors in HDF5 must be stored as KPQ format (for now)
    Do_int = .false.
    if (iDoInt == 1) Do_int = .true.
    call Cho_MOtra(CMO,nTot2,Do_int,ihdf5)
    iOneOnly = 666
  end if
else
  !--------------------------------------------------------------------*
  ! Preliminary step for integral transformation with CD
  if (DoCholesky) then
    call CWTIME(TCR1,TWR1)
    call Cho_X_init(irc,Zero)
    if (irc /= 0) then
      write(u6,*) ' In MoTRA : Cho_X_Init returned non-zero rc = ',irc
      call Abend()
    end if
    call Cho_X_ReoVec(irc) ! get (if not there) CD vecs full stor
    if (irc /= 0) then
      write(u6,*) ' In MoTRA : Cho_X_ReoVec returned non-zero rc = ',irc
      call Abend()
    end if
    call Cho_X_final(irc)
    call CWTIME(TCR2,TWR2)
    tcpu_reo = (TCR2-TCR1)
    write(u6,*)
    write(u6,*) '      Reordering Cholesky vectors to full storage.'
    write(u6,*) '       Elapsed time for the reordering : ',tcpu_reo
    write(u6,*) '       CPU time for the reordering     : ',tcpu_reo
    write(u6,*)
  end if
end if

!----------------------------------------------------------------------*
! Transform the one-electron integrals                                 *
!----------------------------------------------------------------------*
call Tr1Ctl(Ovlp,HOne,Kine,CMO)
!----------------------------------------------------------------------*
! Transform the two-electron integrals                                 *
!----------------------------------------------------------------------*
if (iOneOnly == 0) call Tr2Ctl(CMO)
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
write(u6,*)
call mma_deallocate(CMO)
call mma_deallocate(Ovlp)
call mma_deallocate(Kine)
call mma_deallocate(HOne)

#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  !> close the file ijkl.h5 and turn off HDF5 support.
  call hdf5_close(file_id(1))
  call hdf5_exit()
end if
#endif

call FastIO('STATUS')
ireturn = 0

return

end subroutine Motra

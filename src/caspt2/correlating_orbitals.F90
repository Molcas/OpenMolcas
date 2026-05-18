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

subroutine correlating_orbitals()
!***********************************************************************
!
! Alter the correlating orbital space by various freeze-delete schemes.
!
!***********************************************************************

use InputData, only: Input
use caspt2_global, only: EMP2, LUONEM, NCMO
use caspt2_module, only: BNAME, iAd1m, IfChol, IfQCAN, nAsh, nBas, nBSqT, nDel, nFro, nIsh, nSsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, IDISK, iRC, iSkp, iSym, NDPQ, ntri, nUniqAt
real(kind=wp) :: Dummy(1)
real(kind=wp), allocatable :: CMO_X(:), DPQ(:)

call Get_iScalar('Unique atoms',nUniqAt)

! memory to store MOs
NCMO = NBSQT
call MMA_ALLOCATE(CMO_X,NCMO,Label='CMO_X')

! Read the MOs from the LUONEM file
IDISK = IAD1M(1)
call DDAFILE(LUONEM,2,CMO_X,NCMO,IDISK)

! AFreeze CASPT2 calculation, available only with
! Cholesky or RI type integral representation.
if (Input%AFreeze) then
  if (.not. IfChol) then
    call WarningMessage(2,'AFreeze needs Cholesky/RI.')
    call Quit_OnUserError()
  end if
  write(u6,'(A)') ' Additional orbitals will be frozen or deleted'
  write(u6,'(A,18A4)') ' Selected atoms:  ',(Input%namfro(i),i=1,Input%lnfro)
  write(u6,'(A,8I4)') ' Frozen orbitals before selection:    ',(nfro(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Inactive orbitals before selection:  ',(nish(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Secondary orbitals before selection: ',(nssh(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Deleted orbitals before selection:   ',(ndel(i),i=1,nsym)
  ntri = 0
  do isym=1,nsym
    ntri = ntri+(nbas(isym)+nbas(isym)**2)/2
  end do
  NDPQ = ntri
  call MMA_ALLOCATE(DPQ,NDPQ,Label='DPQ')
  call AFreez(NSYM,NBAS,NFRO,NISH,NASH,NSSH,NDEL,BNAME,size(BNAME),INPUT%NAMFRO,INPUT%LNFRO,DPQ,nDPQ,Input%THRFR,Input%THRDE, &
              IFQCAN,CMO_X,NCMO)
  call MMA_DEALLOCATE(DPQ)
  write(u6,'(A,8I4)') ' Frozen orbitals after selection     ',(nfro(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Inactive orbitals after selection   ',(nish(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Secondary orbitals after selection: ',(nssh(i),i=1,nsym)
  write(u6,'(A,8I4)') ' Deleted orbitals after selection:   ',(ndel(i),i=1,nsym)
end if

! LOV CASPT2 calculation, available only with
! Cholesky or RI type integral representation.
if (Input%LovCASPT2) then
  if (.not. IfChol) then
    call WarningMessage(2,'LOV-CASPT2 needs Cholesky/RI.')
    call Quit_OnUserError()
  end if
  if (IFQCAN == 0) then
    call WarningMessage(2,'LOV-CASPT2 needs Canonical Orbitals.')
    call Quit_OnUserError()
  end if
  if ((Input%thr_atm < Zero) .or. (Input%thr_atm >= One)) then
    write(u6,*) ' Threshold out of range! Must be in [0,1[ '
    call Quit_OnUserError()
  end if

  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' Start LovCASPT2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)') ' Frozen orbitals before selection:   ',(nFro(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Inactive orbitals before selection: ',(nIsh(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

  EMP2 = Zero
  call Lov_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BNAME,size(BNAME),nUniqAt,Input%thr_atm,IFQCAN,Input%DoMP2,Input%DoEnv, &
                  Input%VIRA,EMP2,CMO_X,NCMO)

  if (irc /= 0) then
    write(u6,*) 'LovCASPT2 returned rc= ',irc
    call abend()
  end if
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' End LovCASPT2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  iSkp = 0
  do iSym=1,nSym
    iSkp = max(iSkp,nAsh(iSym))
  end do
  if (iSkp < 1) call xquit(0)
  write(u6,'(A,8I4)') ' Going to perform CASPT2 calculation on the active region only.'
  write(u6,'(A,8I4)')
end if

! Frozen Natural Orbital CASPT2 calculation, available only with
! Cholesky or RI type integral representation.
if (Input%FnoCASPT2) then
  if (.not. IfChol) then
    call WarningMessage(2,'FNO-CASPT2 needs Cholesky/RI.')
    call Quit_OnUserError()
  end if
  if ((Input%vFrac < -One) .or. (Input%vFrac > One)) then
    call WarningMessage(2,'FNO-CASPT2 fraction out of range.')
    write(u6,*) ' Requested fraction of DEcorr or NOs must be'
    write(u6,*) ' between -1.0 and 1.0.'
    call Quit_OnUserError()
  end if

  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' Start FNO-CASPT2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  if (Input%vfrac >= Zero) then
    write(u6,'(A,I3,A)') ' NOs specified as ',int(Input%vfrac*100.0_wp),'% of the total virtual space'
  else
    write(u6,'(A,I3,A)') ' NOs specified as ',100-int(abs(Input%vfrac)*100.0_wp),'% of DEcorr '
  end if
  write(u6,'(A,8I4)') ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

  EMP2 = Zero
  call FNO_CASPT2(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,Input%vfrac,IFQCAN,Input%DoMP2,EMP2,CMO_X,NCMO)

  if (irc /= 0) then
    write(u6,*) 'FNO_CASPT2 returned rc= ',irc
    call abend()
  end if
  write(u6,'(A,8I4)')
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' End FNO-CASPT2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)')
end if

if (Input%GhostDelete) then
  if ((Input%ThrGD < Zero) .or. (Input%ThrGD >= One)) then
    write(u6,*) ' GHOST threshold out of range! Must be in [0,1[ '
    call Quit_OnUserError()
  end if

  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' GHOST virtual space removal'
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)') ' Secondary orbitals before selection:',(nSsh(i),i=1,nSym)
  write(u6,'(A,8I4)') ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

  call Delete_Ghosts(irc,nSym,nBas,nFro,nIsh,nAsh,nSsh,nDel,BNAME,nUniqAt,Input%ThrGD,.true.,CMO_X,Dummy)

  if (irc /= 0) then
    write(u6,*) 'Delete_GHOSTS returned rc= ',irc
    call abend()
  end if
  write(u6,'(A,8I4)')
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)')
end if

! Store the MOs on the LUONEM file
IDISK = IAD1M(1)
call DDAFILE(LUONEM,1,CMO_X,NCMO,IDISK)

call MMA_DEALLOCATE(CMO_X)

! we need to force recanonicalization of the orbitals later
IFQCAN = 0

end subroutine correlating_orbitals

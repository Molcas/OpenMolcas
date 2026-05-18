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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CnstFIFAFIMO(MODE)

use caspt2_global, only: CMOPT2, FIFA, FIFA_all, FIFASA_all, FIMO, FIMO_all, OLag, TraFro
use caspt2_module, only: IfChol, IFDW, IFRMS, IFXMS, NBAS, NBSQT, NFRO, NFROT, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: MODE
integer(kind=iwp) :: iSQ, iSym, iTr, nBasI
real(kind=wp), allocatable :: WRK1(:), WRK2(:)

if (IfChol) then
  !! For DF or CD, we already have FIFA and FIMO in AO,
  !! so just do AO -> MO transformation
  call mma_allocate(WRK1,NBSQT,Label='WRK1')
  call mma_allocate(WRK2,NBSQT,Label='WRK2')
  WRK1(:) = Zero
  WRK2(:) = Zero

  iSQ = 0
  iTR = 0
  do iSym=1,nSym
    ! nOrbI = nOrb(iSym)
    nBasI = nBas(iSym)
    !! FIFA
    if (nFroT == 0) then
      if ((MODE == 0) .and. (IFDW .or. IFRMS)) then
        call SQUARE(FIFA(1+iTr),FIFASA_all(1+iSQ),1,nBasI,nBasI)
      else if (MODE == 1) then
        call SQUARE(FIFA(1+iTr),FIFA_all(1+iSQ),1,nBasI,nBasI)
        !write(u6,*) 'fifa in MO'
        !call sqprt(fifa_all(1+isq),nbasi)
      end if
    else
      call SQUARE(FIFA_all(1+iTr),WRK1,1,nBasI,nBasI)
      !write(u6,*) 'fifasa in AO'
      !call sqprt(wrk1,nbasi)
      if ((MODE == 0) .and. (IFDW .or. IFRMS)) then
        !! with the state-average
        !! FIFASA_all will be natural basis
        call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFASA_all(1+iSQ),WRK1,WRK2)
        !write(u6,*) 'fifasa in MO'
        !call sqprt(fifasa_all(1+isq),nbasi)
      else if (MODE == 1) then
        !! with the state-specific or dynamically weighted
        !! FIFA will be quasi-canonical basis
        call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFA_all(1+iSQ),WRK1,WRK2)
        !write(u6,*) 'fifa in MO'
        !call sqprt(fifa_all(1+isq),nbasi)
        !! canonicalize frozen orbitals
        !! still under investigation, but this is something we
        !! should do to obtain "better" orbital enegies for
        !! methods using state-dependent Fock operators.
        !! We actually need to canonicalize frozen and inactive
        !! orbitals simultaneously?
        if (nFroT /= 0) then
          WRK2(1:nBasI*nBasI) = WRK1(1:nBasI*nBasI)
          call DIAFCK(NBAS(ISYM),FIFA_all,1,NFRO(ISYM),TraFro,NBAS(ISYM),CMOPT2,WRK2)
          CMOPT2(1:NBAS(ISYM)*NFRO(ISYM)) = WRK2(1:NBAS(ISYM)*NFRO(ISYM))
          call OLagTrf(2,iSym,NBSQT,CMOPT2,FIFA_all(1+iSQ),WRK1,WRK2)
        end if
      end if
    end if

    !! FIMO
    !if (MODE == 0) then
    if (nFroT == 0) then
      call SQUARE(FIMO(1+iTr),FIMO_all(1+iSQ),1,nBasI,nBasI)
    else
      call SQUARE(FIMO_all(1+iTr),WRK1,1,nBasI,nBasI)
      call OLagTrf(2,iSym,NBSQT,CMOPT2,FIMO_all(1+iSQ),WRK1,OLag)
      !write(u6,*) 'fimo in MO'
      !call sqprt(fimo_all(1+isq),nbasi)
    end if
    !end if
    iSQ = iSQ+nBasI*nBasI
    iTR = iTR+nBasI*(nBasI+1)/2
  end do
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)

  if (IFXMS .and. (.not. IFDW)) FIFASA_all(1:NBSQT) = FIFA_all(1:NBSQT)
else
  if (nFroT == 0) then
    iSQ = 0
    iTR = 0
    do iSym=1,nSym
      ! nOrbI = nOrb(iSym)
      nBasI = nBas(iSym)
      call SQUARE(FIFA(1+iTr),FIFA_all(1+iSQ),1,nBasI,nBasI)
      call SQUARE(FIMO(1+iTr),FIMO_all(1+iSQ),1,nBasI,nBasI)
      iSQ = iSQ+nBasI*nBasI
      iTR = iTR+nBasI*(nBasI+1)/2
    end do
    if (IFXMS .and. (.not. IFDW)) FIFASA_all(1:NBSQT) = FIFA_all(1:NBSQT)
  end if

  !! XDW or RMS case: call after XDWINI
  !if (MODE == 0) then
  !end if

  !! XMS case: call after GRPINI
  !if (MODE == 1) then
  !end if

  !! SS or MS case: call in dens
  !if (MODE == 2) then
  !  if (IFSADREF) then
  !  else
  !  end if
  !end if
end if

return

end subroutine CnstFIFAFIMO

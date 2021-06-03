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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine RdMBPT(ipCMO,lthCMO,ipEOrb,lthEOr)
!***********************************************************************
!                                                                      *
!     Read the MBPTOUT file genereated by the SCF program              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!     modified by:                                                     *
!     M.G. Schuetz                                                     *
!     University of Lund, Sweden, 1995                                 *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MBPT2_Global, only: nBas, nDsto, nnB
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ipCMO, lthCMO, ipEOrb, lthEOr
integer(kind=iwp) :: i, ipEOrb_t, iStart, iStart_t, iSym
real(kind=wp), allocatable :: CMO_t(:)
logical(kind=iwp), parameter :: Debug = .false.
#include "corbinf.fh"
#include "WrkSpc.fh"

! Read nSym, nBas, nOrb, nOcc, nFro, CMO and orbital energies from COMFILE

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nOrb',nOrb,nSym)
call Get_iArray('nIsh',nOcc,nSym)
call Get_iArray('nFro',nFro,nSym)
if (Debug) then
  write(u6,'(A,8I5)') 'nSym:',nSym
  write(u6,'(A,8I5)') 'nBas:',(nBas(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOrb:',(nOrb(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOcc:',(nOcc(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nFro:',(nFro(i),i=1,nSym)
end if
lthCMO = 0
do iSym=1,nSym
  if (nFro(iSym) /= 0) then
    write(u6,*) 'Some orbitals where frozen in the SCF!'
    call Abend()
  end if
  nDel(iSym) = nBas(iSym)-nOrb(iSym)
  nExt(iSym) = nOrb(iSym)-nOcc(iSym)
  nDsto(iSym) = nDel(iSym)
  lthCMO = lthCMO+nBas(iSym)*nOrb(iSym)
end do

call mma_allocate(CMO_t,lthCMO,Label='CMO_t')
call Get_CMO_(CMO_t,lthCMO)
call GetMem('CMO   ','Allo','Real',ipCMO,lthCMO)

! set MO coefficients of the deleted orbitals to zero
! Observe that these are not included at all in the basis
iStart = ipCMO
iStart_t = 1
do iSym=1,nSym
  call dcopy_(nOrb(iSym)*nBas(iSym),CMO_t(iStart_t),1,Work(iStart),1)
  iStart = iStart+nOrb(iSym)*nBas(iSym)
  iStart_t = iStart_t+nOrb(iSym)*nBas(iSym)
  call dcopy_((nBas(iSym)-nOrb(iSym))*nBas(iSym),[Zero],0,Work(iStart),1)
  iStart = iStart+(nBas(iSym)-nOrb(iSym))*nBas(iSym)
end do
call mma_deallocate(CMO_t)

call Get_OrbE_(ipEOrb_t,lthEOr)
nnB = lthEOr
call GetMem('EOrb  ','Allo','Real',ipEOrb,lthEOr)

! set energies of the deleted orbitals to zero

iStart = ipEOrb
iStart_t = ipEOrb_t
do iSym=1,nSym
  call dcopy_(nOrb(iSym),Work(iStart_t),1,Work(iStart),1)
  iStart = iStart+nOrb(iSym)
  iStart_t = iStart_t+nOrb(iSym)

  call dcopy_(nBas(iSym)-nOrb(iSym),[Zero],0,Work(iStart),1)
  iStart = iStart+nBas(iSym)-nOrb(iSym)
end do
call GetMem('EOrb_t','Free','Real',ipEOrb_t,lthEOr)

return

end subroutine RdMBPT

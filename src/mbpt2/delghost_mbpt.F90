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

subroutine DelGHOST_MBPT(ipCMO,ipCMO_t,lthCMO,ipEOrb,ipEOrb_t,lthEOr)

use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ipCMO, lthCMO, ipEOrb, lthEOr
integer(kind=iwp), intent(inout) :: ipCMO_t, ipEOrb_t
#include "Molcas.fh"
integer(kind=iwp) :: i, irc, iStart, iStart_t, iSym, nUniqAt, nZero(8)
character(len=LenIn8) :: UBName(mxBas)
logical(kind=iwp), parameter :: Debug = .false.
#include "mxdim.fh"
#include "namact.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "WrkSpc.fh"

if (.not. DelGHOST) return

if (Debug) then
  write(u6,'(A,8I5)') 'nSym:',nSym
  write(u6,'(A,8I5)') 'nBas:',(nBas(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOrb:',(nOrb(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nOcc:',(nOcc(i),i=1,nSym)
  write(u6,'(A,8I5)') 'nFro:',(nFro(i),i=1,nSym)
end if
do iSym=1,nSym
  nDel(iSym) = nBas(iSym)-nOrb(iSym)
  nExt(iSym) = nOrb(iSym)-nOcc(iSym)-nFro(iSym)
  nDsto(iSym) = nDel(iSym)
  nZero(iSym) = 0
end do

call GetMem('CMO   ','Allo','Real',ipCMO,lthCMO)

call GetMem('EOrb  ','Allo','Real',ipEOrb,lthEOr)

write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A)') ' GHOST virtual space removal'
write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A,8I4)')
write(u6,'(A,8I4)') ' Secondary orbitals before selection:',(nExt(i),i=1,nSym)
write(u6,'(A,8I4)') ' Deleted orbitals before selection:  ',(nDel(i),i=1,nSym)

call Get_iScalar('Unique atoms',nUniqAt)
call Get_cArray('Unique Basis Names',UBName,LenIn8*nnB)
call Delete_GHOSTS(irc,nSym,nBas,nFro,nOcc,nZero,nExt,nDel,UBName,nUniqAt,thr_ghs,.false.,Work(ipCMO_t),Work(ipEOrb_t))

if (irc /= 0) then
  write(u6,*) 'Delete_GHOSTS returned rc= ',irc
  call abend()
end if
write(u6,'(A,8I4)')
write(u6,'(A)') '-------------------------------------------------------'
write(u6,'(A,8I4)')
write(u6,'(A,8I4)')

! set MO coefficients of the deleted orbitals to zero
! Observe that these are not included at all in the basis
iStart = ipCMO
iStart_t = ipCMO_t
do iSym=1,nSym
  call dcopy_(nOrb(iSym)*nBas(iSym),Work(iStart_t),1,Work(iStart),1)
  iStart = iStart+nOrb(iSym)*nBas(iSym)
  iStart_t = iStart_t+nOrb(iSym)*nBas(iSym)
  call dcopy_((nBas(iSym)-nOrb(iSym))*nBas(iSym),[Zero],0,Work(iStart),1)
  iStart = iStart+(nBas(iSym)-nOrb(iSym))*nBas(iSym)
end do
call GetMem('CMO_t ','Free','Real',ipCMO_t,lthCMO)

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

end subroutine DelGHOST_MBPT

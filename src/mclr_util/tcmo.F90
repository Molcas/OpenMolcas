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

subroutine TCMO(A,isym,ictl)

use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, ipCM, ipMat, nDens
use input_mclr, only: nBas, nOrb, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: A(*)
integer(kind=iwp), intent(in) :: iSym, iCtl
integer(kind=iwp) :: iip(8), ip(8), iRC, iS, jS, nCMOInv, niCMOInv
integer(kind=iwp), allocatable :: iCMOInv(:)
real(kind=wp), allocatable :: CMOInv(:), Temp(:)

call mma_allocate(Temp,nDens,Label='Temp')
call ReLoad(A,isym,norb,nbas)

! irc used in call later must not be uninitialized
! since then MKL library gets upset...
irc = 0

if (ictl == -1) then

  nCMOInv = 0
  niCMOInv = 0
  ip(:) = 0
  iip(:) = 0
  do iS=1,nSym
    if (nbas(is) == 0) cycle
    ip(iS) = 1+nCMOInv
    iip(iS) = 1+nICMOInv
    nCMOInv = nCMOInv+nBas(is)**2
    niCMOInv = niCMOInv+nBas(is)
  end do
  call mma_allocate(CMOInv,nCMOInv,Label='CMOINV')
  call mma_allocate(iCMOInv,niCMOInv,Label='iCMOINV')

  do iS=1,nSym
    if (nbas(is) == 0) cycle
    CMOINV(ip(is):ip(is)+nbas(is)**2-1) = CMO(ipcm(is):ipcm(is)+nbas(is)**2-1)
    call dgetrf_(nBas(is),nBas(is),CMOINV(ip(is)),nBas(is),iCMOINV(iip(is)),irc)
    if (irc /= 0) call SysAbendMsg('tcmo','DGETRF returns non zero',' ')
  end do

  do iS=1,nSym
    js = Mul(is,isym)
    if (nbas(is)*nbas(js) == 0) cycle
    call dgetrs_('T',nbas(is),nbas(js),CMOINV(ip(is)),nBas(is),iCMOINV(iip(is)),A(ipMat(is,js)),nBas(is),irc)
    if (irc /= 0) call SysAbendMsg('tcmo','DGETRS returns non zero',' ')
    call DGETMO(A(ipMat(is,js)),nBas(is),nbas(is),nbas(js),Temp,nbas(js))
    call dgetrs_('T',nbas(js),nbas(is),CMOINV(ip(js)),nBas(js),iCMOINV(iip(js)),Temp,nBas(js),irc)
    if (irc /= 0) call SysAbendMsg('tcmo','DGETRS returns non zero',' ')
    call DGETMO(Temp,nBas(js),nbas(js),nbas(is),A(ipMat(is,js)),nbas(is))
  end do

  call mma_deallocate(CMOInv)
  call mma_deallocate(iCMOInv)

else if (ictl == 1) then

  do iS=1,nSym
    js = Mul(is,isym)
    if (nBas(is)*nBas(js) == 0) cycle
    call DGEMM_('T','N',nOrb(iS),nBas(jS),nBas(iS),One,CMO(ipCM(iS)),nBas(is),A(ipmat(is,js)),nBas(iS),Zero,Temp,nOrb(iS))
    call DGEMM_('N','N',nOrb(is),nOrb(jS),nBas(jS),One,Temp,nOrb(iS),CMO(ipCM(jS)),nBas(jS),Zero,A(ipMat(iS,jS)),nOrb(iS))
  end do

else if (ictl == -2) then

  do iS=1,nSym
    js = Mul(is,isym)
    if (nBas(is)*nBas(js) == 0) cycle
    call DGEMM_('N','N',nBas(iS),nOrb(jS),nOrb(iS),One,CMO(ipCM(iS)),nBas(is),A(ipmat(is,js)),nOrb(iS),Zero,Temp,nBas(iS))
    call DGEMM_('N','T',nBas(is),nBas(jS),nOrb(jS),One,Temp,nBas(iS),CMO(ipCM(jS)),nBas(jS),Zero,A(ipMat(iS,jS)),nBas(iS))
  end do

else

  write(u6,*) 'Oink'
  call SysHalt('tcmo')

end if

call mma_deallocate(Temp)

end subroutine TCMO

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

subroutine OLagFinal(nOLag,nTrf,OLagLoc,Trf)

use caspt2_global, only: CMOPT2, OLagFull, WLag
use caspt2_module, only: IFMSCOUP, iRlxRoot, JSTATE, NBAS, NBAST, NBSQT, NBTRI, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nOLag, nTrf
real(kind=wp), intent(inout) :: OLagLoc(nOLag)
real(kind=wp), intent(in) :: Trf(nTrf)
integer(kind=iwp) :: iBasI, iBasSq, iBasTr, iSym, jBasI, liBasSq, liBasSq2, liBasTr, nBasI
real(kind=wp), allocatable :: WLagLoc(:), WRK(:)

call mma_allocate(WRK,NBSQT,Label='WRK')
call mma_allocate(WLagLoc,NBSQT,Label='WLagLoc')

if (NBSQT /= nOLag) then
  write(u6,*) 'NBSQT /= nOLag in OLagFinal'
  call abend()
end if

WLagLoc(1:NBSQT) = Half*OLagLoc(1:nOLag)
!write(u6,*) 'Wlag square'
!call sqprt(wlag,nbast)

!! W(MO) -> W(AO) using the quasi-canonical orbitals
!! No need to back transform to natural orbital basis
call DGemm_('N','N',nBasT,nBasT,nBasT,One,CMOPT2,nBasT,WLagLoc,nBasT,Zero,WRK,nBasT)
call DGemm_('N','T',nBasT,nBasT,nBasT,One,WRK,nBasT,CMOPT2,nBasT,Zero,WLagLoc,nBasT)

!! square -> triangle for WLag(AO)
WRK(:) = WLagLoc(:)
iBasTr = 1
iBasSq = 1
do iSym=1,nSym
  nBasI = nBas(iSym)
  liBasTr = iBasTr
  liBasSq = iBasSq
  do iBasI=1,nBasI
    do jBasI=1,iBasI
      liBasSq = iBasSq+iBasI-1+nBasI*(jBasI-1)
      if (iBasI == jBasI) then
        WLagLoc(liBasTr) = WRK(liBasSq)
      else
        liBasSq2 = iBasSq+jBasI-1+nBasI*(iBasI-1)
        WLagLoc(liBasTr) = WRK(liBasSq)+WRK(liBasSq2)
      end if
      liBasTr = liBasTr+1
    end do
  end do
  iBasTr = iBasTr+nBasI*(nBasI+1)/2
  iBasSq = iBasSq+nBasI*nBasI
end do
! accumulate W Lagrangian only for MS,XMS,XDW,RMS,
! but not for SS-CASPT2
if ((jState == iRlxRoot) .or. IFMSCOUP) WLag(1:NBTRI) = WLag(1:NBTRI)+WLagLoc(1:NBTRI)
call mma_deallocate(WLagLoc)

!! Transform quasi-canonical -> natural MO basis
!! orbital Lagrangian
call DGemm_('N','N',nBasT,nBasT,nBasT,One,Trf,nBasT,OLagLoc,nBasT,Zero,WRK,nBasT)
call DGemm_('N','T',nBasT,nBasT,nBasT,One,WRK,nBasT,Trf,nBasT,Zero,OLagLoc,nBasT)
!! sufficient only for active
nBasI = nBas(1)
WRK(1:nBasI**2) = OLagLoc(1:nBasI**2)
call DGeSub(WRK,nBas(1),'N',WRK,nBas(1),'T',OLagLoc,nBas(1),nBas(1),nBas(1))
! accumulate orbital Lagrangian only for MS,XMS,XDW,RMS,
! but not for SS-CASPT2
if ((jState == iRlxRoot) .or. IFMSCOUP) OLagFull(1:nOLag) = OLagFull(1:nOLag)+OLagLoc(1:nOLag)

call mma_deallocate(WRK)

end subroutine OLagFinal

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

subroutine Cho_Fock_MoTra(nSym,nBas,nFro,DLT,DSQ,FLT,nFLT,FSQ,ExFac)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nFLT
! TODO: fix intent of these arrays after removing "ip_of_Work"
real(kind=wp), intent(inout) :: DLT(*), DSQ(*), FLT(nFLT), FSQ(*)
real(kind=wp), intent(in) :: ExFac
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ipd, ipDLT, ipDSQ, ipFLT, ipFSQ, ipMOs, ipV, irc, ja, jaa, MOdim, nDen, NScreen, NumV, nXorb(8)
real(kind=wp) :: ChFracMem, dFKmat, dmpk, Thr, Ymax
integer(kind=iwp), external :: ip_of_Work

!****************************************************************
! CALCULATE AND RETURN FMAT DUE TO FROZEN ORBITALS ONLY
!****************************************************************

NScreen = 10
dmpK = 0.1_wp
dFKmat = Zero
nDen = 1
call IZero(nXorb,nSym)

! Initialize Cholesky information

ChFracMem = Zero
call CHO_X_INIT(irc,ChFracMem)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Init returns error code ',irc
  call AbEnd()
end if

ipDSQ = ip_of_Work(DSQ(1))  ! not needed on exit

MOdim = 0
do i=1,nSym
  MOdim = MOdim+nBas(i)**2
end do
call GETMEM('choMOs','allo','real',ipMOs,MOdim)

ipd = ipDSQ
ipV = ipMOs
do i=1,nSym
  if (nBas(i) > 0) then
    Ymax = Zero
    do ja=1,nBas(i)
      jaa = ipd-1+nBas(i)*(ja-1)+ja
      Ymax = max(Ymax,Work(jaa))
    end do
    Thr = 1.0e-8_wp*Ymax
    call CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),NumV,Thr,irc)
    if (irc /= 0) then
      write(u6,*) 'Cho_Fock_Motra: CD_incore returns rc ',irc
      call AbEnd()
    end if

    if (NumV /= nFro(i)) then
      write(u6,'(a,i6,a,i6,a,i6,a,i6,a,i6)') 'Warning! Cho_Fock_Motra: nr of Frozen orbitals from the decomposition of the '// &
                                             'density matrix is ',numV,' in symm. ',i,'; Expected value = ',nFro(i), &
                                             '; Max diagonal of the density in symm. ',i,' is equal to ',Ymax
    end if

  end if

  ipd = ipd+nBas(i)**2
  ipV = ipV+nBas(i)**2
end do

ipDLT = ip_of_Work(DLT(1))
ipFLT = ip_of_Work(FLT(1))
ipFSQ = ip_of_Work(FSQ(1))  ! not needed on exit

call CHO_LK_SCF(irc,nDen,[ipFLT],[ipFSQ],nXorb,nFro,[ipMOs],[ipDLT],Half*ExFac,NScreen,dmpk,dFKmat)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_LK_scf returns error code ',irc
  call AbEnd()
end if

call GADSUM(FLT,nFLT)

call GETMEM('choMOs','free','real',ipMOs,MOdim)

! Finalize Cholesky information

call CHO_X_FINAL(irc)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Final returns error code ',irc
  write(u6,*) 'Try recovery -- continue.'
end if

return

end subroutine Cho_Fock_MoTra

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

subroutine Cho_Fock_MoTra(nSym,nBas,nFro,DLT,DSQ,W_FLT,nFLT,FSQ,ExFac)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6
use Data_Structures, only: DSBA_type, Allocate_DSBA, Deallocate_DSBA

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nFLT
! TODO: fix intent of these arrays after removing "ip_of_Work"
real(kind=wp), intent(inout) :: DLT(*), DSQ(*), W_FLT(nFLT), FSQ(*)
real(kind=wp), intent(in) :: ExFac
integer(kind=iwp) :: i, ipd, ipDLT(1), ipFSQ(1), ipMOs(1), irc, ja, jaa, MOdim, nDen, NScreen, NumV, nXorb(8)
real(kind=wp) :: ChFracMem, dFKmat, dmpk, Thr, Ymax
real(kind=wp), allocatable :: MOs(:)
integer(kind=iwp), external :: ip_of_Work
type (DSBA_Type) FLT

!****************************************************************
! CALCULATE AND RETURN FMAT DUE TO FROZEN ORBITALS ONLY
!****************************************************************

NScreen = 10
dmpK = 0.1_wp
dFKmat = Zero
call IZero(nXorb,nSym)

! Initialize Cholesky information

ChFracMem = Zero
call CHO_X_INIT(irc,ChFracMem)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Init returns error code ',irc
  call AbEnd()
end if

MOdim = 0
do i=1,nSym
  MOdim = MOdim+nBas(i)**2
end do
call mma_allocate(MOs,MOdim,label='chMOs')

ipd = 1
do i=1,nSym
  if (nBas(i) > 0) then
    Ymax = Zero
    do ja=1,nBas(i)
      jaa = ipd-1+nBas(i)*(ja-1)+ja
      Ymax = max(Ymax,DSQ(jaa))
    end do
    Thr = 1.0e-8_wp*Ymax
    call CD_InCore(DSQ(ipd),nBas(i),MOs(ipd),nBas(i),NumV,Thr,irc)
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
end do

nDen = 1
Call Allocate_DSBA(FLT,nBas,nBas,nSym,Case='TRI',Ref=W_FLT)
ipFSQ(1) = ip_of_Work(FSQ(1))  ! not needed on exit
ipMOs(1) = ip_of_Work(MOs(1))
ipDLT(1) = ip_of_Work(DLT(1))
call CHO_LK_SCF(irc,nDen,FLT,ipFSQ,nXorb,nFro,ipMOs,ipDLT,Half*ExFac,NScreen,dmpk,dFKmat)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_LK_scf returns error code ',irc
  call AbEnd()
end if

call Deallocate_DSBA(FLT)

call GADSUM(W_FLT,nFLT)

call mma_deallocate(MOs)

! Finalize Cholesky information

call CHO_X_FINAL(irc)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Final returns error code ',irc
  write(u6,*) 'Try recovery -- continue.'
end if

return

end subroutine Cho_Fock_MoTra

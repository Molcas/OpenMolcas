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

subroutine Cho_Fock_MoTra(nSym,nBas,nFro,W_DLT,W_DSQ,W_FLT,nFLT,W_FSQ,ExFac)

use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nFLT
! TODO: fix intent of these arrays after removing "ip_of_Work"
real(kind=wp), intent(inout) :: W_DLT(*), W_DSQ(*), W_FLT(nFLT), W_FSQ(*)
real(kind=wp), intent(in) :: ExFac
integer(kind=iwp) :: i, irc, ja, nDen, NScreen, NumV, nXorb(8)
real(kind=wp) :: ChFracMem, dFKmat, dmpk, Thr, Ymax
type (DSBA_Type) FLT(1), KLT(1), MOs, DLT, DSQ

!****************************************************************
! CALCULATE AND RETURN FMAT DUE TO FROZEN ORBITALS ONLY
!****************************************************************

NScreen = 10
dmpK = 0.1_wp
dFKmat = Zero
nXorb(:)=0

! Initialize Cholesky information

ChFracMem = Zero
call CHO_X_INIT(irc,ChFracMem)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Init returns error code ',irc
  call AbEnd()
end if

Call Allocate_DT(DLT,nBas,nBas,nSym,aCase='TRI',Ref=W_DLT)
Call Allocate_DT(DSQ,nBas,nBas,nSym,aCase='REC',Ref=W_DSQ)
Call Allocate_DT(MOs,nBas,nBas,nSym)

do i=1,nSym
  if (nBas(i) > 0) then
    Ymax = Zero
    do ja=1,nBas(i)
      Ymax = max(Ymax,DSQ%SB(i)%A2(ja,ja))
    end do
    Thr = 1.0e-8_wp*Ymax
    call CD_InCore(DSQ%SB(i)%A2,nBas(i),MOs%SB(i)%A2,nBas(i),NumV,Thr,irc)
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

end do

nDen = 1
Call Allocate_DT(FLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FLT)
Call Allocate_DT(KLT(1),nBas,nBas,nSym,aCase='TRI',Ref=W_FSQ) ! not needed on exit

call CHO_LK_SCF(irc,nDen,FLT,KLT,nXorb,nFro,[MOs],[DLT],Half*ExFac,NScreen,dmpk,dFKmat)

if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_LK_scf returns error code ',irc
  call AbEnd()
end if

call Deallocate_DT(KLT(1))
call Deallocate_DT(FLT(1))

call GADSUM(W_FLT,nFLT)

call deallocate_DT(MOs)
call deallocate_DT(DSQ)
call deallocate_DT(DLT)

! Finalize Cholesky information

call CHO_X_FINAL(irc)
if (irc /= 0) then
  write(u6,*) 'Cho_Fock_Motra: Cho_X_Final returns error code ',irc
  write(u6,*) 'Try recovery -- continue.'
end if

return

end subroutine Cho_Fock_MoTra

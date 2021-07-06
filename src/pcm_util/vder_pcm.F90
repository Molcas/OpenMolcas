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

subroutine VDer_PCM(nAt,nTs,nS,AtmC,AtmChg,EF_n,EF_e,Tessera,iSphe,DerTes,DerPunt,DerRad,DerCentr,VDer)

use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, nS, iSphe(*)
real(kind=wp), intent(in) :: AtmC(3,nAt), AtmChg(nAt), EF_n(3,*), EF_e(3,*), Tessera(4,*), DerTes(nTs,nAt,3), &
                             DerPunt(nTs,nAt,3,3), DerRad(nS,nAt,3), DerCentr(nS,nAt,3,3)
real(kind=wp), intent(_OUT_) :: VDer(nTs,*)
integer(kind=iwp) :: iAt, iCoord, idx, iTs, L, Lu
real(kind=wp) :: dCoord, DTessNuc, DVNuc, dX, dY, dZ, Fld_e, Fld_n
integer(kind=iwp), external :: IsFreeUnit

! pcm_solvent very temporary! read the potential derivatives from file
Lu = IsFreeUnit(10)
call Molcas_Open(Lu,'DerPot.dat')
!open(1,file='DerPot.dat',status='old',form='formatted')
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    do iTs=1,nTs
      read(1,*) VDer(iTs,idx)
    end do
  end do
end do
close(1)
! pcm_solvent end
! Loop on atoms and coordinates
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord

    ! (Total) derivative of the electronic + nuclear potential

    do iTs=1,nTs
      L = iSphe(iTs)
      ! Derivative of the representative point
      dX = DerPunt(iTs,iAt,iCoord,1)+DerCentr(L,iAt,iCoord,1)
      dY = DerPunt(iTs,iAt,iCoord,2)+DerCentr(L,iAt,iCoord,2)
      dZ = DerPunt(iTs,iAt,iCoord,3)+DerCentr(L,iAt,iCoord,3)
      ! Distance tessera - nucleus
      DTessNuc = sqrt((Tessera(1,iTs)-AtmC(1,iAt))**2+(Tessera(2,iTs)-AtmC(2,iAt))**2+(Tessera(3,iTs)-AtmC(3,iAt))**2)
      dCoord = Tessera(iCoord,iTs)-AtmC(iCoord,iAt)
      ! Deriv. of the nuclear potential (with fixed repres. point)
      DVNuc = -AtmChg(iAt)*dCoord/DTessNuc**3
      ! Electric field (electrons and nuclei) times the derivative of
      ! the repres. point
      Fld_e = EF_e(1,iTs)*dX+EF_e(2,iTs)*dY+EF_e(3,iTs)*dZ
      Fld_n = EF_n(1,iTs)*dX+EF_n(2,iTs)*dY+EF_n(3,iTs)*dZ
      ! Total deriv. of the potential
      VDer(iTs,idx) = VDer(iTs,idx)-Fld_e+DVNuc+Fld_n
      ! pcm_solvent
      ! In MCLR the electron electric field is not computed, probably because of
      ! DoRys set to False in mclr. Setting DoRys to True causes the program
      ! to stop because the abdata.ascii file is not found.
      if ((iat == 1) .and. (icoord == 1) .and. (its == 1)) &
        write(u6,'("In the loop VDer,Fld_e,DVNuc,Fld_n",4f12.6)') VDer(1,idx),Fld_e,DVNuc,Fld_n
      if ((iat == nat) .and. (icoord == 3) .and. (its == 1)) &
        write(u6,'("In the loop VDer,Fld_e,DVNuc,Fld_n",4f12.6)') VDer(1,idx),Fld_e,DVNuc,Fld_n
      ! pcm_solvent end

    end do
  end do
end do
! pcm_solvent
write(u6,'(a)') 'At the end of DerPot,VDer(1,ind),VDer(nTs,ind)'
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    write(u6,'(2f12.6)') VDer(1,idx),VDer(nTs,idx)
  end do
end do
! pcm_solvent end
return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(DerTes)
  call Unused_real_array(DerRad)
end if

end subroutine VDer_PCM

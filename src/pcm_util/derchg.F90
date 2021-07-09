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

subroutine DerChg(nAt,nAt3,nTs,nS,Eps,IAtm,AtmC,AtmChg,DM,DerMat,Tessera,Q,Qtot,QDer,DerTes,DerPunt,DerCentr,DerRad,Der1,Der2, &
                  VDer,Sphere,ISphe)

use Constants, only: Zero, One, Pi
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nAt3, nTs, nS, IAtm(nAt), ISphe(*)
real(kind=wp), intent(in) :: Eps, AtmC(3,nAt), AtmChg(nAt), DM(nTs,*), Tessera(4,*), Q(2,*), DerTes(nTs,nAt,3), &
                             DerPunt(nTs,nAt,3,3), DerCentr(nS,nAt,3,3), DerRad(nS,nAt,3), VDer(nTs,*), Sphere(4,*)
real(kind=wp), intent(_OUT_) :: DerMat(nTs,*), Qtot(*), QDer(3,nAt,*), Der1(*), Der2(*)
integer(kind=iwp) :: iAt, iCoord, idx, iTs
real(kind=wp) :: Diag, Sc_Cond

Diag = -1.0694_wp*sqrt(PI)
Sc_Cond = (Eps-One)/Eps

! Total charges

do iTs=1,nTs
  Qtot(iTs) = Q(1,iTs)+Q(2,iTs)
end do

! Compute the derivative of PCM matrix

! Loop on atoms and coordinates
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    ! Conductor-like case
    call DMat_CPCM(iAt,iCoord,Eps,nTs,nS,nAt,Diag,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

    ! First product: DerDM * Qtot

    call dgemm_('N','N',nTs,1,nTs,-One,DerMat,nTs,Qtot,nTs,Zero,Der1,nTs)
    ! pcm_solvent
    !if ((iat == 5) .and.( icoord == 1)) then
    !  write(u6,'(a)') 'In DerChg first contribution for 5, 1'
    !  do its=1,nts
    !    write(u6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end

    ! Total deriv. of the potential summed up the quantity already
    ! computed (DerDM*Qtot)

    do iTs=1,nTs
      Der1(iTs) = Der1(iTs)+Sc_Cond*VDer(iTs,idx)
    end do

    ! Last product: - DM^-1 * (V^x + DM^x*q)

    ! pcm_solvent
    !if ((iat == 5) .and. (icoord == 1)) then
    !  write(u6,'(a)') 'In DerChg second contribution for 5, 1'
    !  do its=1,nts
    !    write(u6,'(i4,f20.12)') its,der1(its)
    !  end do
    !end if
    ! pcm_solvent end
    call dgemm_('N','N',nTs,1,nTs,-One,DM,nTs,Der1,nTs,Zero,Der2,nTs)
    do iTs=1,nTs
      QDer(iCoord,iAt,iTs) = Der2(iTs)
    end do
  end do
end do

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(nAt3)
  call Unused_integer_array(IAtm)
  call Unused_real_array(AtmC)
  call Unused_real_array(AtmChg)
  call Unused_real_array(DerRad)
  call Unused_real_array(Sphere)
end if

end subroutine DerChg

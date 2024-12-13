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

subroutine CASPT2_BTAMP(LuGAMMA,iS,jS,kS,lS,nBasI,nBasJ,nBasK,nBasL,iOffAO,nBasT,nOcc,CMOPT2,WRK1,WRK2,G_toc)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, LuGAMMA, nBasI, nBasJ, nBasK, nBasL, iOffAO(*), nBasT, nOcc
real(kind=wp), intent(in) :: CMOPT2(*)
real(kind=wp), intent(_OUT_) :: WRK1(*), WRK2(*), G_toc(*)
integer(kind=iwp) :: iBas, iBas0, iRec, jBas, jBas0, kBas, kBas0, lBas, lBas0, Loc
real(kind=wp), parameter :: SCAL = 0.125_wp

! Transform T_{ij}^{rho sigma} to T_{mu nu}^{rho sigma}
! i- and k-shells correspond to rho and sigma (order?)
! The transformed amplitude is stored as D_{j,l,i,k},
! and it will be sorted in integral_util/pget3.f

! It is possible to reduce operations. The third transformation
! can be postponed until j-shell varies.

do kBas0=1,nBasK
  kBas = iOffAO(kS)+kBas0
  do iBas0=1,nBasI
    iBas = iOffAO(iS)+iBas0
    iRec = iBas+nBasT*(kBas-1)
    !! Read the half-transformed-amplitude
    read(unit=LuGamma,rec=iRec) WRK1(1:nOcc*nOcc)
    !! do the remaining (third and fourth) transformation
    call DGemm_('N','N',nBasJ,nOcc,nOcc,One,CMOPT2(1+iOffAO(jS)),nBasT,WRK1,nOcc,Zero,WRK2,nBasJ)
    Loc = nBasJ*nBasL*(iBas0-1+nBasI*(kBas0-1))
    call DGemm_('N','T',nBasJ,nBasL,nOcc,SCAL,WRK2,nBasJ,CMOPT2(1+iOffAO(lS)),nBasT,Zero,G_toc(1+Loc),nBasJ)
  end do
end do

do lBas0=1,nBasL
  lBas = iOffAO(lS)+lBas0
  do jBas0=1,nBasJ
    jBas = iOffAO(jS)+jBas0
    iRec = jBas+nBasT*(lBas-1)
    !! Read the half-transformed-amplitude
    read(unit=LuGamma,rec=iRec) WRK1(1:nOcc*nOcc)
    !! do the remaining (third and fourth) transformation
    call DGemm_('N','N',nBasI,nOcc,nOcc,One,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,Zero,WRK2,nBasI)
    call DGemm_('N','T',nBasI,nBasK,nOcc,SCAL,WRK2,nBasI,CMOPT2(1+iOffAO(kS)),nBasT,Zero,WRK1,nBasI)
    do kBas0=1,nBasK
      kBas = iOffAO(kS)+kBas0
      do iBas0=1,nBasI
        ibas = iOffAO(iS)+iBas0
        Loc = jBas0-1+nBasJ*(lBas0-1+nBasL*(iBas0-1+nBasI*(kBas0-1)))
        G_toc(1+Loc) = G_toc(1+Loc)+WRK1(iBas0+nBasI*(kBas0-1))
      end do
    end do
  end do
end do

do lBas0=1,nBasL
  lBas = iOffAO(lS)+lBas0
  do iBas0=1,nBasI
    iBas = iOffAO(iS)+iBas0
    iRec = iBas+nBasT*(lBas-1)
    !! Read the half-transformed-amplitude
    read(unit=LuGamma,rec=iRec) WRK1(1:nOcc*nOcc)
    !! do the remaining (third and fourth) transformation
    call DGemm_('N','N',nBasJ,nOcc,nOcc,One,CMOPT2(1+iOffAO(jS)),nBasT,WRK1,nOcc,Zero,WRK2,nBasJ)
    call DGemm_('N','T',nBasJ,nBasK,nOcc,SCAL,WRK2,nBasJ,CMOPT2(1+iOffAO(kS)),nBasT,Zero,WRK1,nBasJ)
    do kBas0=1,nBasK
      do jBas0=1,nBasJ
        Loc = jBas0-1+nBasJ*(lBas0-1+nBasL*(iBas0-1+nBasI*(kBas0-1)))
        G_toc(1+Loc) = G_toc(1+Loc)+WRK1(jBas0+nBasJ*(kBas0-1))
      end do
    end do
  end do
end do

do kBas0=1,nBasK
  kBas = iOffAO(kS)+kBas0
  do jBas0=1,nBasJ
    jBas = iOffAO(jS)+jBas0
    if (jBas >= kBas) then
      iRec = jBas+nBasT*(kBas-1)
    else
      iRec = kBas+nBasT*(jBas-1)
    end if
    !! Read the half-transformed-amplitude
    read(unit=LuGamma,rec=iRec) WRK1(1:nOcc*nOcc)
    !! do the remaining (third and fourth) transformation
    if (jbas >= kbas) then
      call DGemm_('N','N',nBasI,nOcc,nOcc,One,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,Zero,WRK2,nBasI)
    else
      call DGemm_('N','T',nBasI,nOcc,nOcc,One,CMOPT2(1+iOffAO(iS)),nBasT,WRK1,nOcc,Zero,WRK2,nBasI)
    end if
    call DGemm_('N','T',nBasI,nBasL,nOcc,SCAL,WRK2,nBasI,CMOPT2(1+iOffAO(lS)),nBasT,Zero,WRK1,nBasI)
    do lBas0=1,nBasL
      do iBas0=1,nBasI
        Loc = jBas0-1+nBasJ*(lBas0-1+nBasL*(iBas0-1+nBasI*(kBas0-1)))
        G_toc(1+Loc) = G_toc(1+Loc)+WRK1(iBas0+nBasI*(lBas0-1))
      end do
    end do
  end do
end do

return

end subroutine CASPT2_BTAMP

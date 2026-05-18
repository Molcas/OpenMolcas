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

subroutine DWDER(OMGDER,HEFF,SLag)
! Computes the derivative of weight factor in XDW-CASPT2

use caspt2_module, only: DWTYPE, NSTATE, ZETA
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: OMGDER(nState,nState), HEFF(nState,nState)
real(kind=wp), intent(inout) :: SLag(nState,nState)
integer(kind=iwp) :: ilStat, jlStat, klStat
real(kind=wp) :: Dag, Dbg, DERAB, DERAC, DEROMG, Ealpha, Ebeta, Egamma, Factor, Hag, Hbg, Scal

do ilStat=1,nState
  Ebeta = HEFF(ilStat,ilStat)
  do jlStat=1,nState
    Ealpha = HEFF(jlStat,jlStat)

    Factor = Zero
    do klStat=1,nState
      Egamma = HEFF(klStat,klStat)
      if (DWType == 1) then
        Factor = Factor+exp(-zeta*(Ealpha-Egamma)**2)
      else if (DWType == 2) then
        Factor = Factor+exp(-zeta*(Ealpha/HEFF(jlStat,klStat))**2)
      else if (DWType == 3) then
        Dag = abs(Ealpha-Egamma)+1.0e-9_wp
        Hag = abs(HEFF(jlStat,klStat))
        Factor = Factor+exp(-zeta*Dag/(sqrt(Hag)+tiny(Hag)))
      end if
    end do

    DEROMG = OMGDER(ilStat,jlStat)

    !! derivative of alpha-beta
    if (DWType == 1) then
      DERAB = exp(-ZETA*(Ealpha-Ebeta)**2)/Factor
      Scal = -Two*ZETA*DERAB*(Ealpha-Ebeta)*DEROMG
      SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
      SLag(ilStat,ilStat) = SLag(ilStat,ilStat)-Scal
    else if (DWType == 2) then
      DERAB = exp(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
      DERAB = DERAB*Two*DEROMG*ZETA*(Ealpha/HEFF(jlStat,ilStat))**2
      Scal = -DERAB/Ealpha
      SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
      Scal = +DERAB/HEFF(jlStat,ilStat)
      SLag(jlStat,ilStat) = SLag(jlStat,ilStat)+Scal
    else if (DWType == 3) then
      Dag = abs(Ealpha-Ebeta)+1.0e-9_wp
      Hag = abs(HEFF(jlStat,ilStat))
      DERAB = exp(-ZETA*Dag/(sqrt(Hag)+tiny(Hag)))/Factor
      DERAB = -ZETA*DERAB*Dag/(sqrt(Hag)+tiny(Hag))*DEROMG
      Scal = DERAB/Dag
      if (Ealpha-Ebeta <= Zero) Scal = -Scal
      SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
      SLag(ilStat,ilStat) = SLag(ilStat,ilStat)-Scal
      Scal = -DERAB/(sqrt(Hag)+tiny(Hag))/sqrt(Hag)*Half
      if (HEFF(jlStat,ilStat) <= Zero) Scal = -Scal
      SLag(jlStat,ilStat) = SLag(jlStat,ilStat)+Scal
    end if

    !! derivative of alpha-gamma
    do klStat=1,nState
      Egamma = HEFF(klStat,klStat)
      if (DWtype == 1) then
        DERAC = exp(-ZETA*(Ealpha-Ebeta)**2)/(Factor*Factor)*exp(-ZETA*(Ealpha-Egamma)**2)
        Scal = Two*ZETA*DERAC*(Ealpha-Egamma)*DEROMG
        SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
        SLag(klStat,klStat) = SLag(klStat,klStat)-Scal
      else if (DWType == 2) then
        DERAC = exp(-ZETA*(Ealpha/HEFF(jlStat,klStat))**2)/Factor*exp(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
        DERAC = -DERAC*Two*DEROMG*ZETA*(Ealpha/HEFF(jlStat,klStat))**2
        Scal = -DERAC/Ealpha
        SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
        Scal = +DERAC/HEFF(jlStat,klStat)
        SLag(jlStat,klStat) = SLag(jlStat,klStat)+Scal
      else if (DWType == 3) then
        Dbg = abs(Ealpha-Egamma)+1.0e-9_wp
        Hbg = abs(HEFF(jlStat,klStat))
        DERAC = exp(-ZETA*Dag/(sqrt(Hag)+tiny(Hag)))/Factor*exp(-ZETA*Dbg/(sqrt(Hbg)+tiny(Hbg)))/Factor
        DERAC = ZETA*DERAC*Dbg/(sqrt(Hbg)+tiny(Hbg))*DEROMG
        Scal = DERAC/Dbg
        if (Ealpha-Egamma <= Zero) Scal = -Scal
        SLag(jlStat,jlStat) = SLag(jlStat,jlStat)+Scal
        SLag(klStat,klStat) = SLag(klStat,klStat)-Scal
        Scal = -DERAC/(sqrt(Hbg)+tiny(Hbg))/sqrt(Hbg)*Half
        if (HEFF(jlStat,klStat) <= Zero) Scal = -Scal
        SLag(jlStat,klStat) = SLag(jlStat,klStat)+Scal
      end if
    end do
  end do
end do

return

end subroutine DWDER

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

      Subroutine DWDER(OMGDER,HEFF,SLag)

      use definitions, only: wp, iwp
      use caspt2_module, only: NSTATE, DWTYPE, ZETA
      use Constants, only: Zero, Two, Half

      implicit none

      real(kind=wp), intent(in) :: OMGDER(nState,nState),               &
     &  HEFF(nState,nState)
      real(kind=wp), intent(inout) :: SLag(nState,nState)

      integer(kind=iwp) :: ilStat, jlStat, klStat
      real(kind=wp) :: Ebeta, Ealpha, Factor, Egamma, Dag, Hag, DEROMG, &
     &  Scal, DERAB, DERAC, Dbg, Hbg
!
!     Computes the derivative of weight factor in XDW-CASPT2
!
      Do ilStat = 1, nState
        Ebeta = HEFF(ilStat,ilStat)
        Do jlStat = 1, nState
          Ealpha = HEFF(jlStat,jlStat)

          Factor = Zero
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWType == 1) Then
              Factor = Factor + exp(-zeta*(Ealpha - Egamma)**2)
            Else If (DWType == 2) Then
              Factor = Factor                                           &
     &          + exp(-zeta*(Ealpha/HEFF(jlStat,klStat))**2)
            Else If (DWType == 3) Then
              Dag = abs(Ealpha - Egamma) + 1.0e-9_wp
              Hag = abs(HEFF(jlStat,klStat))
              Factor = Factor                                           &
     &          + exp(-zeta*Dag/(sqrt(Hag)+tiny(Hag)))
            End If
          End Do

          DEROMG = OMGDER(ilStat,jlStat)

          !! derivative of alpha-beta
          If (DWType == 1) Then
            DERAB = EXP(-ZETA*(Ealpha-Ebeta)**2)/Factor
            Scal = -Two*ZETA*DERAB*(Ealpha-Ebeta)*DEROMG
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
          Else If (DWType == 2) Then
            DERAB = EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
            DERAB = DERAB*Two*DEROMG*ZETA                               &
     &            *(Ealpha/HEFF(jlStat,ilStat))**2
            Scal  = -DERAB/Ealpha
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            Scal  = +DERAB/HEFF(jlStat,ilStat)
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          Else If (DWType == 3) Then
            Dag = abs(Ealpha - Ebeta) + 1.0e-9_wp
            Hag = abs(HEFF(jlStat,ilStat))
            DERAB = EXP(-ZETA*Dag/(sqrt(Hag)+Tiny(Hag)))/Factor
            DERAB = -ZETA*DERAB*Dag/(sqrt(Hag)+tiny(Hag))               &
     &              *DEROMG
            Scal = DERAB/Dag
            If (Ealpha-Ebeta <= Zero) Scal = -Scal
            SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
            SLag(ilStat,ilStat) = SLag(ilStat,ilStat) - Scal
            Scal  = -DERAB/(sqrt(Hag)+tiny(Hag))/sqrt(Hag)*Half
            If (HEFF(jlStat,ilStat) <= Zero) Scal = -Scal
            SLag(jlStat,ilStat) = SLag(jlStat,ilStat) + Scal
          End If

          !! derivative of alpha-gamma
          Do klStat = 1, nState
            Egamma = HEFF(klStat,klStat)
            If (DWtype == 1) Then
              DERAC = EXP(-ZETA*(Ealpha-Ebeta)**2)/(Factor*Factor)      &
     &               *EXP(-ZETA*(Ealpha-Egamma)**2)
              Scal = Two*ZETA*DERAC*(Ealpha-Egamma)*DEROMG
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
            Else If (DWType == 2) Then
              DERAC = EXP(-ZETA*(Ealpha/HEFF(jlStat,klStat))**2)/Factor &
     &              * EXP(-ZETA*(Ealpha/HEFF(jlStat,ilStat))**2)/Factor
              DERAC =-DERAC*Two*DEROMG*ZETA                             &
     &              *(Ealpha/HEFF(jlStat,klStat))**2
              Scal  = -DERAC/Ealpha
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              Scal  = +DERAC/HEFF(jlStat,klStat)
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            Else IF (DWType == 3) Then
              Dbg = abs(Ealpha - Egamma) + 1.0e-9_wp
              Hbg = abs(HEFF(jlStat,klStat))
              DERAC =EXP(-ZETA*Dag/(sqrt(Hag)+tiny(Hag)))/Factor        &
     &              *EXP(-ZETA*Dbg/(sqrt(Hbg)+tiny(Hbg)))/Factor
              DERAC = ZETA*DERAC*Dbg/(sqrt(Hbg)+tiny(Hbg))              &
     &                *DEROMG
              Scal = DERAC/Dbg
              If (Ealpha-Egamma <= Zero) Scal = -Scal
              SLag(jlStat,jlStat) = SLag(jlStat,jlStat) + Scal
              SLag(klStat,klStat) = SLag(klStat,klStat) - Scal
              Scal  = -DERAC/(sqrt(Hbg)+tiny(Hbg))/sqrt(Hbg)*Half
              If (HEFF(jlStat,klStat) <= Zero) Scal = -Scal
              SLag(jlStat,klStat) = SLag(jlStat,klStat) + Scal
            End If
          End Do
        End Do
      End Do

      Return

      End Subroutine DWDER

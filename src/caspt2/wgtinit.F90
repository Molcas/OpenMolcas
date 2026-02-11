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
! Copyright (C) 2019, Stefano Battaglia                                *
!***********************************************************************
subroutine wgtinit(H)

  use definitions,only:wp,iwp,u6
  use caspt2_global,only:iPrGlb
  use PrintLevel, only: DEBUG, VERBOSE
  use caspt2_global, only: DWGT
  use caspt2_module, only: nState, DWType, IfDW, IfXMS, Zeta

  implicit none


  Real(kind=wp),intent(in) :: H(nState,nState)

  Real(kind=wp) :: Ealpha,Ebeta,Egamma,Dab,Dag,xi_ag,xi_ab,Wtot
  Real(kind=wp) :: Hab,Hag

  Integer(kind=iwp) :: I,J,K

  if (IPRGLB >= DEBUG) then
    write (6,*) ' Entered wgtinit.'
  end if

  ! Initialize array of weights with all zeros
  DWGT(:,:)=00D0

  ! Main loop over all states to compute the weights
  do I = 1,nState

    ! If it is an XDW-CASPT2 calculation, the weights are computed
    if (IFDW .and. (zeta >= 0.0_wp)) then
      Ebeta = H(I,I)
      ! Compute normalization factor Wtot, i.e. the sum of all weights
      do J = 1,nState
        Ealpha = H(J,J)
        Wtot = 0.0_wp
        do K = 1,nState
          Egamma = H(K,K)

          ! original XDW-CASPT2, xi = Dab^2
          if (DWType == 1) then
            xi_ag = (Ealpha - Egamma)**2
            ! new XDW-CASPT2, xi = (Haa/Hab)^2
          else if (DWType == 2) then
            xi_ag = (Ealpha/H(J,K))**2
            ! new XDW-CASPT2, xi = Dab/sqrt(Hab), which is DWType == 3
          else
            ! add a small positive constant to numerator to avoid 0/0
            Dag = abs(Ealpha - Egamma) + 1.0e-9_wp
            Hag = abs(H(J,K))
            ! add the smallest value of the same type as Hag to the
            ! denominator to avoid division by 0 which can lead to
            ! segfault with certain compilers
            xi_ag = Dag/(sqrt(Hag)+tiny(Hag))
          end if

          Wtot = Wtot + exp(-zeta*xi_ag)
        end do

        ! original XDW-CASPT2, xi = Dab^2
        if (DWType == 1) then
          xi_ab = (Ealpha - Ebeta)**2
          ! new XDW-CASPT2, xi = (Haa/Hab)^2
        else if (DWType == 2) then
          xi_ab = (Ealpha/H(J,I))**2
          ! new XDW-CASPT2, xi = Dab/sqrt(Hab), which is DWType == 3
        else
          ! add a small positive constant to numerator to avoid 0/0
          Dab = abs(Ealpha - Ebeta) + 1.0e-9_wp
          Hab = abs(H(J,I))
          ! add the smallest value of the same type as Hag to the
          ! denominator to avoid division by 0 which can lead to
          ! segfault with certain compilers
          xi_ab = Dab/(sqrt(Hab)+tiny(Hab))
        end if

        DWGT(I,J) = exp(-zeta*xi_ab)/Wtot

      end do

      ! If it is an XMS-CASPT2 calculation, all the weights are equal,
      ! i.e. they all are 1/nState
    else if (IFXMS .and. (.not. IFDW)) then
      call dcopy_(nState**2, [1.0_wp/nState],0,DWGT,1)

      ! If it is a normal MS-CASPT2, RMS-CASPT2 or a (X)DW-CASPT2 with zeta->infinity
      ! the weight vectors are the standard unit vectors e_1, e_2, ...
    else
      DWGT(I,I) = 1.0d0
    end if

    ! End of loop over states
  end do

  ! In case it is a XDW calculation, print out the weights
  if (IFDW .and. (IPRGLB >= VERBOSE)) then
    write (u6,*) ' Weights calculated with <I|H|I>:'
    call prettyprint(DWGT,nState,nState)
  end if

  return
end

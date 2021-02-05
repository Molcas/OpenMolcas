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

      use definitions, only: wp,iwp,u6
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb

      implicit none

#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      Real(kind=wp),intent(in) :: H(nState,nState)

      Real(kind=wp) :: Ealpha,Ebeta,Egamma,Dab,Dag,xi_ag,xi_ab,Wtot
      Real(kind=wp) :: Hab,Hag

      Integer(kind=iwp) :: I,J,K,IJ


      if (IPRGLB.GE.DEBUG) then
        write(6,*)' Entered wgtinit.'
      end if

      ! Initialize array of weights with all zeros
      call dcopy_(nState**2,[0.0_wp],0,WORK(LDWGT),1)

      ! Main loop over all states to compute the weights
      do I=1,nState

        ! If it is an XDW-CASPT2 calculation, the weights are computed
        if (IFDW .and. (zeta .ge. 0.0_wp)) then
          Ebeta = H(I,I)
          ! Compute normalization factor Wtot, i.e. the sum of all weights
          do J=1,nState
            Ealpha = H(J,J)
            Wtot = 0.0_wp
            ! pref = 0.001D0
            do K=1,nState
              Egamma = H(K,K)

              ! original XDW-CASPT2, xi = Dab^2
              if (DWType == 1) then
                xi_ag = (Ealpha - Egamma)**2
              ! new XDW-CASPT2, xi = (Haa/Hab)^2
              else if (DWType == 2) then
                xi_ag = (Ealpha/H(J,K))**2
              ! new XDW-CASPT2, xi = Dab/sqrt(Hab)
              else if (DWType == 3) then
                ! add a small positive constant to numerator to avoid 0/0
                Dag = abs(Ealpha - Egamma) + 1.0e-9_wp
                Hag = abs(H(J,K))
                ! below 1.0e-9 it is just noise, so set Hag = 0.0
                if (Hag .le. 1.0e-9_wp) then
                  Hag = 0.0_wp
                end if
                xi_ag = Dag/sqrt(Hag)
              else
                ! this should never happen since DWType is checked during input processing
                call Quit_OnUserError
              end if

              Wtot = Wtot + exp(-zeta*xi_ag)
              ! write(u6,*)'exp(-zeta*xi_ag) = ',exp(-zeta*xi_ag)
            end do

            ! original XDW-CASPT2, xi = Dab^2
            if (DWType == 1) then
              xi_ab = (Ealpha - Ebeta)**2
            ! new XDW-CASPT2, xi = (Haa/Hab)^2
            else if (DWType == 2) then
              xi_ab = (Ealpha/H(J,I))**2
            ! new XDW-CASPT2, xi = Dab/sqrt(Hab)
            else if (DWType == 3) then
              ! add a small positive constant to numerator to avoid 0/0
              Dab = abs(Ealpha - Ebeta) + 1.0e-9_wp
              Hab = abs(H(J,I))
              ! below 1.0e-9 it is just noise, so set Hab = 0.0
              if (Hab .le. 1.0e-9_wp) then
                Hab = 0.0_wp
              end if
              xi_ab = Dab/sqrt(Hab)
            else
              ! this should never happen since DWType is checked during input processing
              call Quit_OnUserError
            end if

            IJ = (I-1) + nState*(J-1)
            WORK(LDWGT+IJ) = exp(-zeta*xi_ab)/Wtot

          end do

        ! If it is an XMS-CASPT2 calculation, all the weights are equal,
        ! i.e. they all are 1/nState
        else if (IFXMS.and.(.not.IFDW)) then
          call dcopy_(nState**2,[1.0_wp/nState],0,WORK(LDWGT),1)

        ! If it is a normal MS-CASPT2, RMS-CASPT2 or a (X)DW-CASPT2 with zeta->infinity
        ! the weight vectors are the standard unit vectors e_1, e_2, ...
        else
          WORK(LDWGT + (nState*(I-1)) + (I-1)) = 1.0d0
        end if

      ! End of loop over states
      end do

      ! In case it is a XDW calculation, print out the weights
      if (IFDW.and.(IPRGLB.ge.VERBOSE)) then
        if (IFEFOCK) then
          write(u6,*)' Weights calculated with <I|H0|I>:'
        else
          write(u6,*)' Weights calculated with <I|H|I>:'
        end if
        call prettyprint(WORK(LDWGT),nState,nState)
      end if


      return
      end

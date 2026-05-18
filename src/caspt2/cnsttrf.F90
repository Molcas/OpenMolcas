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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*
      Subroutine CnstTrf(nTrf,Trf0,Trf)

      use caspt2_global, only: TraFro
      use caspt2_module, only: IfChol, NSYM, NFRO, NISH, NRAS1, NRAS2,  &
     &                         NRAS3, NASH, NSSH, NDEL, NBAS
      use definitions, only: wp, iwp
      use Constants, only: One

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: nTrf
      real(kind=wp), intent(in) :: Trf0(nTrf)
      real(kind=wp), intent(_OUT_) :: Trf(nTrf)

      integer(kind=iwp) :: iSQ, iTOrb, ipTrfL, iSym, nBasI, nFroI,      &
     &  nIshI, nAshI, nSshI, nDelI, NR1, NR2, NR3, nCor, nVir, I, J,    &
     &  iIsh, jIsh, IJ, iAsh, jAsh, iSsh, jSsh

      iSQ = 0
      iTOrb = 1 ! LTOrb
      ipTrfL = 0
      Do iSym = 1, nSym
        nBasI = nBas(iSym)
        nFroI = nFro(iSym)
        nIshI = nIsh(iSym)
        nAshI = nAsh(iSym)
        nSshI = nSsh(iSym)
        nDelI = nDel(iSym)
        NR1   = nRAS1(iSym)
        NR2   = nRAS2(iSym)
        NR3   = nRAS3(iSym)
        nCor  = nFroI + nIshI
        nVir  = nSshI + nDelI
        ipTrfL = ipTrfL + iSQ
        !! frozen + inactive
!       Do iIsh = 1, nFroI + nIshI
!         Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = One
!       End Do
        !! frozen
        If (IfChol) Then
          Do I = 1, nFroI
            Do J = 1, nFroI
              Trf(ipTrfL+I+nBasI*(J-1))                                 &
     &          = TraFro(I+nFroI*(J-1))
            End Do
          End Do
        else
          Do iIsh = 1, nFroI
            Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = One
          End Do
        End If
        !! inactive
        Do I = 1, nIshI
          iIsh = nFroI + I
          Do J = 1, nIshI
            jIsh = nFroI + J
            IJ=I-1+nIshI*(J-1)
            Trf(ipTrfL+iIsh+nBasI*(jIsh-1))                             &
     &        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nIshI*nIshI
        !! RAS1 space
        Do I = 1, NR1
          iAsh = nCor + I
          Do J = 1, NR1
            jAsh = nCor + J
            IJ=I-1+NR1*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))                             &
     &        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR1*NR1
        !! RAS2 space
        Do I = 1, NR2
          iAsh = nCor + NR1 + I
          Do J = 1, NR2
            jAsh = nCor + NR1 + J
            IJ=I-1+NR2*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))                             &
     &        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR2*NR2
        !! RAS3 space
        Do I = 1, NR3
          iAsh = nCor + NR1 + NR2 + I
          Do J = 1, NR3
            jAsh = nCor + NR1 + NR2 + J
            IJ=I-1+NR3*(J-1)
            Trf(ipTrfL+iAsh+nBasI*(jAsh-1))                             &
     &        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + NR3*NR3                                         &
!       call sqprt(trf,12)
      ! !! Active
      ! Do iAsh0 = 1, nAshI
      !   iAsh = nCor + iAsh0
      !   Do jAsh0 = 1, nAshI
      !     jAsh = nCor + jAsh0
!     !     Work(ipTrfL+iAsh-1+nBasI*(jAsh-1))
!    *!       = Work(iTOrb+nIshI*nIshI+iAsh0-1+nAshI*(jAsh0-1))
      !     Trf(ipTrfL+iAsh+nBasI*(jAsh-1))
     &!       = Trf0(iTOrb+iAsh0-1+nAshI*(jAsh0-1))
      !   End Do
      ! End Do
!       call sqprt(trf,12)
        !! virtual + deleted (deleted is not needed, though)
!       Do iSsh = nOcc+1, nOcc+nVir
!         Trf(ipTrfL+iSsh+nBasI*(iSsh-1)) = One
!       End Do
        Do I = 1, nVir
          iSsh = nCor + nAshI + I
          Do J = 1, nVir
            jSsh = nCor + nAshI + J
            IJ=I-1+nVir*(J-1)
            Trf(ipTrfL+iSsh+nBasI*(jSsh-1))                             &
     &        = Trf0(iTOrb+IJ)
          End Do
        End Do
        iTOrb = iTOrb + nSshI*nSshI
!       call sqprt(trf,12)
        iSQ = iSQ + nBasI*nBasI

!       n123 = nAshI*nAshI !! just for CAS at present
!       iTOrb = iTOrb + n123 + nSshI*nSshI
!     write(u6,*) 'transformation matrix'
!     call sqprt(trf,nbasi)
      End Do

      Return

      End Subroutine CnstTrf

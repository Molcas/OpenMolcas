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
! Copyright (C) 2010, Jonas Bostrom                                    *
!***********************************************************************
      Real*8 Function Compute_B(irc,kSOk,lSOl,jAOj,nBasFnc,iOpt)
!*****************************************************************
!                                                                *
!     Author Jonas Bostrom, June 2010                            *
!                                                                *
!     Purpose: To do part of MP2 gradient.                       *
!                                                                *
!*****************************************************************
      use ExTerm, only: BMP2
      Implicit Real*8 (a-h,o-z)
      Integer irc,kSOk,lSOl,jAOj,nBasFnc,iOpt
#include "exterm.fh"


      B_mp2 = 0.0d0
      iOff1 = (jAOj)*nBasFnc*nBasFnc + (kSOk-1)*nBasFnc + lSOl
      iOff2 = (jAOj)*nBasFnc*nBasFnc + (lSOl-1)*nBasFnc + kSOk
      B_mp2 = B_mp2 + (Bmp2(iOff1,iOpt)+Bmp2(iOff2,iOpt))/2.0d0
      Compute_B = B_mp2
      irc=0

      End

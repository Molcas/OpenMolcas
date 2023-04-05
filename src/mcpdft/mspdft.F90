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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

module mspdft
  implicit none
  private

#include "mspdft.fh"

  character(len=8) :: mspdftmethod
  logical :: do_rotate = .False.
  integer :: iF1MS,iF2MS,iFxyMS,iFocMS,iIntS,iDIDA,IP2MOt
  integer :: D1AOMS,D1SAOMS

  public :: dogradmspd, mspdftmethod, do_rotate, iF1MS, iF2MS
  public :: iFxyMS, iFocMS, iIntS, iDIDA, IP2MOt, D1AOMS, D1SAOMS

end module mspdft

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
! Copyright (C) 1992, Roland Lindh                                     *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine IniSew_scf(DSCF,EThr,SIntTh,KSDFT)
!***********************************************************************
! input:      EThr,DThr,FThr,DltNTh: div. Threshold values for         *
!               of SCF WF. (only relevant for direct SCF, DSCF=TRUE)   *
! output:     SIntTh: computed cutoff for integrals (prescreening)     *
!                                                                      *
! Note :  the corresponding finalization subroutine is ClsSew          *
!***********************************************************************

use Sizes_of_Seward, only: S
use Gateway_Info, only: ThrInt, CutInt
use OFembed, only: Do_OFemb
use RICD_Info, only: Do_DCCD
use InfSCF, only: nDisc, nCore
use AddCorr, only: Do_Addc, Do_Tw
use Constants, only: One

implicit none
logical DSCF
character(len=*) KSDFT
real*8 EThr, SIntTh
integer nDiff
logical RF_On, Langevin_On, PCM_On, EFP_On
external EFP_On
#include "print.fh"

if (DSCF .or. RF_On() .or. Langevin_On() .or. (KSDFT /= 'SCF') .or. Do_Addc .or. Do_Tw .or. Do_OFemb .or. EFP_On()) then
  nDiff = 0
  if (Langevin_On() .and. (S%iAngMx == 0)) nDiff = 1
  call IniSew(DSCF .or. Langevin_On() .or. PCM_On(),nDiff)
end if

if (Do_DCCD) then
  nCore = 0
  nDisc = 0
end if
if (DSCF) then
  CutInt = EThr*min(1.0D-7,One/dble(S%nDim)**2)
  ThrInt = Cutint
  SIntTh = CutInt
end if

return

end subroutine IniSew_scf

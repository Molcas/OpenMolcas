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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine LovMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,isFNO)
! Purpose: put info in MP2 common blocks.

use ChoMP2, only: C_os, ChkDecoMP2, ChoAlg, Decom_Def, DecoMP2, DoFNO, EOSMP2, ForceBatch, l_Dii, MxQual_Def, MxQualMP2, OED_Thr, &
                  set_cd_thr, SOS_mp2, Span_Def, SpanMP2, ThrMP2, Verbose
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mSym, lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
logical(kind=iwp), intent(in) :: isFNO
integer(kind=iwp) :: iSym
#include "corbinf.fh"

nSym = mSym

do iSym=1,nSym
  nOrb(iSym) = lnOrb(iSym)
  nOcc(iSym) = lnOcc(iSym)
  nFro(iSym) = lnFro(iSym)
  nDel(iSym) = lnDel(iSym)
  nExt(iSym) = lnVir(iSym)
end do

ChoAlg = 2
DecoMP2 = Decom_Def
ThrMP2 = -huge(ThrMP2)
SpanMP2 = Span_Def
MxQualMP2 = MxQual_Def
ChkDecoMP2 = .false.
ForceBatch = .false.
Verbose = .false.
SOS_mp2 = .false.
set_cd_thr = .true.
OED_Thr = 1.0e-8_wp
C_os = 1.3_wp
EOSMP2 = Zero

DoFNO = isFNO
l_Dii = nOcc(1)
do iSym=2,nSym
  l_Dii = l_Dii+nOcc(iSym)
end do

return

end subroutine LovMP2_putInf

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
subroutine LovCASPT2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,isFNO)
! Purpose: put info in MP2 common blocks.

use ChoMP2, only: C_os, ChkDecoMP2, ChoAlg, Decom_Def, DecoMP2, DoFNO, EOSMP2, ForceBatch, l_Dii, MxQual_Def, MxQualMP2, OED_Thr, &
                  set_cd_thr, shf, SOS_mp2, Span_Def, SpanMP2, ThrMP2, Verbose
use cOrbInf, only: nDel, nExt, nFro, nOcc, nOrb, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mSym, lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
logical(kind=iwp), intent(in) :: isFNO

nSym = mSym

nOrb(1:nSym) = lnOrb(1:nSym)
nOcc(1:nSym) = lnOcc(1:nSym)
nFro(1:nSym) = lnFro(1:nSym)
nDel(1:nSym) = lnDel(1:nSym)
nExt(1:nSym) = lnVir(1:nSym)

ChoAlg = 2
DecoMP2 = Decom_Def
ThrMP2 = -9.9e9_wp
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
shf = Zero

DoFNO = isFNO
l_Dii = sum(nOcc(1:nSym))

end subroutine LovCASPT2_putInf

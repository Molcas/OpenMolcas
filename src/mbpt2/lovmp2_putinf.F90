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

subroutine LovMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,X,Y,isFNO)
! Purpose: put info in MP2 common blocks.

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mSym, lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
real(kind=wp), intent(in) :: X(*), Y(*)
logical(kind=iwp), intent(in) :: isFNO
integer(kind=iwp) :: iSym
integer(kind=iwp), external :: ip_of_Work
#include "corbinf.fh"
#include "chomp2_cfg.fh"

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
ip_Dab = ip_of_Work(X(1))
ip_Dii = ip_of_Work(Y(1))
l_Dab = nExt(1)
l_Dii = nOcc(1)
do iSym=2,nSym
  l_Dab = l_Dab+nExt(iSym)**2
  l_Dii = l_Dii+nOcc(iSym)
end do

return

end subroutine LovMP2_putInf

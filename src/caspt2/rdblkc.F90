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
! Copyright (C) 2008, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 2008  PER-AKE MALMQUIST                    *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine RDBLKC(ISYM,ICASE,IVEC,VEC,nVEC)

use definitions, only: iwp, wp
use caspt2_global, only: LUSOLV, IDSCT
use EQSOLV, only: MxSct, ModVec
use caspt2_module, only: NASUP, NISUP, MxCase
#ifdef _DEBUGPRINT_
use caspt2_module, only: CASES
use definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: iSym, iCase, iVec, nVec
real(kind=wp), intent(out) :: VEC(nVec)
integer(kind=iwp) NAS, NIS, NCOEF, MDVEC, IDV, LVEC, IISTA, NCOL, NBLK
#ifdef _DEBUGPRINT_
integer(kind=iwp) I
#endif

! Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUGPRINT_
write(u6,*) ' RDBLKC (Normal repres.)'
write(u6,'(a,i2,a,a,a,i2)') ' Vector nr.',IVEC,'  Case ',CASES(ICASE),' Symm ',ISYM
#endif
NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)
NCOEF = NAS*NIS
if (NCOEF == 0) return
MDVEC = MODVEC(ISYM,ICASE)
!IDV = IDSCT(1,ISYM,ICASE,IVEC)
IDV = IDSCT(1+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
LVEC = 1
do IISTA=1,NIS,MDVEC
  NCOL = min(NIS+1-IISTA,MDVEC)
  NBLK = NAS*NCOL
  call DDAFILE(LUSOLV,2,VEC(LVEC),NBLK,IDV)
  LVEC = LVEC+NBLK
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' First few elements:'
write(u6,'(1x,5f15.6)') (VEC(I),I=1,min(NCOEF,10))
#endif

end subroutine RDBLKC

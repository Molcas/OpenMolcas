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

subroutine RDSCTC(ISCT,ISYM,ICASE,IVEC,VSCT,nVSCT)

use EQSOLV, only: ModVec, MxSCT
use caspt2_global, only: IDSCT, LUSOLV
use caspt2_module, only: MxCASE, NASUP, NISUP
#ifdef _DEBUGPRINT_
use caspt2_module, only: cases
use Definitions, only: u6
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISCT, iSYM, ICASE, IVEC, nVSCT
real(kind=wp), intent(out) :: VSCT(nVSCT)
integer(kind=iwp) :: iDS, MDVEC, NAS, NCOEF, NCOL, NIS, NSCT

! Read coefficient vector from LUSOLV (C repres).
#ifdef _DEBUGPRINT_
write(u6,*) ' RDSCTC (Normal repres.)'
write(u6,'(a,i2,a,a,a,i2,a,i2)') ' Vector nr.',IVEC,'  Case ',CASES(ICASE),' Symm ',ISYM,' Section ',ISCT
#endif
NAS = NASUP(ISYM,ICASE)
NIS = NISUP(ISYM,ICASE)
NCOEF = NAS*NIS
if (NCOEF == 0) return
MDVEC = MODVEC(ISYM,ICASE)
!IDS = IDSCT(ISCT,ISYM,ICASE,IVEC)
IDS = IDSCT(ISCT+MXSCT*(ISYM-1+8*(ICASE-1+MXCASE*(IVEC-1))))
NCOL = min(NIS-MDVEC*(ISCT-1),MDVEC)
NSCT = NAS*NCOL
call DDAFILE(LUSOLV,2,VSCT,NSCT,IDS)
#ifdef _DEBUGPRINT_
write(u6,*) ' First few elements:'
write(u6,'(1x,5f15.6)') VSCT(1:min(10,NSCT))
#endif

end subroutine RDSCTC

!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine getIphInfo(JOBIPH,nconf,LROOTS,IADR15)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: JOBIPH, IADR15(15)
integer(kind=iwp), intent(out) :: nconf, LROOTS
#include "rasdim.fh"
integer(kind=iwp) :: IAD15, NACTEL, ISPIN, NSYM, LSYM, NROOTS, NHOLE1, NELEC3, IPT2
integer(kind=iwp) :: nFro(MxSym), nISh(MxSym), nASh(MxSym), nDel(MxSym), nBas(MxSym), iRoot(MxRoot), nRS1(MxSym), nRS2(MxSym), &
                nRS3(MxSym)
character(len=LenIn8) :: Name(MxOrb), Title(18,MxTit)
character(len=2) :: Header(72)
real(kind=wp) :: POTNUC, Weight(MxRoot)

IAD15=IADR15(1)
call WR_RASSCF_Info(JOBIPH,2,IAD15,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,NAME,LENIN8*mxOrb,NCONF,HEADER,2*72, &
                    TITLE,4*18*mxTit,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)

return

end subroutine getIphInfo

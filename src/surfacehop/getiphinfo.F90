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

use Molcas, only: LenIn, MxOrb, MxRoot, MxSym
use RASDim, only: MxTit
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: JOBIPH, IADR15(15)
integer(kind=iwp), intent(out) :: nconf, LROOTS
integer(kind=iwp) :: IAD15, IPT2, ISPIN, LSYM, NACTEL, nASh(MxSym), nBas(MxSym), nDel(MxSym), NELEC3, nFro(MxSym), NHOLE1, &
                     nISh(MxSym), NROOTS, nRS1(MxSym), nRS2(MxSym), nRS3(MxSym), NSYM
real(kind=wp) :: POTNUC
character(len=LenIn+8) :: Title(18,MxTit)
character(len=2) :: Header(72)
integer(kind=iwp), allocatable :: iRoot(:)
real(kind=wp), allocatable :: Weight(:)
character(len=LenIn+8), allocatable :: AtName(:)

IAD15 = IADR15(1)
call mma_allocate(iRoot,MxRoot,Label='iRoot')
call mma_allocate(Weight,MxRoot,Label='Weight')
call mma_allocate(AtName,MxOrb,Label='AtName')
call WR_RASSCF_Info(JOBIPH,2,IAD15,NACTEL,ISPIN,NSYM,LSYM,NFRO,NISH,NASH,NDEL,NBAS,MxSym,AtName,(LenIn+8)*mxOrb,NCONF,HEADER,2*72, &
                    TITLE,4*18*mxTit,POTNUC,LROOTS,NROOTS,IROOT,MxRoot,NRS1,NRS2,NRS3,NHOLE1,NELEC3,IPT2,WEIGHT)
call mma_deallocate(iRoot)
call mma_deallocate(Weight)
call mma_deallocate(AtName)

end subroutine getIphInfo

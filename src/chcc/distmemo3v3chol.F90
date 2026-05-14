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

subroutine DistMemo3v3chol(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosM1,PosM2,PosM3,PosM4,PosM5,PosK,PosQ,PosT)
! This routine does:
! define initial positions of T,L,M and W arrays,
! used in routine o3v3chol
!
! I/O parameter description:
! maxdim   - maximal dimension of V'
! Posx     - initial positinos of arrays (O-all)
! PosT     - initial and last position (I/O)
!
! requirements of o3v3chol step
! K  - max{v'v'oo}
! Q  - max{v'v'oo}
! V1 - max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
! V2 - max{v'ooo, mv'v', v'v'oo, vo}
! V3 - max{mv'o, v'v'oo, v'ooo}
! V4 - max{v'ooo}
! H1 - max{v'o}
! H2 - max{v'o}
! H3 - max{v'o}
! H4 - max{v'o}
! M1 - max{moo}
! M2 - max{mv'o}
! M3 - max{mv'o}
! M4 - max{moo}
! M5 - max{mv'o}

use Index_Functions, only: nTri_Elem
use chcc_global, only: nc, no, nv, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: maxdim
integer(kind=iwp), intent(out) :: PosV1, PosV2, PosV3, PosV4, PosH1, PosH2, PosH3, PosH4, PosM1, PosM2, PosM3, PosM4, PosM5, PosK, &
                                  PosQ
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: length

!1 Q,K (used also as X,Y)

length = no*no*maxdim*maxdim

PosQ = posT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Q  ',PosQ,length
PosK = posT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM K  ',PosK,length

!2.1 V1 file max{v'oo2, mv'v', v'v'oo, mv'o, moo, vv}
length = maxdim*no*nTri_Elem(no)
if (maxdim*maxdim*nc > length) length = maxdim*maxdim*nc
if (no*no*maxdim*maxdim > length) length = no*no*maxdim*maxdim
if (no*nc*maxdim > length) length = no*nc*maxdim
if (no*nc*no > length) length = no*no*nc
if (nv*nv > length) length = nv*nv

PosV1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V1 ',PosV1,length

!2.2 V2 file - max{v'ooo, mv'v', v'v'oo, vo}

length = no*no*maxdim*maxdim
if (no*no*no*maxdim > length) length = no*no*no*maxdim
if (nc*maxdim*maxdim > length) length = nc*maxdim*maxdim
if (no*nv > length) length = no*nv

PosV2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V2 ',PosV2,length

!2.2 V3 files - max{mv'o, v'v'oo, v'ooo}

length = no*no*maxdim*maxdim
if (no*no*no*maxdim > length) length = no*no*no*maxdim
if (no*nc*maxdim > length) length = nc*no*maxdim

PosV3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V3 ',PosV3,length

!2.4 V4 files - max{v'ooo}

length = no*no*no*maxdim

PosV4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V4 ',PosV4,length

!3 H1-4 files - max{v'o}

length = no*maxdim

PosH1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H1 ',PosH1,length
PosH2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H2 ',PosH2,length
PosH3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H3 ',PosH3,length
PosH4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H4 ',PosH4,length

!4 M1-5 files
!  M1 - max{moo}
!  M2 - max{mv'o}
!  M3 - max{mv'o}
!  M4 - max{moo}
!  M5 - max{mv'o}

length = no*no*nc
PosM1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M1 ',PosM1,length
length = no*nc*maxdim
PosM2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M2 ',PosM2,length
length = no*nc*maxdim
PosM3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M3 ',PosM3,length
length = no*nc*no
PosM4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M4 ',PosM4,length
length = no*nc*maxdim
PosM5 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM M5 ',PosM5,length

if (printkey >= 10) write(u6,*) 'PosT ',PosT

return

end subroutine DistMemo3v3chol

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

subroutine DistMemo3v3t2(maxdim,PosV1,PosV2,PosV3,PosV4,PosH1,PosH2,PosH3,PosH4,PosK,PosQ,PosT)
! This routine does:
! define initial positions of H,V, QK
!
! I/O parameter description:
! maxdim   - maximal dimension of V'
! Posx     - initial positions of arrays (O-all)
! PosT     - initial and last position (I/O)
!
! requirements for o3v3t2:
! H1 - max {v'o}
! H2 - max {v'o}
! H3 - max {v'v',ooo}
! H4 - max {v'o}
! V1 - max {v'ov'o, vv'}
! V2 - max {v'ov'o}
! V3 - max {v'ov'o, o2oo}
! V4 - max {v'ov'o}
! PX - max {v'ov'o}
! QY - max {v'ov'o}

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: maxdim
integer(kind=iwp), intent(out) :: PosV1, PosV2, PosV3, PosV4, PosH1, PosH2, PosH3, PosH4, PosK, PosQ
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: length

!1 Q,K (used also as X,Y)

length = no*no*maxdim*maxdim

PosQ = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM Q  ',PosQ,length
PosK = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM K  ',PosK,length

!2.1 V1 file - max {v'ov'o, vv'}

length = no*no*maxdim*maxdim
if (nv*maxdim > length) length = maxdim*nv

PosV1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V1 ',PosV1,length

!2.2 V2-4 files
!    V2 - max {v'ov'o}
!    V3 - max {v'ov'o, o2oo}
!    V4 - max {v'ov'o}

length = no*no*maxdim*maxdim
if (no*no*no*maxdim > length) length = no*no*no*maxdim

PosV2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V2 ',PosV2,length

length = no*no*maxdim*maxdim
if (nTri_Elem(no)*no*no > length) length = no*no*nTri_Elem(no)

PosV3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V3 ',PosV3,length

length = no*no*maxdim*maxdim
PosV4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM V4 ',PosV4,length

!3 H1,2 files

length = no*maxdim

PosH1 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H1 ',PosH1,length
PosH2 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H2 ',PosH2,length

!3.2 H3 file

length = maxdim*maxdim
if (no*maxdim > length) length = no*maxdim
if (no*no*no > length) then
  length = no*no*no
end if

PosH3 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H3 ',PosH3,length

!3.3 H4 file

length = no*maxdim

PosH4 = PosT
PosT = PosT+length
if (printkey >= 10) write(u6,*) 'DM H4 ',PosH4,length

if (printkey >= 10) write(u6,*) 'PosT ',PosT

return

end subroutine DistMemo3v3t2

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

subroutine DistMemo2v4(NaGrp,NbeGrp,NaSGrp,NbeSgrp,mdGrpa,mdGrpbe,mdSGrpa,mdSGrpbe,PosTau,PosT2n1,PosT2n2,PosT2w,PosL11,PosL12, &
                       PosL21,PosL22,PosL23,PosL24,PosL2W,PosH1,PosH2,PosM1,PosM2,PosW1,PosW2,PosWw,PosWx,PosT,NL2)
! This routine does:
! define initial positions of T,L,M and W arrays,
! described in o2v4ctl routine
!
! I/O parameter description:
! NaGrp    - # of groups in a set (I)
! NbeGrp   - # of groups in be set (I)
! NaSGrp   - # of subgroups in each (a)' group (I)
! NbeSGrp  - # of subgroups in each (be)' group (I)
! mdGrpa   - # maximal dimension od (a)' group (I)
! mdGrpbe  - # maximal dimension od (be)' group (I)
! mdSGrpa  - # maximal dimension od (a)" subgroup (I)
! mdSGrpbe - # maximal dimension od (be)" subgroup (I)
! PosX     - initial positions of arrays (O-all)
! PosT     - next free position (O)
! NL2      - # of L2 vectors (O)

use Index_Functions, only: nTri_Elem
use chcc_global, only: intkey, nc, no, nv, printkey
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NaGrp, NbeGrp, NaSGrp, NbeSgrp, mdGrpa, mdGrpbe, mdSGrpa, mdSGrpbe
integer(kind=iwp), intent(out) :: PosTau, PosT2n1, PosT2n2, PosT2w, PosL11, PosL12, PosL21, PosL22, PosL23, PosL24, PosL2W, PosH1, &
                                  PosH2, PosM1, PosM2, PosW1, PosW2, PosWw, PosWx, NL2
integer(kind=iwp), intent(inout) :: PosT
integer(kind=iwp) :: lenab, lenbega, length

!1 Tau

PosTau = PosT
if (NaGrp == 1) then
  length = no*no*nTri_Elem(nv)
else
  length = no*no*mdGrpa*mdGrpa
end if
PosT = PosT+length
if (printkey >= 10) write(u6,99) 'DM Tau',PosTau

!2 T2n files

posT2n1 = PosT
if ((NbeGrp == 1) .and. (NbeSGrp == 1)) then
  length = nTri_Elem(no)*nTri_Elem(nv)
else
  length = nTri_Elem(no)*mdSGrpbe*mdSGrpbe
end if
PosT = PosT+length

PosT2n2 = PosT
if ((NbeGrp == 1) .and. (NbeSGrp == 1)) then
  length = nTri_Elem(no-1)*nTri_Elem(nv-1)
else
  length = nTri_Elem(no-1)*mdSGrpbe*mdSGrpbe
end if
PosT = PosT+length

PosT2w = PosT
if ((NbeGrp == 1) .and. (NbeSGrp == 1)) then
  length = nTri_Elem(no)*nTri_Elem(nv)
else
  length = nTri_Elem(no)*mdSGrpbe*mdSGrpbe
end if
PosT = PosT+length
if (printkey >= 10) write(u6,99) 'DM T2 ',PosT2n1,PosT2n2,PosT2w

! L1 files

if (intkey == 1) then
  ! integral based
  length = 0
else
  ! cholesky based
  length = nc*mdGrpa*no
end if

if (NaGrp == 1) then
  PosL11 = PosT
  PosL12 = PosT
  PosT = PosT+length
else
  PosL11 = PosT
  PosT = PosT+length
  PosL12 = PosT
  PosT = PosT+length
end if
if (printkey >= 10) write(u6,99) 'DM L1 ',PosL11,PosL12

! L2 files

if (intkey == 1) then
  ! integral based
  length = 0
else
  ! cholesky based
  length = nc*mdGrpa*mdGrpbe
end if

if (NaGrp == 1) then
  if (NbeGrp == 1) then
    ! All L2 are identical
    PosL21 = PosT
    PosL22 = PosT
    PosL23 = PosT
    PosL24 = PosT
    PosT = PosT+length
    NL2 = 1
  else
    ! L21=L23 and L22=L24
    PosL21 = PosT
    PosL23 = PosT
    PosT = PosT+length
    PosL22 = PosT
    PosL24 = PosT
    PosT = PosT+length
    NL2 = 2
  end if
else
  if (NbeGrp == 1) then
    ! L21=L22 and L23=L24
    PosL21 = PosT
    PosL22 = PosT
    PosT = PosT+length
    PosL23 = PosT
    PosL24 = PosT
    PosT = PosT+length
    NL2 = 2
  else
    ! all L2 files are different
    PosL21 = PosT
    PosT = PosT+length
    PosL22 = PosT
    PosT = PosT+length
    PosL23 = PosT
    PosT = PosT+length
    PosL24 = PosT
    PosT = PosT+length
    NL2 = 4
  end if
end if

PosL2W = PosT
if (intkey == 1) then
  ! integral based
  length = 0
else
  ! cholesky based
  if ((NaGrp == 1) .and. (NbeGrp == 1)) then
    length = nc*nTri_Elem(nv)
  else
    length = nc*mdGrpa*mdGrpbe
  end if
  if ((no*mdGrpbe) > length) then
    length = no*mdGrpbe
  end if
  if (nc*no*mdGrpa > length) length = nc*no*mdGrpa
end if

PosT = PosT+length

if (printkey >= 10) write(u6,99) 'DM L2 ',PosL21,PosL22,PosL23,PosL24,PosL2W

! H files

length = no*mdSGrpbe
if (intkey == 1) then
  PosH1 = PosT
  PosT = PosT+length
  PosH2 = PosT
  PosT = PosT+length
else
  PosH1 = PosT
  PosH2 = PosT
end if
if (printkey >= 10) write(u6,99) 'DM H  ',PosH1,PosH2

! M files

length = nc*mdSGrpa*mdSGrpbe

if ((NaGrp == 1) .and. (NaSGrp == 1) .and. (NbeGrp == 1) .and. (NbeSGrp == 1)) then
  ! M1=M2
  PosM1 = PosT
  PosM2 = PosT
  PosT = PosT+length
else
  ! M files are different
  PosM1 = PosT
  PosT = PosT+length
  PosM2 = PosT
  PosT = PosT+length
end if
if (printkey >= 10) write(u6,99) 'DM M  ',PosM1,PosM2

! W1,W2 files

if (NaGrp*NaSGrp == 1) then
  ! W1 and W2 are identical
  length = mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
  if ((no > mdSGrpbe) .and. (intkey == 1)) length = mdSGrpa*mdSGrpbe*mdSGrpa*no
  PosW1 = PosT
  PosW2 = PosT
  PosT = PosT+length
else
  ! W1 and W2 are different
  length = mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
  if ((no > mdSGrpbe) .and. (intkey == 1)) length = mdSGrpa*mdSGrpbe*mdSGrpa*no
  PosW1 = PosT
  PosT = PosT+length
  PosW2 = PosT
  PosT = PosT+length
end if

! Ww file

PosWw = PosT

if ((NaGrp == 1) .and. (NaSGrp == 1)) then
  lenab = nTri_Elem(nv)
else
  lenab = mdSGrpa*mdSGrpa
end if
if ((NbeGrp == 1) .and. (NbeSGrp == 1)) then
  lenbega = nTri_Elem(nv)
else
  lenbega = mdSGrpbe*mdSGrpbe
end if

length = lenab*lenbega
if (intkey == 1) then
  length = mdSGrpa*mdSGrpbe*mdSGrpa*mdSGrpbe
  if (no > mdSGrpbe) length = mdSGrpa*mdSGrpbe*mdSGrpa*no
end if
PosT = PosT+length

! Wx file

PosWx = PosT
length = mdSGrpa*mdSGrpbe*mdSGrpa*no
if (intkey == 0) length = 0
PosT = PosT+length

if (printkey >= 10) then
  write(u6,99) 'DM W  ',PosW1,PosW2,PosWw,PosWx

  write(u6,99) 'PosT ',PosT
end if

return

99 format(a7,10(i10,1x))

end subroutine DistMemo2v4

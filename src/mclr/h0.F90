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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine H0(rdia,MP1,MP2,MQ,isym,nprciv,TimeDep)
! frontend to Jeppe's explicit routines

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: FIMO, H0F, H0S, IDC, Int2, iRefSM, NACOB, NAELCI, NBELCI, NCPCNT, NDPCNT, nGP, NOCOB, NTYP, PSSIGN, SBIDT, &
                     XISPSM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: rdia(*)
integer(kind=iwp), intent(in) :: MP1, MP2, MQ, iSym
integer(kind=iwp), intent(out) :: nprciv
logical(kind=iwp), intent(in) :: TimeDep
integer(kind=iwp) :: i, ieaw, iRC, iSpc, lH0SCR, LH0T, lVec2, MxCSFC, MxDTFC, MXP, nActEl, nDet, nSBDet
real(kind=wp) :: ENA
integer(kind=iwp), allocatable :: SBCNF(:)
real(kind=wp), allocatable :: H0Scr(:), H0T(:), Vec2(:)
real(kind=wp), external :: E2, E2_TD

ispc = 1
NDET = nint(XISPSM(ISYM,ISPC))

NSBDET = MP1+MP2+MQ
MXP = MP1+MP2
LH0T = nTri_Elem(MXP)+MP1*MQ
MXCSFC = max(0,maxval(NCPCNT(1:NTYP)))
MXDTFC = max(0,maxval(NDPCNT(1:NTYP)))

nactel = naelci(1)+nbelci(1)
if (TimeDep) then
  EnA = E2_td(FIMO,Int2,0,-1)
else
  EnA = E2(FIMO,Int2,0,-1)
end if
LH0SCR = max(6*NSBDET,4*NSBDET+4*NOCOB,nTri_Elem(MP1)+MP1**2)
LVEC2 = MXCSFC**2+2*MXDTFC+2*MXDTFC**2+MXDTFC*NACTEL
LVEC2 = max(lvec2,ndet)

if (isym == irefsm) then
  ieaw = 1
else
  ieaw = 2
end if
call mma_allocate(SBIDT,NSBDET,Label='SBIDT')
call mma_allocate(H0S,MXP**2,Label='H0S')
call mma_allocate(H0F,MXP,Label='H0F')
call mma_allocate(H0T,LH0T,Label='H0T')
call mma_allocate(SBCNF,NSBDET,Label='SBCNF')
call mma_allocate(H0SCR,LH0SCR,Label='H0Scr')
call mma_allocate(VEC2,lvec2,Label='Vec2')

call H0MAT_MCLR(H0T,SBIDT,SBCNF,MP1,MP2,MQ,NACOB,NPRCIV,ISYM,IDC,PSSIGN,rDIA,Vec2,H0Scr,ieaw)

do i=1,nprciv
  H0T(nTri_Elem(i)) = H0T(nTri_Elem(i))-ENA
end do
if (NGP) call mkp1(nprciv,SBIDT,H0T,rdia)
!call Triprt('PRECI',' ',H0T,nprciv)
!write(u6,*) (SBIDT(i),i=1,nprciv)
call mma_deallocate(Vec2)
call mma_deallocate(H0Scr)
call mma_deallocate(SBCNF)
call square(H0T,H0S,1,NPRCIV,NPRCIV)
call mma_deallocate(H0T)

irc = 0
call dgetrf_(NPRCIV,NPRCIV,H0S,NPRCIV,H0F,irc)
if (irc /= 0) then
  write(u6,*) 'Sorry but you have an singular ci matrix'
  write(u6,*) 'Set ExpDimension and restart mclr'
  call Abend()
end if

end subroutine H0

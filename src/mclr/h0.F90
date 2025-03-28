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
! frontend to jeppes explicit routines

use Exp, only: H0S, H0F, SBIDT
use negpre, only: nGP
use iso_c_binding
use Arrays, only: Int2, FIMO
use MCLR_Data, only: iRefSM, IDC, PSSIGN
use MCLR_Data, only: NAELCI, NBELCI, XISPSM
use MCLR_Data, only: MXP1, MXP2, MXQ
use MCLR_Data, only: NACOB, NOCOB
use MCLR_Data, only: NTYP, NCPCNT, NDPCNT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: u6

implicit none
real*8 rdia(*)
integer MP1, MP2, MQ, iSym, nprciv
logical TimeDep
integer iSpc, nDet, nSBDet, MXP, LH0T, MxCSFC, MxDTFC, iTyp, nActEl, lH0SCR, ieaw, i, iRC, lVec2
real*8, external :: E2, E2_TD
real*8 ENA
real*8, allocatable :: H0T(:), Vec2(:)
real*8, allocatable, target :: H0Scr(:)
integer, pointer :: iH0Scr(:)
integer, allocatable :: SBCNF(:)

MXP1 = MP1
MXP2 = MP2
MXQ = MQ
ispc = 1
NDET = nint(XISPSM(ISYM,ISPC))

NSBDET = MXP1+MXP2+MXQ
MXP = MXP1+MXP2
LH0T = MXP*(MXP+1)/2+MXP1*MXQ
MXCSFC = 0
MXDTFC = 0
do ITYP=1,NTYP
  MXCSFC = max(MXCSFC,NCPCNT(ITYP))
  MXDTFC = max(MXDTFC,NDPCNT(ITYP))
end do

nactel = naelci(1)+nbelci(1)
if (TimeDep) then
  EnA = E2_td(FIMO,Int2,0,-1)
else
  EnA = E2(FIMO,Int2,0,-1)
end if
LH0SCR = max(6*NSBDET,4*NSBDET+4*NOCOB,MXP1*(MXP1+1)/2+MXP1**2)
LVEC2 = 2*NACTEL+MXCSFC**2+6*MXDTFC+2*MXDTFC**2+max(MXDTFC*NACTEL+2*NACTEL,4*NACOB+2*NACTEL)
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
call c_f_pointer(c_loc(H0SCR),iH0Scr,[LH0SCR])
call mma_allocate(VEC2,lvec2,Label='Vec2')

call H0MAT_MCLR(H0T,SBIDT,SBCNF,MXP1,MXP2,MXQ,NACOB,NPRCIV,ISYM,IDC,PSSIGN,rDIA,Vec2,H0Scr,iH0Scr,ieaw)

do i=1,nprciv
  H0T(i*(i+1)/2) = H0T(i*(i+1)/2)-ENA
end do
if (NGP) call mkp1(nprciv,SBIDT,H0T,rdia)
!call Triprt('PRECI',' ',H0T,nprciv)
!write(u6,*) (SBIDT(i),i=1,nprciv)
call mma_deallocate(Vec2)
nullify(iH0Scr)
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

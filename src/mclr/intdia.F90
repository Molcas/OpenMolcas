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

subroutine INTDIA(DIAG,NSPC,ISPC,ISM,LSPC,IAMCMP,ecore)
! CI diagonal in SD basis for the NCSPC ci spaces defined by
! ISPC,ISM
!
! if IAMCMP /= 0 : then it is assumed that a complex
! hermitian eigenvalued problem is being solved by
! doubling the dimensions. The diagonal is then
! constructed and written out twice

use Str_Info, only: STR, NELEC, NOCTYP
use stdalloc, only: mma_allocate, mma_deallocate
use MCLR_Data, only: IPRDIA
use MCLR_Data, only: iDC, PLSIGN, PSSIGN
use MCLR_Data, only: IASTFI, IBSTFI, ISMOST, MNR1IC, MXR3IC
use MCLR_Data, only: ICISTR
use MCLR_Data, only: NTOOB, NACOB
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use input_mclr, only: nIrrep

implicit none
real*8 DIAG(*)
integer NSPC
! ==============
! Specific Input
! ==============
integer ISPC(NSPC), LSPC(NSPC), ISM(NSPC)
integer IAMCMP
real*8 ECORE
! ==============
! General Input
! ==============
integer idum(1)
real*8, allocatable :: JA(:), KA(:), XA(:), XB(:), SCR(:), H1D(:)
integer, allocatable :: BLTP(:), IOIO(:)
integer LUDIA, MXOCOC, IISPC, NOCTPA, NOCTPB, NLOOP, ILOOP, IATP, IBTP, NAEL, NBEL, MNRS1C, MXRS3C, LLUDIA

! OBS THIS WILL JUST WORK FOR CASSCF/RASSCF RESPONSE
LUDIA = 0

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2(1:8),ntoob,0)
  call dmrg_dim_change_mclr(RGras2(1:8),nacob,0)
end if

! Local memory

call mma_allocate(JA,NTOOB**2,Label='JA')
call mma_allocate(KA,NTOOB**2,Label='KA')
call mma_allocate(XA,NACOB,Label='XA')
call mma_allocate(XB,NACOB,Label='XB')
call mma_allocate(SCR,2*NACOB,Label='SCR')
if (doDMRG) then !yma
  ! wired Call mma_allocate(H1D,NACOB_KLH1D,Label='H1D')
  call mma_allocate(H1D,NACOB,Label='H1D')
else
  call mma_allocate(H1D,NACOB,Label='H1D')
end if
call mma_allocate(BLTP,nIrrep,Label='BLTP')

! Largest NOCTPA*NOCTPB block
MXOCOC = 0
do IISPC=1,NSPC
  NOCTPA = NOCTYP(IASTFI(ISPC(IISPC)))
  NOCTPB = NOCTYP(IBSTFI(ISPC(IISPC)))
  MXOCOC = max(MXOCOC,NOCTPA*NOCTPB)
end do
call mma_allocate(IOIO,MXOCOC,Label='IOIO')
! Diagonal of one-body integrals and coulomb and exchange integrals

call GT1DIA_MCLR(H1D)
call GTJK_MCLR(JA,KA)

! K goes to J - K
KA(:) = JA(:)-KA(:)

! Loop over internal CI spaces

if (IAMCMP == 0) then
  NLOOP = 1
else
  NLOOP = 2
end if
do ILOOP=1,NLOOP
  do IISPC=1,NSPC
    IATP = IASTFI(ISPC(IISPC))
    IBTP = IBSTFI(ISPC(IISPC))
    NAEL = NELEC(IATP)
    NBEL = NELEC(IBTP)
    NOCTPA = NOCTYP(IATP)
    NOCTPB = NOCTYP(IBTP)
    MNRS1C = MNR1IC(ISPC(IISPC))
    MXRS3C = MXR3IC(ISPC(IISPC))

    call ZBLTP(ISMOST(1,ISM(IISPC)),nIrrep,IDC,BLTP,idum)
    call IAIBCM_MCLR(MNRS1C,MXRS3C,NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,IOIO,IPRDIA)

    if (ICISTR <= 1) then
      LLUDIA = 0
    else
      LLUDIA = LUDIA
    end if
    call CIDIA4(NAEL,Str(IATP)%OCSTR,NBEL,Str(IBTP)%OCSTR,NACOB,DIAG,nIrrep,H1D,ISMOST(1,ISM(IISPC)),BLTP,XA,XB,SCR,JA,KA, &
                Str(IATP)%NSTSO,Str(IBTP)%NSTSO,IOIO,NOCTPA,NOCTPB,Str(IATP)%ISTSO,Str(IBTP)%ISTSO,LLUDIA,ECORE,PLSIGN,PSSIGN, &
                IPRDIA,NTOOB,ICISTR)

    if ((ICISTR <= 1) .and. (LUDIA > 0)) then
      ! Each CI space is written in one record
      call ITODS(LSPC(IISPC),1,0,LUDIA)
      call TODSC_MCLR(DIAG,LSPC(IISPC),0,LUDIA)
    end if
  end do
  ! Write end of vector mark
  if (LUDIA > 0) then
    IDUM(1) = -1
    call ITODS(IDUM,1,0,LUDIA)
  end if
end do

call mma_deallocate(IOIO)
call mma_deallocate(BLTP)
call mma_deallocate(H1D)
call mma_deallocate(SCR)
call mma_deallocate(XB)
call mma_deallocate(XA)
call mma_deallocate(KA)
call mma_deallocate(JA)

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(LRras2(1:8),ntoob,0)
  call dmrg_dim_change_mclr(LRras2(1:8),nacob,0)
end if

end subroutine INTDIA

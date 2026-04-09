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

subroutine INTDIA(DIAG,NSPC,ISPC,ISM,IAMCMP,ecore)
! CI diagonal in SD basis for the NCSPC ci spaces defined by
! ISPC,ISM
!
! if IAMCMP /= 0 : then it is assumed that a complex
! hermitian eigenvalued problem is being solved by
! doubling the dimensions. The diagonal is then
! constructed and written out twice

use Str_Info, only: NELEC, NOCTYP, STR
use MCLR_Data, only: IASTFI, IBSTFI, ICISTR, iDC, MNR1IC, MXR3IC, NACOB, NTOOB, PSSIGN
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use input_mclr, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: DIAG(*)
integer(kind=iwp), intent(in) :: NSPC, ISPC(NSPC), ISM(NSPC), IAMCMP
real(kind=wp), intent(in) :: ECORE
integer(kind=iwp) :: IATP, IBTP, idum(1), IISPC, ILOOP, LLUDIA, MNRS1C, MXOCOC, MXRS3C, NAEL, NBEL, NLOOP, NOCTPA, NOCTPB
integer(kind=iwp), allocatable :: BLTP(:), IOIO(:)
real(kind=wp), allocatable :: H1D(:), JA(:), KA(:), SCR(:), XA(:), XB(:)

! OBS THIS WILL JUST WORK FOR CASSCF/RASSCF RESPONSE

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2,ntoob,0)
  call dmrg_dim_change_mclr(RGras2,nacob,0)
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

    call ZBLTP(ISM(IISPC),nIrrep,IDC,BLTP,idum)
    call IAIBCM_MCLR(MNRS1C,MXRS3C,NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,IOIO)

    LLUDIA = 0
    call CIDIA4(NAEL,Str(IATP)%OCSTR,NBEL,Str(IBTP)%OCSTR,NACOB,DIAG,nIrrep,H1D,ISM(IISPC),BLTP,XB,JA,KA,Str(IATP)%NSTSO, &
                Str(IBTP)%NSTSO,IOIO,NOCTPA,NOCTPB,Str(IATP)%ISTSO,Str(IBTP)%ISTSO,LLUDIA,ECORE,PSSIGN,NTOOB,ICISTR)

  end do
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
  call dmrg_dim_change_mclr(LRras2,ntoob,0)
  call dmrg_dim_change_mclr(LRras2,nacob,0)
end if

end subroutine INTDIA

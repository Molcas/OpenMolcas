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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine VVVO_Drv(nSym,nBas,nFro,Keep,iSym,iSymI,iSymJ,iSymK,iSymL,lT2AO,T2AO,vLag,nOcc,nBasT,nBMX,CMO,DPT2AO,DPT2CAO,FPT2AO, &
                    FPT2CAO,DIA,DI,FIFA,FIMO)

use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8), nFro(8), Keep(8), iSym, iSymI, iSymJ, iSymK, iSymL, lT2AO, nOcc, nBasT, nBMX
real(kind=wp), intent(in) :: T2AO(lT2AO), CMO(nBasT**2), DPT2AO(nBasT**2), DPT2CAO(nBasT**2), DIA(nBasT**2), DI(nBasT**2)
real(kind=wp), intent(inout) :: vLag(nBasT**2), FPT2AO(nBasT**2), FPT2CAO(nBasT**2), FIFA(nBasT**2), FIMO(nBasT**2)
integer(kind=iwp) :: LBUF
logical(kind=iwp) :: DoCholesky
real(kind=wp), allocatable :: W1(:), W2(:), WRK(:)

! nAux is the number of occupied orbitals
DoCholesky = .false.
call DecideOnCholesky(DoCholesky)
!if (DoCholesky .and. (ALGO == 2)) then
!  !! I assume ALGO=2 does not exist (it is not documented and under
!  !! debugging, according to rasscf/cho_rasscf_rdinp)
!  call abend()
!end if

call mma_allocate(W2,NBMX*NBMX,Label='W2')
call mma_allocate(WRK,nBasT*nBasT,Label='WRK')

call mma_MaxDBLE(LBUF)
if (DoCholesky) LBUF = NBMX*NBMX+1

! Standard building of the Fock matrix from Two-el integrals

call mma_allocate(W1,LBUF,Label='W1')

if (LBUF < 1+NBMX**2) then
  write(u6,*) ' FockTwo_Drv Error: Too little memory remains for the call to FOCKTWO.'
  write(u6,*) ' Largest allocatable array size LBUF=',LBUF
  write(u6,*) ' Max nr of bf in any symmetry,  NBMX=',NBMX
  write(u6,*) ' Required minimum size     1+NBMX**2=',1+NBMX**2
  write(u6,*) '    (All in Real*8-size words)'
  call ABEND()
end if

if (DoCholesky) then
  call VVVOX2(Keep,iSym,iSymI,iSymJ,iSymK,iSymL,nBasT,vLag,CMO,WRK,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,FIFA,FIMO)
else
  call VVVOX(nSym,nBas,nFro,Keep,iSymI,iSymJ,iSymK,iSymL,NBMX,T2AO,vLag,CMO,nOcc,nBasT,LBUF,W1,W2,WRK,DPT2AO,DPT2CAO,FPT2AO, &
             FPT2CAO,DIA,DI,FIFA,FIMO)
end if
!! vLag must be transposed
!! In VVVOX(2) subroutines, vLag(p,mu) is constructed.
WRK(1:nBasT**2) = vLag(1:nBasT**2)
call trnsps(nBasT,nBasT,WRK,vLag)
call mma_deallocate(WRK)
call mma_deallocate(W1)
call mma_deallocate(W2)

end subroutine VVVO_Drv

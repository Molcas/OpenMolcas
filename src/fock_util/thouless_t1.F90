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

subroutine Thouless_T1(CMO,nSym,nBas,nFro,nOcc,nSsh,T1amp)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nFro(nSym), nOcc(nSym), nSsh(nSym)
real(kind=wp), intent(_OUT_) :: T1amp(*)
integer(kind=iwp) :: i, iDummy(1), iErr, ifr, iOff, iSym, ito, iU, j, jp_C, jp_S, jp_T, jp_X, jU, k, kk, kOff, l_O, l_O2, l_S, &
                     lScr, Lu, nOrb
real(kind=wp) :: Dummy(1), omega
character(len=40) :: OrbTit
real(kind=wp), allocatable :: R(:), S(:), Scr(:), U(:), W(:), X(:), Y(:), Z(:)

! Compute the T1 amplitudes according to Thouless formula
! -------------------------------------------------------

l_S = nBas(1)**2
lScr = nBas(1)*(nFro(1)+nOcc(1))
l_O = nOcc(1)
do iSym=2,nSym
  l_S = l_S+nBas(iSym)**2
  lScr = max(lScr,nBas(iSym)*(nFro(iSym)+nOcc(iSym)))
  l_O = max(l_O,nOcc(iSym))
end do
l_O2 = l_O**2

call mma_allocate(Scr,lScr,Label='Scr')
call mma_allocate(U,lScr,Label='U')

call mma_allocate(W,l_O2,Label='W')
call mma_allocate(Y,l_O2,Label='Y')
call mma_allocate(Z,l_O2,Label='Z')
call mma_allocate(R,l_O2,Label='R')

call mma_allocate(S,l_S,Label='S')
call mma_allocate(X,l_S,Label='X')

call GetOvlp_Localisation(S,'Sqr',nBas,nSym)

Lu = 12
call RdVec('INPORB',Lu,'C',nSym,nBas,nBas,X,Dummy,Dummy,iDummy,OrbTit,1,iErr)

write(u6,*)
write(u6,*) '      Thouless singles amplitudes from: '
write(u6,*) '      '//OrbTit
write(u6,*)

iOff = 0
kOff = 0
do iSym=1,nSym
  jp_S = 1+iOff
  jp_C = 1+iOff+nBas(iSym)*nFro(iSym)
  jp_X = 1+iOff+nBas(iSym)*nFro(iSym)
  nOrb = nOcc(iSym)+nSsh(iSym)

  call GetUmat_T1(U,CMO(jp_C),S(jp_S),X(jp_X),Scr,lScr,nBas(iSym),nOrb,nOcc(iSym))

  iU = 1
  do j=1,nOcc(iSym)
    ifr = 1+nOrb*(j-1)
    ito = 1+nOcc(iSym)*(j-1)
    call dcopy_(nOcc(iSym),U(ifr),1,Scr(ito),1)
    jU = ifr+nOcc(iSym)
    do i=1,nSsh(iSym)
      U(iU) = U(jU)
      iU = iU+1
      jU = jU+1
    end do
  end do

  ! SVD of U in the oo space:   U = Y * w * Z'

  call SVD(nOcc(iSym),nOcc(iSym),nOcc(iSym),Scr,W,.true.,Y,.true.,Z,ierr,R)

  if (ierr /= 0) then
    write(u6,*)
    write(u6,*) ' *** Warning: SVD failed to get singval: ',ierr
    write(u6,*) ' *** Located in Thouless_T1 -- call to SVD .'
    write(u6,*)
    write(u6,*) ' omega= ',(W(k),k=1,nOcc(iSym))
  end if

  call FZero(R,nOcc(iSym)**2)
  do k=1,nOcc(iSym)
    omega = W(k)
    kk = nOcc(iSym)*(k-1)+k
    if (omega > 1.0e-8_wp) then
      R(kk) = One/omega
    end if
  end do

  ! Compute U^-1 = Z * w^-1 * Y'

  call DGEMM_('N','T',nOcc(iSym),nOcc(iSym),nOcc(iSym),One,R,nOcc(iSym),Y,nOcc(iSym),Zero,W,nOcc(iSym))

  call DGEMM_('N','N',nOcc(iSym),nOcc(iSym),nOcc(iSym),One,Z,nOcc(iSym),W,nOcc(iSym),Zero,Scr,nOcc(iSym))

  jp_T = 1+kOff
  call DGEMM_('T','T',nOcc(iSym),nSsh(iSym),nOcc(iSym),One,Scr,nOcc(iSym),U,nSsh(iSym),Zero,T1amp(jp_T),nOcc(iSym))

  iOff = iOff+nBas(iSym)**2
  kOff = kOff+nOcc(iSym)*nSsh(iSym)
end do

call mma_deallocate(X)
call mma_deallocate(S)
call mma_deallocate(R)
call mma_deallocate(Z)
call mma_deallocate(Y)
call mma_deallocate(W)
call mma_deallocate(U)
call mma_deallocate(Scr)

return

end subroutine Thouless_T1

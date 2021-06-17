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
!***********************************************************************
!                                                                      *
! This routine computes the expectation value <s^2> for a UHF          *
! wavefunction.                                                        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
! Ca     -  Alpha orbitals, input.                                     *
! Cb     -  Beta orbitals, input.                                      *
! S      -  Overlap matrix, input.                                     *
! nAlpha -  Number of alpha orbitals occupied per irrep, input.        *
! nBeta  -  Number of beta orbitals occupied per irrep, input.         *
! nBas   -  Number of basis functions per irrep, input.                *
! nOrb   -  Number of orbitals per irrep, input.                       *
! nSym   -  Number of irreps, input.                                   *
! s2     -  <s^2>, output.                                             *
!                                                                      *
!***********************************************************************

subroutine s2calc(Ca,Cb,S,nAlpha,nBeta,nBas,nOrb,nSym,s2)

implicit none
#include "stdalloc.fh"
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 Ca(*)
real*8 Cb(*)
real*8 S(*)
integer nAlpha(*)
integer nBeta(*)
integer nBas(*)
integer nOrb(*)
integer nSym
real*8 s2
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
real*8 sz
real*8 sb
real*8 so
integer iSym
integer npSmat
integer npHalf
integer npTfrm
integer idxCMO
integer idxOvl
integer i
real*8, dimension(:), allocatable :: Smat, Half, Tfrm

!----------------------------------------------------------------------*
! Debug printing stuff                                                 *
!----------------------------------------------------------------------*
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(6,'(a,8i5)') 'nSym  . . . . . .',nSym
write(6,'(a,8i5)') 'nBas  . . . . . .',(nBas(iSym),iSym=1,nSym)
write(6,'(a,8i5)') 'nOrb  . . . . . .',(nOrb(iSym),iSym=1,nSym)
write(6,'(a,8i5)') 'nAlpha  . . . . .',(nAlpha(iSym),iSym=1,nSym)
write(6,'(a,8i5)') 'nBeta . . . . . .',(nBeta(iSym),iSym=1,nSym)
i = 1
do iSym=1,nSym
  call RecPrt('Ca',' ',Ca(i),nBas(iSym),nBas(iSym))
  call RecPrt('Cb',' ',Cb(i),nBas(iSym),nBas(iSym))
  call RecPrt('S ',' ',S(i),nBas(iSym),nBas(iSym))
  i = i+nBas(iSym)**2
end do
#endif
!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
s2 = 0.0d0
!----------------------------------------------------------------------*
! Compute <sz> and N_beta                                              *
!----------------------------------------------------------------------*
sz = 0.0d0
sb = 0.0d0
do iSym=1,nSym
  sz = sz+0.5d0*(nAlpha(iSym)-nBeta(iSym))
  sb = sb+1.0d0*nBeta(iSym)
end do
s2 = s2+sz*(sz+1.0d0)+sb
#ifdef _DEBUGPRINT_
write(6,'(a,f12.6)') 'sz  . . . . . . .',sz
write(6,'(a,f12.6)') 'sb  . . . . . . .',sb
write(6,'(a,f12.6)') 's2  . . . . . . .',s2
#endif
!----------------------------------------------------------------------*
! Compute size of scratch matrices, the allocate.                      *
!----------------------------------------------------------------------*
npSmat = 0
npHalf = 0
npTfrm = 0
do iSym=1,nSym
  npSmat = max(npSmat,nBas(iSym)*nBas(iSym))
  npHalf = max(npHalf,nAlpha(iSym)*nBas(iSym))
  npTfrm = max(npTfrm,nAlpha(iSym)*nBeta(iSym))
end do
#ifdef _DEBUGPRINT_
write(6,'(a,i5)') 'npSmat  . . . . .',npSmat
write(6,'(a,i5)') 'npHalf  . . . . .',npHalf
write(6,'(a,i5)') 'npTfrm  . . . . .',npTfrm
#endif
if (npTfrm == 0) return
call mma_allocate(Smat,npSmat)
call mma_allocate(Half,npHalf)
call mma_allocate(Tfrm,npTfrm)
!----------------------------------------------------------------------*
! Compute sum_a sum_b [ C_alpha* S C_beta ]^2                          *
!----------------------------------------------------------------------*
!call xxdGemul(a,lda,'N/T',b,ldb,'N/T',c,ldc,l,m,n)
! C=A*B, m is contraction index
so = 0.0d0
idxCMO = 1
idxOvl = 1
do iSym=1,nSym
  !write(6,'(a,i1)') 'Irrep ',iSym
  if (nAlpha(iSym)*nBeta(iSym) > 0) then
    call Square(S(idxOvl),Smat,1,nBas(iSym),nBas(iSym))
    call DGEMM_('T','N',nAlpha(iSym),nBas(iSym),nBas(iSym),1.0d0,Ca(idxCMO),nBas(iSym),Smat,nBas(iSym),0.0d0,Half,nAlpha(iSYm))
    !call RecPrt('Half transform','(12f12.6)',Half,nAlpha(iSym),nBas(iSym))
    call DGEMM_('N','N',nAlpha(iSym),nBeta(iSym),nBas(iSym),1.0d0,Half,nAlpha(iSym),Cb(idxCMO),nBas(iSym),0.0d0,Tfrm,nAlpha(iSym))
    !call RecPrt('Transform','(12f12.6)',Tfrm,nAlpha(iSym),nBeta(iSym))
    do i=1,nAlpha(iSym)*nBeta(iSym)
      so = so+Tfrm(i)*Tfrm(i)
    end do
  end if
  idxCMO = idxCMO+nBas(iSym)*nOrb(iSym)
  idxOvl = idxOvl+nBas(iSym)*(nBas(iSym)+1)/2
end do
s2 = s2-so
!write(6,'(a,f12.6)') 'so  . . . . . . .',so
!write(6,'(a,f12.6)') 's2  . . . . . . .',s2
!----------------------------------------------------------------------*
! Deallocate scratch matrix.                                           *
!----------------------------------------------------------------------*
call mma_deallocate(Tfrm)
call mma_deallocate(Half)
call mma_deallocate(Smat)
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine s2calc

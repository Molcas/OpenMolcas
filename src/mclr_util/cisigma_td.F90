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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

#include "compiler_features.h"
#ifdef _WARNING_WORKAROUND_
#  ifdef _LIFETIME_BUG_
#    define _SAVE_TARGET_ , save
#  else
#    define _SAVE_TARGET_
#  endif
#else
#  define _SAVE_TARGET_
#endif

subroutine CISigma_td(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,NT,Have_2_el)
! sigma = H C through the CI interface (CI_Sigma_Op), see CISigma. For
! time-dependent perturbations the vectors have two halves of length
! nConf1: NT = 'N': only the first half; NT = 'S': symmetric operator,
! both halves; otherwise the operator is not symmetric (square
! integrals), and the second half is the transposed operator applied to
! the first half. Int2a is not used.

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use ipPage, only: ipin, ipnout, opout, W
use MCLR_Data, only: ipMat, nConf1, nDens
use CI_Operations, only: CI_Sigma_Op
use Orbital_Spaces, only: nSym
use input_mclr, only: nBas, nCSF, ntAsh, Page, TimeDep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: iState = 1
integer(kind=iwp), intent(in) :: iiSpin, iCSym, iSSym, nInt1, nInt2s, nInt2a, ipCI1, ipCI2
real(kind=wp), target, intent(in) :: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
character, intent(in) :: NT
logical(kind=iwp), intent(in) :: Have_2_el
integer(kind=iwp) :: i, ij, ijkl, iS, j, ji, jilk, jS, k, kl, l, lk, nC, nS
logical(kind=iwp) :: Square
real(kind=wp), allocatable, target :: TI1(:), TI2(:)

#include "macros.fh"
unused_var(Int2a)

! yma: notice the nconf1 if DMRG
if (nconf1 == 0) return

nC = nCSF(iCSym)
nS = nCSF(iSSym)
if (TIMEDEP) then
  Square = (NT == 'T')  ! The operator is not symmetric
  if (page) then
    write(u6,*) 'Page not implemented for Time-dependent perturbations'
    call Abend()
  end if
  call ipin(ipCI1)
  call ipin(ipci2)
  call CI_Sigma_Op(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,Square,nC,W(ipCI1)%A(1:nC),nS,W(ipci2)%A(1:nS))

  if (NT == 'N') return

  if (NT == 'S') then

    ! Symmetric operator, no transpose of integrals needed!
    call CI_Sigma_Op(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,.false.,nC,W(ipCI1)%A(nConf1+1:nConf1+nC), &
                       nS,W(ipci2)%A(nConf1+1:nConf1+nS))

  else

    ! The operator is not symmetric --> transpose the integrals
    call mma_allocate(TI1,nDens,Label='TI1')
    call mma_allocate(TI2,ntash**4,Label='TI2')
    do i=1,ntash
      do j=1,ntash
        ij = i+ntash*(j-1)
        ji = j+ntash*(i-1)
        do k=1,ntash
          do l=1,ntash
            kl = k+ntash*(l-1)
            lk = l+ntash*(k-1)
            if (ij >= kl) then
              ijkl = iTri(ij,kl)
              jilk = iTri(ji,lk)
              TI2(jilk) = int2s(ijkl)
            end if
          end do
        end do
      end do
    end do
    do is=1,nSym
      js = Mul(Mul(icsym,issym),is)
      if (nbas(js)*nbas(is) /= 0) call DGETMO(Int1(ipmat(is,js)),nbas(is),nbas(is),nbas(js),TI1(ipmat(js,is)),nbas(js))
    end do
    call ipin(ipCI1)
    call ipin(ipci2)
    call CI_Sigma_Op(iState,iiSpin,iCSym,iSSym,size(TI1),TI1,size(TI2),TI2,Have_2_el,Square,nC,W(ipCI1)%A(1:nC),nS, &
                       W(ipci2)%A(nConf1+1:nConf1+nS))
    call mma_deallocate(TI1)
    call mma_deallocate(TI2)

  end if

else   ! If not timedep

  if (page) call ipnout(ipci2)
  call ipin(ipCI1)
  call ipin(ipci2)
  call CI_Sigma_Op(iState,iiSpin,iCSym,iSSym,nInt1,Int1,nInt2s,Int2s,Have_2_el,.false.,nC,W(ipCI1)%A(1:nC),nS,W(ipci2)%A(1:nS))
  if (page) call opout(ipci1)

end if

end subroutine CISigma_td

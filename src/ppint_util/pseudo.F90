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

subroutine Pseudo(ai,xi,yi,zi,lit,aj,xj,yj,zj,ljt,gout,intmax,lmn1u,ccr,zcr,nkcrl,nkcru,lcr,ncr,xc,yc,zc,npot)
!***********************************************************************
! This collection of subroutines, adapted from the ARGOS suite of programs,
! calculates integrals of an atomic semi-local one-component (i.e., non-
! relativistic or scalar-relativistic) pseudopotential (PP) between two shells
! of primitive cartesian GTOs. The integrals are contained in field gout(*).
! lmn1u is the maximum l value + 1 of the basis functions
! lproju is the maximum l value + 1 of the pseudopotential
! GTO parameters
!
! The GTOs are given by Ni (x-xi)^li (y-yi)^mi (z-zi)^ni exp(-ai(r-ri)^2)
! and Nj (x-xj)^lj (y-yj)^mj (z-zj)^nj exp(-aj(r-rj)^2). The Ni, Nj are
! chosen such that the integrals are normalized to (2l-1)!!(2m-1)!!(2n-1)!!.
! Input data are
! ai, aj : exponents
! xi,yi,zi, xj,yj,zj : centers
! lit = li+mi+ni+1, ljt = lj+mj+nj+1 : angular-momentum values + 1
! PP parameters
!
! The pseudopotential is of the form V_loc + sum_l V_l P_l, where V_loc and the V_l
! are radial potentials and P_l is the projector onto angular-momentum l.
! Both V_loc and the V_l are given by expansions of the form
! V = sum_i ccr(i) ( r^(ncr(i)-2) exp(-zcr(i)*r^2) ),
! where r is the distance from the PP center (xc, yc, zc).
! The PP parameters are assumed to be taken from a list of PPs for several atoms
! containing the following input data:
! npot : number of radial potentials in the list
! ccr(i) : coefficients
! ncr(i), zcr(i) : exponential parameters
! nkcrl(1,k) : lower index i for V_loc at nucleus k
! nkcru(1,k) : upper index i for V_loc at nucleus k
! nkcrl(l+2,k) : lower index i for V_l at nucleus k
! nkcru(l+2,k) : upper index i for V_l at nucleus k
! lcr(k) : maximum l value + 1 of the PP at nucleus k
! k=kcrs : specifies the PP for which integrals are calculated
! xc,yc,zc : position of nucleus kcrs
!
!***********************************************************************

use ppint_arrays, only: binom, dfac, hpt, hwt, lmf, lml, lmnv, lmx, lmy, lmz, zlm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: OneHalf, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: lproju = 9, imax = 100, kcrs = 1
integer(kind=iwp), intent(in) :: lit, ljt, intmax, lmn1u, nkcrl(lproju+1,kcrs), nkcru(lproju+1,kcrs), lcr(kcrs), ncr(imax), npot
real(kind=wp), intent(in) :: ai, xi, yi, zi, aj, xj, yj, zj, ccr(npot), zcr(npot), xc, yc, zc
real(kind=wp), intent(inout) :: gout(2*intmax)
integer(kind=iwp) :: i, l1max2, lambu, lcru, litot, ljtot, lmax, lmnpwr, lmnvmx, lproju1, ltot1, mproju, ncru, ndfac
real(kind=wp) :: crda(lproju,3), crdb(lproju,3), eps
real(kind=wp), parameter :: facij = Pi**OneHalf

!write(u6,*) 'ncr',(ncr(i),i=1,npot)
!write(u6,*) 'zcr',(zcr(i),i=1,npot)
!write(u6,*) 'ccr',(ccr(i),i=1,npot)
!write(u6,*) 'nkcrl',(nkcrl(i,1),i=1,lcr(1)+1)
!write(u6,*) 'nkcru',(nkcru(i,1),i=1,lcr(1)+1)
!write(u6,*)
!                                                                      *
!***********************************************************************
!                                                                      *
lmnvmx = (lmn1u*(lmn1u+1)*(lmn1u+2))/6
ncru = 0
do i=1,npot
  ncru = max(ncru,ncr(i))
end do
ndfac = max(4*lmn1u+2*lproju-3,6*lproju+3,4*lmn1u-1,2*lmn1u+2*lproju+1,4,ncru+4*lmn1u+2*lproju-1)
lmax = max(1,lmn1u-1+max(lmn1u-1,lproju))
l1max2 = (lmax+1)**2
lmnpwr = (((lmax*(lmax+2)*(lmax+4))/3)*(lmax+3)+(lmax+2)**2*(lmax+4))/16
ltot1 = lit+ljt-1
mproju = 2*lproju+1
lambu = ljt+lproju

call mma_allocate(lmnv,3,lmnvmx,label='lmnv')
call mma_allocate(lmf,l1max2,label='lmf')
call mma_allocate(lml,l1max2,label='lml')
call mma_allocate(lmx,lmnpwr,label='lmx')
call mma_allocate(lmy,lmnpwr,label='lmy')
call mma_allocate(lmz,lmnpwr,label='lmz')
call mma_allocate(binom,lmn1u*(lmn1u+1)/2,label='binom')
call mma_allocate(dfac,ndfac,label='dfac')
call mma_allocate(zlm,lmnpwr,label='zlm')

call lmnvgn(lmn1u,lmnv)
eps = 1.0e-12_wp

call cortab(binom,dfac,eps,lmf,lml,lmx,lmy,lmz,lmax,lmn1u,ndfac,zlm)
!                                                                      *
!***********************************************************************
!                                                                      *
! integrals for V_loc part of PP

lproju1 = lproju+1
call pseud1(ccr,gout,ltot1,ncr,nkcrl,nkcru,zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lproju1,crda,crdb)
!                                                                      *
!***********************************************************************
!                                                                      *
! integrals for V_l P_l parts of PP

lcru = lcr(kcrs)
call pseud2(ccr,gout,lambu,ltot1,mproju,ncr,nkcrl,nkcru,zcr,lit,ljt,ai,aj,xi,yi,zi,xj,yj,zj,xc,yc,zc,kcrs,lcru,lproju1,crda,crdb)

call mma_deallocate(lmnv)
call mma_deallocate(lmf)
call mma_deallocate(lml)
call mma_deallocate(lmx)
call mma_deallocate(lmy)
call mma_deallocate(lmz)
call mma_deallocate(binom)
call mma_deallocate(dfac)
call mma_deallocate(zlm)
call mma_deallocate(hpt)
call mma_deallocate(hwt)

!                                                                      *
!***********************************************************************
!                                                                      *
! final results with normalization factors

!faci = sqrt((Four*ai)**lit/Two*sqrt(Two*ai))/((Four*ai)**(Half*(lit-1))*(Two*ai/Pi)**(Three/Four))
!facj = sqrt((Four*aj)**ljt/Two*sqrt(Two*aj))/((Four*aj)**(Half*(ljt-1))*(Two*aj/Pi)**(Three/Four))
!facij = faci*facj
litot = lit*(lit+1)/2
ljtot = ljt*(ljt+1)/2

call DScal_(litot*ljtot,facij,gout,1)

!do i=1,litot
!  do j=1,ljtot
!    gout(i+(j-1)*litot) = gout(i+(j-1)*litot)*facij
!  end do
!end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Reorder integrals from Argos to Molcas order

call Molcas_Order(GOut,litot,ljtot)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Pseudo

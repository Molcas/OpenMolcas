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
! Copyright (C) 2014, Steven Vancoillie                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!> Faroald, the determinant CI library: the only module the CI interface
!> (ci_lib_util) uses. The data of Faroald stay in its own modules
!> (Faroald_Core, Faroald_RAS, Faroald_CITrans, Faroald_SecQuant); outside
!> they are read-only (the dimensions and flags below). Everything else
!> comes in through the arguments: the setup (Faro_Init, Faro_CITrans_Init,
!> Faro_RAS_Init) and the integrals h(t,u), g(t,u,v,x) over the orbitals of
!> Faroald (chemists' notation), which every routine that needs them takes
!> as arguments; Faroald keeps no copy of them.
module Faroald

use Faroald_Core, only: Faro_Free, Faro_Full, Faro_Init, Faro_nOrb, Faro_ndeta, Faro_ndetb, Faro_nela, Faro_nelb, &
                        Faro_sigma_update, Faro_transition_one_pdm, Faro_transition_two_pdm, Faro_TRACI
use Faroald_RAS, only: Faro_n_ras, Faro_RAS_active, Faro_RAS_count_dets, Faro_RAS_Free, Faro_RAS_Init, Faro_RAS_ndet, &
                       Faro_RAS_sigma, Faro_RAS_traci, Faro_RAS_trans_one_pdm, Faro_RAS_trans_two_pdm
use Faroald_CITrans, only: Faro_conf_dets, Faro_conf_occ, Faro_csf2ras, Faro_csf2sd, Faro_ncsf, Faro_ras2csf, Faro_sd2csf, &
                           Faro_sort, Faro_sym_init
use Faroald_SecQuant, only: Faro_Binom
use Definitions, only: wp, iwp

implicit none
private

public :: Faro_Binom, Faro_CITrans_Free, Faro_CITrans_Init, Faro_conf_dets, Faro_conf_occ, Faro_csf2ras, Faro_csf2sd, &
          Faro_Diag_Cnf_Avg, Faro_Free, Faro_Full, Faro_HDet, Faro_Init, Faro_n_ras, Faro_ncsf, Faro_ndeta, Faro_ndetb, &
          Faro_nela, Faro_nelb, Faro_nOrb, Faro_ras2csf, Faro_RAS_active, Faro_RAS_count_dets, Faro_RAS_Free, Faro_RAS_Init, &
          Faro_RAS_ndet, Faro_RAS_sigma, Faro_RAS_traci, Faro_RAS_trans_one_pdm, Faro_RAS_trans_two_pdm, Faro_sd2csf, &
          Faro_sigma_update, Faro_sort, Faro_sym_init, Faro_transition_one_pdm, Faro_transition_two_pdm, Faro_TRACI

contains

!> The CSF transformation tables of the CSF space: nel electrons in norb
!> orbitals, spin multiplicity mult, the symmetry information nsym, osym
!> (the symmetries of the norb active orbitals) and lsym (the symmetry of
!> the CSF space); C1 is nsym = lsym = 1, osym = 1. RAS: nrsprt partitions,
!> oras the partition of each orbital, nrasel(k) the minimum number of
!> electrons in the partitions 1..k; the orbitals must be the levels of the
!> split graph (by partition, within each by symmetry). CAS is nrsprt = 1,
!> oras = 1, nrasel = nel.
subroutine Faro_CITrans_Init(nel,norb,mult,nsym,osym,lsym,nrsprt,oras,nrasel,nrasmx)

  use Faroald_CITrans, only: Faro_graph_init, Faro_ndet_group, Faro_ncsf_group, Faro_ndo_max, Faro_ndo_min, Faro_ndoc_group, &
                             Faro_nsoc_group, Faro_spintable_create, Faro_spintabs, Faro_spintabs_allocate
  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: nel, norb, mult
  integer(kind=iwp), intent(in) :: nsym, osym(norb), lsym, nrsprt, oras(norb), nrasel(nrsprt), nrasmx(nrsprt)
  integer(kind=iwp) :: ido, iso, nela, nelb

  ! compute alpha/beta subsets
  nela = (nel+(mult-1))/2
  nelb = (nel-(mult-1))/2

  ! determine the range of the configuration groups
  if (nel > norb) then
    Faro_ndo_min = nel-norb
  else
    Faro_ndo_min = 0
  end if
  Faro_ndo_max = nelb

  ! compute the various sizes per group
  call mma_allocate(Faro_ndoc_group,[Faro_ndo_min,Faro_ndo_max],label='ndoc_group')
  call mma_allocate(Faro_nsoc_group,[Faro_ndo_min,Faro_ndo_max],label='nsoc_group')
  call mma_allocate(Faro_ndet_group,[Faro_ndo_min,Faro_ndo_max],label='ndet_group')
  call mma_allocate(Faro_ncsf_group,[Faro_ndo_min,Faro_ndo_max],label='ncsf_group')

  call Faro_spintabs_allocate()
  ! loop over configurations
  do ido=Faro_ndo_min,Faro_ndo_max
    iso = nel-2*ido
    ! compute different block sizes
    Faro_ndoc_group(ido) = Faro_Binom(ido,norb)
    Faro_nsoc_group(ido) = Faro_Binom(iso,norb-ido)
    Faro_ndet_group(ido) = Faro_Binom(nela-ido,iso)
    Faro_ncsf_group(ido) = Faro_ndet_group(ido)-Faro_Binom(nela-ido+1,iso)
    ! compute+store spin table for this configuration
    Faro_spintabs(ido)%ndet = Faro_ndet_group(ido)
    Faro_spintabs(ido)%ncsf = Faro_ncsf_group(ido)
    call Faro_spintable_create(iso,nelb-ido,Faro_spintabs(ido))
  end do

  ! the file order of the CSFs (Faro_sort)
  call Faro_sym_init(nel,norb,nsym,osym,lsym,nrsprt,oras,nrasel,nrasmx)
  call Faro_graph_init(nel,norb,mult)

end subroutine Faro_CITrans_Init

!> Free the CSF transformation tables (Faro_CITrans_Init).
subroutine Faro_CITrans_Free()

  use Faroald_CITrans, only: Faro_graph_free, Faro_ndet_group, Faro_ncsf_group, Faro_ndoc_group, Faro_nsoc_group, &
                             Faro_spintabs_free, Faro_sym_free
  use stdalloc, only: mma_deallocate

  call mma_deallocate(Faro_ndoc_group,safe='*')
  call mma_deallocate(Faro_nsoc_group,safe='*')
  call mma_deallocate(Faro_ndet_group,safe='*')
  call mma_deallocate(Faro_ncsf_group,safe='*')
  call Faro_spintabs_free()
  call Faro_graph_free()
  call Faro_sym_free()

end subroutine Faro_CITrans_Free

!> The configuration average of the CSF diagonal of the CSF space set up
!> by Faro_CITrans_Init (in its configuration order): the average of the
!> determinant energies of each configuration, from its occupations, with
!> the integrals h, g over the norb orbitals of Faroald.
subroutine Faro_Diag_Cnf_Avg(norb,h,g,n,Diag)

  use stdalloc, only: mma_deallocate
  use Constants, only: Zero, One, Two

  integer(kind=iwp), intent(in) :: norb, n
  real(kind=wp), intent(in) :: h(norb,norb), g(norb,norb,norb,norb)
  real(kind=wp), intent(out) :: Diag(n)
  integer(kind=iwp) :: i, ic, j, na, nb, ncf, nd2, no
  integer(kind=iwp), allocatable :: dco(:), kf(:), ncs(:), sco(:)
  real(kind=wp) :: E, ni, nj, P

  call Faro_conf_occ(ncf,kf,ncs,dco,sco)
  do ic=1,ncf
    nd2 = popcnt(dco(ic))
    no = popcnt(sco(ic))
    na = Faro_nela-nd2
    nb = Faro_nelb-nd2
    P = Zero
    if (no > 1) P = real(na*(na-1)+nb*(nb-1),kind=wp)/real(no*(no-1),kind=wp)
    E = Zero
    do i=1,norb
      if (btest(dco(ic),i-1)) then
        ni = Two
        E = E+g(i,i,i,i)
      else if (btest(sco(ic),i-1)) then
        ni = One
      else
        cycle
      end if
      E = E+ni*h(i,i)
      do j=i+1,norb
        if (btest(dco(ic),j-1)) then
          nj = Two
        else if (btest(sco(ic),j-1)) then
          nj = One
        else
          cycle
        end if
        E = E+ni*nj*g(i,i,j,j)
        if ((ni == Two) .and. (nj == Two)) then
          E = E-Two*g(i,j,j,i)
        else if ((ni == Two) .or. (nj == Two)) then
          E = E-g(i,j,j,i)
        else
          E = E-P*g(i,j,j,i)
        end if
      end do
    end do
    Diag(kf(ic):kf(ic)+ncs(ic)-1) = E
  end do
  call mma_deallocate(kf)
  call mma_deallocate(ncs)
  call mma_deallocate(dco)
  call mma_deallocate(sco)

end subroutine Faro_Diag_Cnf_Avg

!> <aL bL|H|aR bR> between Faroald determinants (alpha and beta strings as
!> bit strings, orbital t in bit t-1, alpha string before the beta string)
!> by the Slater-Condon rules, with the integrals h, g over the norb
!> orbitals (chemists' notation), H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx
!> (E_tu E_vx - delta_uv E_tx), as in Faro_sigma_update.
function Faro_HDet(norb,h,g,aL,bL,aR,bR) result(E)

  use Constants, only: Zero, Half

  real(kind=wp) :: E
  integer(kind=iwp), intent(in) :: norb, aL, bL, aR, bR
  real(kind=wp), intent(in) :: h(norb,norb), g(norb,norb,norb,norb)
  integer(kind=iwp) :: dA, dB, nA, nB, i, j, p, q, r, s

  dA = ieor(aL,aR)
  dB = ieor(bL,bR)
  nA = popcnt(dA)/2
  nB = popcnt(dB)/2
  E = Zero
  if (nA+nB > 2) return

  if (nA+nB == 0) then
    do i=1,norb
      if (btest(aL,i-1)) then
        E = E+h(i,i)
        do j=1,norb
          if (btest(aL,j-1)) E = E+Half*(g(i,i,j,j)-g(i,j,j,i))
          if (btest(bL,j-1)) E = E+g(i,i,j,j)
        end do
      end if
      if (btest(bL,i-1)) then
        E = E+h(i,i)
        do j=1,norb
          if (btest(bL,j-1)) E = E+Half*(g(i,i,j,j)-g(i,j,j,i))
        end do
      end if
    end do
  else if ((nA == 1) .and. (nB == 0)) then
    E = HDet_Single(norb,h,g,aL,aR,bR)
  else if ((nA == 0) .and. (nB == 1)) then
    E = HDet_Single(norb,h,g,bL,bR,aR)
  else if (nA == 2) then
    E = HDet_Double(norb,g,aL,aR)
  else if (nB == 2) then
    E = HDet_Double(norb,g,bL,bR)
  else
    ! one alpha (p <- q) and one beta (r <- s) excitation
    p = Bit_Orb(norb,iand(dA,aL),1)
    q = Bit_Orb(norb,iand(dA,aR),1)
    r = Bit_Orb(norb,iand(dB,bL),1)
    s = Bit_Orb(norb,iand(dB,bR),1)
    E = Bit_Par(aR,p,q)*Bit_Par(bR,r,s)*g(p,q,r,s)
  end if

end function Faro_HDet

! Faro_HDet: single excitation p <- q in the string xL/xR of one spin, yR
! the string of the other spin
function HDet_Single(norb,h,g,xL,xR,yR) result(E)

  real(kind=wp) :: E
  integer(kind=iwp), intent(in) :: norb, xL, xR, yR
  real(kind=wp), intent(in) :: h(norb,norb), g(norb,norb,norb,norb)
  integer(kind=iwp) :: dX, iC, j, p, q

  dX = ieor(xL,xR)
  p = Bit_Orb(norb,iand(dX,xL),1)
  q = Bit_Orb(norb,iand(dX,xR),1)
  iC = iand(xL,xR)
  E = h(p,q)
  do j=1,norb
    if (btest(iC,j-1)) E = E+g(p,q,j,j)-g(p,j,j,q)
    if (btest(yR,j-1)) E = E+g(p,q,j,j)
  end do
  E = Bit_Par(xR,p,q)*E

end function HDet_Single

! Faro_HDet: double excitation p,r <- q,s in the string xL/xR of one spin
function HDet_Double(norb,g,xL,xR) result(E)

  real(kind=wp) :: E
  integer(kind=iwp), intent(in) :: norb, xL, xR
  real(kind=wp), intent(in) :: g(norb,norb,norb,norb)
  integer(kind=iwp) :: aI, dX, p, q, r, s

  dX = ieor(xL,xR)
  p = Bit_Orb(norb,iand(dX,xL),1)
  r = Bit_Orb(norb,iand(dX,xL),2)
  q = Bit_Orb(norb,iand(dX,xR),1)
  s = Bit_Orb(norb,iand(dX,xR),2)
  aI = ibset(ibclr(xR,q-1),p-1)
  E = Bit_Par(xR,p,q)*Bit_Par(aI,r,s)*(g(p,q,r,s)-g(p,s,r,q))

end function HDet_Double

! the n-th orbital (bit) set in x
function Bit_Orb(norb,x,n) result(k)

  integer(kind=iwp) :: k
  integer(kind=iwp), intent(in) :: norb, x, n
  integer(kind=iwp) :: m

  m = 0
  do k=1,norb
    if (btest(x,k-1)) then
      m = m+1
      if (m == n) return
    end if
  end do
  k = 0

end function Bit_Orb

! sign of a+_t a_u acting on the string x: the parity of the number of
! orbitals occupied in x strictly between t and u
function Bit_Par(x,t,u) result(iSgn)

  integer(kind=iwp) :: iSgn
  integer(kind=iwp), intent(in) :: x, t, u
  integer(kind=iwp) :: lo, hi

  lo = min(t,u)
  hi = max(t,u)
  iSgn = 1-2*poppar(ibits(x,lo,hi-lo-1))

end function Bit_Par

end module Faroald

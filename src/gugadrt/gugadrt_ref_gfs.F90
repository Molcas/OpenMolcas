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

subroutine gugadrt_ref_gfs(nel,ndj,locu,nm)

use gugadrt_global, only: lsm_inn, max_ref, norb_dz, norb_inn, nstart_act, spin
use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nel, nm
integer(kind=iwp), intent(out) :: ndj, locu(8,max_ref)
integer(kind=iwp) :: i, l1, l2, l3, l4, l5, l6, l7, l8, ldj, lh, lhe, lhs, lhsm(8), lm, lpsum, m, m1, m2, m3, m4, m5, m6, m7, m8, &
                     mdj, mys, ne_act, ne_s, nes, npair, nre
integer(kind=iwp), allocatable :: lscu(:,:)

call mma_allocate(lscu,[0,8],[1,max_ref],label='lscu')

ne_act = nel-2*norb_dz
ne_s = nint(spin*2)
lhs = nstart_act
lhe = norb_inn
lhsm(1:8) = 0
do lh=lhs,lhe
  lm = lsm_inn(lh)
  lhsm(lm) = lhsm(lm)+1
end do
mdj = 0
do nes=ne_s,ne_act,2
  do l1=0,lhsm(1)
    do l2=0,lhsm(2)
      do l3=0,lhsm(3)
        do l4=0,lhsm(4)
          do l5=0,lhsm(5)
            do l6=0,lhsm(6)
              do l7=0,lhsm(7)
                do l8=0,lhsm(8)
                  lpsum = l1+l2+l3+l4+l5+l6+l7+l8
                  if (lpsum /= nes) cycle
                  mys = 1
                  if (mod(l1,2) == 1) mys = Mul(mys,1)
                  if (mod(l2,2) == 1) mys = Mul(mys,2)
                  if (mod(l3,2) == 1) mys = Mul(mys,3)
                  if (mod(l4,2) == 1) mys = Mul(mys,4)
                  if (mod(l5,2) == 1) mys = Mul(mys,5)
                  if (mod(l6,2) == 1) mys = Mul(mys,6)
                  if (mod(l7,2) == 1) mys = Mul(mys,7)
                  if (mod(l8,2) == 1) mys = Mul(mys,8)
                  if (mys /= nm) cycle
                  mdj = mdj+1
                  lscu(0,mdj) = lpsum
                  lscu(1,mdj) = l1
                  lscu(2,mdj) = l2
                  lscu(3,mdj) = l3
                  lscu(4,mdj) = l4
                  lscu(5,mdj) = l5
                  lscu(6,mdj) = l6
                  lscu(7,mdj) = l7
                  lscu(8,mdj) = l8
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do
ndj = 0
do m=1,mdj
  npair = (ne_act-lscu(0,m))/2
  do l1=0,lhsm(1)-lscu(1,m)
    do l2=0,lhsm(2)-lscu(2,m)
      do l3=0,lhsm(3)-lscu(3,m)
        do l4=0,lhsm(4)-lscu(4,m)
          do l5=0,lhsm(5)-lscu(5,m)
            do l6=0,lhsm(6)-lscu(6,m)
              do l7=0,lhsm(7)-lscu(7,m)
                outer: do l8=0,lhsm(8)-lscu(8,m)
                  lpsum = l1+l2+l3+l4+l5+l6+l7+l8
                  if (lpsum == npair) then
                    m1 = l1*2+lscu(1,m)
                    m2 = l2*2+lscu(2,m)
                    m3 = l3*2+lscu(3,m)
                    m4 = l4*2+lscu(4,m)
                    m5 = l5*2+lscu(5,m)
                    m6 = l6*2+lscu(6,m)
                    m7 = l7*2+lscu(7,m)
                    m8 = l8*2+lscu(8,m)
                    do ldj=1,ndj
                      if (all([m1,m2,m3,m4,m5,m6,m7,m8] == locu(:,ldj))) then
                        cycle outer
                      end if
                    end do
                    ndj = ndj+1
                    locu(1,ndj) = m1
                    locu(2,ndj) = m2
                    locu(3,ndj) = m3
                    locu(4,ndj) = m4
                    locu(5,ndj) = m5
                    locu(6,ndj) = m6
                    locu(7,ndj) = m7
                    locu(8,ndj) = m8
                  end if
                end do outer
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

call mma_deallocate(lscu)

do nre=1,ndj
  write(u6,'(5x,i6,8i3)') nre,(locu(i,nre),i=1,8)
end do

return

end subroutine gugadrt_ref_gfs

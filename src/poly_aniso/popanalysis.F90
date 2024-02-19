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

subroutine PopAnalysis(nneq,neq,exch,nexch,nmax,lmax,NmaxPop,Z)
! computes the density matrix of each interacting site in the coupled eigenstate Zi
! N - (Integer) - total number of the exchange states ( total number of basis functions )
! Z - (Complex*16 array, size N) - contains the coupled eigenvector (a given column of the
!      entire exchange Z array)
! Nsites - total number of interacting magnetic sites
! ibas(N,Nsites) - Integer array denoting the coupled basis
! nneq -- number of non equivalent sites
! neq(Nneq) number of equivalent sites of type i

use Constants, only: cZero
use Definitions, only: u6

implicit none
#include "stdalloc.fh"
! main input variables
integer, intent(in) :: nneq, exch, nmax, lmax
integer, intent(in) :: neq(nneq), nexch(nneq)
integer, intent(in) :: NmaxPop
complex(kind=8), intent(in) :: Z(exch,exch)
! local variables
integer :: i, j, l, isite, i1, i2, nb1, nb2, nb3, tmp, il, nb
integer :: nind(lmax,2), intc(lmax)
integer :: ibas(exch,lmax)
complex(kind=8), allocatable :: pop(:,:,:,:) ! pop(exch,lmax,nmax,nmax)
character(len=50) :: fmtline

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'enter POPULATION ANALYSIS Subroutine'
write(u6,'(A, i3)') '      nmax = ',nmax
write(u6,'(A, i3)') '      exch = ',exch
write(u6,'(A, i3)') '   NmaxPop = ',NmaxPop
write(u6,'(A, i3)') '      lmax = ',lmax
write(u6,'(A, i3)') '      nneq = ',nneq
write(u6,'(A,8i3)') '    neq(i) = ',(neq(i),i=1,nneq)
write(u6,'(A,8i3)') '  nexch(i) = ',(nexch(i),i=1,nneq)
#endif
call mma_allocate(pop,exch,lmax,nmax,nmax,'pop')
call zcopy_(exch*lmax*nmax*nmax,[cZero],0,pop,1)
! fill some general arrays:
! generate the tables:
nind(:,:) = 0
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do
intc(:) = 0
ibas(:,:) = 0
intc(1) = 1
if (lmax > 1) then
  do i=2,lmax
    i1 = nind(i-1,1)
    intc(i) = intc(i-1)*nexch(i1)
  end do
end if
do nb=1,exch
  nb1 = nb-1
  do i=1,lmax
    ibas(nb,lmax-i+1) = nb1/intc(lmax-i+1)
    nb1 = nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
  end do
end do
isite = 0

pop(:,:,:,:) = cZero
i = int(24+13*lmax)
write(fmtline,'(A,i3,A,i2,A)') '(',i,'A)'
write(u6,fmtline) ('-',j=1,i)
i = int((24+13*lmax-19)/2)
write(fmtline,'(A,i3,A,i3,A)') '(',i,'X,A,',i,'X)'
write(u6,fmtline) 'POPULATION ANALYSIS'
i = int((8+13*lmax-19)/2)
write(fmtline,'(A,i3,A,i3,A)') '(',i,'X,A,',i,'X)'
write(u6,fmtline) '(i.e. diagonal value of density matrices of the interacting sites)'
i = int(24+13*lmax)
write(fmtline,'(A,i3,A)') '(',i,'A)'
write(u6,fmtline) ('-',j=1,i)
write(fmtline,'(A,i3,A)') '(A,',lmax,'A)'
write(u6,fmtline) '--------|-----|',('------------|',i=1,lmax)
write(u6,fmtline) 'Exchange|Basis|',('   center   |',i=1,lmax)
write(fmtline,'(A,i3,A)') '(A,',lmax,'(A,i2,A))'
write(u6,fmtline) ' state  | set |',('     ',i,'     |',i=1,lmax)
write(fmtline,'(A,i3,A)') '(A,',lmax,'A)'
write(u6,fmtline) '--------|-----|',('------------|',i=1,lmax)

! loop over all states for which we want density matrices and
! expectation values to be calculated
do nb1=1,NmaxPop
  do l=1,lmax
    isite = nind(l,1)
    do i1=1,nexch(isite)
      do i2=1,nexch(isite)
        ! sum over all other components of other sites
        do nb2=1,exch
          do nb3=1,exch
            if (ibas(nb2,l)+1 /= i1) cycle
            if (ibas(nb3,l)+1 /= i2) cycle
            tmp = 0
            do il=1,lmax
              if (il /= l) tmp = tmp+(ibas(nb2,il)-ibas(nb3,il))**2
            end do

            if (tmp == 0) pop(nb1,l,i1,i2) = pop(nb1,l,i1,i2)+conjg(Z(nb2,nb1))*Z(nb3,nb1)

          end do
        end do
      end do
    end do
  end do

  do i=1,nmax ! maxim local basis
    write(fmtline,'(A,i2,A)') '(2x,i4,2x,A,1x,i2,2x,A,',lmax,'(1x,F10.8,1x,A))'
    write(u6,fmtline) nb1,'|',i,'|',(real(pop(nb1,l,i,i)),'|',l=1,lmax)
  end do
  write(fmtline,'(A,i2,A)') '(A,',lmax,'A)'
  write(u6,fmtline) '--------|-----|',('------------|',i=1,lmax)
end do ! nb1

call mma_deallocate(pop)

! compute and print the calculated expectation values:
!write(u6,*)
!write(u6,'(A)') 'EXPECTATION VALUES'
!write(u6,'(5A)') '--------|----|',('------------------------------|',i=1,4)
!write(u6,'(5A)') 'Exchange|Site|','   MAGNETIC MOMENT (M=-L-2S)  |','        SPIN MOMENT (S)       |', &
!                 '      ORBITAL MOMENT (L)      |','     TOTAL MOMENT (J=L+S)     |'
!write(u6,'(5A)') ' state  | Nr.|',('     X         Y         Z    |',i=1,4)
!write(u6,'(5A)') '--------|----|',('------------------------------|',i=1,4)
!do nb1=1,NmaxPop
!  ! we proceed to compute expectation values for this nb1 exchange state
!  Mx(:) = cZero
!  My(:) = cZero
!  Mz(:) = cZero
!  Sx(:) = cZero
!  Sy(:) = cZero
!  Sz(:) = cZero
!  Jx(:) = cZero
!  Jy(:) = cZero
!  Jz(:) = cZero
!  Lx(:) = cZero
!  Ly(:) = cZero
!  Lz(:) = cZero
!
!  do l=1,lmax
!    isite = nind(ind_exch(l),1)
!    do i1=1,nexch(isite)
!      do i2=1,nexch(isite)
!        Mx(l) = Mx(l)+pop(nb1,l,i1,i2)*M(isite,1,i1,i2)
!        My(l) = My(l)+pop(nb1,l,i1,i2)*M(isite,2,i1,i2)
!        Mz(l) = Mz(l)+pop(nb1,l,i1,i2)*M(isite,3,i1,i2)
!        Sx(l) = Sx(l)+pop(nb1,l,i1,i2)*S(isite,1,i1,i2)
!        Sy(l) = Sy(l)+pop(nb1,l,i1,i2)*S(isite,2,i1,i2)
!        Sz(l) = Sz(l)+pop(nb1,l,i1,i2)*S(isite,3,i1,i2)
!      end do
!    end do
!    Lx(l) = -Mx(l)-g_e*Sx(l)
!    Ly(l) = -My(l)-g_e*Sy(l)
!    Lz(l) = -Mz(l)-g_e*Sz(l)
!    Jx(l) =  Lx(l)+Sx(l)
!    Jy(l) =  Ly(l)+Sy(l)
!    Jz(l) =  Lz(l)+Sz(l)
!  end do
!
!  do l=1,lmax
!    if (l == int((lmax+1)/2)) then
!      write(u6,'(i5,3x,A,1x,i2,1x,A,4(3(F9.5,1x),A))') nb1,'|',l,'|',real(Mx(l)),real(My(l)),real(Mz(l)),'|',real(Sx(l)), &
!                                                       real(Sy(l)),real(Sz(l)),'|',real(Lx(l)),real(Ly(l)),real(Lz(l)),'|', &
!                                                       real(Jx(l)),real(Jy(l)),real(Jz(l)),'|'
!    else
!      write(u6,'(8x,A,1x,i2,1x,A,4(3(F9.5,1x),A))') '|',l,'|',real(Mx(l)),real(My(l)),real(Mz(l)),'|',real(Sx(l)),real(Sy(l)), &
!                                                    real(Sz(l)),'|',real(Lx(l)),real(Ly(l)),real(Lz(l)),'|',real(Jx(l)), &
!                                                    real(Jy(l)),real(Jz(l)),'|'
!    end if
!  end do
!  write(u6,'(5A)') '--------|----|',('------------------------------|',i=1,4)
!end do

return

end subroutine PopAnalysis

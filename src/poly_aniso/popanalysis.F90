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
      Subroutine PopAnalysis(nneq,neq,exch,nexch,nmax,lmax,             &
     & NmaxPop,Z)
!  computes the density matrix of each interacting site in the coupled eigenstate Zi
!  N - (Integer) - total number of the exchange states ( total number of basis functions )
!  Z - (Complex*16 array, size N) - contains the coupled eigenvector (a given column of the
!       entire exchange Z array)
!  Nsites - total number of interacting magnetic sites
!  ibas(N,Nsites) - Integer array denoting the coupled basis
!  nneq -- number of non equivalent sites
!  neq(Nneq) number of equivalent sites of type i
!
      Implicit None
#include "stdalloc.fh"
      Integer, parameter        :: wp=kind(0.d0)
! main input variables
      Integer, intent(in)          :: nneq, exch, nmax, lmax
      Integer, intent(in)          :: neq(nneq), nexch(nneq)
      Integer, intent(in)          :: NmaxPop
      Complex(kind=8), intent(in) :: Z(exch,exch)
! local variables
      Integer          :: i,j,l,isite,i1,i2,nb1,nb2,nb3,tmp,il,nb
      Integer          :: nind(lmax,2),intc(lmax)
      Integer          :: ibas(exch,lmax)
!     pop(exch,lmax,nmax,nmax)
      Complex(kind=8), allocatable :: pop(:,:,:,:)
      Character(len=50):: fmtline
      Logical          :: DBG
      DBG=.false.

      If(DBG) Then
         Write(6,'(A)') 'enter POPULATION ANALYSIS Subroutine'
         Write(6,'(A, i3)') '      nmax = ',nmax
         Write(6,'(A, i3)') '      exch = ',exch
         Write(6,'(A, i3)') '   NmaxPop = ',NmaxPop
         Write(6,'(A, i3)') '      lmax = ',lmax
         Write(6,'(A, i3)') '      nneq = ',nneq
         Write(6,'(A,8i3)') '    neq(i) = ',(neq(i),i=1,nneq)
         Write(6,'(A,8i3)') '  nexch(i) = ',(nexch(i),i=1,nneq)
      End If
      Call mma_allocate(pop,exch,lmax,nmax,nmax,'pop')
      Call zcopy_(exch*lmax*nmax*nmax,[(0.0_wp,0.0_wp)],0,pop,1)
! fill some general arrays:
! generate the tables:
      nind(:,:)=0
      l=0
      Do i=1,nneq
         Do j=1,neq(i)
         l=l+1
         nind(l,1)=i
         nind(l,2)=j
         End Do
      End Do
      intc(:)=0
      ibas(:,:)=0
      intc(1)=1
      If (lmax.gt.1) Then
        Do i=2,lmax
        i1=nind(i-1, 1)
        intc(i)=intc(i-1)*nexch(i1)
        End Do
      End If
      Do nb=1,exch
      nb1=nb-1
         Do i=1,lmax
         ibas(nb, lmax-i+1)= nb1 / intc(lmax-i+1)
         nb1=nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
         End Do
      End Do
      isite=0
!
      pop=(0.0_wp,0.0_wp)
      i=int(24+13*lmax)
      Write(fmtline,'(A,i3,A,i2,A)') '(',i,'A)'
      Write(6,fmtline) ('-',j=1,i)
      i=int((24+13*lmax-19)/2)
      Write(fmtline,'(A,i3,A,i3,A)') '(',i,'X,A,',i,'X)'
      Write(6,fmtline)   'POPULATION ANALYSIS'
      i=int((8+13*lmax-19)/2)
      Write(fmtline,'(A,i3,A,i3,A)') '(',i,'X,A,',i,'X)'
      Write(6,fmtline)   '(i.e. diagonal value of density'//            &
     & ' matrices of the interacting sites)'
      i=int(24+13*lmax)
      Write(fmtline,'(A,i3,A)') '(',i,'A)'
      Write(6,fmtline) ('-',j=1,i)
      Write(fmtline,'(A,i3,A)') '(A,',lmax,'A)'
      Write(6,fmtline) '--------|-----|',                               &
     & ('------------|',i=1,lmax)
      Write(6,fmtline) 'Exchange|Basis|',                               &
     & ('   center   |', i=1,lmax)
      Write(fmtline,'(A,i3,A)') '(A,',lmax,'(A,i2,A))'
      Write(6,fmtline) ' state  | set |',                               &
     & ('     ',i,'     |',i=1,lmax)
      Write(fmtline,'(A,i3,A)') '(A,',lmax,'A)'
      Write(6,fmtline) '--------|-----|',                               &
     & ('------------|',i=1,lmax)

!     loop over all states for which we want density matrices and
!     expectation values to be calculated
      Do nb1=1,NmaxPop
         Do l=1,lmax
         isite=nind(l,1)
            Do i1=1,nexch(isite)
               Do i2=1,nexch(isite)
!     sum over all other components of other sites
                  Do nb2=1,exch
                     Do nb3=1,exch
                        If(ibas(nb2,l)+1 .ne. i1) Go To 11
                        If(ibas(nb3,l)+1 .ne. i2) Go To 11
                           tmp=0
                           Do il=1,lmax
                           If (il.eq.l) Go To 10
                   tmp=tmp+(ibas(nb2,il)-ibas(nb3,il))**2
 10                        continue
                           End Do

                        If(tmp.gt.0 ) Go To 11

      pop(nb1,l,i1,i2)=pop(nb1,l,i1,i2)+conjg(Z(nb2,nb1))*Z(nb3,nb1)

 11                     continue
                     End Do
                  End Do
               End Do
            End Do
         End Do

         Do i=1,nmax ! maxim local basis
         Write(fmtline,'(A,i2,A)') '(2x,i4,2x,A,1x,i2,2x,A,',           &
     &                              lmax,'(1x,F10.8,1x,A))'
         Write(6,fmtline) nb1,'|',i,'|',                                &
     & ( dble(pop(nb1,l,i,i)),'|',l=1,lmax)
         End Do
      Write(fmtline,'(A,i2,A)') '(A,',lmax,'A)'
      Write(6,fmtline) '--------|-----|',                               &
     & ('------------|',i=1,lmax)
      End Do ! nb1


      call mma_deallocate(pop)

      Return
      End


! compute and print the calculated expectation values:
!      Write(6,*)
!      Write(6,'(A)') 'EXPECTATION VALUES'
!      Write(6,'(5A)') '--------|----|',
!     & ('------------------------------|',i=1,4)
!      Write(6,'(5A)') 'Exchange|Site|',
!     & '   MAGNETIC MOMENT (M=-L-2S)  |',
!     & '        SPIN MOMENT (S)       |',
!     & '      ORBITAL MOMENT (L)      |',
!     & '     TOTAL MOMENT (J=L+S)     |'
!      Write(6,'(5A)') ' state  | Nr.|',
!     & ('     X         Y         Z    |',i=1,4)
!      Write(6,'(5A)') '--------|----|',
!     & ('------------------------------|',i=1,4)
!      Do nb1=1,NmaxPop
!c  we proceed to compute expectation values for this nb1 exchange state
!      Mx(:) = cZero
!      My(:) = cZero
!      Mz(:) = cZero
!      Sx(:) = cZero
!      Sy(:) = cZero
!      Sz(:) = cZero
!      Jx(:) = cZero
!      Jy(:) = cZero
!      Jz(:) = cZero
!      Lx(:) = cZero
!      Ly(:) = cZero
!      Lz(:) = cZero

!         Do l=1,lmax
!         isite = nind(ind_exch(l),1)
!           Do i1=1,nexch(isite)
!              Do i2=1,nexch(isite)
!       Mx(l)=Mx(l)+pop(nb1,l,i1,i2)*M(isite,1,i1,i2)
!       My(l)=My(l)+pop(nb1,l,i1,i2)*M(isite,2,i1,i2)
!       Mz(l)=Mz(l)+pop(nb1,l,i1,i2)*M(isite,3,i1,i2)
!       Sx(l)=Sx(l)+pop(nb1,l,i1,i2)*S(isite,1,i1,i2)
!       Sy(l)=Sy(l)+pop(nb1,l,i1,i2)*S(isite,2,i1,i2)
!       Sz(l)=Sz(l)+pop(nb1,l,i1,i2)*S(isite,3,i1,i2)
!              End Do
!            End Do
!       Lx(l)=-Mx(l)-g_e*Sx(l)
!       Ly(l)=-My(l)-g_e*Sy(l)
!       Lz(l)=-Mz(l)-g_e*Sz(l)
!       Jx(l)= Lx(l)+Sx(l)
!       Jy(l)= Ly(l)+Sy(l)
!       Jz(l)= Lz(l)+Sz(l)
!         End Do
!
!         Do l=1,lmax
!         If(l.eq.int((lmax+1)/2)) Then
!      Write(6,'(i5,3x,A,1x,i2,1x,A,4(3(F9.5,1x),A))')
!     & nb1,'|',l,'|',
!     & dble(Mx(l)),dble(My(l)),dble(Mz(l)),'|',
!     & dble(Sx(l)),dble(Sy(l)),dble(Sz(l)),'|',
!     & dble(Lx(l)),dble(Ly(l)),dble(Lz(l)),'|',
!     & dble(Jx(l)),dble(Jy(l)),dble(Jz(l)),'|'
!         Else
!      Write(6,'(8x,A,1x,i2,1x,A,4(3(F9.5,1x),A))')
!     & '|',l,'|',
!     & dble(Mx(l)),dble(My(l)),dble(Mz(l)),'|',
!     & dble(Sx(l)),dble(Sy(l)),dble(Sz(l)),'|',
!     & dble(Lx(l)),dble(Ly(l)),dble(Lz(l)),'|',
!     & dble(Jx(l)),dble(Jy(l)),dble(Jz(l)),'|'
!         End If
!         End Do
!      Write(6,'(5A)') '--------|----|',
!     & ('------------------------------|',i=1,4)
!      End Do
!

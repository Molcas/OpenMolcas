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
! Copyright (C) 2021, Vladislav Kochetov                               *
!***********************************************************************
subroutine kab
  use rhodyn_data
  use rhodyn_utils, only: transform, dashes
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none
!***********************************************************************
!
! Purpose :  calculate dissipation rates k_ab
!
!***********************************************************************

  integer :: max_i, max_j, iii, jjj, n_sf
  real(8) ::  max_k
  real(8), dimension(:,:), allocatable :: G,G_SF,omega_ab,n_w,J_w,r_ab_SO
  real(8), dimension(:), allocatable :: freq, temp_gk
  complex(8), dimension(:,:),allocatable ::k_ab,gamma_pd,gamma_pd_basis
  complex(8), dimension(:,:,:),allocatable :: G_SO
  character(len=256):: format1

  format1='(2(I8),4(G15.8,X))'
  n_sf = sum(nconf)

  write(*,*) 'N_SF=', n_sf
  write(*,*) 'Nmode=', Nmode
  write(*,*) 'Nconftot=', nconftot
  write(*,*) 'Nstate=', Nstate

  call mma_allocate(G,Nmode,n_sf)
  call mma_allocate(G_SF,Nmode,Nstate)
  call mma_allocate(G_SO,Nmode,Nstate,Nstate)
  call mma_allocate(freq,Nmode)
	call mma_allocate(temp_gk,Nstate)
  call mma_allocate(r_ab_SO,Nstate,Nstate)

! reading quintets and triplets

  call dashes()
  write(*,*) ' Begin reading HR-factors data '
  call dashes()

  jj=0
  iii=0

  open(11,file='HRFACT',status='old',iostat=error)
  if (error/=0) write(*,*) 'reading file HRFACT failed!'
  read(11,*)
  do k=1,Nmode
    read(11,*) freq(k)
  enddo
  read(11,*)
  if (HRSO) then
    do i=1,NState
      read(11,*)
      do k=1,Nmode
        write(*,*) 'k=', k
			  do j=1,Nstate
          read(11,'(E16.8)',advance='no') temp_gk(j)
			    write(*,*) temp_gk(j)
			    G_SO(k,i,j) = cmplx(temp_gk(j))
			  enddo
			  read(11,*)
      enddo
    enddo
  else
    do k=1,Nmode
      read(11,*)(G(k,j),j=2,n_sf)
      G(k,1)=0d0
      write(*,*) k, (G(k,j),j=2,n_sf)
    enddo
!!!!
! the gradient for each SF states is calculated using RASSCF while
! the other data are calculated in DFT (Gaussian), so there is a shift
! from the ground to ground state G(I,1)!=0, here we shift
! all the SF states with respect to ground states
    do k=1,n
      if(k==1)then
        iii=0
      else
        iii=iii+nconf(k-1)*ispin(k-1)
      endif
      jj=jj+nconf(k)
      do i=1,Nmode
        jjj=iii
        do j=jj-nconf(k)+1,jj
          do ii=1,ispin(k)
            jjj=jjj+1
            G_SF(i,jjj)=G(i,j)
          enddo
        enddo
      enddo
    enddo
    call dashes(72)
    do i=1,Nmode
      write(*,*)(G_SF(i,j),j=1,nconftot)
    enddo
! contruct the G_SO matrix
    do i=1,NState
      do j=1,NState
        do k=1,Nmode
          do ii=1,Nstate
            G_SO(k,i,j)=G_SO(k,i,j)+&
                conjg(SO_CI(ii,i))*G_SF(k,ii)*SO_CI(ii,j)
          enddo
        enddo
      enddo
    enddo
	endif !HRSO
  close(11)

! change the unit to au
  Freq=Freq*cmtoau

  call dashes()
  write(*,*)' End read data '
  call dashes()

! Spectrum density

  call mma_allocate(k_ab,Nstate,Nstate)
  call mma_allocate(n_w,Nstate,Nstate)
  call mma_allocate(J_w,Nstate,Nstate)
  call mma_allocate(omega_ab,Nstate,Nstate)

  write(*,*) 'Gamma=',gamma, 'Hartree',gamma/cmtoau, 'cm-1'

  open (20,file='kab_out.dat',status='replace')
  J_w=0d0
  do i=1,Nstate
! let us cut the spectral density (no), if yes then j=1,i
    do j=1,Nstate
      omega_ab(i,j)=E_SO(i)-E_SO(j)
      do k=1,Nmode
        J_w(i,j)=J_w(i,j)+abs(G_SO(k,i,j))**2*Freq(k)**2*&
                  (abs(omega_ab(i,j))*Freq(k)*gamma/&
                  ((abs(omega_ab(i,j))**2-Freq(k)**2)**2+&
                  omega_ab(i,j)**2*gamma**2))
      enddo
      if (i==j) then
        n_w(i,j)=0d0
      else
        n_w(i,j)=1d0/(exp(abs(omega_ab(i,j))/(k_b*T))-1d0)
      endif
! note that the K_ab matrix should be real matrix, but here we define
! it to be complex for the convernient of the following transform
      if (omega_ab(i,j)>=0d0) then
        k_ab(i,j)=4*(1+n_w(i,j))*J_w(i,j)
      else
        k_ab(i,j)=4*n_w(i,j)*J_w(i,j)
      endif
      if(real(k_ab(i,j))>=(0.01/autoeV))then
        if(omega_ab(i,j)>=0d0)then
          write (20,format1) i,j,dble(k_ab(i,j))*autoev, &
                omega_ab(i,j)**2,(1.0+n_w(i,j)),J_w(i,j)
        else
            write (20,format1) i,j,dble(k_ab(i,j))*autoev, &
                 omega_ab(i,j)**2,n_w(i,j),J_w(i,j)
          endif
        endif

      enddo
    enddo

    close(20)

! for this  spectral density the pure decay matrix gamma_pd equal
! to 0.0d0 always because when the frequency equal
! to 0d0, J_w(0d0)=0d0 so gamma_pd=0
        call mma_allocate(gamma_pd,Nstate,Nstate)
        gamma_pd=0d0

        if (ipglob>3) then
          call dashes()
          write(*,*)' Print matrix K_ab '
          call dashes()

          open(13,file='Kab_matrix_eV.dat',status='replace')
!!vk!! write procedure for printing matrices
    do i=1,Nstate
      write(13,*)(dble(k_ab(i,j))*autoev,j=1,Nstate)
    enddo
    max_k=0d0
    do i=1,Nstate
      do j=1,Nstate
        if (real(k_ab(i,j))>=max_k) then
          max_k=k_ab(i,j)
          max_i=i
          max_j=j
        endif
      enddo
    enddo
    write(*,*) Max_I,Max_J,Max_K*autoev,' eV', &
      omega_ab(max_i,max_j),1+n_w(Max_I,Max_J),J_w(Max_I,Max_J)
    close(13)
  endif


! transform the k_ab and gamma_pd matrix to the required basis

  call mma_allocate(gamma_pd_basis,d,d)

  select case (basis)
  case ('CSF')
    ! k_ab SO->CSF
    call transform(k_ab,CSF2SO,kab_basis,.False.)
    ! gamma_pd SO->CSF
    call transform(gamma_pd,CSF2SO,gamma_pd_basis,.False.)
  case ('SF')
    ! K_ab SO->SF
    call transform(k_ab,SO_CI,kab_basis,.False.)
    ! gamma_pd SO->SF
    call transform(gamma_pd,SO_CI,gamma_pd_basis,.False.)
    ! Kab_basis = dble(Kab_basis)
  case ('SO')
    kab_basis=k_ab
    gamma_pd_basis=gamma_pd
  end select

! print out the bigger kab_basis

  open (22,file='max_Kab_basis.dat',status='replace')

  max_k=0d0
  do i=1,d
    do j=1,d
      if (abs(kab_basis(i,k))>=max_k) then
        max_k=abs(kab_basis(i,j))
        max_i=i
        max_j=j
      endif
    enddo
  enddo

  write(22,*)'the maximum of Kab in ', basis
  write(22,'(2(I8),G15.8,A)') max_i,max_j,max_k*autoev,' eV'

  do i=1,d
    do j=1,d
      if(abs(kab_basis(i,j))>=(0.01/autoeV))then
        write(22,'(2(I8),3(G15.8,X))')i,j,abs(kab_basis(i,j)), &
        dble(kab_basis(I,J))*autoev,aimag(kab_basis(I,J))*autoev
      endif
    enddo
  enddo


!        if (basis=='SO') then
!          Max_K=0d0
!          do i=161,Nstate
!            do j=161,Nstate
!              if (abs(Kab_basis(i,j))>=Max_K) then
!                Max_K=abs(Kab_basis(i,j))
!                Max_I=i
!                Max_J=j
!              endif
!            enddo
!          enddo

!          write(22,*)'Maximum of Kab in core states in basis of ', Basis
!          write(22,'(2(I8),G15.8,A)')Max_I,Max_J,Max_K*autoev,' eV'

!          do i=1,Nstate
!            do j=1,Nstate
!              if(abs(Kab_basis(i,j))<=(0.01/autoeV).and.&
!               abs(Kab_basis(i,j))>=(0.005/autoeV))then
!              write(22,'(2(I8),3(G15.8,X))')I,J,abs(Kab_basis(i,j)),&
!              dble(Kab_basis(i,j))*autoev,aimag(Kab_basis(i,j))*autoev
!              endif
!            enddo
!          enddo

!      elseif (basis=='SF') then

!        Max_k=0d0
!        do i=26,175
!          do j=26,175
!            if(abs(Kab_basis(I,J))>=Max_K)then
!              Max_K=abs(Kab_basis(i,j))
!              Max_i=i
!              Max_j=j
!            endif
!          Enddo
!        Enddo

!        write(22,*) 'Maximum of Kab in quintet core states'//&
!                    ' in basis of ', basis
!        write(22,'(2(I8),G15.8,A)') Max_I,Max_J,Max_K*autoev,' eV'

!        Max_k=0d0
!        do i=311,Nstate
!          do j=311,Nstate
!            if(abs(Kab_basis(I,J))>=Max_K)then
!              Max_K=abs(Kab_basis(I,J))
!              Max_i=i
!              Max_j=j
!            endif
!          Enddo
!        Enddo

!        write(22,*)'Maximum of Kab in triplets core states'//&
!                   ' in basis of ', Basis
!        write(22,'(2(I8),G15.8,A)')Max_i,Max_j,Max_k*autoev,' eV'

!        do i=1,Nstate
!          do j=1,Nstate
!            if(abs(Kab_basis(i,J))<=(0.01/autoeV).and.&
!               abs(Kab_basis(I,J))>=(0.005/autoeV)) then
!              write(22,'(2(I8),3(G15.8,X))')I,J,abs(Kab_basis(I,J)),&
!              dble(Kab_basis(I,J))*autoev,aimag(Kab_basis(I,J))*autoev
!            endif
!          enddo
!        enddo
!      endif

  close(22)

! contruct the matrix (k_bar)_ij=0.5*sum_k[(kab_basis)_ik+(kab_basis)_jk]

  do i=1,d
    do j=1,d
      do k=1,d
        k_bar_basis(j,i) = abs(k_bar_basis(j,i)+0.5d0* &
                         (kab_basis(j,k)+kab_basis(i,k)))
      enddo
    enddo
  enddo

  K_bar_basis=K_bar_basis+gamma_pd_basis
  open (23,file='r_ab_SO.dat',status='replace')

! check whether coherence dephasing r_ab=0.5*[sum_e(K_ae+K_be)]>0 in SO basis
  do i=1,Nstate
    do j=1,Nstate
      do k=1,Nstate
        r_ab_SO(j,i) = r_ab_SO(j,i)+0.5d0* &
                       (dble(K_ab(j,k))+dble(K_ab(i,k)))
      enddo
    enddo
  enddo

  do i=1,Nstate
    do j=1,Nstate
      if (r_ab_SO(i,j)<=0d0) then
        write(23,*)'R_ab_SO<0', i, j,r_ab_SO(i,j)
      endif
    enddo
  enddo
  close(23)
  write(*,*) 'End k_ab'

  call mma_deallocate(G)
  call mma_deallocate(G_SF)
  call mma_deallocate(G_SO)
  call mma_deallocate(Freq)
  call mma_deallocate(temp_gk)
  call mma_deallocate(r_ab_SO)
  call mma_deallocate(k_ab)
  call mma_deallocate(n_w)
  call mma_deallocate(J_w)
  call mma_deallocate(omega_ab)
  call mma_deallocate(gamma_pd)
  call mma_deallocate(gamma_pd_basis)

end

subroutine cut_matrices()
  use rhodyn_data
  use rhodyn_utils, only: dashes, removeLineAndColumn, removeColumn
  implicit none
!
! cut dynamics matrices from size (lrootstot,lrootstot) to size (d,d)
! not checked yet
  complex(8), dimension(:,:), allocatable :: d1, d2, d3 

  call dashes()
  write(*,*) 'Begin cut_matrices'
  call dashes()

  ! cut dynamics matrices
  call removeLineAndColumn(hamiltonian, istates)
  call removeLineAndColumn(density0, istates)
  if (flag_dyson) call removeLineAndColumn(dysamp_bas,istates)
  ! cut dipole moment
  allocate(d1(lrootstot,lrootstot))
  allocate(d2(lrootstot,lrootstot))
  allocate(d3(lrootstot,lrootstot))
  d1 = dipole_basis(:,:,1)
  d2 = dipole_basis(:,:,2)
  d3 = dipole_basis(:,:,3)
  call removeLineAndColumn(d1,istates)
  call removeLineAndColumn(d2,istates)
  call removeLineAndColumn(d3,istates)
  deallocate(dipole_basis)
  allocate(dipole_basis(d,d,3))
  dipole_basis(:,:,1) = d1
  dipole_basis(:,:,2) = d2
  dipole_basis(:,:,3) = d3
  deallocate(d1)
  deallocate(d2)
  deallocate(d3)

  ! cut transform marices
  call removeColumn(U_CI,istates)
  call removeColumn(CSF2SO,istates)
  call removeLineAndColumn(SO_CI,istates)
end subroutine cut_matrices
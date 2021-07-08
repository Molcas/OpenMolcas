subroutine sortci(N1,A,WR,C)
  use rhodyn_data
  use rhodyn_utils, only: dashes, transform
  implicit none
  integer::N1,INFO,LWORK
  real(8), dimension(N1,N1) :: A, C, diag, B
  real(8), dimension(N1) :: WR
  real(8), dimension (2*N1)::WORK
  B=A
  LWORK=2*N1
  call dsyev('V', 'U', N1, A, N1, WR, WORK, LWORK, INFO)
  if (INFO/=0) then
    write(*,*) 'ERROR in sortci'
    call Abend()
  endif
  call dsyev('V', 'U', N1, A, N1, WR, WORK, LWORK, INFO)
  C=A
  if (ipglob>3) then
    call transform(B,C,diag)
    call dashes(72)
    write(*,*) 'Printout the diagonalized H(RASSCF) matrix'
    call dashes(72)
    do k=1,10
      write(*,*) (diag(k,l),l=1,10)
    enddo
  endif
end subroutine sortci

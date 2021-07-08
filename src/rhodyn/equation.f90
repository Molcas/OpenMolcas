subroutine equation (time,rhot,res)
  use rhodyn_data
  implicit none
!
! Liouville equation is solved here
!
  real(8) :: time
  complex(8),dimension(:,:) :: rhot, res
  procedure(pulse_func) :: pulse

  if (flag_pulse) call pulse(hamiltonian,hamiltoniant,time)
  call ZGEMM('N','N',d,d,d,-onei,hamiltoniant,d,rhot,d,zero,res,d)
  call ZGEMM('N','N',d,d,d, onei,rhot,d,hamiltoniant,d, one,res,d)

! auger decay part
  if (flag_decay.or.ion_diss/=0d0) then
    call ZGEMM('N','N',d,d,d,one,decay,d,rhot,d,one,res,d)
  endif

  if (flag_diss) then
    do i=1,d
      do j=1,d
        if (i/=j) then
          res(i,j) = res(i,j) - K_bar_basis(i,j) * rhot(i,j)
        endif
        res(i,i) = res(i,i) - Kab_basis(i,j) * rhot(i,i) + &
                              Kab_basis(j,i) * rhot(j,j)
      enddo
    enddo
  endif

end

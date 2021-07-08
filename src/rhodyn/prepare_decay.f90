subroutine prepare_decay
  use rhodyn_data
  use rhodyn_utils, only: transform, mult, dashes
  implicit none
! energy and time should fullfill relation delta_E*delta_t=h
! construct the decay in SOC states basis sets

  if (ipglob>3) write(*,*) 'Begin of prepare_decay'

  decay=zero
! Auger decay
  if (flag_decay) then
    do i=Nval+1,Nval+N_L3
      decay(i,i)=-tau_L3/2/pi
    enddo
    do i=Nval+N_L3+1,Nstate
      decay(i,i)=-tau_L2/2/pi
    enddo
    if (basis=='CSF') then
      call mult(CSF2SO,decay,tmp)
      call mult(tmp,CSF2SO,decay,.False.,.True.)
    elseif (basis=='SF') then
      call mult(SO_CI,decay,tmp)
      call mult(tmp,SO_CI,decay,.False.,.True.)
    endif
  endif

! ionization
  if (flag_dyson.and.ion_diss/=0d0) then
    ii=1
    do k=1,N
      do i=ii,(ii+nconf(k)*ispin(k)-1)
        if (ion_blocks(k)) decay(i,i) = decay(i,i) - ion_diss
      enddo
      ii=ii+nconf(k)*ispin(k)
    enddo
    if (basis=='CSF') then
      call mult(U_CI_compl,decay,tmp)
      call mult(tmp,U_CI_compl,decay,.False.,.True.)
    elseif (basis=='SO') then
      call mult(SO_CI,decay,tmp,.True.,.False.)
      call mult(tmp,SO_CI,decay)
    endif
  endif

!!!!!!!!!!
  if (ipglob>4) then
    call dashes()
    write(*,*) 'Decay matrix'
    do i=1,Nstate
      write(*,*)(decay(i,j),j=1,Nstate)
    enddo
    call dashes()
  endif
!!!!!!!!!!

  if (ipglob>3) write(*,*) 'End of prepare_decay'

end

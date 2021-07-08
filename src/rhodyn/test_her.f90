subroutine test (densityt_time, time)
  use rhodyn_data
  implicit none

  complex(8), dimension(Nstate,Nstate) :: densityt_time
  real(8) :: time,abserror

  abserror=0d0
  do i=1,Nstate
    do j=(i+1),Nstate
      if ((abs(dble(densityt_time(i,j))-dble(densityt_time(j,i)))>=&
                      threshold).and.(abs(dble(densityt_time(i,j))-&
                      dble(densityt_time(j,i)))>=abserror)) then
        abserror = abs(dble(densityt_time(i,j))-&
                       dble(densityt_time(j,i)))
      endif
      if ((abs(aimag(densityt_time(i,j))+aimag(densityt_time(j,i)))>=&
          threshold).and.(abs(aimag(densityt_time(i,j))+&
          aimag(densityt_time(j,i)))>=abserror)) then
        abserror = abs(aimag(densityt_time(i,j)) + aimag(densityt_time(j,i)))
      endif
    enddo
  enddo
  if (abserror>=threshold) then
    write(6,'(2(A,X,G28.16,X))')'time=',time/fstoau,'error=',abserror
  endif
end

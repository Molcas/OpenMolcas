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

subroutine JKQPar(N1,N2,HEXCH,Jpar)

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: N1, N2
complex(kind=8), intent(in) :: HEXCH(N1,N1,N2,N2)
complex(kind=8), intent(out) :: Jpar((N1-1),-(N1-1):(N1-1), (N2-1),-(N2-1):(N2-1))
! local variables
integer :: k1, k2, q1, q2, ipr, i1, i2, j1, j2, i
real(kind=8) :: F, THRS, R, knm(12,0:12)
complex(kind=8) :: F1, F2, F12, CI
complex(kind=8) :: trace_exch, trace, fact
complex(kind=8), allocatable :: O1(:,:), W1(:,:)
complex(kind=8), allocatable :: O2(:,:), W2(:,:)
complex(kind=8), allocatable :: O1_O2(:,:,:,:) !(N1,N1,N2,N2)
complex(kind=8), allocatable :: O1_W2(:,:,:,:) !(N1,N1,N2,N2)
complex(kind=8), allocatable :: W1_O2(:,:,:,:) !(N1,N1,N2,N2)
complex(kind=8), allocatable :: W1_W2(:,:,:,:) !(N1,N1,N2,N2)
complex(kind=8) :: Jcc((N1-1),0:(N1-1),(N2-1),0:(N2-1))
complex(kind=8) :: Jcs((N1-1),0:(N1-1),(N2-1),0:(N2-1))
complex(kind=8) :: Jsc((N1-1),0:(N1-1),(N2-1),0:(N2-1))
complex(kind=8) :: Jss((N1-1),0:(N1-1),(N2-1),0:(N2-1))
external :: trace_exch, trace
logical :: DBG
!real(kind=8) :: cm_to_MHz
!-----------------------------------------------------------------------

knm = 0.0_wp
call Set_knm(knm)
call mma_allocate(O1,N1,N1,'O1')
call mma_allocate(O2,N2,N2,'O2')
call mma_allocate(W1,N1,N1,'W1')
call mma_allocate(W2,N2,N2,'W2')
call mma_allocate(O1_O2,N1,N1,N2,N2,'O1_O2')
call mma_allocate(O1_W2,N1,N1,N2,N2,'O1_W2')
call mma_allocate(W1_O2,N1,N1,N2,N2,'W1_O2')
call mma_allocate(W1_W2,N1,N1,N2,N2,'W1_W2')

!-----------------------------------------------------------------------

!cm_to_MHz = 29979.2458_wp
DBG = .false.
ipr = 1

! we need to project now the HEXCH: in products of ITOs
!  HEXCH = SUM(rank1,proj1,rank2,proj2)=
!         { B(rank1,proj1,rank2,proj2)* O1(rank1,proj1) * O2(rank2,proj2) }
! Naoya definition
! eq.40 in DoI:10.1103/PhysRevB.91.174438
Jpar = (0.0_wp,0.0_wp)
do k1=1,N1-1
  do q1=0,k1

    do k2=1,N2-1
      do q2=0,k2

        call zcopy_(N1*N1,[(0.0_wp,0.0_wp)],0,O1,1)
        call zcopy_(N2*N2,[(0.0_wp,0.0_wp)],0,O2,1)
        call zcopy_(N1*N1,[(0.0_wp,0.0_wp)],0,W1,1)
        call zcopy_(N2*N2,[(0.0_wp,0.0_wp)],0,W2,1)
        ! get the ITOs for each site:
        call Stewens_matrixel(k1,q1,N1,O1,W1,ipr)
        call Stewens_matrixel(k2,q2,N2,O2,W2,ipr)

        if ((k1 == 1) .and. (q1 == 1)) then
          call pa_prmat('JKQPar:  O1(k1=1,q1=1): ',O1,n1)
          call pa_prmat('JKQPar:  W1(k1=1,q1=1): ',W1,n1)
        else if ((k1 == 1) .and. (q1 == 0)) then
          call pa_prmat('JKQPar:  O1(k1=1,q1=0): ',O1,n1)
          call pa_prmat('JKQPar:  W1(k1=1,q1=0): ',W1,n1)
        end if

        !if (DBG) then
        !  write(6,*)
        !  write(6,'(A)') '--------------------------------'
        !  write(6,*)
        !  write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O1  k1 = ',k1,'q1 =',q1
        !  write(6,*)
        !  do i1=1,N1
        !    write(6,'(20(2ES14.7,1x))') ((0.5_wp, 0.0_wp)*(cmplx((-1)**q1)*W1(i1,j1)+O1(i1,j1)),j1=1,N1)
        !  end do
        !  write(6,*)
        !  write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W1  k1 = ',k1,'q1 =',q1
        !  write(6,*)
        !  do i1=1,N1
        !    write(6,'(20(2ES14.7,1x))') ((0.0_wp,-0.5_wp)*(cmplx((-1)**q1)*W1(i1,j1)-O1(i1,j1)),j1=1,N1)
        !  end do
        !  write(6,*)
        !  write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O2  k2 = ',k2,'q2 =',q2
        !  write(6,*)
        !  do i2=1,N2
        !    write(6,'(20(2ES14.7,1x))') ((0.5_wp, 0.0_wp)*(cmplx((-1)**q2)*W2(i2,j2)+O2(i2,j2)),j2=1,N2)
        !  end do
        !  write(6,*)
        !  write(6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W2  k2 = ',k2,'q2 =',q2
        !  write(6,*)
        !  do i2=1,N2
        !    write(6,'(20(2ES14.7,1x))') ((0.0_wp,-0.5_wp)*(cmplx((-1)**q2)*W2(i2,j2)-O2(i2,j2)),j2=1,N2)
        !  end do
        !end if

        ! Build 4 coupled tensor products:
        ! O1-O2
        ! O1-W2
        ! W1-O2
        ! W1-W2
        call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,O1_O2,1)
        call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,O1_W2,1)
        call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,W1_O2,1)
        call zcopy_(N1*N1*N2*N2,[(0.0_wp,0.0_wp)],0,W1_W2,1)
        do i1=1,N1
          do j1=1,N1
            do i2=1,N2
              do j2=1,N2
                O1_O2(i1,j1,i2,j2) = O1(i1,j1)*O2(i2,j2)
                O1_W2(i1,j1,i2,j2) = O1(i1,j1)*W2(i2,j2)
                W1_O2(i1,j1,i2,j2) = W1(i1,j1)*O2(i2,j2)
                W1_W2(i1,j1,i2,j2) = W1(i1,j1)*W2(i2,j2)
              end do
            end do
          end do
        end do

        ! SP_HZFSO=trace(nDIMcf,HCF,DIP_O)
        ! SP_HZFSW=trace(nDIMcf,HCF,DIP_W)
        ! SP_MOW  =trace(nDIMcf,DIP_O,DIP_W)

        ! B(N,-M)=SP_HZFSO/SP_MOW
        ! B(N, M)=SP_HZFSW/SP_MOW

        FACT = (0.0_wp,0.0_wp)
        FACT = trace(N1,O1,W1)*trace(N2,O2,W2)
        !write(6,'(A,4I3,2F20.13)') 'k1,q1,k2,q2,FACT=',k1,q1,k2,q2,FACT

        Jpar(k1,-q1,k2,-q2) = trace_exch(N1,N2,HEXCH,O1_O2)/FACT
        Jpar(k1,q1,k2,-q2) = trace_exch(N1,N2,HEXCH,W1_O2)/FACT
        Jpar(k1,-q1,k2,q2) = trace_exch(N1,N2,HEXCH,O1_W2)/FACT
        Jpar(k1,q1,k2,q2) = trace_exch(N1,N2,HEXCH,W1_W2)/FACT

        ! Jpar( k1,-q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, O1_O2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1, q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, W1_O2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1,-q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, O1_W2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1, q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, W1_W2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT

        if (DBG) then
          if ((q1 == 0) .and. (q2 == 0)) then
            write(6,'(A,2ES18.10)') 'FACT = ',FACT
            if (abs(Jpar(k1,-q1,k2,-q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,0,k2,0)
          else if ((q1 == 0) .and. (q2 /= 0)) then
            write(6,'(A,2ES18.10)') 'FACT = ',FACT
            if (abs(Jpar(k1,0,k2,-q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',-q2,')=',Jpar(k1,0,k2,-q2)
            if (abs(Jpar(k1,0,k2,q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,0,k2,q2)
          else if ((q1 /= 0) .and. (q2 == 0)) then
            write(6,'(A,2ES18.10)') 'FACT = ',FACT
            if (abs(Jpar(k1,0,k2,-q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',q2,')=',Jpar(k1,-q1,k2,0)
            if (abs(Jpar(k1,0,k2,q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,q1,k2,0)
          else if ((q1 /= 0) .and. (q2 /= 0)) then
            if (abs(Jpar(k1,-q1,k2,-q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',-q2,')=',Jpar(k1,-q1,k2,-q2)
            if (abs(Jpar(k1,q1,k2,-q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',-q2,')=',Jpar(k1,q1,k2,-q2)
            if (abs(Jpar(k1,-q1,k2,q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',q2,')=',Jpar(k1,-q1,k2,q2)
            if (abs(Jpar(k1,q1,k2,q2)) > 0.5d-13) &
              write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,q1,k2,q2)
          end if
        end if ! DBG

      end do
    end do
  end do
end do

Jcc((N1-1),0:(N1-1),(N2-1),0:(N2-1)) = (0.0_wp,0.0_wp)
Jcs((N1-1),0:(N1-1),(N2-1),0:(N2-1)) = (0.0_wp,0.0_wp)
Jsc((N1-1),0:(N1-1),(N2-1),0:(N2-1)) = (0.0_wp,0.0_wp)
Jss((N1-1),0:(N1-1),(N2-1),0:(N2-1)) = (0.0_wp,0.0_wp)

do k1=1,N1-1
  do k2=1,N2-1
    do q1=0,k1
      do q2=0,k2

        F1 = cmplx((-1)**(-q1),0,wp)
        F2 = cmplx((-1)**(-q2),0,wp)
        F12 = cmplx((-1)**(-q1-q2),0,wp)
        CI = (0.0_wp,-1.0_wp)

        Jcc(k1,q1,k2,q2) = Jpar(k1,q1,k2,q2)+F2*Jpar(k1,q1,k2,-q2)+F1*Jpar(k1,-q1,k2,q2)+F12*Jpar(k1,-q1,k2,-q2)

        Jss(k1,q1,k2,q2) = Jpar(k1,q1,k2,q2)-F2*Jpar(k1,q1,k2,-q2)-F1*Jpar(k1,-q1,k2,q2)+F12*Jpar(k1,-q1,k2,-q2)

        Jcs(k1,q1,k2,q2) = (Jpar(k1,q1,k2,q2)-F2*Jpar(k1,q1,k2,-q2)+F1*Jpar(k1,-q1,k2,q2)-F12*Jpar(k1,-q1,k2,-q2))*CI

        Jsc(k1,q1,k2,q2) = (Jpar(k1,q1,k2,q2)+F2*Jpar(k1,q1,k2,-q2)-F1*Jpar(k1,-q1,k2,q2)-F12*Jpar(k1,-q1,k2,-q2))*CI

        if (DBG) then
          if (abs(Jcc(k1,q1,k2,q2)) > 0.5d-13) &
            write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jcc(',k1,',',q1,',',k2,',',q2,')=',Jcc(k1,q1,k2,q2)
          if (abs(Jcs(k1,q1,k2,q2)) > 0.5d-13) &
            write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jcs(',k1,',',q1,',',k2,',',q2,')=',Jcs(k1,q1,k2,q2)
          if (abs(Jsc(k1,q1,k2,q2)) > 0.5d-13) &
            write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jsc(',k1,',',q1,',',k2,',',q2,')=',Jsc(k1,q1,k2,q2)
          if (abs(Jss(k1,q1,k2,q2)) > 0.5d-13) &
            write(6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jss(',k1,',',q1,',',k2,',',q2,')=',Jss(k1,q1,k2,q2)
        end if

      end do
    end do
  end do
end do

if (DBG) then
  THRS = 1.0E-16_wp
  write(6,'(A,ES18.7)') 'Real Exchange parameters with values larger than: ',THRS
  write(6,'(120A)') ('-',i=1,119),'|'
  write(6,'(A)') ' k1 |  q1 || k2 |  q2 |-------- O1-O2 --------|-------- O1-W2 --------|-------- W1-O2 --------|'// &
                 '-------- W1-W2 --------|'
  do k1=1,N1-1
    do q1=-k1,k1
      write(6,'(A)') '----|-----||----|-----|-----------------------|-----------------------|-----------------------|'// &
                     '-----------------------|'
      do k2=1,N2-1
        do q2=-k2,k2

          F = knm(k1,abs(q1))*knm(k2,abs(q2))
          R = abs(Jcc(k1,abs(q1),k2,abs(q2)))+abs(Jcs(k1,abs(q1),k2,abs(q2)))+abs(Jsc(k1,abs(q1),k2,abs(q2)))+ &
              abs(Jss(k1,abs(q1),k2,abs(q2)))

          if (R > THRS) then

            if ((q1 < 0) .and. (q2 < 0)) then

              write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),4(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                         dble(Jcc(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                         dble(Jcs(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                         dble(Jsc(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                         dble(Jss(k1,abs(q1),k2,abs(q2)))*F,'|'

            else if ((q1 >= 0) .and. (q2 < 0)) then

              write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                         dble(Jcs(k1,abs(q1),k2,abs(q2)))*F,'|'

            else if ((q1 < 0) .and. (q2 >= 0)) then

              write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                         dble(Jsc(k1,abs(q1),k2,abs(q2)))*F,'|'

            else ! (q1 > 0) .and. (q2 > 0)

              write(6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                         dble(Jcc(k1,abs(q1),k2,abs(q2)))*F,'|'

            end if
          end if

        end do
      end do
    end do
  end do
  write(6,'(120A)') ('-',i=1,119),'|'

end if

!-----------------------------------------------------------------------
call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)
call mma_deallocate(O1_O2)
call mma_deallocate(O1_W2)
call mma_deallocate(W1_O2)
call mma_deallocate(W1_W2)

return

end subroutine JKQPar

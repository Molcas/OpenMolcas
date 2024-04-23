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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne, Onei
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: N1, N2
complex(kind=wp), intent(in) :: HEXCH(N1,N1,N2,N2)
complex(kind=wp), intent(out) :: Jpar(N1-1,-(N1-1):N1-1,N2-1,-(N2-1):N2-1)
integer(kind=iwp) :: i2, ipr, j2, k1, k2, q1, q2
real(kind=wp) :: knm(12,0:12)
complex(kind=wp) :: F1, F12, F2, fact
#ifdef _DEBUGPRINT_
real(kind=wp) :: F, R
#endif
complex(kind=wp), allocatable :: Jcc(:,:,:,:), Jcs(:,:,:,:), Jsc(:,:,:,:), Jss(:,:,:,:), O1(:,:), O1_O2(:,:,:,:), O1_W2(:,:,:,:), &
                                 O2(:,:), W1(:,:), W1_O2(:,:,:,:), W1_W2(:,:,:,:), W2(:,:)
real(kind=wp), parameter :: THRS = 1.0e-16_wp !, cm_to_MHz = cm_s*1.0e-6_wp
complex(kind=wp), external :: trace, trace_exch
!-----------------------------------------------------------------------

call Set_knm(knm)
call mma_allocate(O1,N1,N1,'O1')
call mma_allocate(O2,N2,N2,'O2')
call mma_allocate(W1,N1,N1,'W1')
call mma_allocate(W2,N2,N2,'W2')
call mma_allocate(O1_O2,N1,N1,N2,N2,'O1_O2')
call mma_allocate(O1_W2,N1,N1,N2,N2,'O1_W2')
call mma_allocate(W1_O2,N1,N1,N2,N2,'W1_O2')
call mma_allocate(W1_W2,N1,N1,N2,N2,'W1_W2')
call mma_allocate(Jcc,[1,N1-1],[0,N1-1],[1,N2-1],[0,N2-1],label='Jcc')
call mma_allocate(Jcs,[1,N1-1],[0,N1-1],[1,N2-1],[0,N2-1],label='Jcs')
call mma_allocate(Jsc,[1,N1-1],[0,N1-1],[1,N2-1],[0,N2-1],label='Jsc')
call mma_allocate(Jss,[1,N1-1],[0,N1-1],[1,N2-1],[0,N2-1],label='Jss')

!-----------------------------------------------------------------------

ipr = 1

! we need to project now the HEXCH: in products of ITOs
!  HEXCH = SUM(rank1,proj1,rank2,proj2)=
!         { B(rank1,proj1,rank2,proj2)* O1(rank1,proj1) * O2(rank2,proj2) }
! Naoya definition
! eq.40 in Doi:10.1103/PhysRevB.91.174438
Jpar(:,:,:,:) = cZero
do k1=1,N1-1
  do q1=0,k1

    do k2=1,N2-1
      do q2=0,k2

        O1(:,:) = cZero
        O2(:,:) = cZero
        W1(:,:) = cZero
        W2(:,:) = cZero
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

#       ifdef _DEBUGPRINT_
        !write(u6,*)
        !write(u6,'(A)') '--------------------------------'
        !write(u6,*)
        !write(u6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O1  k1 = ',k1,'q1 =',q1
        !write(u6,*)
        !do i1=1,N1
        !  write(6,'(20(2ES14.7,1x))') (Half*((-One)**q1*W1(i1,j1)+O1(i1,j1)),j1=1,N1)
        !end do
        !write(u6,*)
        !write(u6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W1  k1 = ',k1,'q1 =',q1
        !write(u6,*)
        !do i1=1,N1
        !  write(u6,'(20(2ES14.7,1x))') (-Half*Onei*((-One)**q1*W1(i1,j1)-O1(i1,j1)),j1=1,N1)
        !end do
        !write(u6,*)
        !write(u6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   O2  k2 = ',k2,'q2 =',q2
        !write(u6,*)
        !do i2=1,N2
        !  write(u6,'(20(2ES14.7,1x))') (Half*((-One)**q2*W2(i2,j2)+O2(i2,j2)),j2=1,N2)
        !end do
        !write(u6,*)
        !write(u6,'(5x,a,i3,3x,A,I3)') 'JKQPAR:   W2  k2 = ',k2,'q2 =',q2
        !write(u6,*)
        !do i2=1,N2
        !  write(u6,'(20(2ES14.7,1x))') (-Half*Onei*((-One)**q2*W2(i2,j2)-O2(i2,j2)),j2=1,N2)
        !end do
#       endif

        ! Build 4 coupled tensor products:
        ! O1-O2
        ! O1-W2
        ! W1-O2
        ! W1-W2
        do i2=1,N2
          do j2=1,N2
            O1_O2(:,:,i2,j2) = O1(:,:)*O2(i2,j2)
            O1_W2(:,:,i2,j2) = O1(:,:)*W2(i2,j2)
            W1_O2(:,:,i2,j2) = W1(:,:)*O2(i2,j2)
            W1_W2(:,:,i2,j2) = W1(:,:)*W2(i2,j2)
          end do
        end do

        ! SP_HZFSO=trace(nDIMcf,HCF,DIP_O)
        ! SP_HZFSW=trace(nDIMcf,HCF,DIP_W)
        ! SP_MOW  =trace(nDIMcf,DIP_O,DIP_W)

        ! B(N,-M)=SP_HZFSO/SP_MOW
        ! B(N, M)=SP_HZFSW/SP_MOW

        FACT = cZero
        FACT = trace(N1,O1,W1)*trace(N2,O2,W2)
        !write(u6,'(A,4I3,2F20.13)') 'k1,q1,k2,q2,FACT=',k1,q1,k2,q2,FACT

        Jpar(k1,-q1,k2,-q2) = trace_exch(N1,N2,HEXCH,O1_O2)/FACT
        Jpar(k1,q1,k2,-q2) = trace_exch(N1,N2,HEXCH,W1_O2)/FACT
        Jpar(k1,-q1,k2,q2) = trace_exch(N1,N2,HEXCH,O1_W2)/FACT
        Jpar(k1,q1,k2,q2) = trace_exch(N1,N2,HEXCH,W1_W2)/FACT

        ! Jpar( k1,-q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, O1_O2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1, q1, k2,-q2 ) = trace_exch(N1,N2,HEXCH, W1_O2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1,-q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, O1_W2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT
        ! Jpar( k1, q1, k2, q2 ) = trace_exch(N1,N2,HEXCH, W1_W2) * Knm(k1,ABS(q1)) * Knm(k2,ABS(q2)) / FACT

#       ifdef _DEBUGPRINT_
        if ((q1 == 0) .and. (q2 == 0)) then
          write(u6,'(A,2ES18.10)') 'FACT = ',FACT
          if (abs(Jpar(k1,-q1,k2,-q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,0,k2,0)
        else if ((q1 == 0) .and. (q2 /= 0)) then
          write(u6,'(A,2ES18.10)') 'FACT = ',FACT
          if (abs(Jpar(k1,0,k2,-q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',-q2,')=',Jpar(k1,0,k2,-q2)
          if (abs(Jpar(k1,0,k2,q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,0,k2,q2)
        else if ((q1 /= 0) .and. (q2 == 0)) then
          write(u6,'(A,2ES18.10)') 'FACT = ',FACT
          if (abs(Jpar(k1,0,k2,-q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',q2,')=',Jpar(k1,-q1,k2,0)
          if (abs(Jpar(k1,0,k2,q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,q1,k2,0)
        else if ((q1 /= 0) .and. (q2 /= 0)) then
          if (abs(Jpar(k1,-q1,k2,-q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',-q2,')=',Jpar(k1,-q1,k2,-q2)
          if (abs(Jpar(k1,q1,k2,-q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',-q2,')=',Jpar(k1,q1,k2,-q2)
          if (abs(Jpar(k1,-q1,k2,q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',-q1,',',k2,',',q2,')=',Jpar(k1,-q1,k2,q2)
          if (abs(Jpar(k1,q1,k2,q2)) > 0.5e-13_wp) &
            write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jpar(',k1,',',q1,',',k2,',',q2,')=',Jpar(k1,q1,k2,q2)
        end if
#       endif

      end do
    end do
  end do
end do

Jcc(:,:,:,:) = cZero
Jcs(:,:,:,:) = cZero
Jsc(:,:,:,:) = cZero
Jss(:,:,:,:) = cZero

do k1=1,N1-1
  do k2=1,N2-1
    do q1=0,k1
      do q2=0,k2

        F1 = (-cOne)**(-q1)
        F2 = (-cOne)**(-q2)
        F12 = (-cOne)**(-q1-q2)

        Jcc(k1,q1,k2,q2) = Jpar(k1,q1,k2,q2)+F2*Jpar(k1,q1,k2,-q2)+F1*Jpar(k1,-q1,k2,q2)+F12*Jpar(k1,-q1,k2,-q2)

        Jss(k1,q1,k2,q2) = Jpar(k1,q1,k2,q2)-F2*Jpar(k1,q1,k2,-q2)-F1*Jpar(k1,-q1,k2,q2)+F12*Jpar(k1,-q1,k2,-q2)

        Jcs(k1,q1,k2,q2) = -Onei*(Jpar(k1,q1,k2,q2)-F2*Jpar(k1,q1,k2,-q2)+F1*Jpar(k1,-q1,k2,q2)-F12*Jpar(k1,-q1,k2,-q2))

        Jsc(k1,q1,k2,q2) = -Onei*(Jpar(k1,q1,k2,q2)+F2*Jpar(k1,q1,k2,-q2)-F1*Jpar(k1,-q1,k2,q2)-F12*Jpar(k1,-q1,k2,-q2))

#       ifdef _DEBUGPRINT_
        if (abs(Jcc(k1,q1,k2,q2)) > 0.5e-13_wp) &
          write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jcc(',k1,',',q1,',',k2,',',q2,')=',Jcc(k1,q1,k2,q2)
        if (abs(Jcs(k1,q1,k2,q2)) > 0.5e-13_wp) &
          write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jcs(',k1,',',q1,',',k2,',',q2,')=',Jcs(k1,q1,k2,q2)
        if (abs(Jsc(k1,q1,k2,q2)) > 0.5e-13_wp) &
          write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jsc(',k1,',',q1,',',k2,',',q2,')=',Jsc(k1,q1,k2,q2)
        if (abs(Jss(k1,q1,k2,q2)) > 0.5e-13_wp) &
          write(u6,'(A,2(i2,A,i3,A),2ES18.10)') 'Jss(',k1,',',q1,',',k2,',',q2,')=',Jss(k1,q1,k2,q2)
#       endif

      end do
    end do
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,'(A,ES18.7)') 'Real Exchange parameters with values larger than: ',THRS
write(u6,'(2A)') repeat('-',119),'|'
write(u6,'(A)') ' k1 |  q1 || k2 |  q2 |-------- O1-O2 --------|-------- O1-W2 --------|-------- W1-O2 --------|'// &
                '-------- W1-W2 --------|'
do k1=1,N1-1
  do q1=-k1,k1
    write(u6,'(A)') '----|-----||----|-----|-----------------------|-----------------------|-----------------------|'// &
                    '-----------------------|'
    do k2=1,N2-1
      do q2=-k2,k2

        F = knm(k1,abs(q1))*knm(k2,abs(q2))
        R = abs(Jcc(k1,abs(q1),k2,abs(q2)))+abs(Jcs(k1,abs(q1),k2,abs(q2)))+abs(Jsc(k1,abs(q1),k2,abs(q2)))+ &
            abs(Jss(k1,abs(q1),k2,abs(q2)))

        if (R > THRS) then

          if ((q1 < 0) .and. (q2 < 0)) then

            write(u6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),4(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                        real(Jcc(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                        real(Jcs(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                        real(Jsc(k1,abs(q1),k2,abs(q2)))*F,'| ', &
                                                                        real(Jss(k1,abs(q1),k2,abs(q2)))*F,'|'

          else if ((q1 >= 0) .and. (q2 < 0)) then

            write(u6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                        real(Jcs(k1,abs(q1),k2,abs(q2)))*F,'|'

          else if ((q1 < 0) .and. (q2 >= 0)) then

            write(u6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                        real(Jsc(k1,abs(q1),k2,abs(q2)))*F,'|'

          else ! (q1 > 0) .and. (q2 > 0)

            write(u6,'( 2((1x,I2,1x,A),(1x,I3,1x,A)),2(ES21.14,1x,A))') k1,'|',q1,'||',k2,'|',q2,'| ', &
                                                                        real(Jcc(k1,abs(q1),k2,abs(q2)))*F,'|'

          end if
        end if

      end do
    end do
  end do
end do
write(u6,'(2A)') repeat('-',119),'|'
#endif

!-----------------------------------------------------------------------
call mma_deallocate(O1)
call mma_deallocate(O2)
call mma_deallocate(W1)
call mma_deallocate(W2)
call mma_deallocate(O1_O2)
call mma_deallocate(O1_W2)
call mma_deallocate(W1_O2)
call mma_deallocate(W1_W2)
call mma_deallocate(Jcc)
call mma_deallocate(Jcs)
call mma_deallocate(Jsc)
call mma_deallocate(Jss)

return

end subroutine JKQPar

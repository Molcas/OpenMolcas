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
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                                                                      C
!  THIS ROUTINE PERFORMS A SPLINE FIT TO A SET OF POINTS (XFIT,YFIT).  C
!  IT INTERPOLATES A SET OF POINTS (XINT,YINT), AND LOCATES ALL        C
!  EXTREME VALUES (XEXT,YEXT) AND DETERMINES THE TYPE (IEXT).          C
!  VARIOUS TYPES OF BOUNDARY CONDITIONS CAN BE CHOOSEN (IBOUND),       C
!  AT PRESENT IBOUND=1 AND 2 ARE IMPLEMENTED.                          C
!                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!IFG IBOUND=2 implemented, but never used or tested, so disabling it
#undef _IBOUND2_

subroutine SPLINE(XFIT,YFIT,NFIT,XINT,YINT,NIN,XEXT,YEXT,IEXT,NEXT,IBOUND)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three
use Definitions, only: wp, iwp, u6

implicit none
! NFIT, NIN ARE INPUT ONLY
integer(kind=iwp), intent(in) :: NFIT, NIN, IBOUND
! NEXT IS INPUT (DIMENSION BOUND) AND OUTPUT (NR OF EXTREMA).
integer(kind=iwp), intent(inout) :: NEXT
real(kind=wp), intent(in) :: XFIT(NFIT), YFIT(NFIT)
real(kind=wp), intent(out) :: XINT(NIN), YINT(NIN), XEXT(NEXT), YEXT(NEXT)
integer(kind=iwp), intent(out) :: IEXT(NEXT)
real(kind=wp) :: VSUM, T0(2), A0, A1, A2, B0, B1, P, Q, S, SE, T, T1, X, Y
real(kind=wp), allocatable :: MatA(:,:), VecB(:), VecC(:), VecD(:), VecH(:)
integer(kind=iwp) :: NOS, NS, I, J, IT, IP
#ifdef _IBOUND2_
real(kind=wp) :: C0, C1, C2, S1, S2
#endif

NOS = NFIT-1
call mma_allocate(MatA,NFIT,NFIT,label='MatA')
call mma_allocate(VecB,NFIT,label='VecB')
call mma_allocate(VecC,NFIT,label='VecC')
call mma_allocate(VecD,NOS,label='VecD')
call mma_allocate(VecH,NOS,label='VecH')
MatA(:,:) = Zero
do I=1,NOS
  VecH(I) = XFIT(I+1)-XFIT(I)
  VecD(I) = (YFIT(I+1)-YFIT(I))/VecH(I)
end do
select case(IBOUND)
  case(1)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  IBOUND=1:                                                       C
    !  BOTH END POINTS ARE EXTRAPOLATED WITH STRAIGHT LINES.           C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    MatA(1,1) = Two*VecH(1)
    MatA(1,2) = VecH(1)
    VecB(1) = Three*VecH(1)*VecD(1)
    MatA(NFIT,NFIT-1) = VecH(NOS)
    MatA(NFIT,NFIT) = Two*VecH(NOS)
    VecB(NFIT) = Three*VecH(NOS)*VecD(NOS)
    do I=2,NOS
      MatA(I,I-1) = VecH(I)
      MatA(I,I) = Two*(VecH(I)+VecH(I-1))
      MatA(I,I+1) = VecH(I-1)
      VecB(I) = Three*(VecH(I-1)*VecD(I)+VecH(I)*VecD(I-1))
    end do
    call dminv(NFIT,NFIT,MatA)
    do I=1,NFIT
      VSUM = Zero
      do J=1,NFIT
        VSUM = VSUM+MatA(I,J)*VecB(J)
      end do
      VecC(I) = VSUM
    end do
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  PERFORM INTERPOLATION (EXTRAPOLATION).                          C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    do IP=1,NIN
      if (XINT(IP) < XFIT(1)) then
        YINT(IP) = YFIT(1)+VecC(1)*(XINT(IP)-XFIT(1))
      else if (XINT(IP) > XFIT(NFIT)) then
        YINT(IP) = YFIT(NFIT)+VecC(NFIT)*(XINT(IP)-XFIT(NFIT))
      else
        NS = 1
        do while (XINT(IP) > XFIT(NS+1))
          NS = NS+1
        end do
        T = (XINT(IP)-XFIT(NS))/VecH(NS)
        T1 = One-T
        YINT(IP) = T*YFIT(NS+1)+T1*YFIT(NS)+VecH(NS)*T*T1*(T1*(VecC(NS)-VecD(NS))-T*(VecC(NS+1)-VecD(NS)))
      end if
    end do
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  LOCATE ALL EXTREMUM POINTS.                                     C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    NEXT = 0
    do NS=1,NOS
      A2 = Three*(VecC(NS)+VecC(NS+1)-Two*VecD(NS))
      A1 = Two*(Three*VecD(NS)-Two*VecC(NS)-VecC(NS+1))
      A0 = VecC(NS)
      B1 = Two*A2
      B0 = A1
      if (abs(A2) > (abs(A1)+abs(A0))*1.0e-3_wp) then
      !***** A2 LARGE *****
        P = -A1/A2/Two
        Q = A0/A2
        S = P*P-Q
        if (abs(S) < 1.0e-10_wp) then
        !***** DOUBLE ROOT *****
          IT = 2
          T = P
          if ((T >= Zero).and.(T <= One)) then
            X = XFIT(NS)+T*VecH(NS)
            Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
            !write(u6,2001) TEXT(IT),X,Y
            NEXT = NEXT+1
            XEXT(NEXT) = X
            YEXT(NEXT) = Y
            IEXT(NEXT) = IT
          end if
        else if (S > Zero) then
        !***** TWO ROOTS *****
          T0(1) = P+sqrt(S)
          T0(2) = P-sqrt(S)
          do I=1,2
            T = T0(I)
            if ((T >= Zero).and.(T <= One)) then
              SE = B0+T*B1
              X = XFIT(NS)+T*VecH(NS)
              Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
              IT = 1
              if (SE > Zero) IT = 3
              !write(u6,2001) TEXT(IT),X,Y
              NEXT = NEXT+1
              XEXT(NEXT) = X
              YEXT(NEXT) = Y
              IEXT(NEXT) = IT
            end if
          end do
        end if
      else if (abs(A1) >= 0.9_wp*abs(A0)) then
      !***** A2 SMALL *****
        T = -A0/A1
        T = -(A0+A2*T*T)/A1
        T = -(A0+A2*T*T)/A1
        T = -(A0+A2*T*T)/A1
        if ((T >= Zero).and.(T <= One)) then
          IT = 1
          if (B0 > Zero) IT = 3
          X = XFIT(NS)+T*VecH(NS)
          Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
          !write(u6,2001) TEXT(IT),X,Y
          NEXT = NEXT+1
          XEXT(NEXT) = X
          YEXT(NEXT) = Y
          IEXT(NEXT) = IT
        end if
      end if
    end do

#ifdef _IBOUND2_
  case(2)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  IBOUND=2:                                                       C
    !  EXTRAPOLATION TO PLUS INFINITY IS DONE WITH A STRAIGHT          C
    !  HORIZONTAL LINE.                                                C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    MatA(NFIT,NFIT) = ONe
    VecB(NFIT) = Zero
    MatA(NOS,NOS) = Two
    MatA(NOS,NFIT) = One
    VecB(NOS) = Three*VecD(NOS)
    do I=2,NOS
      MatA(I-1,I-1) = VecH(I)
      MatA(I-1,I) = Two*(VecH(I)+VecH(I-1))
      MatA(I-1,I+1) = VecH(I-1)
      VecB(I-1) = Three*(VecH(I-1)*VecD(I)+VecH(I)*VecD(I-1))
    end do
    call dminv(NFIT,NFIT,MatA)
    do I=1,NFIT
      VSUM = Zero
      do J=1,NFIT
        VSUM = VSUM+MatA(I,J)*VecB(J)
      end do
      VecC(I) = VSUM
    end do
    S1 = VecC(1)
    S2 = (VecC(1)+Two*VecC(2)-Three*VecD(1))/VecH(1)
    C2 = -S2/S1
    C1 = S2/C2/C2/EXP(-C2*XFIT(1))
    C0 = YFIT(1)-C1*EXP(-C2*XFIT(1))
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  PERFORM INTERPOLATION (EXTRAPOLATION).                          C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    do IP=1,NIN
      if (XINT(IP) < XFIT(1)) then
        YINT(IP) = C0+C1*EXP(-C2*XINT(IP))
      else if (XINT(IP) > XFIT(NFIT)) then
        YINT(IP) = YFIT(NFIT)+VecC(NFIT)*(XINT(IP)-XFIT(NFIT))
      else
        NS = 1
        do while (XINT(IP) > XFIT(NS+1))
          NS = NS+1
        end do
        T = (XINT(IP)-XFIT(NS))/VecH(NS)
        T1 = One-T
        YINT(IP) = T*YFIT(NS+1)+T1*YFIT(NS)+VecH(NS)*T*T1*(T1*(VecC(NS)-VecD(NS))-T*(VecC(NS+1)-VecD(NS)))
      end if
    end do
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !                                                                  C
    !  LOCATE ALL EXTREMUM POINTS.                                     C
    !                                                                  C
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    NEXT = 0
    do NS=1,NOS
      A2 = Three*(VecC(NS)+VecC(NS+1)-Two*VecD(NS))
      A1 = Two*(Three*VecD(NS)-Two*VecC(NS)-VecC(NS+1))
      A0 = VecC(NS)
      B1 = Two*A2
      B0 = A1
      if (abs(A2) > (abs(A1)+abs(A0))*1.0e-3_wp) then
      !***** A2 LARGE *****
        P = -A1/A2/Two
        Q = A0/A2
        S = P*P-Q
        if (abs(S) < 1.0e-10_wp) then
        !***** DOUBLE ROOT *****
          IT = 2
          T = P
          if ((T >= Zero).and.(T <= One)) then
            X = XFIT(NS)+T*VecH(NS)
            Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
            !write(u6,2001) TEXT(IT),X,Y
            NEXT = NEXT+1
            XEXT(NEXT) = X
            YEXT(NEXT) = Y
            IEXT(NEXT) = IT
          end if
        else if (S > Zero) then
        !***** TWO ROOTS *****
          T0(1) = P+sqrt(S)
          T0(2) = P-sqrt(S)
          do I=1,2
            T = T0(I)
            if ((T >= Zero).and.(T <= One)) then
              SE = B0+T*B1
              X = XFIT(NS)+T*VecH(NS)
              Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
              IT = 1
              if (SE > Zero) IT = 3
              !write(u6,2001) TEXT(IT),X,Y
              NEXT = NEXT+1
              XEXT(NEXT) = X
              YEXT(NEXT) = Y
              IEXT(NEXT) = IT
            end if
          end do
        end if
      else if (abs(A1) >= 0.9_wp*abs(A0)) then
      !***** A2 SMALL *****
        T = -A0/A1
        T = -(A0+A2*T*T)/A1
        T = -(A0+A2*T*T)/A1
        T = -(A0+A2*T*T)/A1
        if ((T >= Zero).and.(T <= One)) then
          IT = 1
          if (B0 > Zero) IT = 3
          X = XFIT(NS)+T*VecH(NS)
          Y = YFIT(NS)+VecH(NS)*T*(A0+T*(A1/Two+T*A2/Three))
          !write(u6,2001) TEXT(IT),X,Y
          NEXT = NEXT+1
          XEXT(NEXT) = X
          YEXT(NEXT) = Y
          IEXT(NEXT) = IT
        end if
      end if
    end do
#endif

  case default
    write(u6,*) 'SPLINE Error: IBOUND should be 1 or 2.'
    write(u6,'(1x,a,2i8)') 'But IBOUND=',IBOUND
    call Abend()
end select

call mma_deallocate(MatA)
call mma_deallocate(VecB)
call mma_deallocate(VecC)
call mma_deallocate(VecD)
call mma_deallocate(VecH)

return

!2001 format(1x,a,'  X=',E12.5,' ,  Y=',E12.5)

end subroutine SPLINE

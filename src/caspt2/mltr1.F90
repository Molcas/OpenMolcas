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

subroutine MLTR1(IMLTOP,LST1,X,nX,F,nF,Y,nY)
! Given a lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L1,L2,L3,L4 for short, for a given
! item, and with V=VAL1(L4),
! compute the Rank-1 updates, if IMLTOP=0,
! X(L1,p,q):= Add V*F(L2,p)*Y(L3,q), p=1..LEN1, q=1..LEN2
! else, the conjugate expression
! Y(L3,q):= Add V*F(L2,p)*X(L1,p,q), p=1..LEN1, q=1..LEN2
! or for IMLTOP=2, compute
!     F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
! However, strides etc can vary from case to case: The indices
! may appear in 1st, 2nd or 3rd position, Fortran-address-wise,
! thus an independent stride is to be given for each of the
! indices in the above 'formal' index ordering, as follows:
! The formal Y(p,q,r) is accessed as
! Y(1+INCX1*(p-1)+INCX2*(q-1)+INCX3*(r-1)), etc.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, MyRank, nProcs
#endif
use Sigma_data, only: INCF1, INCF2, INCX1, INCX2, INCX3, INCY1, INCY2, LEN1, LEN2, NFR1, NLST1, VAL1
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, LST1(4,NLST1), nX, nF, nY
real(kind=wp), intent(inout) :: X(nX), F(nF), Y(nY)
integer(kind=iwp) :: I, I_F, ILST, ILST1_IOFF, ILST1_SKIP, IX, IY, L1, L2, L3, L4
real(kind=wp) :: A, V
real(kind=wp), external :: DDot_

!SVC: determine outer loop properties
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  ILST1_IOFF = MYRANK+1
  ILST1_SKIP = NPROCS
else
#endif
  ILST1_IOFF = 1
  ILST1_SKIP = 1
#ifdef _MOLCAS_MPP_
end if
#endif

select case (IMLTOP)
  case (0)
    do ILST=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      I_F = INCF1*(L2-1)+1
      IY = INCY1*(L3-1)+1
      do I=1,LEN1
        A = V*F(I_F)
        call DAXPY_(LEN2,A,Y(IY),INCY2,X(IX),INCX3)
        IX = IX+INCX2
        I_F = I_F+INCF2
      end do
    end do
  case (1)
    do ILST=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      I_F = INCF1*(L2-1)+1
      IY = INCY1*(L3-1)+1
      do I=1,LEN2
        Y(IY) = Y(IY)+V*DDOT_(LEN1,F(I_F),INCF2,X(IX),INCX2)
        IX = IX+INCX3
        IY = IY+INCY2
      end do
    end do
  case default
    do ILST=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST)
      L2 = LST1(2,ILST)
      L3 = LST1(3,ILST)
      L4 = LST1(4,ILST)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      I_F = INCF1*(L2-1)+1
      IY = INCY1*(L3-1)+1
      ! F(L2,p) := Add V*X(L1,p,q)*Y(L3,q)
      do I=1,LEN2
        A = V*Y(IY)
        call DAXPY_(LEN1,A,X(IX),INCX2,F(I_F),INCF2)
        IX = IX+INCX3
        IY = IY+INCY2
      end do
    end do
end select

NFR1 = NFR1+2*NLST1*LEN1*LEN2

end subroutine MLTR1

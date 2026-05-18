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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MLTMV(IMLTOP,LST1,nLST1,X,nX,F,nF,Y,nY)
! Given a lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L1,L2,L3,L4 for short, for a given
! item, and with V=VAL1(L4),
! compute, for IMLTOP=0 or 1, respectively.
!    X(L1,i) := Add V*F(L2,a)*Y(L3,i,a), i=1..LEN1, a=1..LEN2
! or Y(L3,i,a):= Add V*F(L2,a)*X(L1,i), i=1..LEN1, a=1..LEN2
! or for IMLTOP=2, compute
!     F(L2,a) := Add V*X(L1,i)*Y(L3,i,a)
! However, strides etc can vary from case to case: The indices
! may appear in 1st, 2nd or 3rd position, Fortran-address-wise,
! thus an independent stride is to be given for each of the
! indices in the above 'formal' index ordering, as follows:
! The formal X(p,q,r) is accessed as
! X(1+INCX1*(p-1)+INCX2*(q-1)+INCX3*(r-1)), etc.

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs, Is_Real_Par
#endif
use Sigma_data, only: INCF1, INCF2, INCX1, INCY2, INCY3, LEN1, LEN2, NFMV, VAL1, INCY1, INCX2
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, nLST1, nX, nF, nY
real(kind=wp), intent(inout) :: X(nX), F(nF), Y(nY)
integer(kind=iwp), intent(in) :: LST1(4,NLST1)
integer(kind=iwp) ILST1_IOFF, ILST1_SKIP, ILST1, L1, L2, L3, L4, I, IF, IX, IY
real(kind=wp) A, V
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
    do ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST1)
      L2 = LST1(2,ILST1)
      L3 = LST1(3,ILST1)
      L4 = LST1(4,ILST1)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      IF = INCF1*(L2-1)+1
      if ((IF < 1) .or. (IF > nF)) cycle
      IY = INCY1*(L3-1)+1
      ! X(L1,i) := Add V*F(L2,a)*Y(L3,i,a), i=1..LEN1, a=1..LEN2
      do I=1,LEN1
        X(IX) = X(IX)+V*DDOT_(LEN2,F(IF),INCF2,Y(IY),INCY3)
        IX = IX+INCX2
        IY = IY+INCY2
      end do
    end do
  case (1)
    do ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST1)
      L2 = LST1(2,ILST1)
      L3 = LST1(3,ILST1)
      L4 = LST1(4,ILST1)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      IF = INCF1*(L2-1)+1
      if ((IF < 1) .or. (IF > nF)) cycle
      IY = INCY1*(L3-1)+1
      ! or Y(L3,i,a):= Add V*F(L2,a)*X(L1,i), i=1..LEN1, a=1..LEN2
      do I=1,LEN2
        A = V*F(IF)
        call DAXPY_(LEN1,A,X(IX),INCX2,Y(IY),INCY2)
        IY = IY+INCY3
        IF = if+INCF2
      end do
    end do
  case default
    do ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
      L1 = LST1(1,ILST1)
      L2 = LST1(2,ILST1)
      L3 = LST1(3,ILST1)
      L4 = LST1(4,ILST1)
      V = VAL1(L4)
      IX = INCX1*(L1-1)+1
      IF = INCF1*(L2-1)+1
      if ((If < 1) .or. (If > nF)) cycle
      IY = INCY1*(L3-1)+1
      ! F(L2,a) := Add V*X(L1,i)*Y(L3,i,a)
      do I=1,LEN1
        A = V*X(IX)
        call DAXPY_(LEN2,A,Y(IY),INCY3,F(If),INCF2)
        IX = IX+INCX2
        IY = IY+INCY2
      end do
    end do
end select

NFMV = NFMV+2*NLST1*LEN1*LEN2

end subroutine MLTMV

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

subroutine MLTDXP(IMLTOP,LST1,LST2,X,nX,F,nF,Y,nY)
! Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L11,L12,L13,L14 for short, for a given
! item, and with V1=VAL1(L14), and similar for the other list,
! compute, for IMLTOP=0 or 1 resp.,
!     X(L11,L21,a) := Add V1*V2*F(L12,L22)*Y(L13,L23,a)
! or  Y(L13,L23,a) := Add V1*V2*F(L12,L22)*X(L11,L21,a), a=1,LEN1
! or for IMLTOP=2, compute
!     F(L12,L22) := Add V1*V2*X(L11,L21,a)*Y(L13,L23,a)
! However, strides etc can vary from case to case: The indices
! may appear in 1st, 2nd or 3rd position, Fortran-address-wise,
! thus an independent stride is to be given for each of the
! indices in the above 'formal' index ordering, as follows:
! The formal X(p,q,r) is accessed as
! X(1+INCX1*(p-1)+INCX2*(q-1)+INCX3*(r-1)), etc.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, MyRank, nProcs
#endif
use Sigma_data, only: INCF1, INCF2, INCX1, INCX2, INCX3, INCY1, INCY2, INCY3, LEN1, NFDXP, NLST1, NLST2, VAL1, VAL2
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, LST1(4,NLST1), LST2(4,NLST2), nX, nF, nY
real(kind=wp), intent(inout) :: X(nX), F(nF), Y(nY)
integer(kind=iwp) :: I_F, ILST1, ILST1_IOFF, ILST1_SKIP, ILST2, IX, IY, L11, L12, L13, L14, L21, L22, L23, L24
real(kind=wp) :: A, V, V1, V2
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
      L11 = LST1(1,ILST1)
      L12 = LST1(2,ILST1)
      L13 = LST1(3,ILST1)
      L14 = LST1(4,ILST1)
      V1 = VAL1(L14)
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        IX = 1+INCX1*(L11-1)+INCX2*(L21-1)
        IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
        I_F = 1+INCF1*(L12-1)+INCF2*(L22-1)
        A = V1*V2*F(I_F)
        call DAXPY_(LEN1,A,Y(IY),INCY3,X(IX),INCX3)
      end do
    end do
  case (1)
    do ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
      L11 = LST1(1,ILST1)
      L12 = LST1(2,ILST1)
      L13 = LST1(3,ILST1)
      L14 = LST1(4,ILST1)
      V1 = VAL1(L14)
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        IX = 1+INCX1*(L11-1)+INCX2*(L21-1)
        IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
        I_F = 1+INCF1*(L12-1)+INCF2*(L22-1)
        A = V1*V2*F(I_F)
        call DAXPY_(LEN1,A,X(IX),INCX3,Y(IY),INCY3)
      end do
    end do
  case default
    do ILST1=ILST1_IOFF,NLST1,ILST1_SKIP
      L11 = LST1(1,ILST1)
      L12 = LST1(2,ILST1)
      L13 = LST1(3,ILST1)
      L14 = LST1(4,ILST1)
      V1 = VAL1(L14)
      do ILST2=1,NLST2
        L21 = LST2(1,ILST2)
        L22 = LST2(2,ILST2)
        L23 = LST2(3,ILST2)
        L24 = LST2(4,ILST2)
        V2 = VAL2(L24)
        IX = 1+INCX1*(L11-1)+INCX2*(L21-1)
        IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
        I_F = 1+INCF1*(L12-1)+INCF2*(L22-1)
        A = V1*V2*F(I_F)
        V = V1*V2
        F(I_F) = F(I_F)+V*DDOT_(LEN1,X(IX),INCX3,Y(IY),INCY3)
      end do
    end do
end select
NFDXP = NFDXP+2*NLST1*NLST2*LEN1

end subroutine MLTDXP

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

subroutine MLTSCA(IMLTOP,LST1,NLST1,LST2,NLST2,X,nX,F,nF,Y,nY)
! Given two lists with entries LST1(4,ITEM), ITEM=1,NLST1, the
! four entries called L11,L12,L13,L14 for short, for a given
! item, and with V1=VAL1(L14), and similar for the other list,
! compute, for IMLTOP=0 or 1 respectively,
!     X(L11,L21) := Add V1*V2*F(L12,L22)*Y(L13,L23)
!  or Y(L13,L23) := Add V1*V2*F(L12,L22)*X(L11,L21)
! or for IMLTOP=2, compute
!     F(L12,L22) := Add V1*V2*X(L11,L21)*Y(L13,L23)

#ifdef _MOLCAS_MPP_
use Para_Info, only: MyRank, nProcs, Is_Real_Par
#endif
use Sigma_data, only: INCF1, INCF2, INCX1, INCX2, INCY1, INCY2, NFSCA, VAL1, VAL2
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: IMLTOP, NLST1, NLST2, nX, nF, nY
real(kind=wp), intent(inout) :: X(nX), F(nF), Y(nY)
integer(kind=iwp), intent(in) :: LST1(4,NLST1), LST2(4,NLST2)
integer(kind=iwp) if, ILST1, ILST1_IOFF, ILST1_SKIP, ILST2, IX, IY, L11, L12, L13, L14, L21, L22, L23, L24
real(kind=wp) V1, V2

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

if (IMLTOP == 0) then
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
      if = 1+INCF1*(L12-1)+INCF2*(L22-1)
      IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
      X(IX) = X(IX)+V1*V2*F(if)*Y(IY)
    end do
  end do
else if (IMLTOP == 1) then
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
      if = 1+INCF1*(L12-1)+INCF2*(L22-1)
      IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
      Y(IY) = Y(IY)+V1*V2*F(if)*X(IX)
    end do
  end do
else
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
      if = 1+INCF1*(L12-1)+INCF2*(L22-1)
      IY = 1+INCY1*(L13-1)+INCY2*(L23-1)
      F(if) = F(if)+V1*V2*X(IX)*Y(IY)
    end do
  end do
end if

NFSCA = NFSCA+4*NLST1*NLST2

end subroutine MLTSCA

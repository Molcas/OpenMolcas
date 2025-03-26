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

!#define _DEBUGPRINT_
subroutine GETINT_MCLR(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL,ieaw)
! Outer routine for accessing integral block

use Arrays, only: pInt2, KINT2, KINT2A
use MCLR_Data, only: Square
#ifdef _DEBUGPRINT_
use MCLR_Data, only: NOBPTS
#endif
use input_mclr, only: nsMOB

implicit none
real*8 XINT(*)
integer ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, ICOUL, ieaw
#ifdef _DEBUGPRINT_
integer nI, nK, nIK, nJ, nL, nJL, nIJ, nKL
#endif

if (.not. square) then
  if (ieaw /= 0) then
    call GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,KINT2a,pINT2,NSMOB,ICOUL,ieaw)
  else
    call GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,KINT2,pINT2,NSMOB,ICOUL,ieaw)
  end if
else
  call GETINC_ABS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,KINT2,pINT2,NSMOB,ICOUL)
end if

#ifdef _DEBUGPRINT_
NI = NOBPTS(ITP,ISM)
NK = NOBPTS(KTP,KSM)
if (IKSM == 0) then
  NIK = NI*NK
else
  NIK = NI*(NI+1)/2
end if
NJ = NOBPTS(JTP,JSM)
NL = NOBPTS(LTP,LSM)
if (JLSM == 0) then
  NJL = NJ*NL
else
  NJL = NJ*(NJ+1)/2
end if
if (ICOUL == 0) then
  write(6,*) ' 2 electron integral block for TS blocks'
  write(6,*) ' Ixchng :',IXCHNG
  write(6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
  call WRTMAT(XINT,NIK,NJL,NIK,NJL)
else
  write(6,*) ' Integrals in Coulomb form'
  write(6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
  NIJ = NI*NJ
  NKL = NK*NL
  call WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
end if
#endif

end subroutine GETINT_MCLR

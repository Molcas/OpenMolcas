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
subroutine GETINT_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IKSM,JLSM,ICTL,ieaw)
! Outer routine for accessing integral block

use MCLR_Data, only: KINT2, KINT2a, Square
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Index_Functions, only: nTri_Elem
use MCLR_Data, only: NOBPTS
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: XINT(*)
integer(kind=iwp), intent(in) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IKSM, JLSM, ICTL, ieaw
integer(kind=iwp) :: iXChng, iCoul
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nI, nK, nIK, nJ, nL, nJL, nIJ, nKL
#endif

!write(u6,*) 'square in getint_td',square
if (.not. square) then

  IXCHNG = 0
  ICOUL = 0
  if (ictl == 2) ICOUL = 1
  if (ictl == 3) ICOUL = 1
  if (ictl == 4) IXCHNG = 1
  if (ieaw /= 0) then
    call GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,KINT2a,ICOUL,ieaw)
  else
    call GETINC_ABT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,KINT2,ICOUL,ieaw)
  end if
else
  call GETINC_ABS_td(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IKSM,JLSM,KINT2,ICTL)

end if

#ifdef _DEBUGPRINT_
NI = NOBPTS(ITP,ISM)
NK = NOBPTS(KTP,KSM)
if (IKSM == 0) then
  NIK = NI*NK
else
  NIK = nTri_Elem(NI)
end if
NJ = NOBPTS(JTP,JSM)
NL = NOBPTS(LTP,LSM)
if (JLSM == 0) then
  NJL = NJ*NL
else
  NJL = nTri_Elem(NJ)
end if
if (ICOUL == 0) then
  write(u6,*) ' 2 electron integral block for TS blocks'
  write(u6,*) ' Ixchng :',IXCHNG
  write(u6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
  call WRTMAT(XINT,NIK,NJL,NIK,NJL)
else
  write(u6,*) ' Integrals in Coulomb form'
  write(u6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
  NIJ = NI*NJ
  NKL = NK*NL
  call WRTMAT(XINT,NIJ,NKL,NIJ,NKL)
end if
#endif

!stop ' Jeppe forced me to stop in GETINT'

end subroutine GETINT_td

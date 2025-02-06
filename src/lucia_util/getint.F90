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

subroutine GETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
! Outer routine for accessing integral block

use lucia_data, only: NOBPTS, NTOOBS
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: XINT(*)
integer(kind=iwp) :: ITP, ISM, JTP, JSM, KTP, KSM, LTP, LSM, IXCHNG, IKSM, JLSM, ICOUL
integer(kind=iwp) :: NI, NIK, NJ, NJL, NK, NL, NTEST

NTEST = 0

if (NTEST >= 1) then
  !write(u6,*) ' I_USE_SIMTRH in GETINT =',I_USE_SIMTRH
  write(u6,*) ' GETINT : ICOUL = ',ICOUL
  write(u6,*) 'ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM :'
  write(u6,'(8I4)') ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM
end if
! Read integrals in in RASSCF format
call GETINCN_RASSCFS(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,IXCHNG,IKSM,JLSM,ICOUL)

if (NTEST /= 0) then
  if (ITP == 0) then
    NI = NTOOBS(ISM)
  else
    NI = NOBPTS(ITP,ISM)
  end if

  if (KTP == 0) then
    NK = NTOOBS(KSM)
  else
    NK = NOBPTS(KTP,KSM)
  end if

  if (IKSM == 0) then
    NIK = NI*NK
  else
    NIK = NI*(NI+1)/2
  end if

  if (JTP == 0) then
    NJ = NTOOBS(JSM)
  else
    NJ = NOBPTS(JTP,JSM)
  end if

  if (LTP == 0) then
    NL = NTOOBS(LSM)
  else
    NL = NOBPTS(LTP,LSM)
  end if

  if (JLSM == 0) then
    NJL = NJ*NL
  else
    NJL = NJ*(NJ+1)/2
  end if
  write(u6,*) ' 2 electron integral block for TS blocks'
  write(u6,*) ' Ixchng :',IXCHNG
  write(u6,*) ' After GETINC'
  write(u6,'(1X,4(A,I2,A,I2,A))') '(',ITP,',',ISM,')','(',JTP,',',JSM,')','(',KTP,',',KSM,')','(',LTP,',',LSM,')'
  call WRTMAT(XINT,NIK,NJL,NIK,NJL)
end if

end subroutine GETINT

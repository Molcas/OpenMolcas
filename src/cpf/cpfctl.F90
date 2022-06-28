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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine CPFCTL(C,S,W,TPQ,ENP,EPP,BST,EPB,AP,BIJ,CN,TEMP)

use cpf_global, only: DETOT, ETHRE, ETOT, ICASE, ICONV, ICPF, IDIIS, INCPF, INDX, IPRINT, IRC, ISDCI, ITER, ITPUL, JSY, MAXIT, &
                      MAXITP, NTMAX, POTNUC
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: C(*), S(*), EPP(*), BST((MAXIT+1)**2)
real(kind=wp), intent(in) :: W(*), ENP(*)
real(kind=wp), intent(_OUT_) :: TPQ(*), EPB(*), AP(*)
real(kind=wp), intent(out) :: BIJ((MAXIT+1)**2), CN(MAXIT+1), TEMP(NTMAX)
integer(kind=iwp) :: I, IP, ITP
real(kind=wp) :: C0, DECORR, DELE, EENP, ENER, ETOTT

call EPSBIS(JSY,INDX,C,W,EPB)
call EPSPRIM(JSY,INDX,C,S,EPP)
IP = IRC(4)
call VECSUM_CPFMCPF(EPP,ENER,IP)
ETOTT = ENER+POTNUC
DELE = ETOTT-ETOT
ETOT = ETOTT
if (ITER == 1) then
  write(u6,'(1X,A)') ' ITER      TOTAL ENERGY          CORR ENERGY           DECREASE'
end if
write(u6,'(1X,I3,3(5X,F16.8))') ITER,ETOT,ENER,DELE
if ((abs(DELE) < ETHRE) .and. (ITPUL /= 1)) ICONV = 1
if ((ICONV == 0) .and. (ITER /= MAXIT)) then
  ! If more iterations should be done.
  IDIIS = 0
  if (ITPUL == MAXITP) IDIIS = 1
  call APPRIM(EPP,EPB,TPQ,AP,ENP,TEMP,ICASE)
  call CUPDATE(JSY,INDX,C,S,AP,BST,ENP)
  ITP = ITPUL+1
  call DIIS_CPF(C,S,BST,MAXIT,BIJ,ITP,CN)
else
  if (ICONV == 1) write(u6,37)
  if (ICONV == 0) write(u6,38)
  if (ISDCI == 1) write(u6,30) ETOT
  if (ICPF == 1) write(u6,35) ETOT
  if (INCPF == 1) write(u6,39) ETOT
  if ((ISDCI == 0) .and. (ICPF == 0) .and. (INCPF == 0)) write(u6,36) ETOT
  write(u6,31) ENER,POTNUC
  if (ISDCI == 1) call Add_Info('E_SDCI',[ETOT],1,8)
  if (ICPF == 1) call Add_Info('E_CPF',[ETOT],1,8)
  if (INCPF == 1) call Add_Info('E_ACPF',[ETOT],1,8)
  if ((ISDCI == 0) .and. (ICPF == 0) .and. (INCPF == 0)) call Add_Info('E_MCPF',[ETOT],1,8)
  if (ISDCI /= 0) then
    EENP = ENP(IRC(4))
    C0 = One/sqrt(EENP)
    DECORR = ENER*(EENP-One)
    DETOT = ETOT+DECORR
    write(u6,32) DETOT
    write(u6,33) DECORR,C0
  end if

  if (IPRINT > 5) then
    if (IPRINT > 5) write(u6,34) (ENP(I),I=1,IRC(4))
  end if
end if

return

30 format(/,5X,'FINAL CI ENERGY',6X,F17.8)
31 format(5X,'FINAL CORRELATION ENERGY',F14.8,'  REFERENCE ENERGY',F17.8)
32 format(5X,'DAVIDSON CORR. ENERGY',F17.8)
33 format(5X,'DAVIDSON CORRECTION',F19.8,'  C0 = ',F12.6)
34 format(/,(5X,'ENP',5F10.6))
35 format(/,5X,'FINAL CPF ENERGY',5X,F17.8)
36 format(/,5X,'FINAL MCPF ENERGY',5X,F17.8)
37 format(/,5X,'CALCULATION CONVERGED')
38 format(/,5X,'CALCULATION NOT COMPLETELY CONVERGED')
39 format(/,5X,'FINAL ACPF ENERGY',4X,F17.8)

end subroutine CPFCTL

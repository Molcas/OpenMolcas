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
! Copyright (C) 2011, Jeppe Olsen                                      *
!               2011, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CHECK_BLOCKS_FOR_BK_APPROX(IATP,IBTP,JATP,JBTP,IASM,IBSM,JASM,JBSM,IOCTPA,IOCTPB,I_DO_EXACT_BLOCK)
! Check whether block <IATP, IBTP! H! JATP, JBTP> should
! be calculated exactly or by BK approx
! Input
! ======
! IATP IBTP JATP JBTP : Supergroups (relative numbers)
! IOCTPA, IOBTPB : Offset for type
! Output
! ======
! I_DO_EXACT_BLOCK = 1 => Do exact block
!                  = 0 => Set block to zero
!                  =-1 => Use diagonal aproximation
! Giovanni +Jeppe Olsen, Sept 2011, on a bench at Torre Normanna, Sicily

use lucia_data, only: IBSPGPFTP, MXPNGAS, NELFSPGP, NGAS
use spinfo, only: IOCCPSPC, NGASBK
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: IATP, IBTP, JATP, JBTP, IASM, IBSM, JASM, JBSM
integer(kind=iwp), intent(out) :: IOCTPA, IOCTPB, I_DO_EXACT_BLOCK
integer(kind=iwp) :: ICHECK_OCC_IN_ACCSPC, IOCC(MXPNGAS), IOCC_IN, JOCC(MXPNGAS), JOCC_IN
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IGAS
#endif

IOCTPA = IBSPGPFTP(1)
IOCTPB = IBSPGPFTP(2)
#ifdef _DEBUGPRINT_
write(u6,*) 'IOCTPA, IOCTPB',IOCTPA,IOCTPB
write(u6,*) 'IATP,IBTP,JATP,JBTP'
write(u6,'(5x,I4,5x,I4,5x,I4,5X,I4)') IATP,IBTP,JATP,JBTP
write(u6,*) 'NELFSPGP (IA), (IB), (JA), (JB)'
write(u6,*) (NELFSPGP(igas,IOCTPA-1+IATP),igas=1,NGAS)
write(u6,*) (NELFSPGP(igas,IOCTPB-1+IBTP),igas=1,NGAS)
write(u6,*) (NELFSPGP(igas,IOCTPA-1+JATP),igas=1,NGAS)
write(u6,*) (NELFSPGP(igas,IOCTPB-1+JBTP),igas=1,NGAS)
#endif
IOCC(1:NGAS) = NELFSPGP(1:NGAS,IOCTPA-1+IATP)+NELFSPGP(1:NGAS,IOCTPB-1+IBTP)
JOCC(1:NGAS) = NELFSPGP(1:NGAS,IOCTPA-1+JATP)+NELFSPGP(1:NGAS,IOCTPB-1+JBTP)
#ifdef _DEBUGPRINT_
write(u6,*) ' Routine CHECK_BLOCKS_FOR_BK_APPROX is speaking!'
write(u6,*) ' I am doing BK-type of approximation'
write(u6,*) ' Min and Max for subspace with exact Hamiltonian'
write(u6,*) ' ==============================================='
write(u6,*) 'NGASBK : ',NGASBK
write(u6,*) '              Min. Occ.      Max. Occ.'
do IGAS=1,NGASBK
  write(u6,'(A,I2,10X,I3,9X,I3)') '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
end do
#endif

IOCC_IN = ICHECK_OCC_IN_ACCSPC(IOCC,IOCCPSPC,NGASBK,20)
JOCC_IN = ICHECK_OCC_IN_ACCSPC(JOCC,IOCCPSPC,NGASBK,20)
! = 1 if the Occupation is IN
! = 0 of the Occupation is OUT
! If both occupation are outside of IOCCPSPC, we make approximations
if ((IOCC_IN == 0) .and. (JOCC_IN == 0)) then
  ! If the blocks are identical use diagonal approximation
  if ((IATP == JATP) .and. (IASM == JASM) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) then
    ! Compute as diagonal
    I_DO_EXACT_BLOCK = -1
  else
    ! Neglect
    I_DO_EXACT_BLOCK = 0
  end if
else
  ! At least one block is in PSPC, so calculate exactly
  I_DO_EXACT_BLOCK = 1
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' CHECK_BLOCKS_FOR_BK_APPROX is speaking'
write(u6,'(A, 4I4)') ' Input blocks IA, IB, JA, JB = ',IATP,IBTP,JATP,JBTP
write(u6,'(A, 4I4)') ' Input blocks IASM, IBSM, JASM, JBSM = ',IASM,IBSM,JASM,JBSM
write(u6,'(A,I4)') ' I_DO_EXACT_BLOCK = ',I_DO_EXACT_BLOCK
#endif

end subroutine CHECK_BLOCKS_FOR_BK_APPROX

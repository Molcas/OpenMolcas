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
! Copyright (C) 2015, Roland Lindh                                     *
!***********************************************************************

subroutine MK_PROP(PROP,IPROP,ISTATE_,JSTATE_,LABEL,ITYPE,BUFF,NBUFF,DENS,NDENS,MASK,ISY12,IOFF)
!***********************************************************************
!     Objective: to compute the transition property between state      *
!                ISTATE and JSTATE of property IPROP.                  *
!                                                                      *
!     This routine will be generalized to a direct routine later.      *
!                                                                      *
!     Author: Roland Lindh, Uppsala University, 23 Dec. 2015           *
!***********************************************************************

use OneDat, only: sOpSiz
use Cntrl, only: ICOMP, IPUSED, NPROP, NSTATE, PNAME, PNUC, PORIG
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IPROP, ISTATE_, JSTATE_, ITYPE, NBUFF, NDENS, MASK, ISY12, IOFF(8)
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), BUFF(NBUFF), DENS(NDENS,4)
character(len=8) :: LABEL
integer(kind=iwp) :: I12, IC, IDUM(1), IINT, IOPT, IPOS, IRC, ISCHK, ISTATE, ISY1, ISY2, JSTATE, NB1, NB12, NB2, NSIZ
real(kind=wp) :: PSUM
real(kind=wp), external :: DDot_

ISTATE = max(ISTATE_,JSTATE_)
JSTATE = min(ISTATE_,JSTATE_)

if (LABEL(1:4) == 'ASD ') LABEL(1:5) = 'MAGXP'
IC = ICOMP(IPROP)
!write(u6,*) 'Mk_Prop: Label=',Label
!write(u6,*) 'Mk_Prop:    IC=',IC
IOPT = ibset(0,sOpSiz)
NSIZ = 0
call iRDONE(IRC,IOPT,LABEL,IC,IDUM,ISCHK)
if (IRC == 0) NSIZ = IDUM(1)
if (mod(ISCHK/MASK,2) == 0) return
IOPT = 0
! Rulin: The 'spin-dependent' part of hyperfine contribution
if (LABEL(1:5) == 'MAGXP') then
  call HFCSD(LABEL,IC,BUFF,NBUFF,NSIZ,ISCHK)
  LABEL(1:5) = 'ASD  '
else
  call RDONE(IRC,IOPT,LABEL,IC,BUFF,ISCHK)
end if
!write(u6,*) 'NBUFF,NSIZ=',NBUFF,NSIZ
if ((IRC /= 0) .and. (LABEL(1:4) /= 'TMOM')) then
  write(u6,*)
  write(u6,'(6X,A)') '*** ERROR IN SUBROUTINE MK_PROP ***'
  write(u6,'(6X,A)') '  FAILED IN READING FROM  ONEINT'
  write(u6,'(6X,A,A)') '  LABEL     = ',LABEL
  write(u6,'(6X,A,I2)') '  COMPONENT = ',IC
  write(u6,*)
  return
end if
IPUSED(IPROP) = 1
! IF THIS IS THE FIRST CALL TO THE SUBROUTINE, PICK UP SOME DATA:
!if (ICALL == 0) then
!-SVC: for safety reasons, always pick up the data, since it is not
!      necessarily done on the first call!
! PICK UP THE ORIGIN COORDINATES:
PORIG(1,IPROP) = BUFF(NSIZ+1)
PORIG(2,IPROP) = BUFF(NSIZ+2)
PORIG(3,IPROP) = BUFF(NSIZ+3)
! PICK UP THE NUCLEAR CONTRIBUTION FROM INTEGRAL BUFFER
if ((PNAME(IPROP)(1:3) == 'ASD') .or. (PNAME(IPROP)(1:3) == 'PSO')) then
  write(u6,*) 'Removing nuclear contrib from ASD and PSO:'
  PNUC(IPROP) = Zero
else if ((ITYPE == 2) .or. (ITYPE == 4)) then
  PNUC(IPROP) = Zero
else
  PNUC(IPROP) = BUFF(NSIZ+4)
end if
IINT = 1
PSUM = Zero
do ISY1=1,nIrrep
  NB1 = NBASF(ISY1)
  if (NB1 == 0) cycle
  do ISY2=1,ISY1
    I12 = MUL(ISY1,ISY2)
    if (iand(2**(I12-1),ISCHK) == 0) cycle
    NB2 = NBASF(ISY2)
    if (NB2 == 0) cycle
    NB12 = NB1*NB2
    if (ISY1 == ISY2) NB12 = (NB12+NB1)/2
    if (I12 == ISY12) then
      IPOS = IOFF(ISY1)+1
      PSUM = PSUM+DDOT_(NB12,BUFF(IINT),1,DENS(IPOS,ITYPE),1)
    end if
    IINT = IINT+NB12
  end do
end do
!write(u6,*) 'PSUM=',PSUM,LABEL,IC
! IN THE CASE OF MULTIPOLES, CHANGE SIGN TO ACCOUNT FOR THE NEGATIVE
! ELECTRONIC CHARGE AS COMPARED TO THE NUCLEAR CONTRIBUTION.
if (LABEL(1:5) == 'MLTPL') PSUM = -PSUM
! In the case of AMFI integrals, they should be multiplied by 2.
! The reason for the factor of two is that this program uses spin
! (uses Clebsch-Gordan coefficients to define Wigner-Eckart
!  reduced matrix elements of spin-tensor properties)
! while the AMFI authors used Pauli matrices.
if (LABEL(1:4) == 'AMFI') PSUM = Two*PSUM
PROP(ISTATE,JSTATE,IPROP) = PSUM
if ((ITYPE == 1) .or. (ITYPE == 3)) then
  PROP(JSTATE,ISTATE,IPROP) = PSUM
else
  PROP(JSTATE,ISTATE,IPROP) = -PSUM
end if

end subroutine MK_PROP

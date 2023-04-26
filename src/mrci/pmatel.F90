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

subroutine PMATEL(ISTATE,JSTATE,PROP,PINT,SMAT,CNO,OCC,SFOLD,AFOLD,TDAO)

use mrci_global, only: BNAME, IPCOMP, NBAS, NBAST, NBTRI, NCMO, NPROP, NRROOT, NSYM, PNAME, PNUC, PORIG, PTYPE
use Symmetry_Info, only: Mul
use OneDat, only: sNoNuc, sNoOri, sOpSiz
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ISTATE, JSTATE
real(kind=wp), intent(_OUT_) :: PROP(NRROOT,NRROOT,NPROP), SMAT(*)
real(kind=wp), intent(out) :: PINT(NBTRI+4), SFOLD(NBTRI), AFOLD(NBTRI)
real(kind=wp), intent(in) :: CNO(NCMO), OCC(NBAST), TDAO(NBAST,NBAST)
integer(kind=iwp) :: I, ICALL = 0, ICOMP, IDUM(1), IDUMMY, IEND, IFROM, IJ, IOPT, IPROP, IRTC, ISTA, ISY, ISY1, ISY12, ISY2, &
                     ISYMLB, ITO, J, NB1, NB12, NB2, NSIZ
real(kind=wp) :: SGN, X
character(len=8) :: LABEL
real(kind=wp), external :: DDOT_

if (ISTATE == JSTATE) then
  ! READ OVERLAP INTEGRALS FROM TRAONE.
  IOPT = ibset(ibset(0,sNoOri),sNoNuc)
  LABEL = 'MLTPL  0'
  ICOMP = 1
  call RDONE(IRTC,IOPT,LABEL,ICOMP,SMAT,IDUMMY)
  ! CALCULATE AND WRITE MULLIKEN CHARGES.
  write(u6,*)
  write(u6,'(A,I2)') ' MULLIKEN CHARGES FOR STATE NR ',ISTATE
  call CHARGE(NSYM,NBAS,BNAME,CNO,OCC,SMAT,2,.true.,.true.)
  write(u6,*) ' ',repeat('*',70)
end if
! FOLD TDAO SYMMETRICALLY (ANTI-SYMM) INTO SFOLD (AFOLD):
! MOLCAS2 UPDATE: SYMMETRY-BLOCKED STORAGE.
IJ = 0
IEND = 0
do ISY=1,NSYM
  ISTA = IEND+1
  IEND = IEND+NBAS(ISY)
  do I=ISTA,IEND
    do J=ISTA,I-1
      IJ = IJ+1
      SFOLD(IJ) = TDAO(I,J)+TDAO(J,I)
      AFOLD(IJ) = TDAO(I,J)-TDAO(J,I)
    end do
    IJ = IJ+1
    SFOLD(IJ) = TDAO(I,I)
    AFOLD(IJ) = Zero
  end do
end do
NSIZ = 0
IOPT = 0
do IPROP=1,NPROP
  ! PICK UP MATRIX ELEMENTS FROM ONE-ELECTRON FILE:
  call iRDONE(IRTC,ibset(IOPT,sOpSiz),PNAME(IPROP),IPCOMP(IPROP),IDUM,ISYMLB)
  if (IRTC == 0) NSIZ = IDUM(1)
  call RDONE(IRTC,IOPT,PNAME(IPROP),IPCOMP(IPROP),PINT,ISYMLB)
  ! SEPARATE OUT THE OPERATOR GAUGE ORIGIN, AND NUCLEAR CONTRIBUTION:
  if (ICALL == 0) then
    PORIG(1,IPROP) = PINT(NSIZ+1)
    PORIG(2,IPROP) = PINT(NSIZ+2)
    PORIG(3,IPROP) = PINT(NSIZ+3)
    PNUC(IPROP) = PINT(NSIZ+4)
  end if
  if (ISYMLB /= 1) then
    ! NON-DIAGONAL SYMMETRY BLOCKS MUST BE COMPRESSED AWAY:
    IFROM = 1
    ITO = 1
    do ISY1=1,NSYM
      NB1 = NBAS(ISY1)
      if (NB1 == 0) cycle
      do ISY2=1,ISY1
        NB2 = NBAS(ISY2)
        if (NB2 == 0) cycle
        ISY12 = MUL(ISY1,ISY2)
        if (mod(ISYMLB,2**(ISY12)) == 0) cycle
        NB12 = NB1*NB2
        if (ISY12 == 1) then
          NB12 = (NB12+NB1)/2
          if (IFROM > ITO) call DCOPY_(NB12,PINT(IFROM),1,PINT(ITO),1)
          ITO = ITO+NB12
        end if
        IFROM = IFROM+NB12
      end do
    end do
    NSIZ = ITO
  end if
  ! PUT DDOT OF TR DENS MATRIX AND INTEGRALS INTO PROPER MATRIX ELEMENT
  ! FOR MULTIPOLES, USE NEGATIVE SIGN OF ELECTRONIC PART.
  SGN = One
  if (PNAME(IPROP)(1:5) == 'MLTPL') SGN = -SGN
  if (PTYPE(IPROP) == 'HERM') then
    X = SGN*DDOT_(NSIZ,SFOLD,1,PINT,1)
    PROP(ISTATE,JSTATE,IPROP) = X
    PROP(JSTATE,ISTATE,IPROP) = X
  else
    X = SGN*DDOT_(NSIZ,AFOLD,1,PINT,1)
    PROP(ISTATE,JSTATE,IPROP) = X
    PROP(JSTATE,ISTATE,IPROP) = -X
  end if
end do
ICALL = 1

return

end subroutine PMATEL

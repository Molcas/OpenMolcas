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

subroutine NATORB_CPF(D,CM,CMO,DSYM,CAO,OCC,M)

use cpf_global, only: ICH, IPRINT, IROW, NASH, NBAS, NFRO, NISH, NORB, NVIR
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: D(*), CM(*)
real(kind=wp), intent(_OUT_) :: CMO(*), DSYM(*), CAO(*), OCC(*)
integer(kind=iwp), intent(in) :: M
integer(kind=iwp) :: I, II, IJ, IJ0, ILAS, IP, ISTA, J, JP, kp, M1, NBM, NBP, NFR, NI, NIJ, NJ, NORBM, NORBM2
real(kind=wp) :: TERM

! FORM DENSITY MATRIX BY SYMMETRY
NBM = NBAS(M)
if (NBM == 0) return
if (NORB(M) == 0) return
if (IPRINT >= 15) write(u6,1) M
NBP = 0
M1 = M-1
do I=1,M1
  NBP = NBP+NORB(I)
end do
NFR = NFRO(M)
NORBM = NFR+NISH(M)+NASH(M)+NVIR(M)
if (NORBM == 0) return
NORBM2 = IROW(NORBM+1)
DSYM(1:NORBM2) = Zero
II = 0
do I=1,NFR
  II = II+I
  DSYM(II) = Two
end do
! REST OF DENSITY MATRIX
IJ = 0
do I=1,NORBM
  do J=1,I
    IJ = IJ+1
    NI = ICH(NBP+I)
    if (NI < 0) cycle
    NJ = ICH(NBP+J)
    if (NJ < 0) cycle
    NIJ = IROW(NI)+NJ
    if (NJ > NI) NIJ = IROW(NJ)+NI
    DSYM(IJ) = D(NIJ)
  end do
end do
! DIAGONALIZE
call JACSCF(DSYM,CMO,OCC,NORBM,-1,1.0e-11_wp)
OCC(1:NORBM) = -OCC(1:NORBM)
call ORDER(CMO,OCC,NORBM)
if (IPRINT >= 15) write(u6,30)
ILAS = 0
OCC(1:NORBM) = -OCC(1:NORBM)
do I=1,NORBM
  ISTA = ILAS+1
  ILAS = ILAS+NORBM
  if (IPRINT >= 15) write(u6,40) I,OCC(I),(CMO(J),J=ISTA,ILAS)
end do
! TRANSFORM TO AO-BASIS
if (IPRINT >= 15) write(u6,45)
IJ0 = -NORBM
kp = 1
do I=1,NORBM
  IJ0 = IJ0+NORBM

  do IP=1,NBM
    TERM = Zero
    IJ = IJ0
    JP = IP-NBM
    do J=1,NORBM
      IJ = IJ+1
      JP = JP+NBM
      TERM = TERM+CMO(IJ)*CM(JP)
    end do
    CAO(kp+IP-1) = TERM
  end do

  if (IPRINT >= 15) write(u6,40) I,OCC(I),(CAO(IP),IP=kp,kp+nbm-1)
  kp = kp+NBM
end do

return

1 format(///,5X,'SYMMETRY NUMBER',I3)
30 format(//,5X,'NATURAL ORBITALS IN MO-BASIS',//,7X,'OCCUPATION NUMBER',5X,'COEFFICIENTS')
40 format(/,5X,I4,F10.6,5F10.6,/(19X,5F10.6))
45 format(//,5X,'NATURAL ORBITALS IN AO-BASIS',//,11X,'OCCUPATION NUMBER',5X,'COEFFICIENTS')

end subroutine NATORB_CPF

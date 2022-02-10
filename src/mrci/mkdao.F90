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

subroutine MKDAO(CNO,OCC,DAO)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension CNO(NCMO), OCC(NBAST), DAO(NBAST,NBAST)

call FZERO(DAO,NBAST**2)
IB = 1
ICNO = 1
do ISYM=1,NSYM
  IB1 = IB
  NB = NBAS(ISYM)
  do I=1,NB
    X = OCC(IB)
    call DGER(NB,NB,X,CNO(ICNO),1,CNO(ICNO),1,DAO(IB1,IB1),NBAST)
    IB = IB+1
    ICNO = ICNO+NB
  end do
end do

return

end subroutine MKDAO

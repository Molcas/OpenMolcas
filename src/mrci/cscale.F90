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

subroutine CSCALE(INDX,INTSYM,C,X)

implicit real*8(A-H,O-Z)
dimension C(*), INDX(*), INTSYM(*)
#include "SysDef.fh"
#include "mrci.fh"
!PAM97 external UNPACK
!PAM97 integer UNPACK
!Statement function
!PAM96 JSYM(L) = UNPACK(INTSYM((L+9)/10),3*mod(L-1,10)+1,3)+1
!JSYM(L) = JSUNP(INTSYM,L)

!do II1 = IRC(3)+1,IRC(4)
II1 = IRC(3)+1
30 if (II1 > IRC(4)) goto 10
!if (JSYM(II1) /= LSYM) goto 40
if (JSUNP(INTSYM,II1) /= LSYM) goto 40
NA = INDX(II1)
MA = 1
if (NVIRT < 1) goto 620
!do MA=1,NVIRT
720 C(NA+NDIAG(MA)) = X*C(NA+NDIAG(MA))
!end do
MA = MA+1
if (MA <= NVIRT) goto 720
620 continue
40 II1 = II1+1
goto 30
10 continue
!end do

return

end subroutine CSCALE

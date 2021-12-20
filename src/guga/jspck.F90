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
!***********************************************************************
      SUBROUTINE JSPCK(INTSYM,L,ISYM)
      DIMENSION INTSYM(*)

!      INTW=INTSYM((L+9)/10)
!      IPOW=2**(27-3*MOD(L-1,10))
!      INTW=INTW+(ISYM-1)*IPOW
!      INTSYM((L+9)/10)=INTW
       My=(L+9)/10
      IPOW=2**(27-3*MOD(L-1,10))
      INTSYM(My)=INTSYM(My)+(ISYM-1)*IPOW
      RETURN
      END

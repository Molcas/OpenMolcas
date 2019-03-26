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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

      Subroutine Loop_Kriging(qIntLast)
        use globvar
        Real*8 qIntLast
        write (6,*) "qIntlast", qIntLast
        Write (6,*) 'Kriging values in Loop Kriging'
        Write (6,*) 'iter:', iter
        Write (6,*) 'nInter', nInter
        Write (6,*) 'npxAI', npxAI
      end
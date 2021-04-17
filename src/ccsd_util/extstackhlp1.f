************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2006, Pavel Neogrady                                   *
************************************************************************
        subroutine extstackhlp1 (a,b,dimij,dimb,bb)
        integer dimij,dimb,bb,ij
        real*8 a(1:dimij),b(1:dimij,1:dimb)
        do ij=1,dimij
        a(ij)=b(ij,bb)
        end do
        return
        end

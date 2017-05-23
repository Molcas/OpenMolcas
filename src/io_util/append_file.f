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
* Copyright (C) 2015,2016, Valera Veryazov                             *
************************************************************************
            Subroutine Append_file(iUnit)
            iset=0
            rewind(iUnit)
10          read(iUnit,*,err=20,end=20)
            iset=iset+1
            goto 10
20          rewind(iUnit)
            do i=1,iset
            read(iUnit,*)
            enddo
            return
            end
c same for unformatted files
            Subroutine Append_file_u(iUnit)
            iset=0
            rewind(iUnit)
10          read(iUnit,err=20,end=20)
            iset=iset+1
            goto 10
20          rewind(iUnit)
            do i=1,iset
            read(iUnit)
            enddo
            return
            end

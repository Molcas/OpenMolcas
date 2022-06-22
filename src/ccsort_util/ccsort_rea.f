************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine ccsort_rea (lun,length,vector)
c
c     this routine read length-R8 numbers from opened unformatted file
c     with number lun form the given possition as one record
c
c     lun    - Logical unit number of file, where mediate is stored (Input)
c     length - # of R8 numbers to be read  (Input)
c     vector - space, where numbers are stored after reading  (Output)

c
       integer lun,length,i
       real*8 vector(1:length)
c
       read (lun) (vector(i),i=1,length)
c
       return
       end

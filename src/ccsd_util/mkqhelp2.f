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
       subroutine mkqhelp2 (vector,dimv,length,factor)
c
c     this routine do vector = vector*factot
c     vector - multilyied vector (I/O)
c     dimv   - dimension of vecrot
c     length - length of vector to be multiplyied
c     factor - scaling factor
c
c     $N.B. this routine should be substitued by mv0s3v
c
       integer dimv,length
       real*8 vector(1:dimv)
       real*8 factor
c
c     help variable
c
       integer n
c
       if (length.gt.0) then
       do 10 n=1,length
       vector(n)=vector(n)*factor
 10     continue
       end if
c
       return
       end

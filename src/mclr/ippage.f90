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
Module ipPage
Integer, Parameter:: Max_CI_Vectors=40
Integer, Parameter:: On_Disk=0, In_Memory=1, Null_Vector=2
Integer, Parameter:: dWrite=0, Write=1, Read=2

!         ip_Mem : memory pointer
!         n  : Length of CI-vector
!         ida: disk address

Integer:: ip_Mem(0:Max_CI_Vectors)
Integer::      n(0:Max_CI_Vectors)
Integer::    ida(0:Max_CI_Vectors)
Integer:: Status(0:Max_CI_Vectors)

Integer:: n_CI_Vectors=0
Integer:: iDisk_Addr_End=0
Integer:: Lu_ip=-99
Logical:: DiskBased=.False.
End Module ipPage

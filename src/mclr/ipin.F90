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
! Copyright (C) Anders Bernhardsson                                    *
!***********************************************************************
       Integer Function ipin(ii)
!
!      Object: return pointer to vector ii with a length of n(ii) and
!              make the vector available in memory as W(ii)%Vec
!
!
       use ipPage
       Implicit Integer (a-h,o-z)
!
       nn=n(ii)
       ipin = ipin1(ii,nn)
!
       Return
       End

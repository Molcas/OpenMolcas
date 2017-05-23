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
      subroutine par_range(n,i,j)
      ! Distribute 'n' evenly over processes and return
      ! the range (i,j) of this particular process.
      ! If there is no valid range, then j<i will be returned,
      ! so it is possible to use it consistently for looping.
#include "para_info.fh"
      nqot = n / nprocs
      nrem = n - nqot * nprocs
      if (myrank .lt. nrem) then
        i = myrank * (nqot + 1) + 1
        j = i + nqot
      else
        i = nrem * (nqot + 1) + (myrank - nrem) * nqot + 1
        j = i + nqot - 1
      end if
      end

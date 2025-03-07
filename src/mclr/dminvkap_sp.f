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
      SubRoutine DMInvKap_sp(rin,rout,isym)
*
*     _____     -1
*     Kappa  = M  Kappa
*          ip   pq     iq
*
*
*     In: rMFact        Factorized preconditioner (diagonal part
*                         of the electronic hessian that couple
*                         rotations with one common index)
*     In,Out rOut       Orbital rotaotion
*
*     iSym              Symmetry of rotation
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use MCLR_Data, only: nDens2
      Implicit None
      Integer iSym
      Real*8 rin(*),rout(*)

      Real*8, Allocatable:: Temp(:)
*
      Call mma_allocate(Temp,ndens2,Label='Temp')
      Call Uncompress(rin,Temp,isym)
*
      Call Compress(Temp,rout,isym)
      Call mma_deallocate(Temp)

c Avoid unused argument warnings
      end SubRoutine DMInvKap_sp

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
      Subroutine GR_DArray(Array,nArray)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
*
      integer(kind=iwp), intent(in):: nArray
      real(kind=wp), intent(inout):: Array(nArray)
*
#ifdef _MOLCAS_MPP_
      real(kind=wp) TCpu1,TWall1,TCpu2,TWall2

      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
      Call CWTime(TCpu1,TWall1)
      Call GADGOP(Array,nArray,'+')
      Call CWTime(TCpu2,TWall2)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Array)
#endif
*
      End Subroutine GR_DArray

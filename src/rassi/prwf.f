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
      SUBROUTINE PRWF(SGS,CIS,ISYCI,CI,CITHR)
      use definitions, only: iwp, wp
      use gugx, only: SGStruct, CIStruct
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      Type (SGStruct), intent(in):: SGS
      Type (CIStruct), intent(in):: CIS
      integer(kind=iwp) ISYCI
      real(kind=wp), intent(in):: CI(*), CITHR

      Integer(kind=iwp), Allocatable:: ICS(:)
      Integer(kind=iwp) NLEV,NMIDV

      NLEV  =SGS%nLev
      NMIDV =CIS%nMidV

      CALL mma_allocate(ICS,NLEV,Label='ICS')
      CALL PRWF1(SGS,CIS,NLEV,NMIDV,SGS%ISM,ICS,CIS%NOCSF,
     &           CIS%IOCSF,CIS%NOW,CIS%IOW,ISYCI,CI,CITHR)
      CALL mma_deallocate(ICS)

      END SUBROUTINE PRWF

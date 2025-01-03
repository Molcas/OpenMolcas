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
Module CSM_data
use lucia_data, only: mxpobs
Private mxpobs


!stuff from csm.fh
!        NSMSX        :        Number of symmetries single ex  (nIrrep)
!        NSMDX        :        Number of symmetries double ex  (nIrrep)
!        NSMST        :         Number of symmetries string        (nIrrep)
!        NSMCI        :        Number of symmetries CI Space   (nIrrep)
!        NSMXT        :        Not in use                         (nIrrep)
!        ITSSX        :        Total symmetrix single excitation (1)
!        ITSDX        :        Total symmetrix double excitation (1)
!        ITSXT        :         Not in use                          (1)
Integer    NSMSX,NSMDX,NSMST,NSMCI,NSMXT,ITSSX,ITSDX,ITSXT

!stuff from csmprd.fh
!        ADASX        : symmetry operator         initialized in syminf
!        ASXAD        :                           initialized in syminf
!        ADSXA        :                           initialized in syminf
!        SXSXDX       :                           initialized in syminf
!        SXDXSX       :                           initialized in syminf
!
integer ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),              &
              ADSXA(MXPOBS,2*MXPOBS),                             &
              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

End Module CSM_data

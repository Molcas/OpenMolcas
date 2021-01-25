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
*     nSh2BF   - field with dim. nIrrep,nShlls
*                # SO functions in irrep for shell iShell
*     iShOff   - field with dim. nIrrep,nShlls
*                position of 1st component of iShell in irrep
*     iSh2Sh   - field with dim. nShlls*nIrrep
*                pseudo shell index of iShell in irrep
*                (not any iShell contributes to all irreps)
*     nShIrp(0:nIrrep-1)
*              - # of shells contributing to each irrep
*     nShBFMx  - largest shell size (of all irrep, all shells)
*     iSO2Sh   - field with dim. nDim
*                shell the SO index iSO in irp belongs to
*     icntr  - field holding center number for each shell
*
      Module Index_arrays
      Integer nShIrp(0:7),nShBFMx
      Integer, Allocatable:: nShBF(:,:), iShOff(:,:), iSh2Sh(:,:),
     &                       iSO2Sh(:), iCntr(:)
      End Module Index_arrays

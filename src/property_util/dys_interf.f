************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2018, Jesper Norell                                    *
*               2018, Joel Creutzberg                                  *
************************************************************************
      Subroutine Dys_Interf(SO,i_root,i_file,NZ,CMO,ENE,OCC)
************************************************************************
!     Subroutine to generate .DysOrb and .molden files for Dyson orbitals
!     heavily based on the interf subroutine.

!     The UHF logical parameter to molden_interface could be used to
!     produce separate alpha and beta orbitals.
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "casvb.fh"
#include "WrkSpc.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='INTERF  ')
      Character*30 Filename
      Character*30 orbfile
      Character*80 Note
      INTEGER SO ! 0=SF, 1=SO-Real, 2=SO-Imaginary

      Dimension ENE(*)
      Dimension CMO(*)
      Dimension OCC(*)
*                                                                      *

************************************************************************
*                                                                      *
*     Make the .DysOrb file (for use with Molden_interface)
*
      IF (SO.EQ.0) THEN
       Write(orbfile,'(A10,I0,A1,I0)') 'DYSORB.SF.',i_root,".",i_file
      ELSE IF (SO.EQ.1) THEN
       Write(orbfile,'(A10,I0,A1,I0,A3)') 'DYSORB.SO.',i_root,
     &  ".",i_file,".Re"
      ELSE IF (SO.EQ.2) THEN
       Write(orbfile,'(A10,I0,A1,I0,A3)') 'DYSORB.SO.',i_root,
     &  ".",i_file,".Im"
      ENDIF
      Note='Temporary orbital file for the MOLDEN interface.'
      LuTmp=50
      LuTmp=IsFreeUnit(LuTmp)
      Call WrVec(orbfile,LuTmp,'COE',1,NZ,NZ,
     &     CMO,OCC,ENE,5,Note)

* Parameter number 4 should in principle give the symmetry of the
* orbital. For now we will put them all in the same symmetry 1.
* This could be improved later.

*                                                                      *
************************************************************************
*                                                                      *
*     Call the generic MOLDEN interface

      IF (SO.EQ.0) THEN
       Write(filename,'(A10,I0,A1,I0)') 'MD_DYS.SF.',i_root,
     &  ".",i_file
      ELSE IF (SO.EQ.1) THEN
       Write(filename,'(A10,I0,A1,I0,A3)') 'MD_DYS.SO.',i_root,
     &  ".",i_file,".Re"
      ELSE IF (SO.EQ.2) THEN
       Write(filename,'(A10,I0,A1,I0,A3)') 'MD_DYS.SO.',i_root,
     &  ".",i_file,".Im"
      ELSE
      ENDIF
      Call Molden_Interface(0,orbfile,filename,.False.)
*
* The first parameter is a logical for UHF which could be used to
* separate alpha and beta orbital in the future.
* Functionality enabled by the last paramter is not understood, but
* probably not useful for us.
************************************************************************
*                                                                      *
      Return
      End

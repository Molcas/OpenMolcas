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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine moscow_cvb()
      write(6,*)' Casvb dummy routine called : MOSCOW'
      return
      entry service_cvb()
      write(6,*)' Casvb dummy routine called : SERV'
      return
      entry rtransf_plc()
      write(6,*)' Molint dummy routine called : rtransf_plc'
      return
      entry perfloc_plc()
      write(6,*)' Molint dummy routine called : perfloc_plc'
      return
      entry plcconst_plc()
      write(6,*)' Molint dummy routine called : plcconst_plc'
      return
      entry rconstr_plc()
      write(6,*)' Molint dummy routine called : rconstr_plc'
      return
      entry getr_plc()
      write(6,*)' Molint dummy routine called : getr_plc'
      return
      entry qget_plc()
      write(6,*)' Molint dummy routine called : qget_plc'
      return
      end

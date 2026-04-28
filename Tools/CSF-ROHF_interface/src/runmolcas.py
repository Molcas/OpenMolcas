#!/usr/bin/env python3
"""
MOLCAS runner script - Python version
Runs MOLCAS computation with specified input file and monitors RASSCF iterations

Creates output files:
- .pylog: Python execution log with detailed information
- .iterdata: RASSCF iteration data (iteration number and RDM energy) for convergence analysis
"""

import os
import sys
import subprocess
import shutil
import argparse
import importlib.util
import time
import re
import threading
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Any


def _load_settings_value(name: str) -> Any:
    settings_path = Path(__file__).with_name("settings.py")
    spec = importlib.util.spec_from_file_location("settings", settings_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load settings from {settings_path}")

    settings_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(settings_module)

    return getattr(settings_module, name)


MOLCAS_PATH = _load_settings_value("MOLCAS_PATH")


class TeeOutput:
    """Class to redirect output to both console and log file"""
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log = open(log_file, 'w', buffering=1)  # Line buffered
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()
        
    def close(self):
        self.log.close()


class TeeError:
    """Class to redirect error output to both console and log file"""
    def __init__(self, log_file):
        self.terminal = sys.stderr
        self.log = open(log_file, 'a', buffering=1)  # Line buffered, append mode
        
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        
    def flush(self):
        self.terminal.flush()
        self.log.flush()
        
    def close(self):
        self.log.close()


def setup_environment(scratch_dir=None):
    """Set up environment variables for MOLCAS
    
    Parameters
    ----------
    scratch_dir : str, optional
        Path to scratch directory. If None, uses "./scratch" (default: None)
    """
    # Module loading equivalent (these would need to be handled by the system)
    # In Python, we assume the modules are already loaded or handle via environment

    # Set MOLCAS environment variables
    if not MOLCAS_PATH:
        print("Error: MOLCAS_PATH is not set in src/settings.py.", flush=True)
        print("Please edit src/settings.py and set MOLCAS_PATH to your OpenMolcas installation directory.", flush=True)
        raise SystemExit(1)

    molcas_path = MOLCAS_PATH
    os.environ["MOLCAS"] = molcas_path
    os.environ["MOLCASEXE"] = os.path.join(molcas_path, "pymolcas")

    # Set temporary directory
    tmpdir = scratch_dir if scratch_dir else "./scratch"
    os.environ["TMPDIR"] = tmpdir
    os.environ["WorkDir"] = tmpdir

    return tmpdir, os.environ["MOLCASEXE"]


def parse_molcas_iteration_data(log_file, iteration_number):
    """Parse MOLCAS log file to extract iteration data for a specific iteration
    
    The iteration data appears after "I read the following energies" lines, but there may be
    additional diagnostic output in between.
    Format: "        1   0    5    1  -107.43700441    0.00E+00  -5.59E-03    7  10 1 -1.49E-02*  1.00   0.00    SX    NO    0:00:01"
    """
    if not os.path.exists(log_file):
        return None
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        # Find the iteration data that comes after "I read the following energies" 
        # for the specific iteration number
        for i, line in enumerate(lines):
            # Look for "I read the following energies" line
            if "I read the following energies" in line:
                # Check subsequent lines until we find our iteration data
                # Increased search range to handle additional diagnostic output
                for j in range(1, 50):  # Search up to 50 lines after
                    if i + j < len(lines):
                        next_line = lines[i + j].strip()
                        # Look for line that starts with our iteration number
                        import re
                        # Pattern to match iteration data line:
                        # "        1   0    5    1  -107.43700441    0.00E+00  -5.59E-03    7  10 1 -1.49E-02*  1.00   0.00    SX    NO    0:00:01"
                        pattern = rf'^\s*{iteration_number}\s+(\d+)\s+(\d+)\s+(\d+)\s+([-\d\.E\+]+)\s+([-\d\.E\+\*]+)\s+([-\d\.E\+\*]+)\s+(\d+)(\s+(\d+))?\s+(\d+)\s+([-\d\.E\+\*]+)\s+([\d\.]+)\s+([-\d\.]+)\s+(\w+)\s+(\w+)\s+([\d:]+)'
                        
                        match = re.search(pattern, next_line)
                        if match:
                            ci_iter = int(match.group(1))
                            sx_iter = int(match.group(2))
                            ci_root = int(match.group(3))
                            rasscf_energy = float(match.group(4))
                            energy_change = match.group(5).replace('*', '')  # Remove * marker
                            max_rot_param = match.group(6).replace('*', '')  # Remove * marker
                            # max BLB fields could be merged with > 999 orbitals
                            if match.group(8) is None:
                                max_blb_elem = int(match.group(7)[:-4])
                                max_blb_val = int(match.group(7)[-4:])
                            else:
                                max_blb_elem = int(match.group(7))
                                max_blb_val = int(match.group(9))
                            max_blb_idx = int(match.group(10))
                            max_blb_value = match.group(11).replace('*', '')  # Remove * marker
                            level_shift = float(match.group(12))
                            ln_srch_min = float(match.group(13))
                            step_type = match.group(14)
                            qn_update = match.group(15)
                            walltime = match.group(16)
                            
                            return {
                                'ci_iter': ci_iter,
                                'sx_iter': sx_iter,
                                'ci_root': ci_root,
                                'rasscf_energy': rasscf_energy,
                                'energy_change': energy_change,
                                'max_rot_param': max_rot_param,
                                'max_blb_elem': max_blb_elem,
                                'max_blb_val': max_blb_val,
                                'max_blb_idx': max_blb_idx,
                                'max_blb_value': max_blb_value,
                                'level_shift': level_shift,
                                'ln_srch_min': ln_srch_min,
                                'step_type': step_type,
                                'qn_update': qn_update,
                                'walltime': walltime
                            }
                            
    except Exception as e:
        print(f"Warning: Could not parse MOLCAS log file for iteration {iteration_number}: {e}", flush=True)
    
    return None


def check_for_neci_message(log_file, last_position):
    """Check if 'Run spin-free GUGA NECI externally.' appears in log file after last_position
    
    Returns:
        tuple: (found: bool, new_position: int)
            found: True if the message was found after last_position
            new_position: Current file size to use for next check
    """
    if not os.path.exists(log_file):
        return False, last_position
    
    try:
        with open(log_file, 'r') as f:
            # Seek to last known position
            f.seek(last_position)
            # Read new content
            new_content = f.read()
            # Get current position
            new_position = f.tell()
            
        # Check if the message appears in new content
        if "Run spin-free GUGA NECI externally." in new_content:
            return True, new_position
        else:
            return False, new_position
            
    except Exception as e:
        print(f"Warning: Could not read log file: {e}", flush=True)
        return False, last_position


def extract_rdm_energy_from_fciqmc_output(output_file):
    """Extract RDM energy from FCIQMC output file
    
    Looks for line containing "TOTAL" like:
    *TOTAL ENERGY* CALCULATED USING THE *REDUCED DENSITY MATRICES*: -5.0927404720233E+03
    
    Parameters
    ----------
    output_file : str
        Path to FCIQMC output file
    
    Returns
    -------
    float
        RDM energy extracted from output
    
    Raises
    ------
    RuntimeError
        If RDM energy line not found in output
    """
    if not os.path.exists(output_file):
        raise RuntimeError(f"FCIQMC output file not found: {output_file}")
    
    try:
        with open(output_file, 'r') as f:
            for line in f:
                if "TOTAL" in line:
                    # Extract the last column (energy value)
                    parts = line.strip().split()
                    if parts:
                        energy_str = parts[-1]
                        return float(energy_str)
        
        raise RuntimeError("RDM energy line not found in FCIQMC output")
    
    except Exception as e:
        raise RuntimeError(f"Failed to extract RDM energy from FCIQMC output: {e}")


def run_external_fciqmc(fciqmc_dir, fcidump_path, neci_command, workdir):
    """Run external FCIQMC calculation
    
    Parameters
    ----------
    fciqmc_dir : str
        Path to FCIQMC directory containing input file
    fcidump_path : str
        Path to FCIDUMP file to copy
    neci_command : str
        Command to execute NECI (e.g., "mpirun -np 10 /path/to/neci.exe")
    workdir : str
        Scratch directory where RDM files should be copied back
    
    Returns
    -------
    float
        RDM energy extracted from FCIQMC output
    
    Raises
    ------
    RuntimeError
        If FCIQMC execution fails or RDM files not created
    """
    print(f"Running external FCIQMC in directory: {fciqmc_dir}", flush=True)
    
    # Get absolute path of fciqmc_dir to ensure proper working directory
    fciqmc_dir_abs = os.path.abspath(fciqmc_dir)
    
    # Copy FCIDUMP to fciqmc_dir with name "FCIDUMP"
    fcidump_dest = os.path.join(fciqmc_dir_abs, "FCIDUMP")
    shutil.copy(fcidump_path, fcidump_dest)
    print(f"Copied {fcidump_path} to {fcidump_dest}", flush=True)
    
    # Execute FCIQMC - output files will be created in fciqmc_dir
    fciqmc_output = os.path.join(fciqmc_dir_abs, "out")
    fciqmc_error = os.path.join(fciqmc_dir_abs, "err")
    
    # Parse command to handle shell execution properly
    cmd = f"cd {fciqmc_dir_abs} && {neci_command} input > out 2> err"
    print(f"Executing FCIQMC with command: {cmd}", flush=True)
    print(f"Working directory: {fciqmc_dir_abs}", flush=True)
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=False,
            text=True
        )
        
        if result.returncode != 0:
            # Read error file if it exists
            error_msg = ""
            if os.path.exists(fciqmc_error):
                with open(fciqmc_error, 'r') as f:
                    error_msg = f.read()
            raise RuntimeError(f"FCIQMC execution failed with return code {result.returncode}\n{error_msg}")
        
        print("FCIQMC execution completed", flush=True)
        
    except Exception as e:
        raise RuntimeError(f"Failed to execute FCIQMC: {e}")
    
    # Check for RDM files
    rdm_files = ['DMAT.1', 'PSMAT.1', 'PAMAT.1']
    missing_files = []
    
    for rdm_file in rdm_files:
        rdm_src = os.path.join(fciqmc_dir_abs, rdm_file)
        if not os.path.exists(rdm_src):
            missing_files.append(rdm_file)
    
    if missing_files:
        raise RuntimeError(f"FCIQMC did not create required RDM files: {missing_files}")
    
    # Copy RDM files to scratch directory
    for rdm_file in rdm_files:
        rdm_src = os.path.join(fciqmc_dir_abs, rdm_file)
        rdm_dest = os.path.join(workdir, rdm_file)
        shutil.copy(rdm_src, rdm_dest)
        print(f"Copied {rdm_src} to {rdm_dest}", flush=True)
    
    # Extract RDM energy from output
    rdm_energy = extract_rdm_energy_from_fciqmc_output(fciqmc_output)
    print(f"Extracted RDM energy from FCIQMC output: {rdm_energy}", flush=True)
    
    return rdm_energy


def monitor_status_file(filename, workdir, CSF_stepvec, stop_event=None, debug=False, sleep_interval=0.5, manual_mode=False, fciqmc_dir=None, neci_command=None):
    """Monitor the status file for RASSCF iterations and write to NEWCYCLE when needed
    
    Creates a .iterdata file containing RASSCF iteration numbers, MOLCAS iteration data,
    and RDM energies for easy analysis of convergence behavior.
    
    Parameters
    ----------
    filename : str
        Base filename for MOLCAS calculation
    workdir : str
        Scratch/working directory
    CSF_stepvec : array-like
        CSF step vector
    stop_event : threading.Event, optional
        Event to signal when monitoring should stop
    debug : bool, optional
        Enable debug mode
    sleep_interval : float, optional
        Sleep interval for monitoring
    manual_mode : bool, optional
        Enable manual FCIQMC execution mode
    fciqmc_dir : str, optional
        Path to FCIQMC directory (required if manual_mode=True)
    neci_command : str, optional
        NECI execution command (required if manual_mode=True)
    """
    import IntegralClass
    import GUGA_diag
    import gen_spinfree_rdm as gen_rdm
    
    status_file = f"{filename}.status"
    molcas_log_file = f"{filename}.log"
    last_iteration = 0
    max_wait_time = 3600 * 100000 # Maximum wait time in seconds
    start_time = time.time()
    IntegralClass_instance = None
    rdm_energy = 0.0
    no_status_file_count = 0  # Track how long status file has been missing
    log_file_position = 0  # Track position in log file to detect new messages
    
    # Initialize iteration data file
    iterdata_file = f"{filename}.iterdata"
    with open(iterdata_file, 'w') as f:
        f.write("# RASSCF Iteration Data from MOLCAS Log\n")
        f.write("# Columns: Iter, CI_iter, SX_iter, CI_root, RASSCF_energy, Energy_change, max_ROT_param, max_BLB_elem, max_BLB_val, max_BLB_idx, max_BLB_value, Level_shift, Ln_srch_min, Step_type, QN_update, Walltime, RDM_Energy\n")
        f.write("# Iter  CI_iter  SX_iter  CI_root     RASSCF_energy        Energy_change       max_ROT_param       max_BLB_elem    max_BLB_val     max_BLB_idx        max_BLB_value     Level_shift     Ln_srch_min      Step_type   QN_update      Walltime           RDM_Energy\n")
    print(f"Created iteration data file: {iterdata_file}", flush=True)
    
    print(f"Monitoring status file: {status_file}", flush=True)
    
    while time.time() - start_time < max_wait_time:
        # Check if stop event is set (MOLCAS process finished)
        if stop_event and stop_event.is_set():
            print("Stop event detected, monitor exiting", flush=True)
            break
            
        if os.path.exists(status_file):
            no_status_file_count = 0  # Reset counter
            try:
                with open(status_file, 'r') as f:
                    content = f.read()
                
                # Look for RASSCF iteration pattern
                match = re.search(r'RASSCF:\s+Iteration\s+(\d+)', content)

                if match:
                    current_iteration = int(match.group(1))
                    current_iteration = int(match.group(1))
                    
                    if current_iteration > last_iteration:
                        print(f"Detected NEW RASSCF iteration {current_iteration}", flush=True)
                        
                        # Move RDM files to scratch directory on first iteration only
                        if IntegralClass_instance is None:
                            print("Moving RDM files to scratch directory...", flush=True)
                            rdm_files = ['DMAT.1', 'PSMAT.1', 'PAMAT.1']
                            missing_files = []
                            
                            for rdm_file in rdm_files:
                                if os.path.exists(rdm_file):
                                    dest_path = os.path.join(workdir, rdm_file)
                                    shutil.copy(rdm_file, dest_path)
                                    print(f"Moved {rdm_file} to {dest_path}", flush=True)
                                else:
                                    missing_files.append(rdm_file)
                            
                            if missing_files:
                                print(f"ERROR: Required RDM files not found: {missing_files}", flush=True)
                                print("Cannot proceed without RDM files - killing entire job...", flush=True)
                                os._exit(1)
                        
                        # First wait for "Run spin-free GUGA NECI externally." message in log file
                        print(f"Waiting for 'Run spin-free GUGA NECI externally.' message in log file for iteration {current_iteration}...", flush=True)
                        neci_message_found = False
                        neci_wait_start = time.time()
                        max_wait_for_neci = max_wait_time
                        
                        while time.time() - neci_wait_start < max_wait_for_neci:
                            message_found, log_file_position = check_for_neci_message(molcas_log_file, log_file_position)
                            if message_found:
                                neci_message_found = True
                                print(f"Detected 'Run spin-free GUGA NECI externally.' message for iteration {current_iteration}", flush=True)
                                break
                            time.sleep(1)  # Check every second
                        
                        if not neci_message_found:
                            print("ERROR: NECI message not found in log file after waiting", flush=True)
                            print("Killing entire job...", flush=True)
                            os._exit(1)
                        
                        # Now wait for integral file to be created for this iteration
                        print(f"Waiting for integral file for iteration {current_iteration}...", flush=True)
                        fcidmp_path = os.path.join(workdir, f"{filename}.FciDmp")
                        max_wait_for_file = max_wait_time
                        wait_start = time.time()
                        file_found = False
                        
                        # Since we delete/rename the file after each iteration, just wait for it to exist
                        while time.time() - wait_start < max_wait_for_file:
                            if os.path.exists(fcidmp_path):
                                file_found = True
                                break
                            time.sleep(1)  # Check every second
                        
                        if not file_found:
                            print("ERROR: Integral file not found after waiting", flush=True)
                            print("Killing entire job...", flush=True)
                            os._exit(1)
                        
                        # First, rename/move the FCIDUMP file before reading
                        # This ensures we wait for a fresh file in the next iteration
                        if debug:
                            # In debug mode, print the status file content before renaming
                            print(f"=== Debug: MOLCAS status file content before processing iteration {current_iteration} ===", flush=True)
                            try:
                                with open(status_file, 'r') as f:
                                    status_content = f.read()
                                    print(status_content, flush=True)
                            except Exception as e:
                                print(f"Debug: Could not read status file: {e}", flush=True)
                            print("=== End of status file content ===", flush=True)
                            
                            # In debug mode, rename to preserve the file
                            final_fcidmp_path = os.path.join(workdir, f"{filename}.FciDmp.iter{current_iteration}")
                            shutil.move(fcidmp_path, final_fcidmp_path)
                            print(f"Debug mode: Moved FCIDUMP to {final_fcidmp_path}", flush=True)
                        else:
                            # In normal mode, rename to a temporary name, will delete after reading
                            final_fcidmp_path = os.path.join(workdir, f"{filename}.FciDmp.processing")
                            shutil.move(fcidmp_path, final_fcidmp_path)
                            print(f"Moved FCIDUMP to {final_fcidmp_path} for processing", flush=True)
                        
                        print(f"Reading integral file from {final_fcidmp_path}...", flush=True)
                        
                        try:
                            if manual_mode:
                                # Manual FCIQMC mode: run external FCIQMC
                                if not fciqmc_dir or not neci_command:
                                    print("ERROR: Manual mode requires fciqmc_dir and neci_command", flush=True)
                                    os._exit(1)
                                
                                rdm_energy = run_external_fciqmc(fciqmc_dir, final_fcidmp_path, neci_command, workdir)
                                
                                # In manual mode, we still clean up the FCIDUMP file unless in debug mode
                                if not debug:
                                    os.remove(final_fcidmp_path)
                                    print(f"Deleted temporary FCIDUMP file", flush=True)
                            else:
                                # Normal mode: calculate RDM energy using IntegralClass
                                # Re-initialize integral class every iteration to read updated integrals
                                IntegralClass_instance = IntegralClass.FCIDUMPReader(final_fcidmp_path)
                                # Calculate RDM energy here using your specific method
                                DiagElement_instance = GUGA_diag.DiagElement(len(CSF_stepvec), IntegralClass_instance)
                                rdm_energy = DiagElement_instance.calc_diag_elem(CSF_stepvec, add_core=True)
                                print(f"RDM energy calculated: {rdm_energy}", flush=True)
                                
                                # In normal mode, delete the temporary file after processing
                                if not debug:
                                    os.remove(final_fcidmp_path)
                                    print(f"Deleted temporary FCIDUMP file", flush=True)
                                    
                        except Exception as e:
                            print(f"ERROR: Failed to process integral file or calculate energy: {e}", flush=True)
                            print("Killing entire job...", flush=True)
                            os._exit(1)
                        
                        # Write RDM energy to NEWCYCLE file
                        newcycle_file = os.path.join(workdir, "NEWCYCLE")
                        with open(newcycle_file, 'w') as f:
                            f.write(f"{rdm_energy}\n")
                        print(f"Written RDM energy ({rdm_energy}) to {newcycle_file}", flush=True)
                        
                        # Write iteration data to .iterdata file
                        # Monitor the log file until MOLCAS iteration data appears
                        print(f"Waiting for MOLCAS iteration data in log file for iteration {current_iteration}...", flush=True)
                        molcas_data = None
                        max_wait_for_molcas_data = max_wait_time
                        molcas_data_wait_start = time.time()
                        
                        while time.time() - molcas_data_wait_start < max_wait_for_molcas_data:
                            molcas_data = parse_molcas_iteration_data(molcas_log_file, current_iteration)
                            if molcas_data:
                                print(f"MOLCAS iteration data found for iteration {current_iteration}", flush=True)
                                break
                            time.sleep(1)  # Check every second
                        
                        if not molcas_data:
                            print(f"Warning: MOLCAS iteration data not found for iteration {current_iteration} after waiting {max_wait_for_molcas_data} seconds", flush=True)
                        
                        with open(iterdata_file, 'a') as f:
                            if molcas_data:
                                f.write(f"{current_iteration:6d}  {molcas_data['ci_iter']:8d}  {molcas_data['sx_iter']:8d}  {molcas_data['ci_root']:8d}  {molcas_data['rasscf_energy']:16.8f}  {molcas_data['energy_change']:>16s}  {molcas_data['max_rot_param']:>16s}  {molcas_data['max_blb_elem']:13d}  {molcas_data['max_blb_val']:11d}  {molcas_data['max_blb_idx']:11d}  {molcas_data['max_blb_value']:>16s}  {molcas_data['level_shift']:12.6f}  {molcas_data['ln_srch_min']:12.2f}  {molcas_data['step_type']:>10s}  {molcas_data['qn_update']:>8s}  {molcas_data['walltime']:>12s}  {rdm_energy:16.10f}\n")
                                print(f"Written iteration data with MOLCAS info (iter={current_iteration}, molcas_energy={molcas_data['rasscf_energy']}, rdm_energy={rdm_energy}) to {iterdata_file}", flush=True)
                            else:
                                # Fallback if MOLCAS data not available after timeout
                                f.write(f"{current_iteration:6d}  {'N/A':>8s}  {'N/A':>8s}  {'N/A':>8s}  {'N/A':>16s}  {'N/A':>16s}  {'N/A':>16s}  {'N/A':>13s}  {'N/A':>11s}  {'N/A':>11s}  {'N/A':>16s}  {'N/A':>12s}  {'N/A':>12s}  {'N/A':>10s}  {'N/A':>8s}  {'N/A':>12s}  {rdm_energy:16.10f}\n")
                                print(f"Warning: Written N/A values for MOLCAS iteration data for iteration {current_iteration}", flush=True)
                        
                        last_iteration = current_iteration
                else:
                    # Check if we're no longer in RASSCF iterations
                    if last_iteration > 0:
                        # Look for completion indicators
                        if any(indicator in content for indicator in ['RASSCF ENDS', 'SEWARD', 'MCKINLEY', 'End of calculation']):
                            print("RASSCF iterations completed", flush=True)
                            break
            except (IOError, OSError):
                # File might be being written to, wait a bit
                pass
        else:
            # Status file doesn't exist
            no_status_file_count += 1
            # If we've processed iterations and status file is gone for 10 seconds, assume completion
            if last_iteration > 0 and no_status_file_count > 20:  # 20 * 0.5s = 10 seconds
                print("Status file disappeared after iterations, assuming completion", flush=True)
                break
        
        time.sleep(sleep_interval)  # Check interval configurable via command line
    
    if time.time() - start_time >= max_wait_time:
        print("Monitoring timeout reached", flush=True)
    else:
        print("Monitor thread exiting normally", flush=True)


def run_molcas(filename, CSF_stepvec, debug=False, sleep_interval=0.5, manual_mode=False, fciqmc_dir=None, neci_command=None, scratch_dir=None):
    """Run MOLCAS computation with the given filename and monitor for RASSCF iterations
    
    Parameters
    ----------
    scratch_dir : str, optional
        Path to scratch directory. If None, uses "./scratch" (default: None)
    """
    tmpdir, molcas_exe = setup_environment(scratch_dir)

    try:
        # Clean up and create temporary directory
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.makedirs(tmpdir, exist_ok=True)

        # Delete status file if it exists from previous run
        status_file = f"{filename}.status"
        if os.path.exists(status_file):
            os.remove(status_file)
            print(f"Removed existing status file: {status_file}", flush=True)

        # Prepare the input filename
        input_file = f"{filename}.inp"

        # Run MOLCAS command
        cmd = [molcas_exe, "-nt", "1", "-np", "1", "-b", "1", "-f", input_file]

        print(f"Running MOLCAS with input file: {input_file}", flush=True)
        print(f"Command: {' '.join(cmd)}", flush=True)

        # Start MOLCAS process
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Create stop event for monitor thread
        stop_event = threading.Event()
        
        # Start monitoring in a separate thread
        monitor_thread = threading.Thread(
            target=monitor_status_file, 
            args=(filename, tmpdir, CSF_stepvec, stop_event, debug, sleep_interval, manual_mode, fciqmc_dir, neci_command)
        )
        monitor_thread.daemon = False  # Changed to non-daemon so it completes properly
        monitor_thread.start()

        # Wait for MOLCAS to complete
        stdout, stderr = process.communicate()
        
        # Signal monitor thread to stop
        stop_event.set()

        if process.returncode == 0:
            print("MOLCAS computation completed successfully!", flush=True)
        else:
            print(f"MOLCAS computation failed with return code: {process.returncode}", flush=True)

        if stdout:
            print("STDOUT:", stdout, flush=True)
        if stderr:
            print("STDERR:", stderr, flush=True)

        # Wait for monitor thread to complete (increased timeout)
        print("Waiting for monitor thread to complete...", flush=True)
        monitor_thread.join(timeout=60)  # Increased from 5 to 60 seconds
        if monitor_thread.is_alive():
            print("WARNING: Monitor thread did not complete in time", flush=True)
        else:
            print("Monitor thread completed successfully", flush=True)


    except subprocess.CalledProcessError as e:
        print(f"Error running MOLCAS: {e}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)
    finally:
        # Clean up temporary directory (unless in debug mode)
        if debug:
            print(f"Debug mode: Preserving scratch directory: {tmpdir}")
        else:
            if os.path.exists(tmpdir):
                shutil.rmtree(tmpdir)
                print(f"Cleaned up temporary directory: {tmpdir}")


def run_molcas_with_csf(filename, csf_stepvec, debug=False, sleep_interval=0.5, user_rdm=False, manual_mode=False, fciqmc_dir=None, neci_command=None, scratch_dir=None):
    """
    Library function to run MOLCAS computation with specified filename and CSF.
    
    This function can be imported and used from other Python scripts.
    
    Creates output files:
    - {filename}.pylog: Python execution log
    - {filename}.iterdata: RASSCF iteration data (iteration number and RDM energy)
    
    Parameters
    ----------
    filename : str
        Base filename for the MOLCAS input file (without .inp extension)
    csf_stepvec : array-like
        CSF step vector values (will be converted to numpy array)
    debug : bool, optional
        Enable debug mode to preserve scratch directory and FCIDUMP files (default: False)
    sleep_interval : float, optional
        Sleep interval in seconds for status file monitoring (default: 0.5)
    user_rdm : bool, optional
        If True, skip RDM file generation and use user-provided files (default: False)
    manual_mode : bool, optional
        If True, use external FCIQMC execution instead of internal RDM energy calculation (default: False)
    fciqmc_dir : str, optional
        Path to FCIQMC directory containing input file (required if manual_mode=True)
    neci_command : str, optional
        NECI execution command (required if manual_mode=True)
    scratch_dir : str, optional
        Path to scratch directory. If None, uses "./scratch" (default: None)
    
    Returns
    -------
    None
    
    Raises
    ------
    FileNotFoundError
        If the input file does not exist
    RuntimeError
        If MOLCAS computation fails
    
    Examples
    --------
    >>> import numpy as np
    >>> from runmolcas import run_molcas_with_csf
    >>> csf = np.array([1, 2, 3, 4])
    >>> run_molcas_with_csf('benzene', csf)
    """
    import gen_spinfree_rdm as gen_rdm
    
    # Set up logging to .pylog file
    log_filename = f"{filename}.pylog"
    tee_stdout = TeeOutput(log_filename)
    tee_stderr = TeeError(log_filename)
    
    # Save original stdout/stderr
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    # Redirect output
    sys.stdout = tee_stdout
    sys.stderr = tee_stderr
    
    try:
        # Print header with calculation variables
        print("="*80)
        print(f"MOLCAS Python Runner - Execution Log")
        print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*80)
        print()
        print("CALCULATION VARIABLES:")
        print(f"  Basename:        {filename}")
        print(f"  Input file:      {filename}.inp")
        print(f"  Log file:        {log_filename}")
        print(f"  Iteration data:  {filename}.iterdata")
        print(f"  Debug mode:      {debug}")
        print(f"  Sleep interval:  {sleep_interval} seconds")
        print(f"  Scratch dir:     {scratch_dir}")
        print(f"  User RDM mode:   {user_rdm}")
        print(f"  Manual mode:     {manual_mode}")
        if manual_mode:
            print(f"  FCIQMC dir:      {fciqmc_dir}")
            print(f"  NECI command:    {neci_command}")
        print(f"  CSF step vector: {csf_stepvec}")
        print(f"  CSF length:      {len(csf_stepvec)}")
        print("="*80)
        print()
        
        # Convert to numpy array if needed
        CSF_stepvec = np.array(csf_stepvec)
        print(f"CSF step vector: {CSF_stepvec}")
        
        # Check if input file exists
        input_file = f"{filename}.inp"
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file '{input_file}' not found!")
        
        # Validate manual mode parameters
        if manual_mode:
            if not fciqmc_dir:
                raise ValueError("Manual mode requires fciqmc_dir parameter")
            if not neci_command:
                raise ValueError("Manual mode requires neci_command parameter")
            
            # Convert to absolute path if relative
            if not os.path.isabs(fciqmc_dir):
                fciqmc_dir = os.path.abspath(fciqmc_dir)
            
            # Check if fciqmc_dir exists
            if not os.path.isdir(fciqmc_dir):
                raise FileNotFoundError(f"FCIQMC directory not found: {fciqmc_dir}")
            
            # Check if input file exists in fciqmc_dir
            fciqmc_input = os.path.join(fciqmc_dir, "input")
            if not os.path.exists(fciqmc_input):
                raise FileNotFoundError(f"FCIQMC input file not found: {fciqmc_input}")
            
            print(f"Manual mode validated:")
            print(f"  FCIQMC directory: {fciqmc_dir}")
            print(f"  FCIQMC input file: {fciqmc_input}")
            print(f"  NECI command: {neci_command}")
            print()
        
        # Generate or check RDM files
        if user_rdm:
            # Check if user-provided RDM files exist
            print("User RDM mode enabled - checking for user-provided RDM files...")
            rdm_files = ['DMAT.1', 'PSMAT.1', 'PAMAT.1']
            missing_files = []
            for rdm_file in rdm_files:
                if os.path.exists(rdm_file):
                    print(f"  Found: {rdm_file}")
                else:
                    missing_files.append(rdm_file)
            
            if missing_files:
                raise FileNotFoundError(f"User RDM mode enabled but required files not found: {missing_files}")
            
            print("All required RDM files found.")
        else:
            # Generate RDM files in current directory (will be moved to scratch later)
            print("Generating RDM files...")
            gen_rdm.write_spinfree_1rdm(CSF_stepvec, filename='DMAT.1', thr=1e-12)
            gen_rdm.write_spinfree_2rdm(CSF_stepvec, thr=1e-12)
        
        # Run MOLCAS computation
        run_molcas(filename, CSF_stepvec, debug=debug, sleep_interval=sleep_interval, 
                   manual_mode=manual_mode, fciqmc_dir=fciqmc_dir, neci_command=neci_command,
                   scratch_dir=scratch_dir)
        
        # Print footer
        print()
        print("="*80)
        print(f"Execution completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("="*80)
        
    finally:
        # Restore original stdout/stderr
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        
        # Close log files
        tee_stdout.close()
        tee_stderr.close()


def main():
    """Main function to handle command line arguments and run MOLCAS"""
    parser = argparse.ArgumentParser(description="Run MOLCAS computation with specified input file")
    parser.add_argument("filename", help="Base filename for the MOLCAS input file (without .inp extension)")
    parser.add_argument("stepvec", nargs='+', type=int, help="CSF step vector values")
    parser.add_argument("-d", "--debug", action="store_true", 
                        help="Enable debug mode (preserve scratch directory and FCIDUMP files)")
    parser.add_argument("-s", "--sleep", type=float, default=0.5, dest="sleep_interval",
                        help="Sleep interval in seconds for status file monitoring (default: 0.5)")
    parser.add_argument("-u", "--user-rdm", action="store_true", dest="user_rdm",
                        help="Use user-provided RDM files (DMAT.1, PSMAT.1, PAMAT.1) instead of generating them")
    parser.add_argument("-m", "--manual", type=str, dest="fciqmc_dir", metavar="FCIQMC_DIR",
                        help="Enable manual FCIQMC mode with path to FCIQMC directory (relative or absolute)")
    parser.add_argument("--scratch", type=str, dest="scratch_dir", metavar="SCRATCH_DIR",
                        help="Path to scratch directory (default: ./scratch)")

    args = parser.parse_args()

    # If manual mode is enabled, prompt for NECI command
    neci_command = None
    if args.fciqmc_dir:
        print("Manual FCIQMC mode enabled")
        print(f"FCIQMC directory: {args.fciqmc_dir}")
        
        # Prompt user for NECI execution command
        print("\nPlease enter the NECI execution command.")
        print("Example: mpirun -np 10 /home/ab/cd/ef/neci.exe")
        neci_command = input("NECI command: ").strip()
        
        if not neci_command:
            print("Error: NECI command cannot be empty")
            sys.exit(1)
        
        print(f"Using NECI command: {neci_command}")
        print()

    try:
        # Use the library function
        run_molcas_with_csf(
            args.filename, 
            args.stepvec, 
            debug=args.debug, 
            sleep_interval=args.sleep_interval, 
            user_rdm=args.user_rdm,
            manual_mode=(args.fciqmc_dir is not None),
            fciqmc_dir=args.fciqmc_dir,
            neci_command=neci_command,
            scratch_dir=args.scratch_dir
        )
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()

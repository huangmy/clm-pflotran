#!/usr/bin/env python

"""
Program to manage and run CLM-PFLOTRAN regression tests

Author: Ben Andre <bandre@lbl.gov>

"""

from __future__ import print_function
from __future__ import division

import argparse
import copy
from collections import deque
import csv
import datetime
import gzip
import math
import os
import platform
import pprint
import re
import shutil
import string
import subprocess
import sys
import textwrap
import time
import traceback

if sys.version_info[0] == 2:
    import ConfigParser as config_parser
else:
    import configparser as config_parser


class TestMachine(object):
    """Claass to hold machine specific information required to run clm-pflotran tests.


    Input file for machine information is a standard config file with
    the following format. The host name should be the output of the
    command "hostname -s".

    [hostname]
    machine = userdefined
    os = Darwin
    compiler = gnu
    mpi_vendor = mpich
    max_np = 
    gmake =
    cc = __path_to_cc__
    cxx = __path_to_cxx__
    fc = __path_to_fc__
    mpicc = __path_to_mpicc__
    mpicxx = __path_to_mpicxx__
    mpifc = __path_to_mpifc__
    mpiexec = __path_to_mpiexec__
    fflags = -fno-range-check
    slibs = -L/opt/local/lib -lnetcdff -L/opt/local/lib -lnetcdf
    netcdf_path = /opt/local
    blas_flags = -framework Accelerate
    ldflags = -all_load

    """

    def __init__(self):
        self._builder_id = None
        self._userdefined = False
        self._create_case_info = {}
        self._env_build_info = {}
        self._env_run_info = {}
        self._env_mach_pes_info = {}
        self._macros_info = {}


    def __str__(self):
        txt = "  Machine : {0}\n".format(self._builder_id)
        txt += "    create_case :\n"
        for p in self._create_case_info:
            txt += "      {0} : {1}\n".format(p, self._create_case_info[p])

        txt += "    env_build :\n"
        for p in self._env_build_info:
            txt += "      {0} : {1}\n".format(p, self._env_build_info[p])

        txt += "    macros :\n"
        for p in self._macros_info:
            txt += "      {0} : {1}\n".format(p, self._macros_info[p])

        txt += "    env_run :\n"
        for p in self._env_run_info:
            txt += "      {0} : {1}\n".format(p, self._env_run_info[p])

        txt += "    env_mach_pes :\n"
        for p in self._env_mach_pes_info:
            txt += "      {0} : {1}\n".format(p, self._env_mach_pes_info[p])

        return txt

    def setup(self, config_filename):
        self._set_builder_id()
        self._read_config_file(config_filename)

    def _set_builder_id(self):
        # import socket; hostname = socket.gethostname()
        hostname = platform.node()
        index = hostname.find(".")
        if index > 0:
            self._builder_id = hostname[0:index]
        else:
            self._builder_id = hostname

    def _read_config_file(self, filename):
        # NOTE(bja, 2013-10-12) reading the config file in the class
        # makes it harder to unit test the class!
        self._config_filename = filename
        config = config_parser.SafeConfigParser()
        config.read(self._config_filename)
        if not config.has_section(self._builder_id):
            raise Exception("ERROR: machine config does not contain a section for"
                            "current host: {0}".format(self._builder_id))
        host_info = list_to_dict(config.items(self._builder_id))
        
        machine = self._set_host_property(host_info, "machine")
        
        if machine == "userdefined":
            self._userdefined = True
        else:
            self._userdefined = False

        # create_case info
        self._create_case_info["machine"] = self._set_host_property(host_info, "machine")
        self._create_case_info["compiler"] = self._set_host_property(host_info, "compiler")

        # env_build info
        if self._userdefined is True:
            self._env_build_info["os"] = self._set_host_property(host_info, "os")
            self._env_build_info["mpilib"] = self._set_host_property(host_info, "mpi_vendor")
            self._env_build_info["compiler"] = self._set_host_property(host_info, "compiler")
            self._env_build_info["gmake"] = self._set_host_property(host_info, "gmake")

        # env_mach_pes info
        if self._userdefined is True:
            self._env_mach_pes_info["max_np"] = self._set_host_property(host_info, "max_np")

        # macros info
        if self._userdefined is True:
            self._macros_info["os"] = self._set_host_property(host_info, "os")
            self._macros_info["cc"] = self._set_host_property(host_info, "cc")
            self._macros_info["cxx"] = self._set_host_property(host_info, "cxx")
            self._macros_info["fc"] = self._set_host_property(host_info, "fc")
            self._macros_info["mpicc"] = self._set_host_property(host_info, "mpicc")
            self._macros_info["mpicxx"] = self._set_host_property(host_info, "mpicxx")
            self._macros_info["mpifc"] = self._set_host_property(host_info, "mpifc")
            self._macros_info["mpiexec"] = self._set_host_property(host_info, "mpiexec")
            self._macros_info["mpi_vendor"] = self._set_host_property(host_info, "mpi_vendor")
            self._macros_info["fflags"] = self._set_host_property(host_info, "fflags")
            self._macros_info["slibs"] = self._set_host_property(host_info, "slibs")
            self._macros_info["netcdf_path"] = self._set_host_property(host_info, "netcdf_path")
            self._macros_info["blas_flags"] = self._set_host_property(host_info, "blas_flags")
            self._macros_info["ldflags"] = self._set_host_property(host_info, "ldflags", required=False)
            if host_info.has_key("ice_hack"):
                self._macros_info["ice_hack"] = self._set_host_property(host_info, "ice_hack")
#        self._ = self._set_host_property(host_info, "")

    def _set_host_property(self, host_info, prop_name, required=True):
        value = None
        if host_info.has_key(prop_name):
            value = host_info[prop_name]
        else:
            if required:
                raise Exception("ERROR: host info for '{0}' must include field '{1}'".format(
                    self._builder_id, prop_name))
            else:
                value = ""
        return value

class RegressionTest(object):
    """Class to setup and run a coupled CLM-PFLOTRAN problem then compare
    the results.

    It also generates a stand alone shell script that should exactly
    recreate the test case using only the CESM scripts tools
    (create_newcase, cesm_setup, xmlchange, etc). Does not rely on
    external tools like runCLM.py

    Test information is read from a standard configuration file with
    format below. Some key things to note:

    * The case name is the file name with cfg removed

    * cases are placed in the root of the cesm directory (same level as scripts and models)

    * the build and run directories are placed in ${CASE_DIR}/bld and
      ${CASE_DIR}/run respectively

    * executable is either 'build' or the name of the case which contains a reusable executable

    * user_nl_XXX sections keword value pairs for user_nl_XXX sections
      should be exactly what would go into the user_nl file if you
      were creating it yourself.

        * Strings in user_nl_XXX sections should be properly quoted.

        * The only exception is paths, which must be unquoted and will
          have the correct absolute path for the machine and test
          environment prepended.

    * The clm_config_opts and clm_bldnl_opts are APPENDED to the
      existing values rather than overwritting.

    * Strings in xml sections should NOT be quoted


    [case]
    resolution = CLM_USRDAT
    compset = I1850CLM45CN
    np = 1
    executable = build
    
    [site_data]
    site = US-Brw
    data_dir = scripts/ccsm_utils/Tools/lnd/clm/PTCLM/PTCLM_sitedata/
    site_group = AmeriFlux
    x_pts = 1
    y_pts = 1
    
    [user_nl_clm]
    use_pflotran = .true.
    finidat = lnd/clm2/initdata/US-Brw_I1850CLM45CN.clm2.r.0601-01-01-00000.nc
    fsurdat = lnd/clm2/surfdata_map/surfdata_1x1pt_US-Brw_simyr1850_ugrid_c131015.nc
    paramfile = lnd/clm2/paramdata/clm_params.c130821.US-Brw.nc
    
    [user_nl_datm] 
    &shr_strdata_nml taxmode = 'cycle, extend'
    
    [env_build]
    clm_config_opts = -clm4me on -vsoilc_centbgc no-cent -nofire
    
    [env_run]
    run_refdate = 0001-05-01
    stop_option = ndays
    stop_n = 2
    clm_bldnml_opts = -mask navy
    datm_mode = CLM1PT
    datm_clmncep_yr_start = 1998
    datm_clmncep_yr_end = 2006
    datm_clmncep_yr_align = 1
    
    [pflotran]
    mesh_dir = sgrid-1x1

    """

    def __init__(self):
        self._logfile = None
        self._case_script = None
        self._case_script_filename = None
        self._name = None
        self._userdefined = False
        self._cesm_root_dir = None
        self._case_info = {}
        self._regression_info = {}
        self._site_data = {}
        self._clm_namelist_info = {}
        self._datm_namelist_info = {}
        self._datm_streams_info = {}
        self._env_run_info = {}
        self._env_build_info = {}
        self._env_mach_pes_info = {}
        self._use_pflotran = False
        self._user_data_name = None
        self._pflotran_info = {}
        self._local_config = {}
        self._macros_info = {}
        self._error_logs = []
        self._timeout = 600.0
        self._txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
        self._status = {"create" : False,
                        "configure" : False,
                        "build" : False,
                        "run" : False,
                        "test" : False,}


    def __str__(self):
        txt = "  Test : {0}\n".format(self._name)
        txt += "    cesm_root_dir :\n"
        txt += "      {0}\n".format(self._cesm_root_dir)
        txt += "    case :\n"
        for p in self._case_info:
            txt += "      {0} : {1}\n".format(p, self._case_info[p])

        txt += "    site data :\n"
        for p in self._site_data:
            txt += "      {0} : {1}\n".format(p, self._site_data[p])

        txt += "    clm namelist :\n"
        for p in self._clm_namelist_info:
            txt += "      {0} : {1}\n".format(p, self._clm_namelist_info[p])

        txt += "    datm namelist :\n"
        for p in self._datm_namelist_info:
            txt += "      {0} : {1}\n".format(p, self._datm_namelist_info[p])

        txt += "    env_run :\n"
        for p in self._env_run_info:
            txt += "      {0} : {1}\n".format(p, self._env_run_info[p])

        txt += "    env_build :\n"
        for p in self._env_build_info:
            txt += "      {0} : {1}\n".format(p, self._env_build_info[p])

        txt += "    env_mach_pes :\n"
        for p in self._env_mach_pes_info:
            txt += "      {0} : {1}\n".format(p, self._env_mach_pes_info[p])

        txt += "    pflotran :\n"
        for p in self._pflotran_info:
            txt += "      {0} : {1}\n".format(p, self._pflotran_info[p])

        txt += "    local info :\n"
        for p in self._local_config:
            txt += "      {0} : {1}\n".format(p, self._local_config[p])

        txt += "    testing info :\n"
        for p in self._regression_info:
            txt += "      {0} : {1}\n".format(p, self._regression_info[p])

        return txt

    def name(self):
        return self._name

    def setup(self, config_filename, machine, local_config, cesm_root_dir, suitelog):
        """Setup the test object by reading the config files and combining
        machine and config info into the appropriate data
        structures. Also does as much error checking as possible along the way.

        """
        filename = os.path.basename(config_filename)
        self._name = filename[: filename.rfind(".cfg")]
        print("CLM-PFLOTRAN Regression Test : {0}".format(self.name()), file=suitelog)
        self._cesm_root_dir = cesm_root_dir
        self._case_root_dir = "{0}/test_cases".format(os.path.abspath(self._cesm_root_dir))
        try:
            os.makedirs(self._case_root_dir)
        except:
            pass
        self._case_dir = "{0}/{1}".format(self._case_root_dir, self._name)

        self._read_test_config_file(config_filename, suitelog)
        self._create_logfile()
        self._add_machine_info_to_test(machine)
        self._check_local_config(local_config, suitelog)
        self._initialize_case_script()
        # FIXME(bja, 2013-10) need to link_dirtree if necessary....

    def run_test(self, suitelog):
        status = self._create_case(suitelog)
        if status != 0:
            raise RuntimeError("{0} : create case failed.".format(self.name()))
        status = self._configure_case(suitelog)
        if status != 0:
            raise RuntimeError("{0} : configure case failed.".format(self.name()))
        status = self._build_case(suitelog)
        if status != 0:
            raise RuntimeError("{0} : build case failed.".format(self.name()))
        status = self._run_case(suitelog)
        if status != 0:
            raise RuntimeError("{0} : run case failed.".format(self.name()))
        status = self._test_case(suitelog)
        if status != 0:
            raise RuntimeError("{0} : test case failed.".format(self.name()))
        self._finalize_case()

    def _create_case(self, suitelog):
        """Create the case by running the create_newcase command
        """
        print("Creating case...", end='', file=suitelog)
        if os.path.isdir(self._case_dir):
            old_dir = "{0}.old".format(self._case_dir)
            if os.path.isdir(old_dir):
                shutil.rmtree(old_dir)
            os.rename(self._case_dir, old_dir)

        scripts_dir = "{0}/{1}".format(self._cesm_root_dir, "scripts")
        os.chdir(scripts_dir)
        print("cd {0}".format(scripts_dir), file=self._case_script)
        cmd = []
        cmd.append("./create_newcase")
        cmd.append("-case")
        cmd.append("{0}".format(self._case_dir))
        cmd.append("-res")
        cmd.append(self._case_info["resolution"])
        cmd.append("-compset")
        cmd.append(self._case_info["compset"])
        cmd.append("-mach")
        cmd.append(self._case_info["machine"])
        cmd.append("-compiler")
        cmd.append(self._case_info["compiler"])
        print("# Creating case with command :", file=self._case_script)
        status = self._run_command(cmd)
        if status == 0:
            self._status["create"] = True
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)
        return status


    def _configure_case(self, suitelog):
        """Configure the case by creating a machine specific macros file,
        modifies the env_*.xml files, user_nl_* files, and calls
        cesm_setup.

        """
        print("Configuring case :", file=suitelog)
        print("# Configuring case :", file=self._case_script)
        os.chdir(self._case_dir)
        print("cd {0}".format(self._case_dir), file=self._case_script)
        status = 0
        if self._userdefined:
            status += self._create_macros(suitelog)
            status += self._modify_env_mach_pes(suitelog)
        status += self._modify_env_build(suitelog)
        status += self._modify_env_run(suitelog)
        status += self._run_cesm_setup(suitelog)
        if self._userdefined:
            status += self._modify_run_script(suitelog)
        # NOTE: these must come *after* cesm_setup is called:
        status += self._create_clm_user_namelist(suitelog)
        status += self._create_datm_user_namelist(suitelog)
        status += self._modify_datm_streams(suitelog)
        if self._use_pflotran is True:
            status += self._copy_pflotran_files(suitelog)
        self._case_script.flush()
        # FIXME(bja, 2013-10) search for error messages in configure
        if status == 0:
            self._status["configure"] = True
        return status

    def _build_case(self, suitelog):
        """Run the case.build command
        """
        print("Building case...", end='', file=suitelog)
        print("# building case :", file=self._case_script)
        build_error_re = re.compile("^[\s]*ERROR: cat (/.+[\d]{6}-[\d]{6})[\s]*$")
        status = 0
        if self._case_info["executable"] == "build":
            # build the case
            cmd = ["./{0}.build".format(self._name)]
            status += self._run_command(cmd)
            # FIXME(bja, 2013-10) search for error messages in build
            with open(self._logfile, 'r') as testlog:
                for line in testlog:
                    build_error = build_error_re.match(line)
                    if build_error:
                        print("\n\n  Found build error log :", file=suitelog)
                        print(build_error.group(1), file=suitelog)
                        self._error_logs.append(build_error.group(1))
            if status == 0:
                self._status["build"] = True
        else:
            # assume the executable keyword points to a case name that has an executable
            exe_case_dir = "{0}/{1}".format(self._case_root_dir, self._case_info["executable"])
            if not os.path.isdir(exe_case_dir):
                raise Exception("ERROR: could not find case dir for executable reuse:\n{0}".format(exe_case_dir))

            exe_path = "{0}/bld/cesm.exe".format(exe_case_dir)
            if not os.path.isfile(exe_path):
                raise Exception("ERROR: could not find executable for reuse:\n{0}".format(exe_path))
            status += self._run_xml_change("env_build.xml", "EXEROOT", "{0}/bld".format(exe_case_dir))
            status += self._run_xml_change("env_build.xml", "BUILD_COMPLETE", "TRUE")
            print("#./{0}.build".format(self._name), file=self._case_script)
        if status == 0:
            self._status["build"] = True
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)
        return status

    def _run_case(self, suitelog):
        """Run the case.run command. The sleep command should have been
        commented out during the configure step.

        """
        print("Running case...", end='', file=suitelog)
        print("# Running case :", file=self._case_script)
        run_script = "{0}.run".format(self._name)
        cmd = ["./{0}".format(run_script)]
        status = self._run_command(cmd)
        status += self._check_test_runtime_error(suitelog)
        if status == 0:
            self._status["run"] = True
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)
        return status

    def _test_case(self, suitelog):
        """Compare the test results....

        Using quick and dirty comparison tool:
        ./quick-compare-nc.py -b baseline.nc -c compare.nc -t "1.0e-7 absolute" -f FSH

        """
        print("Testing case...", end='', file=suitelog)
        baseline_filename = "{0}/{1}".format(self._local_config["regression_dir"],
                                             self._regression_info["file"])
        current_filename = "{0}/run/{1}".format(self._case_dir, self._regression_info["file"])
        status = 0
        for field in self._regression_info:
            if field != "file":
                cmd = ["{0}/clm4-pf-tools/regression_tests/quick-compare-nc.py".format(self._cesm_root_dir),
                       "-b", baseline_filename,
                       "-c", current_filename,
                       "-t", "{0}".format(self._regression_info[field]), # need to quote tolerance info!
                       "-f", field,
                   ]
                #print("\n\t{0}".format(' '.join(cmd)), file=suitelog)
                status += self._run_command(cmd)
                
        if status == 0:
            self._status["test"] = True
            print(" passed.", file=suitelog)
        else:
            print(" failed.", file=suitelog)
        return status

    def _finalize_case(self):
        """Do any cleanup necessary (e.g. close the case script file)
        """
        self._finalize_case_script()

    def status(self, suitelog):
        """Report the overall test status.
        """
        print("Test '{0}' status :\n    {1}".format(self.name(), self._status), file=suitelog)
        stat = 0
        for s in self._status:
            if self._status[s] is False:
                stat += 1
    
        if stat > 0:
            print("F", end='')
            # dump case info as present in the current object
            print(self, file=suitelog)
            # dump the case generation script
            print("\n\n", file=suitelog)
            print("*** Warning: this script may be incorrect/incomplete!", file=suitelog)
            with open(self._case_script_filename, 'r') as tmp:
                for line in tmp:
                    print(line, end='', file=suitelog)
            # dump the case generation log
            print("\n\n", file=suitelog)
            with open("{0}/test_cases/{1}.log".format(self._cesm_root_dir, self.name()), 'r') as tmp:
                for line in tmp:
                    print(line, end='', file=suitelog)
            # dump any error logs identified during previous steps
            print("\n\n", file=suitelog)
            for log in self._error_logs:
                if not os.path.isfile(log):
                    print("Error log identified in previous stage could not be found: ", file=suitelog)
                    print(log, file=suitelog)
                else:
                    print("Dumping error log file: ", file=suitelog)
                    print(log, file=suitelog)
                    with open(log, 'r') as tmp:
                        for line in tmp:
                            print(line, end='', file=suitelog)

        else:
            print(".", end='')
        sys.stdout.flush()

        # only want to report a single pass/fail for each test:
        summary = 0
        if stat > 0:
            summary = 1
        return summary

    def _check_test_runtime_error(self, suitelog):
        """Check for runtime error messages in standard out
        """
        status = 0
        status += self._check_test_did_not_complete(suitelog)
        status += self._check_test_petsc_error(suitelog)
        return status

    def _check_test_did_not_complete(self, suitelog):
        """Check for "Model did not complete" error in the run standard output (run logfile)
        """
        status = 0
        did_not_complete_re = re.compile("^Model did not complete - see (.+)$")
        with open(self._logfile, 'r') as run_stdout:
            for line in run_stdout.readlines():
                match = did_not_complete_re.match(line)
                if match:
                    status += 1
                    print("ERROR: check_test_did_not_complete(): model did not complete regex found a match!", file=suitelog)
                    print("    {0}".format(match.group(0)), file=suitelog)
        return status

    def _check_test_petsc_error(self, suitelog):
        """Check the cesm log file for petsc error messages.
        """
        status = 0
        petsc_error_re = re.compile("PETSC ERROR: ")
        petsc_error = False
        run_filenames = os.listdir("{0}/run".format(self._case_dir))
        cesm_log = ''
        for f in run_filenames:
            if f.startswith("cesm.log"):
                cesm_log = f
        if cesm_log is '':
            print("check_test_petsc_error(): Could not find cesm.log file in run directory!", file=suitelog)
            petsc_error = True
        cesm_log = "{0}/run/{1}".format(self._case_dir, cesm_log)
        log_data = None
        if cesm_log.endswith(".gz"):
            with gzip.open(cesm_log, "rb") as f:
                log_data = f.readlines()
        else:
            with open(cesm_log, 'r') as f:
                log_data = f.readlines()
        if log_data is not None:
            for line in log_data:
                match = petsc_error_re.search(line)
                if match:
                    status += 1
                    petsc_error = True

        if petsc_error == True:
            print("ERROR: check_test_petsc_error(): petsc error regex found a match in : {0}!".format(cesm_log), file=suitelog)
        return status

    def _create_logfile(self):
        # TODO(bja, 2013-10) add repo info to top of log file
        self._logfile = "{0}/{1}.log".format(self._case_root_dir, self._name)
        if os.path.isfile(self._logfile):
            os.rename(self._logfile, self._logfile + ".old")

    def _initialize_case_script(self):
        # TODO(bja, 2013-10) add repo info to top of script
        now = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
        self._case_script_filename = "{0}/{1}-{2}.sh".format(self._case_root_dir, self._name, now)
        self._case_script = open(self._case_script_filename, 'w')
        print("#/usr/bin/env bash", file=self._case_script)
        print("#", file=self._case_script)
        print("# Script automatically generate by clm-pflotran-tests.py", file=self._case_script)

    def _finalize_case_script(self):
        self._case_script.close()

    def _add_machine_info_to_test(self, machine):
        """Merge machine specific information into the test information
        """
        self._userdefined = machine._userdefined
        self._case_info.update(machine._create_case_info)
        self._env_build_info.update(machine._env_build_info)
        self._env_run_info.update(machine._env_run_info)
        self._macros_info.update(machine._macros_info)
        self._env_mach_pes_info.update(machine._env_mach_pes_info)

        # check for required parameters.
        txt = ""

        if txt is not "":
            raise Exception(txt)

    def _run_command(self, command):
        """Generic function to run a shell command, with timout limit, and
        append output to a log file.

        """
        print(" ".join(command), file=self._case_script)
        try:
            with open(self._logfile, 'a') as run_stdout:
                start = time.time()
                proc = subprocess.Popen(command,
                                        shell=False,
                                        stdout=run_stdout,
                                        stderr=subprocess.STDOUT)
                while proc.poll() is None:
                    time.sleep(0.1)
                    if time.time() - start > self._timeout:
                        proc.kill()
                        time.sleep(0.1)
                        message = self._txtwrap.fill(
                            "ERROR: job '{0}' has exceeded timeout limit of "
                            "{1} seconds.".format(self.name(), self._timeout))
                        print(''.join(['\n', message, '\n']), file=testlog)
                cmd_status = abs(proc.returncode)
        except Exception as e:
            print("ERROR: Running command :\n    '{0}'".format(" ".join(command)))
            print(e)
            cmd_status = 1
        return cmd_status

    def _read_test_config_file(self, filename, suitelog):
        """Read the test specific configuration file and extract the various
        sections.

        Required values should be checked in these sections and return
        errors if they are not found.

        """
        if filename is None:
            raise Exception("ERROR: RegressionTest.read_config_file: must provide a config filename")
        self._config_filename = filename
        config = config_parser.SafeConfigParser()
        config.read(self._config_filename)

        self._check_case_info(config, suitelog)
        self._check_site_data(config, suitelog)
        self._check_clm_namelist_info(config, suitelog)
        self._check_datm_namelist_info(config, suitelog)
        self._check_datm_streams_info(config, suitelog)
        self._check_env_build_info(config, suitelog)
        self._check_env_run_info(config, suitelog)
        self._check_pflotran_input(config, suitelog)
        self._check_regression_input(config, suitelog)

    def _check_for_param(self, param_dict, name, section, default=None):
        """Generic routine to check for parameters in a configuration file
        section dictionary.

        Missing required parameters are reported as errors.

        Defaults are added for optional params.

        """
        txt = ""
        if not param_dict.has_key(name):
            if default is not None:
                param_dict[name] = default
            else:
                txt += ("ERROR: test configuration '{0}': '{1}' section "
                        "must have a '{2}' field\n".format(
                            self._config_filename, section, name))
        return txt

    def _check_case_info(self, config, suitelog):
        section = "case"
        if not config.has_section(section):
            raise Exception("ERROR: test config file '{0}' does not contain a"
                            "'{1}' section.".format(self._config_filename, section))
        # NOTE(bja, 2014-02) must set raw=True becaues the compset can
        # contain '%' characters, and ConfigParser will try to
        # interpolate them.
        self._case_info = list_to_dict(config.items(section, raw=True))
        self._case_info["name"] = self._name
        # check for required parameters.
        txt = ""
        txt += self._check_for_param(self._case_info, "compset", section)
        txt += self._check_for_param(self._case_info, "resolution", section)
        txt += self._check_for_param(self._case_info, "np", section, default="1")
        txt += self._check_for_param(self._case_info, "executable", section)

        if txt is not "":
            raise Exception(txt)

    def _check_site_data(self, config, suitelog):
        section = "site_data"
        if not config.has_section(section):
            print("NOTE: test config file '{0}' does not contain a"
                            "'{1}' section.".format(self._config_filename, section), file=suitelog)
        else:
            self._site_data = list_to_dict(config.items(section))
            # check for required parameters.
            txt = ""
            txt += self._check_for_param(self._site_data, "site", section)
            txt += self._check_for_param(self._site_data, "data_dir", section)
            txt += self._check_for_param(self._site_data, "site_group", section)
            txt += self._check_for_param(self._site_data, "x_pts", section)
            txt += self._check_for_param(self._site_data, "y_pts", section)
    
            if txt is not "":
                raise RuntimeError(txt)

            self._user_data_name = "{0}x{1}pt_{2}".format(self._site_data['x_pts'],
                                                          self._site_data['y_pts'],
                                                          self._site_data['site'])

            #data_sets = ["sitedata", "pdtdata", "soildata"]
            data_sets = ["sitedata"]
            # NOTE (bja, 20131015) it looks like all the site groups use the same csv formats
            for d in data_sets:
                filename = "{0}/{1}/{2}_{3}.txt".format(self._cesm_root_dir, self._site_data["data_dir"], self._site_data["site_group"], d)
                if not os.path.isfile(filename):
                    raise RuntimeError("Could not find site data file : {0}".format(filename))
                with open(filename, 'r') as csv_file:
                    csvdata = csv.reader(csv_file, delimiter=',', quotechar='"')
                    for line in csvdata:
                        if line[0] == self._site_data["site"]:
                            if d == "sitedata":
                                self._site_data["name"] = line[1]
                                self._site_data["state"] = line[2]
                                self._site_data["lon"] = line[3]
                                self._site_data["lat"] = line[4]
                                self._site_data["elev"] = line[5]
                                self._site_data["startyear"] = line[6]
                                self._site_data["endyear"] = line[7]
                                self._site_data["alignyear"] = line[8]

    def _check_clm_namelist_info(self, config, suitelog):
        section = "user_nl_clm"
        if not config.has_section(section):
            print("NOTE: test config file '{0}' does not contain a "
                  "'{1}' section.".format(self._config_filename, section), file=suitelog)
        self._clm_namelist_info = list_to_dict(config.items(section))
        # check for required parameters.
        txt = ""
        txt += self._check_for_param(
            self._clm_namelist_info, "use_pflotran", section, default=".true.")

        if txt is not "":
            raise Exception(txt)

    def _check_datm_namelist_info(self, config, suitelog):
        section = "user_nl_datm"
        if not config.has_section(section):
            print("NOTE: test config file '{0}' does not contain a "
                  "'{1}' section.".format(self._config_filename, section), file=suitelog)
        else:
            self._datm_namelist_info = list_to_dict(config.items(section))
            # check for required parameters.
            txt = ""
    
            if txt is not "":
                raise Exception(txt)

    def _check_datm_streams_info(self, config, suitelog):
        section = "datm_streams"
        if not config.has_section(section):
            print("NOTE: test config file '{0}' does not contain a "
                  "'{1}' section.".format(self._config_filename, section), file=suitelog)
        else:
            self._datm_streams_info = list_to_dict(config.items(section))
            # check for required parameters.
            txt = ""
    
            if txt is not "":
                raise Exception(txt)

    def _check_env_build_info(self, config, suitelog):
        section = "env_build"
        if not config.has_section(section):
            print("NOTE: test configuration '{0}' does not contain an "
                  "'{1}' section.".format(self._name, section), file=suitelog)
            return
        self._env_build_info = list_to_dict(config.items(section))
        # check for required parameters.
        txt = ""

        if txt is not "":
            raise Exception(txt)

    def _check_env_run_info(self, config, suitelog):
        section = "env_run"
        if not config.has_section(section):
            print("NOTE: test configuration '{0}' does not contain an "
                  "'{1}' section.".format(self._name, section), file=suitelog)
            return
        self._env_run_info = list_to_dict(config.items(section))
        # check for required parameters.
        txt = ""

        if txt is not "":
            raise Exception(txt)

    def _check_pflotran_input(self, config, suitelog):
        section = "pflotran"
        if self._clm_namelist_info["use_pflotran"] == ".true.":
            #print("Setting use_pflotran = True")
            self._use_pflotran = True
            if not config.has_section(section):
                msg = ("ERROR: {0} : clm_namelist sets 'use_pflotran = .true.' "
                       "but does not have a pflotran_files".format(self._name))
                raise Exception(msg)
            else:
                self._pflotran_info = list_to_dict(config.items(section))
                if not self._pflotran_info.has_key("mesh_dir"):
                    raise Exception("ERROR: {0} : pflotran mesh_dir missing".format(self._name))

                # maping files not needed because they are in the
                # pflotran input file. Probably easier to scrape that
                # file and then copy the files to the run directory
                # instead of listing here....

#XXX        else:
#XXX            print("Leaving use_pflotran == False")
#XXX            print("use_pflotran = '{0}'".format(self._clm_namelist_info["use_pflotran"]))

    def _check_regression_input(self, config, suitelog):
        """Setup data needed for regression testing physics output.
        """
        section = "regression"
        if not config.has_section(section):
            raise Exception("ERROR: test config file '{0}' does not contain a"
                            "'{1}' section.".format(self._config_filename, section))
        self._regression_info = list_to_dict(config.items(section))
        # check for required parameters.
        txt = ""
        txt += self._check_for_param(self._regression_info, "file", section)
        # must be at least one field to check
        if len(self._regression_info) < 2:
            txt += ("ERROR: test config file '{0}' '{1}' section must contain "
                    "at least one parameter!".format(self._config_filename, section))

        if txt is not "":
            raise Exception(txt)


    def _check_local_config(self, local_config, suitelog):
        self._local_config = local_config
        # check for required parameters
        txt = ""
        txt += self._check_for_param(self._local_config, "PETSC_DIR", "command-line")
        txt += self._check_for_param(self._local_config, "PETSC_ARCH", "command-line")
        txt += self._check_for_param(self._local_config, "PFLOTRAN_DIR", "command-line")
        txt += self._check_for_param(self._local_config, "data_dir", "command-line")
        txt += self._check_for_param(self._local_config, "regression_dir", "command-line")

        if txt is not "":
            raise Exception(txt)

        txt = ""
        petsc_dir = self._local_config["PETSC_DIR"]
        if not os.path.isabs(petsc_dir):
            # if we received a relative path, make it absolute
            petsc_dir = os.path.abspath(petsc_dir)
            self._local_config["PETSC_DIR"] = petsc_dir
        if not os.path.isdir(petsc_dir):
            txt += "ERROR: specified PETSC_DIR is not a directory : {0}".format(petsc_dir)
        
        petsc_arch = "{0}/{1}".format(petsc_dir, self._local_config["PETSC_ARCH"])
        if not os.path.isdir(petsc_arch):
            txt += "ERROR: specified PETSC_ARCH is not a directory : {0}".format(petsc_arch)
        
        libpetsc_static = "{0}/lib/libpetsc.a".format(petsc_arch)
        libpetsc_dynamic = "{0}/lib/libpetsc.dylib".format(petsc_arch)

        if not os.path.isfile(libpetsc_static) and not os.path.isfile(libpetsc_dynamic):
            txt += "ERROR: could not find libpetsc in PETSC_DIR/PETSC_ARCH/lib!"

        pflotran_dir = self._local_config["PFLOTRAN_DIR"]
        if not os.path.isabs(pflotran_dir):
            # if we received a relative path, make it absolute
            pflotran_dir = os.path.abspath(pflotran_dir)
            self._local_config["PFLOTRAN_DIR"] = pflotran_dir
        if not os.path.isdir(pflotran_dir):
            txt += "ERROR: specified PFLOTRAN_DIR is not a directory : {0}".format(pflotran_dir)
        libpflotran = "{0}/src/clm-pflotran/libpflotran.a".format(pflotran_dir)
        if not os.path.isfile(libpflotran):
            txt += "ERROR: could not find libpflotran in PFLOTRAN_DIR. Expected : {0}".format(libpflotran)

        if not os.path.isdir(self._local_config["data_dir"]):
            txt += "ERROR: specified input data directory does not exist : {0}".format(self._local_config["data_dir"])
        if not os.path.isdir(self._local_config["regression_dir"]):
            txt += "ERROR: specified regression data directory does not exist : {0}".format(self._local_config["regression_dir"])

        if txt is not "":
            raise Exception(txt)

    def _run_cesm_setup(self, suitelog):
        print("  Running cesm_setup...", end='', file=suitelog)
        status = 0
        cmd = ["./cesm_setup"]
        stat = self._run_command(cmd)
        if stat is not 0:
            status = 1

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _modify_run_script(self, suitelog):
        # note, must come after call to cesm_setup
        print("  Modifying run script...", end='', file=suitelog)
        status = 0
        try:
            run_script = "{0}.run".format(self._name)
            run_script_orig = "{0}.orig".format(run_script)
            shutil.copy(run_script, run_script_orig)
            sleep_re = re.compile("^sleep 25")
            mpiexec_re = re.compile("^#mpiexec")
            with open(run_script, 'w') as new_script, open(run_script_orig, 'r') as orig_script:
                for line in orig_script:
                    new_line = line
                    if sleep_re.match(line):
                        new_line = "#{0}".format(line)
                    if mpiexec_re.match(line):
                        new_line = "{0} {1}".format(self._macros_info["mpiexec"], line[8:])
                    new_script.write(new_line)

            # now write the same changes to the shell script
            print("# Modify run script", file=self._case_script)
            mpirun_cmd="s@#mpiexec@{0}@".format(self._macros_info["mpiexec"])
            sleep_cmd="s@^sleep@#sleep@"
            print("perl -w -i -p -e '{0}' {1}.run".format(mpirun_cmd, self._name), file=self._case_script)
            print("perl -w -i -p -e '{0}' {1}.run".format(sleep_cmd, self._name), file=self._case_script)
            print(" done.", file=suitelog)
        except Exception as error:
            status += 1
            print(error, file=suitelog)
        return status

    def _run_xml_change(self, xmlfile, identifier, value, append=False):
        ident = identifier.upper()
        cmd = ["./xmlchange", "-file", xmlfile, "-id", ident, "-val", value]
        if append is True:
            cmd.append("-append")
        #print(" ".join(cmd))
        status = self._run_command(cmd)
        if status is not 0:
            print("ERROR: xmlchange failed on command: {0}".format(" ".join(cmd)))
        return status

    def _create_macros(self, suitelog):
        print("  Creating Macros...", end='', file=suitelog)
        print("# Create Macros", file=self._case_script)
        status = 0
        try:
            macros_dict = {}
            macros_dict.update(self._macros_info)
            macros_dict.update(self._env_build_info)
            macros_dict.update(self._local_config)
            macros_dict["case_dir"] = self._case_dir
            macros_dict["petsc_include"] = "include {PETSC_DIR}/conf/variables".format(
                PETSC_DIR=self._local_config["PETSC_DIR"])
            macros_dict["user_include"] = ("-I{PETSC_DIR}/{PETSC_ARCH}/include "
                                           "-I{PETSC_DIR}/include "
                                           "-I{PFLOTRAN_DIR}/src/clm-pflotran".format(
                                               PETSC_DIR=self._local_config["PETSC_DIR"],
                                               PETSC_ARCH=self._local_config["PETSC_ARCH"],
                                               PFLOTRAN_DIR=self._local_config["PFLOTRAN_DIR"]))
            macros_dict["user_fflags"] = "-DCLM_PFLOTRAN"
            macros_dict["user_ldflags"] = ("{PFLOTRAN_DIR}/src/clm-pflotran/libpflotran.a "
                                           "-L{PETSC_DIR}/{PETSC_ARCH}/lib "
                                           "$(PETSC_LIB) -lnetcdff -lnetcdf".format(
                                               PETSC_DIR=self._local_config["PETSC_DIR"],
                                               PETSC_ARCH=self._local_config["PETSC_ARCH"],
                                               PFLOTRAN_DIR=self._local_config["PFLOTRAN_DIR"]))
            macros_dict["compiler_uppercase"] = self._case_info["compiler"].upper()
            macros_dict["ice_hack"] = ""
            if self._macros_info.has_key("ice_hack"):
                macros_dict["ice_hack"] = "{0}/{1}".format(self._case_dir, self._macros_info["ice_hack"])


            macros_str = macros_template.substitute(macros_dict)
            with open("Macros", 'w') as macros:
                macros.write(macros_str)
            print("cat > Macros <<\EOF", file=self._case_script)
            print(macros_str, file=self._case_script)
            print("EOF", file=self._case_script)
            print(" done.", file=suitelog)
        except Exception as error:
            status = 1
            print(error, file=suitelog)
        return status

    def _modify_env_build(self, suitelog):
        xmlfile = "env_build.xml"
        print("  Modifying {0}...".format(xmlfile), end='', file=suitelog)
        print("# Modifying : {0}".format(xmlfile), file=self._case_script)
        status = 0
        # set any machine/problem specific flags 
        for i in self._env_build_info:
            if i.lower() == "clm_config_opts":
                status += self._run_xml_change(xmlfile, i, self._env_build_info[i], append=True)
            else:
                status += self._run_xml_change(xmlfile, i, self._env_build_info[i])

        # set universal flags
        status += self._run_xml_change(xmlfile, "DEBUG", "TRUE")
        status += self._run_xml_change(xmlfile, "SUPPORTED_BY", "'clm-pflotran test case'")
        bld_dir = "{0}/{1}".format(os.getcwd(), "bld")
        os.mkdir(bld_dir)
        status += self._run_xml_change(xmlfile, "EXEROOT", bld_dir)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _modify_env_run(self, suitelog):
        xmlfile = "env_run.xml"
        print("  Modifying {0}...".format(xmlfile), end='', file=suitelog)
        print("# Modifying : {0}".format(xmlfile), file=self._case_script)
        status = 0
        # set any machine/problem specific flags 
        for i in self._env_run_info:
            if i.lower() == "datm_mode":
                if (self._env_run_info[i] == "CLM1PT" and not self._site_data.has_key('site')):
                    status += 1
                    print("ERROR: modify_env_run : datm_mode = CLM1PT requires site data.", file=self._logfile)
                status += self._run_xml_change(xmlfile, i.upper(), self._env_run_info[i])
                # if the user didn't override the start/stop/align
                # site data in the env_run block, we have to write it
                # here!
                if not self._env_run_info["datm_clmncep_yr_start"]:
                    status += self._run_xml_change(xmlfile, "DATM_CLMNCEP_YR_START", self._site_data["startyear"])
                if not self._env_run_info["datm_clmncep_yr_end"]:
                    status += self._run_xml_change(xmlfile, "DATM_CLMNCEP_YR_END", self._site_data["endyear"])
                if not self._env_run_info["datm_clmncep_yr_align"]:
                    status += self._run_xml_change(xmlfile, "DATM_CLMNCEP_YR_ALIGN", self._site_data["alignyear"])
            elif i.lower() == "clm_bldnml_opts":
                status += self._run_xml_change(xmlfile, i.upper(), '"{0}"'.format(self._env_run_info[i]), append=True)
            else:
                status += self._run_xml_change(xmlfile, i.upper(), self._env_run_info[i])

        # set universal flags
        din_loc = "{0}/cesm-inputdata".format(self._local_config["data_dir"])
        status += self._run_xml_change(xmlfile, "DIN_LOC_ROOT", din_loc)
        status += self._run_xml_change(xmlfile, "DIN_LOC_ROOT_CLMFORC", "'$DIN_LOC_ROOT'")

        run_dir = "{0}/{1}".format(os.getcwd(), "run")
        os.mkdir(run_dir)
        status += self._run_xml_change(xmlfile, "rundir", run_dir)
        if self._user_data_name is not None:
            status += self._run_xml_change(xmlfile, "CLM_USRDAT_NAME", self._user_data_name)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _modify_env_mach_pes(self, suitelog):
        xmlfile = "env_mach_pes.xml"
        print("  Modifying {0}...".format(xmlfile), end='', file=suitelog)
        print("# Modifying : {0}".format(xmlfile), file=self._case_script)
        status = 0
        np = "{0}".format(min(int(self._case_info["np"]),
                              int(self._env_mach_pes_info["max_np"])))
        status += self._run_xml_change(xmlfile, "NTASKS_ATM", np)
        status += self._run_xml_change(xmlfile, "NTASKS_LND", np)
        status += self._run_xml_change(xmlfile, "NTASKS_ICE", np)
        status += self._run_xml_change(xmlfile, "NTASKS_OCN", np)
        status += self._run_xml_change(xmlfile, "NTASKS_CPL", np)
        status += self._run_xml_change(xmlfile, "NTASKS_GLC", np)
        status += self._run_xml_change(xmlfile, "NTASKS_ROF", np)
        status += self._run_xml_change(xmlfile, "NTASKS_WAV", np)
        status += self._run_xml_change(xmlfile, "MAX_TASKS_PER_NODE", np)
        status += self._run_xml_change(xmlfile, "TOTALPES", np)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _create_clm_user_namelist(self, suitelog):
        print("  Creating user_nl_clm...", end='', file=suitelog)
        print("# Modify user_nl_clm", file=self._case_script)
        status = 0
        clm_user_nl_filename = os.path.abspath("user_nl_clm")
        if not os.path.isfile(clm_user_nl_filename):
            print("ERROR: Could not find user_nl_clm file!", file=self._logfile)
            status += 1
        else:
            print("cat >> user_nl_clm << EOF", file=self._case_script)
            need_cesm_input_path = ['finidat', 'fsurdat', 'paramfile']
            with open(clm_user_nl_filename, 'a') as user_nl:
                for i in self._clm_namelist_info:
                    if i in need_cesm_input_path:
                        val = "'{0}/cesm-inputdata/{1}'".format(self._local_config["data_dir"],
                                                 self._clm_namelist_info[i])
                    else:
                        val = self._clm_namelist_info[i]
                    user_nl.write("{0} = {1}\n".format(i, val))
                    print("{0} = {1}".format(i, val), file=self._case_script)
                if self._use_pflotran:
                    user_nl.write("&clm_pflotran_inparm\n")
                    user_nl.write("    pflotran_prefix = '{0}'\n".format(self._name))
                    user_nl.write("/\n")
                    self._case_script.write("&clm_pflotran_inparm\n")
                    self._case_script.write("    pflotran_prefix = '{0}'\n".format(self._name))
                    self._case_script.write("/\n")

            print("EOF", file=self._case_script)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _create_datm_user_namelist(self, suitelog):
        print("  Creating user_nl_datm...", end='', file=suitelog)
        print("# Modify user_nl_datm", file=self._case_script)
        status = 0
        datm_user_nl_filename = os.path.abspath("user_nl_datm")
        if not os.path.isfile(datm_user_nl_filename):
            print("ERROR: Could not find user_nl_datm file!", file=self._logfile)
            status += 1
        else:
            print("cat >> user_nl_datm << EOF", file=self._case_script)
            need_cesm_input_path = []
            with open(datm_user_nl_filename, 'a') as user_nl:
                for i in self._datm_namelist_info:
                    if i in need_cesm_input_path:
                        val = "'{0}/cesm-inputdata/{1}'".format(self._local_config["data_dir"],
                                                 self._datm_namelist_info[i])
                    else:
                        val = self._datm_namelist_info[i]
                    user_nl.write("{0} = {1}\n".format(i, val))
                    print("{0} = {1}".format(i, val), file=self._case_script)

            print("EOF", file=self._case_script)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _modify_datm_streams(self, suitelog):
        print("  Modify datm streams...", end='', file=suitelog)
        status = 0
        if not self._datm_streams_info:
            print(" (N/A)", end='', file=suitelog)
        else:
            print("# Modify datm streams", file=self._case_script)
            for f in self._datm_streams_info:
                if f.find("clm1pt.clm_usrdat") != -1:
                    stream_name = "datm.streams.txt.{0}".format(f.upper())
                elif f.find("presaero") != -1:
                    stream_name = "datm.streams.txt.{0}".format(f)
                else:
                    raise RuntimeError("Unsupported datm.streams.txt type : '{0}'".format(f))

                stream_orig = "{0}/CaseDocs/{1}".format(self._case_dir, stream_name)
                stream_copy = "{0}/user_{1}".format(self._case_dir, stream_name)
                # don't want to overwrite the file if we already modified it....
                if not os.path.isfile(stream_copy):
                    shutil.copyfile(stream_orig, stream_copy)
                    print("cp {0} {1}".format(stream_orig, stream_copy), file=self._case_script)

                stream_changes = self._datm_streams_info[f].split(":")
                for i, s in enumerate(stream_changes):
                    stream_changes[i] = s.strip()

                # these files aren't valid xml according to python, so
                # we hav to modify manually....
                # assume each ';' delimited is a substitution delimited by an '@'
                for c in stream_changes:
                    change = c.split("@")
                    pattern = change[0].strip()
                    substitution = change[1].strip()
                    shell_str = "perl -w -i -p -e 's@{0}@{1}@' {2}".format(pattern, substitution, stream_copy)
                    print(shell_str, file=self._case_script)
                    with open(stream_copy, 'r') as stream_data, open("{0}.i".format(stream_copy), 'w') as stream_edit:
                        for line in stream_data:
                            line = re.sub(pattern, substitution, line)
                            stream_edit.write(line)

                shutil.copyfile("{0}.i".format(stream_copy), stream_copy)
                #print("cp {0}.i {1}".format(stream_copy, stream_copy), file=self._case_script)


        if status is 0:
            print(" done.", file=suitelog)
        else:
            print(" failed.", file=suitelog)

        return status

    def _copy_pflotran_files(self, suitelog):
        print("  Copying pflotran input files...", end='', file=suitelog)
        print("# Copying pflotran input files.", file=self._case_script)
        status = 0
        txt = ""
        data_root = self._local_config["data_dir"]
        pflotran_in = "{0}/pflotran/{1}.in".format(data_root, self._name)
        if not os.path.isfile(pflotran_in):
            txt += "ERROR: Could not find pflotran input file : {0}\n".format(pflotran_in)
            status += 1
        copy_path = "run/{0}.in".format(self._name)
        shutil.copyfile(pflotran_in, copy_path)
        print("cp {0} {1}".format(pflotran_in, copy_path), file=self._case_script)

        mesh_dir = "{0}/pflotran/meshes-maps/{1}".format(data_root, self._pflotran_info["mesh_dir"])
        if not os.path.isdir(mesh_dir):
            txt += "ERROR: Could not find pflotran mesh directory : {0}\n".format(mesh_dir)
            status += 1
        mesh_files = [ f for f in os.listdir(mesh_dir) if os.path.isfile(os.path.join(mesh_dir, f)) ]

        for f in mesh_files:
            mesh_file = "{0}/{1}".format(mesh_dir, f)
            case_copy = "run/{0}".format(f)
            shutil.copyfile(mesh_file, case_copy)
            print("cp {0} {1}".format(mesh_file, case_copy), file=self._case_script)

        if status is 0:
            print(" done.", file=suitelog)
        else:
            if txt is not "":
                print(txt, file=suitelog)
            print(" failed.", file=suitelog)

        return status

macros_template = string.Template("""#
# Makefile Macros generated from clm4_5_04/scripts/ccsm_utils/Machines/config_compilers.xml using
# COMPILER=gnu
# OS=Darwin
# MACH=userdefined
#

${petsc_include}

CPPDEFS+= -DFORTRANUNDERSCORE -DNO_R16 -D${os} -DCPR${compiler_uppercase}
#CPPDEFS+= -DFORTRANUNDERSCORE -DNO_R16 -DSYS${os}  -D${os} -DCPR${compiler}

#SLIBS+=# USERDEFINED $$(shell $$(NETCDF_PATH)/bin/nc-config --flibs)
SLIBS+=${slibs}

CONFIG_ARGS:=

CXX_LINKER:=FORTRAN

ESMF_LIBDIR:=

FC_AUTO_R8:= -fdefault-real-8

FFLAGS:= -O -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none
FFLAGS+=${fflags}

FFLAGS_NOOPT:= -O0

FIXEDFLAGS:=  -ffixed-form

FREEFLAGS:= -ffree-form

MPICC:=${mpicc}

MPICXX:=${mpicxx}

MPIFC:=${mpifc}

MPI_LIB_NAME:=

MPI_PATH:=

NETCDF_PATH:=${netcdf_path}

PNETCDF_PATH:=

SCC:=${cc}

SCXX:=${cxx}

SFC:=${fc}

SUPPORTS_CXX:=TRUE

# linking to external libraries...
USER_INCLDIR:=${user_include}
FFLAGS+=${user_fflags}
LDFLAGS +=

ifeq ($$(DEBUG), TRUE)
   FFLAGS += -g -Wall
endif

ifeq ($$(compile_threaded), true)
   LDFLAGS += -fopenmp
   CFLAGS += -fopenmp
   FFLAGS += -fopenmp
endif

ifeq ($$(MODEL), cism)
   CMAKE_OPTS += -D CISM_GNU=ON
endif

ifeq ($$(MODEL), driver)
   LDFLAGS += ${user_ldflags} ${ldflags}
   # mac os blas/lapack
   LDFLAGS += ${blas_flags}
   # NOTE(bandre): ugly hack to get around linking error on some machines...
   LDFLAGS += ${ice_hack}
endif

""")

def list_to_dict(input_list, upper_case=False):
    output_dict = {}
    for item in input_list:
        key = item[0]
        value = item[1]
        if upper_case is True:
            key = key.upper()
        output_dict[key] = value
    return output_dict

def read_local_config(filename):
    """The local configuration file contains information specific to the
    current testing environment (in contrast to the general machine
    information and test information which will generally be constant)

    The local configuration file is a standard config file with:

    [petsc]
    petsc_dir = __absolute_path_to_petsc_dir__
    petsc_arch = $PETSC_ARCH
    
    [pflotran]
    pflotran_dir = __absolute_path_to_pflotran_dir__
    
    [data]
    data_dir = __absolute_path_to_data_dir__


    """
    config = config_parser.SafeConfigParser()
    config.read(filename)
    local_config = {}

    def _check_section(conf, section):
        if not conf.has_section(section):
            raise RuntimeError("ERROR: local config file must contain a '{0}' section".format(section))

    def _get_section_option(conf, section, option, upper_case=False):
        sect = list_to_dict(conf.items(section))
        if not sect.has_key(option) and not sect.has_key(option.upper()):
            raise RuntimeError("ERROR: local config section '{0}' must contain"
                            "a '{1}' keyword.".format(section, option))
        if option.rfind("dir") != -1:
            # make sure directory paths are absolute paths
            if not os.path.isabs(sect[option]):
                raise RuntimeError("ERROR: local config file : '{0}' directory paths must be absolute paths.".format(option))

        opt = {option : sect[option]}
        if upper_case is True:
            opt = {option.upper() : sect[option]}

        return opt

    section = "petsc"
    _check_section(config, section)
    keys = ["petsc_dir", "petsc_arch"]
    for k in keys:
        local_config.update(_get_section_option(config, section, k, upper_case=True))

    section = "pflotran"
    _check_section(config, section)
    key = "pflotran_dir"
    local_config.update(_get_section_option(config, section, key, upper_case=True))

    section = "data"
    _check_section(config, section)
    keys = ["data_dir"]
    for k in keys:
        local_config.update(_get_section_option(config, section, k))

    section = "regression"
    _check_section(config, section)
    keys = ["regression_dir"]
    for k in keys:
        local_config.update(_get_section_option(config, section, k))

    return local_config

def clear_pflotran_environment(local_config):
    """We always want use the petsc and pflotran environment variables
    specified by the user local config and never use variables from
    the environment. 

    NOTE: the build systems require the environment variables so we
    overwrite them instead of nullifying them.

    """
    env = ["PETSC_DIR", "PETSC_ARCH", "PFLOTRAN_DIR"]
    for i in env:
        os.environ[i] = local_config[i]
        
def append_command_to_log(command, testlog, tempfile):
    """
    Append the results of a shell command to the test log
    """
    print("$ {0}".format(" ".join(command)), file=testlog)
    testlog.flush()
    with open(tempfile, "w") as tempinfo:
        subprocess.call(command, shell=False,
                        stdout=tempinfo,
                        stderr=subprocess.STDOUT)
        # NOTE(bja) 2013-06 : need a short sleep to ensure the
        # contents get written...?
        time.sleep(0.01)
    with open(tempfile, 'r') as tempinfo:
        shutil.copyfileobj(tempinfo, testlog)


def setup_test_suite_log(txtwrap, local_config):
    """
    Create the test suite log and try to add some useful information about
    the environment, petsc and pflotran.
    """
    now = datetime.datetime.today().strftime("%Y-%m-%d_%H-%M-%S")
    filename = "clm-pflotran-tests-{0}.testlog".format(now)
    testlog = open(filename, 'w')
    print("  Test log file : {0}".format(filename))

    print_section_seperator("CLM-PFLOTRAN Regression Test Log", testlog)
    print("Date : {0}".format(now), file=testlog)
    print("System Info :", file=testlog)
    print("    platform : {0}".format(sys.platform), file=testlog)
    test_dir = os.getcwd()
    print("Test directory : ", file=testlog)
    print("    {0}".format(test_dir), file=testlog)

    tempfile = "{0}/tmp-pflotran-regression-test-info.txt".format(test_dir)

    print_section_seperator("Provenance :", testlog)

    print("\nCLM-PFLOTRAN repository status :", file=testlog)
    print("--------------------------------", file=testlog)
    if os.path.isdir("{0}/../../.hg".format(test_dir)):
        cmd = ["hg", "parent"]
        append_command_to_log(cmd, testlog, tempfile)
        cmd = ["hg", "status", "-q"]
        append_command_to_log(cmd, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        print("    unknown", file=testlog)

    print("\nCLM-PFLOTRAN data repository status :", file=testlog)
    print("-------------------------------------", file=testlog)
    repo = local_config["data_dir"]
    if os.path.isdir("{0}/.hg".format(repo)):
        os.chdir(repo)
        cmd = ["hg", "parent"]
        append_command_to_log(cmd, testlog, tempfile)
        cmd = ["hg", "status", "-q"]
        append_command_to_log(cmd, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        print("    unknown", file=testlog)

    print("\nCLM-PFLOTRAN regression repository status :", file=testlog)
    print("-------------------------------------------", file=testlog)
    repo = local_config["regression_dir"]
    if os.path.isdir("{0}/.hg".format(repo)):
        os.chdir(repo)
        cmd = ["hg", "parent"]
        append_command_to_log(cmd, testlog, tempfile)
        cmd = ["hg", "status", "-q"]
        append_command_to_log(cmd, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        print("    unknown", file=testlog)

    print("\nPFLOTRAN repository status :", file=testlog)
    print("----------------------------", file=testlog)
    repo = local_config["PFLOTRAN_DIR"]
    if os.path.isdir("{0}/.hg".format(repo)):
        os.chdir(repo)
        cmd = ["hg", "parent"]
        append_command_to_log(cmd, testlog, tempfile)
        cmd = ["hg", "status", "-q"]
        append_command_to_log(cmd, testlog, tempfile)
        print("\n\n", file=testlog)
    else:
        print("    unknown", file=testlog)

    print("PETSc information :", file=testlog)
    print("-------------------", file=testlog)
    petsc_dir = local_config["PETSC_DIR"]
    if petsc_dir:
        print("    PETSC_DIR : {0}".format(petsc_dir), file=testlog)
        petsc_arch = local_config["PETSC_ARCH"]
        if petsc_arch:
            print("    PETSC_ARCH : {0}".format(petsc_arch), file=testlog)

        os.chdir(petsc_dir)
        print("    petsc repository status :", file=testlog)
        if os.path.isdir("{0}/.git".format(petsc_dir)):
            cmd = ["git", "log", "-1", "HEAD"]
            append_command_to_log(cmd, testlog, tempfile)
            cmd = ["git", "status", "-u", "no"]
            append_command_to_log(cmd, testlog, tempfile)
        elif os.path.isdir("{0}/.hg".format(petsc_dir)):
            cmd = ["hg", "parent"]
            append_command_to_log(cmd, testlog, tempfile)
            cmd = ["hg", "status", "-q"]
            append_command_to_log(cmd, testlog, tempfile)
        else:
            print("    No git or hg directory was found in your PETSC_DIR",
                  file=testlog)
        os.chdir(test_dir)
        print("\n\n", file=testlog)
    else:
        print("    PETSC_DIR was not defined.", file=testlog)

    os.remove(tempfile)
    print(80 * '-', file=testlog)
    return testlog

def print_section_seperator(section, outfile):
    print(80 * '-', file=outfile)
    print(section, file=outfile)

def summary_report(run_time, report, outfile):
    """
    Overall summary of test results
    """
    print(70 * '-', file=outfile)
    print("Regression test summary:", file=outfile)
    print("    Total run time: {0:4g} [s]".format(run_time), file=outfile)
    test_count = 0
    num_failures = 0
    for filename in report:
        test_count += 1
        num_failures += report[filename]

    print("    Tests run : {0}".format(test_count), file=outfile)

    success = True
    if num_failures > 0:
        print("    Failed : {0}".format(num_failures), file=outfile)
        success = False

    if success:
        print("    All tests passed.", file=outfile)

    print("\n", file=outfile)
    return num_failures


def print_no_baseline_warning(testlog):
    print(80 * '*')
    print("*")
    print("* WARNING: only minimal testing against a baseline is implemented!")
    print("*")
    print(80 * '*')

    print(80 * '*', file=testlog)
    print("*", file=testlog)
    print("* WARNING: only minimal testing against a baseline is implemented!", file=testlog)
    print("*", file=testlog)
    print(80 * '*', file=testlog)


def commandline_options():
    parser = argparse.ArgumentParser(description="Run a CLM-PFLOTRAN regression "
                                     "tests.")

    parser.add_argument("--backtrace", action='store_true', default=False,
                        help='show exception backtraces as extra debugging '
                        'output')

    parser.add_argument('-d', "--debug", action='store_true', default=False,
                        help='extra debugging output')

    parser.add_argument("-l", "--local-config-file",
                        nargs=1, default=None, required=True,
                        help="local configuration file to use (location of petsc, pflotran for this test)")

    parser.add_argument("-m", "--machines-config-file",
                        nargs=1, default=None, required=True,
                        help="machines configuration file to use")

    parser.add_argument("-t", "--test-config-files",
                        nargs='+', default=None, required=True,
                        help="test configuration file(s) to use")

    parser.add_argument("-r", "--cesm-root-dir",
                        nargs=1, default=None, required=True,
                        help="path to the root cesm directory to use")

    options = parser.parse_args()
    return options



def main(options):
    # must be done before setting up the suite log to get the correct pflotran/petsc info
    local_config = read_local_config(options.local_config_file[0])
    clear_pflotran_environment(local_config)

    txtwrap = textwrap.TextWrapper(width=78, subsequent_indent=4*" ")
    suitelog = setup_test_suite_log(txtwrap, local_config)

    working_dir = os.getcwd()

    print_no_baseline_warning(suitelog)

    status = {}

    machine = TestMachine()
    machine.setup(options.machines_config_file[0])
    print_section_seperator("Machine Info :", suitelog)
    print(machine, file=suitelog)

    cesm_root_dir = os.path.abspath(options.cesm_root_dir[0])

    print_section_seperator("Running tests :", suitelog)
    start_time = time.time()
    for test_config in options.test_config_files:
        os.chdir(working_dir)
        print_section_seperator('', suitelog)
        test = RegressionTest()
        try:
            test.setup(test_config, machine, local_config, cesm_root_dir, suitelog)
            if options.debug:
                print(test, file=suitelog)
            test.run_test(suitelog)
        except Exception as e:
            print(str(e), file=suitelog)
        status[test.name()] = test.status(suitelog)

    stop_time = time.time()
    print("")

    run_time = stop_time - start_time
    summary_report(run_time, status, suitelog)
    run_status = summary_report(run_time, status, sys.stdout)

    suitelog.close()

    return run_status

if __name__ == "__main__":
    options = commandline_options()
    try:
        status = main(options)
        sys.exit(status)
    except Exception as e:
        print(str(e))
        if options.backtrace:
            traceback.print_exc()
        sys.exit(1)

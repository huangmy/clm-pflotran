#!/usr/bin/env python

import os, sys, csv, time, math
from optparse import OptionParser
import Scientific.IO.NetCDF
from Scientific.IO import NetCDF
#from Numeric import *


#DMR 4/16/13
#python ./runCLM.py does the following:
#  1. Call routines to create point data (makepointdata.py, makemetdata.py)
#  2. Set point and case-specific namelist options
#  2. configure case
#  3. build (compile) CESM with clean_build first if requested
#  4. apply patch for transient CO2 if transient run
#  6. apply user-specified PBS and submit information
#  7. submit job to PBS queue if requested.
#
#  For reproducibility, a copy of the current call_PTCLM.py is saved
#  to the newly created case directory.  This is for informational
#  purposes only - the script should not be executed from within
#  the case directory.
#
# FMY 6/6/2013
# modified to work for CLM4-pf (CLM4.5.06, with PFLOTRAN interface) version used by NGEE-Arctic
#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--caseidprefix", dest="mycaseid", default="", \
                  help="Unique identifier to include as a prefix to the case name")
parser.add_option("--caseroot", dest="caseroot", default='./', \
                  help = "case root directory (default = ./, i.e., under scripts/)")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--runroot", dest="runroot", default="../runs", \
                  help="Directory where the run would be created")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  help = "input data directory for CESM (required)")
parser.add_option("--cesmdir", dest="cesmdir", default='..', \
                  help = "cesm directory (default = .., i.e., upper of scripts/)")
parser.add_option("--compset", dest="compset", default='I1850CLM45CN', \
                  help = "component set to use (required)")
parser.add_option("--machine", dest="machine", default = 'userdefined', \
                  help = "machine to use ---- \n"  
                         "default = userdefined \n"
                         "options = checking by ./create_newcase -list machines \n "
                         "NOTE: make sure the option you chose well-defined in config_machines.xml")
parser.add_option("--compiler", dest="compiler", default = 'gnu', \
                  help = "compiler to use on machine ---- \n"
                         "   default = gnu) \n "
                         "   options = intel,ibm, pgi,pathscale,gnu,cray,lahey,userdefined \n "
                         "NOTE: make sure the option you chose well-defined in config_compilers.xml")
parser.add_option("--mpilib", dest="mpilib", default = 'mpi-serial', \
                  help = "mpi library to use (default = mpi-serial)"
                         "options=openmpi,mpich,mpt,ibm,mpi-serial, BUT upon your system")
parser.add_option("--no_fire", dest="nofire", action="store_true", \
                  default=False, help="Turn off fire algorightms")
#parser.add_option("--C13", dest="C13", default=False, \
#                  help = 'Switch to turn on C13', action="store_true")
parser.add_option("--vertsoilc", dest="vsoilc", default=False, \
                  help = 'To turn on CN with multiple soil layers, excluding CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--centbgc", dest="centbgc", default=False, \
                  help = 'To turn on CN with multiple soil layers, CENTURY C module (CLM4ME on as well)', action="store_true")
parser.add_option("--CH4", dest="CH4", default=False, \
                  help = 'To turn on CN with CLM4me', action="store_true")
#parser.add_option("--maxpft", dest="maxpft", default=17, \
#                  help = 'user-defined max no of pft, excluding crop pfts. Must provide --parm_file')
parser.add_option("--parm_file", dest="parm_file", default="", \
                  help = 'CLM user-defined physiological parameter file')
parser.add_option("--co2_file", dest="co2_file", default="fco2_datm_1765-2007_c100614.nc", \
                  help = 'CLM transient CO2 file for diagnostic option')
parser.add_option("--ad_spinup", action="store_true", \
                  dest="ad_spinup", default=False, \
                  help = 'Run accelerated decomposition spinup (note: exit-ad will run in the end as well')
parser.add_option("--nyears_ad_spinup", dest="ny_ad", default=600, \
                  help = 'number of years to run ad_spinup')
parser.add_option("--branch", dest="branch", default=False, \
          help = 'Switch for branch run', action="store_true")
parser.add_option("--coldstart", dest="coldstart", default=False, \
                  help = "set cold start (mutually exclusive w/finidat)", \
                  action="store_true")
parser.add_option("--finidat_case", dest="finidat_case", default='', \
                  help = "case containing initial data file to use" \
                  +" (should be in your runroot directory)")
parser.add_option("--finidat", dest="finidat", default='', \
                  help = "initial data file to use" \
                  +" (should be in the runroot directory)")
parser.add_option("--finidat_year", dest="finidat_year", default=-1, \
                  help = "model year of initial data file (default is" \
                  +" last available)")
parser.add_option("--run_units", dest="run_units", default='nyears', \
                  help = "run length units (ndays, nyears)")
parser.add_option("--run_n", dest="run_n", default=600, \
                  help = "run length (in run units)")
parser.add_option("--align_year", dest="align_year", default="1850", \
                  help = 'Alignment year (transient run only)')
parser.add_option("--rmold", dest="rmold", default=False, action="store_true", \
                  help = 'Remove old case directory with same name' \
                  +" before create/update new one")
parser.add_option("--srcmods_loc", dest="srcmods_loc", default='', \
                  help = 'Copy sourcemods from this location')
parser.add_option("--hist_userdefined", dest="hist_file", default='', \
                  help = 'user=defined hist file')
parser.add_option("--clean_build", dest="clean_build", default=False, \
                  help = 'Perform clean build before building', \
                  action="store_true")
parser.add_option("--debug_build", dest="debug_build", default=False, \
                  help = 'Perform debug build', \
                  action="store_true")
parser.add_option("--clean_config", dest="clean_config", default=False, \
                  help = 'Perform clean setup before setting-up', \
                  action="store_true")
parser.add_option("--no_config", dest="no_config", default=False, \
                  help = 'do NOT configure case', action="store_true")
parser.add_option("--no_build", dest="no_build", default=False, \
                  help = 'do NOT build CESM', action="store_true")
parser.add_option("--no_submit", dest="no_submit", default=False, \
                  help = 'do NOT submit CESM to queue', action="store_true")
parser.add_option("--cleanlogs",dest="cleanlogs", help=\
                   "Removes temporary and log files that are created",\
                   default=False,action="store_true")
parser.add_option("--queue", dest="queue", default='essg08q', \
                  help = 'PBS submission queue')
parser.add_option("--np", dest="np", default=1, \
                  help = 'number of processors')
parser.add_option("--ninst", dest="ninst", default=1, \
                  help = 'number of land model instances')
parser.add_option("--tstep", dest="tstep", default=0.5, \
                  help = 'CLM timestep (hours)')
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--makemetdata", dest="makemet", default=False, \
		  help = 'Generate meteorology', action="store_true")
parser.add_option("--nopointdata", action="store_true", \
                  dest="nopointdata", help="Do NOT make point data (use data already created)", \
                  default=False)
parser.add_option("--ugriddir", dest="ugriddir", default="none", \
                  help = "subdirectory under 'ccsm_input' for the following 0.5x0.5 datasets to make point data: \n"                
                  "(1)domain.360x720_ORCHIDEE0to360.100409.nc; \n" 
                  "(2)surfdata_360x720cru_simyr1850_c130415.nc (for I1850CLM45CN); \n" 
                  "(3)surfdata.pftdyn_0.5x0.5_simyr1850-2010.nc (for I20TRCLM45CN)")
parser.add_option("--soilgrid", dest="soilgrid", default=False, \
                  help = 'Use gridded soil data', action="store_true")
parser.add_option("--regional", dest="regional", default=False, \
                   help="flag for regional", action="store_true")
parser.add_option("--refcase", dest="refcase" , default='none', \
                  help = 'Use already compiled CLM case')

(options, args) = parser.parse_args()

#------------------- arguments ------------------------------------------------------------

# current directory
PTCLMdir = os.getcwd()

# cesm model directory
csmdir = os.path.abspath(options.cesmdir)
scriptsdir = csmdir+'/scripts'

# case directory
if (options.machine == 'userdefined' or \
options.caseroot == '' or \
(os.path.exists(options.caseroot) == False)):
    caseroot = csmdir+'/cases'
else:
    caseroot = os.path.abspath(options.caseroot)
print('CASE root directory: '+options.caseroot)

# case run root directory
if (options.machine == 'userdefined' or \
options.runroot == '' or \
(os.path.exists(options.runroot) == False)):
    runroot = csmdir+"/runs"
else:
    runroot = os.path.abspath(options.runroot)
print('CASE RUN root directory: '+runroot)

#check for valid input data directory
if (options.machine == 'userdefined' or \
options.ccsm_input == '' or \
(os.path.exists(options.ccsm_input) == False)):
    print('Error:  invalid input data directory')
    sys.exit()
else:
    ccsm_input = os.path.abspath(options.ccsm_input)

#check for valid compset
compset = options.compset
if (compset != 'I1850CLM45CN' and compset != 'ICLM45CN' and compset != 'I20TRCLM45CN'):
    print('Error:  please enter one of following valid options for compset')
    print('        (I1850CLM45CN, ICLM45CN, I20TRCLM45CN)')
    sys.exit()

#check consistency of options   
if (compset == 'I20TRCLM45CN'):
    #ignore spinup option if transient compset
    if (options.ad_spinup):
        print('Spinup options not available for transient compset.')
        sys.exit()        
    
#get full path of finidat file
finidat      = options.finidat
finidat_year = int(options.finidat_year)

if (finidat == '' and options.finidat_case == ''):   # not user-defined
    if (options.coldstart==False and compset == "I1850CLM45CN"):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+options.site+ \
                                  '_I1850CLM45CN_ad_spinup'
        else:
            options.finidat_case = options.site+'_I1850CLM45CN_ad_spinup'
    
        if (options.finidat_year == -1):
            finidat_year = int(options.ny_ad)+1
    
    if (compset == "I20TRCLM45CN"):
        if (options.mycaseid != ''):
            options.finidat_case = options.mycaseid+'_'+options.site+ \
                                  '_I1850CLM45CN'
        else:
            options.finidat_case = options.site+'_I1850CLM45CN'
        
        if (options.finidat_year == -1):
            finidat_year = 1850
        
        #finidat is required for transient compset
        if (os.path.exists(runroot+'/'+options.finidat_case) == False):
            print('Error:  must provide initial data file for I20TRCLM45CN compset, OR, '+ \
                  runroot+'/'+options.finidat_case+' existed as refcase')
            sys.exit()


if (options.finidat_case != ''):
    finidat_yst = str(finidat_year)
    if (finidat_year >= 100 and finidat_year < 1000):
        finidat_yst = '0'+str(finidat_year)
    if (finidat_year >= 10 and finidat_year < 100):
        finidat_yst = '00'+str(finidat_year)
    if (finidat_year < 10):
        finidat_yst = '000'+str(finidat_year)

    if (finidat == ''):
        finidat = runroot+'/'+options.finidat_case+'/run/'+ \
                  options.finidat_case+'.clm2.r.'+finidat_yst+ \
                  '-01-01-00000.nc'

#construct default casename
casename    = options.site+"_"+compset
if (options.mycaseid != ""):
    casename = options.mycaseid+'_'+casename
if (options.ad_spinup):
    casename = casename+'_ad_spinup'

#get site year information
PTCLMfiledir = PTCLMdir+'/PTCLM_files'

sitedatadir = os.path.abspath(PTCLMfiledir+'/PTCLM_sitedata')
os.chdir(sitedatadir)
AFdatareader = csv.reader(open(options.sitegroup+'_sitedata.txt',"rb"))
for row in AFdatareader:
    if row[0] == options.site:
        startyear=int(row[6])
        endyear=int(row[7])
        alignyear = int(row[8])
        if (options.regional == True):
            numxpts=int(row[9])
            numypts=int(row[10])
        else:
            numxpts=1
            numypts=1
os.chdir(PTCLMdir)

#get simyr
mysimyr=1850
if (options.compset == 'ICLM45CN'):
    mysimyr=2000


#construct case dir (clm4.5 case has two parts: directory+casename
if (caseroot != "./"):
    casedir=caseroot+"/"+casename
else:
    casedir=casename

#Check for existing case directory
if (os.path.exists(casedir)):
    print('Warning:  Case directory exists and --rmold not specified')
    var = raw_input('proceed (p), remove old (r), or exit (x)? ')
    if var[0] == 'r':
        os.system('rm -rf '+casedir)
    if var[0] == 'x':
        sys.exit()    
print("CASE directory is: "+casedir+"\n")

#construct case build and run directory
blddir=runroot+"/"+casename
print ("CASE exeroot is: "+blddir+"\n")
rundir=runroot+"/"+casename+"/run"
print ("CASE rundir is: "+rundir+"\n")

#Environment variable hacks
os.putenv("CLM_USRDAT_NAME", str(numxpts)+"x"+str(numypts)+"pt_"+options.site)    
os.putenv("DOMAINPATH", options.ccsm_input+'/share/domains/domain.clm')

#------------------- make point data for site -------------------------------

os.chdir(PTCLMdir)
ptcmd = 'python makepointdata.py --compset '+compset+ \
              ' --site '+options.site+' --sitegroup '+options.sitegroup+ \
              ' --csmdir '+csmdir+' --ccsm_input '+options.ccsm_input

ptcmd = ptcmd + ' --grid_input /'+options.ugriddir
if (options.metdir != 'none'):
    ptcmd = ptcmd + ' --metdir '+options.metdir
if (options.makemet):
    ptcmd = ptcmd + ' --makemetdata'
if (options.soilgrid):
    ptcmd = ptcmd + ' --soilgrid'
if (options.regional):
    ptcmd = ptcmd + ' --regional'
os.system(ptcmd)

# default pft parameter file
pftfile = ccsm_input+'/lnd/clm2/pftdata/pft-physiology.c130503.'+options.site+'.nc'
os.system('cp -f '+ccsm_input+'/lnd/clm2/pftdata/pft-physiology.c130503.nc ' \
              + pftfile)
# new or user-defined pft-phys file if desired
if (options.parm_file != ''):
    pftfile = ccsm_input+'/lnd/clm2/pftdata/' + \
                  options.parm_file + '.' + options.site + '.nc'
    os.system('cp -f '+ccsm_input+'/lnd/clm2/pftdata/'+ options.parm_file + \
               pftfile)

#------------------IF no refcase, create, configure and build -----------------------------------------
#set number of run years for ad-spinup cases
if (options.ny_ad != options.run_n and options.ad_spinup):
    options.run_n = options.ny_ad

if (options.run_n != 600 and options.compset== "I20TRCLM45CN"):
    options.run_n = endyear - mysimyr +1

# set 'debug' mode
debugoption = ''
if (options.debug_build):
    debugoption = ' -confopts _D'
    print ("case build will be configured for debugging \n")

# set machine and compiler information, new addtions in clm45
compileroption = ' -compiler gnu'
if (options.compiler != "gnu"):
    compileroption = ' -compiler '+options.compiler

mpiliboption = ' -mpilib mpi-serial'
if (options.mpilib != "mpi-serial"):
    mpiliboption = ' -mpilib '+options.mpilib

os.chdir(scriptsdir)
if (options.refcase == 'none'):
    #create new case
    print ('./create_newcase -case '+casedir+' -mach '+options.machine + \
                 ' -compset '+ options.compset +' -res CLM_USRDAT ' + \
                 compileroption + mpiliboption + debugoption)
    os.system('./create_newcase -case '+casedir+' -mach '+options.machine + \
                 ' -compset '+ options.compset +' -res CLM_USRDAT ' + \
                 compileroption + mpiliboption + debugoption + \
                  ' > create_newcase.log')
    if (os.path.isdir(casedir)):
        print(casename+' created.  See create_newcase.log for details')
        os.system('mv create_newcase.log ' + casedir +"/"+casename+"_case.log")
    else:
        print('failed to create case.  See create_newcase.log for details')

    # go to newly created case directory
    os.chdir(casedir)

#----------- env_build.xml modification ---------------------------
    
    # for some reasons, user-defined a few critical root direcotories are not easily modified in machine files
    if (options.runroot != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_build.xml -id ' \
                  +'EXEROOT -val '+runroot+'/'+casename) 

    # turn off rof module
    os.system('./xmlchange -file env_build.xml -id ' \
                  +'RTM_MODE -val NULL') 

    # clm4_5 cn config options
    clmcn_opts = "-phys clm4_5 -bgc cn"
    # Koven's multiple soil layer, nitri-denitri, and bsw consistency (excluding century C model)
    # clm4me has to be turned on as well      
    if (options.vsoilc):
        if (options.centbgc):
            clmcn_opts += " -clm4me on -vsoilc_centbgc on"
        else:
            clmcn_opts += " -clm4me on -vsoilc_centbgc no-cent"
    else:
        # CLM4ME option only
        if (options.CH4):
            clmcn_opts += " -clm4me on"
   
    #turn off fire module 
    if (options.nofire):
            clmcn_opts += " -nofire"

    os.system('./xmlchange -file env_build.xml -id ' \
                  +'CLM_CONFIG_OPTS -val "'+clmcn_opts+'"')
    print ("CLM module options: " + clmcn_opts +"\n")

#---------- env_run.xml modification ------------------------------------
    # input/run/output directory
    if (options.runroot != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'RUNDIR -val '+rundir) 
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S -val TRUE') 
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DOUT_S_ROOT -val '+runroot+'/archives/'+casename) 
    if (options.ccsm_input != '' or options.machine == "userdefined"):
        os.system('./xmlchange -file env_run.xml -id ' \
                  +'DIN_LOC_ROOT -val '+ccsm_input) 
    
    # datm options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_MODE -val CLM1PT') 
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_START -val '+str(startyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_END -val '+str(endyear))
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'DATM_CLMNCEP_YR_ALIGN -val '+str(alignyear))

    # run timestep
    if (options.tstep != 0.5):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'ATM_NCPL -val '+str(int(24/float(options.tstep))))
    
    # run-type adjusting
    if (options.branch):
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_TYPE -val branch')
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFDATE -val '+finidat_yst+'-01-01')
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFCASE -val '+options.finidat_case)
    else:
        if (options.ad_spinup==False and options.coldstart==False):
            os.system('./xmlchange -file env_run.xml -id ' \
                      +'RUN_REFDATE -val '+finidat_yst+'-01-01')
    
    #adds capability to run with transient CO2
    if (compset == 'I20TRCLM45CN'):
        os.system('./xmlchange -file env_run.xml -id ' \
                          +'CCSM_BGC -val CO2A')
        os.system('./xmlchange -file env_run.xml -id ' \
                          +'CLM_CO2_TYPE -val diagnostic')       

    # user-defined running stop options
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_OPTION -val '+options.run_units)
    os.system('./xmlchange -file env_run.xml -id ' \
                  +'STOP_N -val '+str(options.run_n))

    #User-defined resolution
    os.system('./xmlchange -file env_run.xml -id CLM_BLDNML_OPTS ' \
                  +'-val "-mask navy"')
    os.system('./xmlchange -file env_run.xml -id CLM_USRDAT_NAME ' \
                  +'-val '+str(numxpts)+'x'+str(numypts)+'pt_'+options.site)

    # extra build namelist options for CLM
    stdout  = os.popen("./xmlquery -valonly -silent CLM_BLDNML_OPTS")
    env_val = stdout.read().rstrip( )   
    # ad-spinup with exit-spinup included
    if (options.ad_spinup):
        env_val += " -bgc_spinup on"
        os.system('./xmlchange -file env_run.xml -id ' \
                      +'CLM_BLDNML_OPTS -val "'+env_val+'"')
    
    # not yet figured out how to set this option, obviously not here
    #if (options.C13):
    #        clmcn_opts += " -use_c13 on"

    # not yet figured out how to set this option, obviously not here
    #modification of maxpft
    #if (options.maxpft != 17):
    #        clmcn_opts += " -maxpft "+str(options.maxpft)
    
        
#--------- env_mach_pes.xml modification ------------------------------------
    # normal pt mode
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val 1')
    os.system('./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val 1')
    
    #if number of land instances > 1
    if (int(options.ninst) > 1):
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NINST_LND -val '+options.ninst)
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_LND -val '+options.ninst)

    #if running with > 1 processor
    if (int(options.np) > 1):
        os.system('./xmlchange -file env_mach_pes.xml -id ' \
                      +'NTASKS_LND -val '+options.np)

#------- cesm setup -------------------------------------------
    #clean configure if requested prior to configure
    if (options.clean_config):
        os.system('./cesm_setup -clean')
        os.system('rm -f Macro')
        os.system('rm -f user-nl-*')

    if (options.no_config == False):
        os.system('./cesm_setup > configure.log')
        
        # geneated *.run lacks of csh command
        myinput  = open("./"+casename+".run")  
        myoutput = open("./"+casename+"temp.run",'w')
        myoutput.write("#! /bin/csh -f \n")
        for s in myinput:
            myoutput.write(s)
        myoutput.close()
        myinput.close()
        os.system("mv "+casename+"temp.run "+casename+".run")  
        os.system("chmod a+x ./"+casename+".run")

    else:
        print("Warning:  No case configure performed.  PTCLM will not " \
                  +"make any requested modifications to env_*.xml files.  Exiting.")
        sys.exit()    

    #clm user-defined namelist from a file 
    output = open("user_nl_clm",'w')
    output.write('&clm_inparm\n')
        
    # clm output hist user-defined file
    if (options.hist_file != ''):
        hvars_file = open(options.hist_file)
        for s2 in hvars_file:
            myline = s2.strip()
            output.write(myline+"\n")
            hvars_file.close()
            
    if (finidat != ''):
        output.write(" finidat = '"+finidat+"'\n")
    
    output.write(" fsurdat = '"+ccsm_input+"/lnd/clm2/surfdata/" + \
                 "surfdata_"+str(numxpts)+"x"+str(numypts)+"pt_"+options.site+ \
                 "_simyr"+str(mysimyr)+".nc'\n")
    
    if (compset == 'I20TRCLM45CN'):
        output.write(" fpftdyn = '"+ccsm_input+"/lnd/clm2/surfdata/" + \
                          "surfdata.pftdyn_"+str(numxpts)+"x"+str(numypts)+"pt_"+ \
                          options.site+".nc'\n")
    
    if (pftfile != ''):
        output.write(" fpftcon=   '" + pftfile + "'\n")

    output.close()

    #copy sourcemods
    if (options.srcmods_loc != ''):
        if (os.path.exists(options.srcmods_loc) == False):
            print('Invalid srcmods directory.  Exiting')
            sys.exit()
        options.srcmods_loc = os.path.abspath(options.srcmods_loc)
        os.system('cp -r '+options.srcmods_loc+'/* ./SourceMods')   
   
#------- build clm45 within cesm  
    #clean build if requested prior to build
    if (options.clean_build):
        os.system('./'+casename+'.clean_build')
        os.system('rm -rf '+runroot+'/'+casename+'/run/*.*log.*')
        os.system('rm -rf '+runroot+'/'+casename+'/run/*.nc')
        os.system('rm -rf '+runroot+'/'+casename+'/*.*log.*')
        os.system('rm -rf '+runroot+'/'+casename+'/*.exe.*')
    
    #compile cesm
    if (options.no_build == False):
        os.system('./'+casename+'.build')
        
    #transient CO2 patch for transient run (datm buildnml mods)
    #if (compset == "I20TRCLM45CN"):

    
#--------------- make necessary modificaitons to run script for OIC
    if (options.machine == "linux_pgi"):
        myinput  = open("./"+casename+".run")
        myoutput = open("./"+casename+"temp.run",'w')
        for s in myinput:
            if s[6:8]  == '-N':
                myoutput.write("#PBS -N "+casename+"\n")
            elif s[9:14] == 'batch':
                myoutput.write("#PBS -q "+options.queue+"\n")
            #elif s[0:4] == 'cd /':                 
                #output.write("cd "+csmdir+'/scripts/'+casename+"\n")
                #os.chdir(casedir)
            elif s[0:7] == './Tools':
                myoutput.write("cd "+casedir+"\n")
                myoutput.write(s)
            elif s[0:14] =="##PBS -l nodes":
                myoutput.write("#PBS -l nodes="+str((int(options.np)-1)/8+1)+ \
                                 ":ppn="+str(min(int(options.np),8))+"\n")  
            elif s[9:17] == 'walltime':
                myoutput.write("#PBS -l walltime=48:00:00\n") 
            elif s[0:5] == '##PBS':
                myoutput.write(s.replace("##PBS","#PBS"))
            elif s[0:7] == '   exit':
                myoutput.write('   #exit 2')
            elif s[0:10] == '   #mpirun':
                myoutput.write("   mpirun -np "+str(options.np)+" --hostfile $PBS_NODEFILE ./ccsm.exe >&! ccsm.log.$LID\n")
            elif s[0:5] == 'sleep':
                myoutput.write("sleep 5\n")
            else:
                myoutput.write(s)
        myoutput.close()
        myinput.close()
        os.system("mv "+casename+"temp.run "+casename+".run")  
    

#----------------------------Reference case set ------------------------------------------
# the following has not yet checked for CLM4.5
else:  

    os.chdir(rundir)
    incasename  = options.refcase+'_REFCASE_'+options.compset
    if (options.ad_spinup):
        incasename = incasename + '_ad_spinup'
    os.system('mkdir -p '+casename+'/run')
    os.chdir(casename+'/run')
    print 'Copying files from '+incasename+' to '+casename
    os.system('cp ../../'+incasename+'/run/*_in* .')
    os.system('cp ../../'+incasename+'/run/ccsm.exe .')
    os.system('cp ../../'+incasename+'/run/*.nml .')
    os.system('cp ../../'+incasename+'/run/*eam* .')
    os.system('cp ../../'+incasename+'/run/*.rc .')
   

    #Change generic site/case name to actual site/case name in namelst files
    os.system('chmod u+w *')
    os.system('sed -e s/'+incasename+'/'+casename+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    os.system('sed -e s/REFCASE/'+options.site+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    ptstr = str(numxpts)+'x'+str(numypts)+'pt'
    os.system('sed -e s/1x1pt/'+ptstr+'/ig  lnd_in > lnd_in_tmp')
    os.system('mv lnd_in_tmp lnd_in')
    os.system('sed -e s/'+incasename+'/'+casename+'/ig  datm_atm_in > datm_atm_in_tmp')
    os.system('mv datm_atm_in_tmp datm_atm_in')
    os.system('sed -e s/REFCASE/'+options.site+'/ig  datm_atm_in > datm_atm_in_tmp')
    os.system('mv datm_atm_in_tmp datm_atm_in')
    os.system('sed -e s/1x1pt/'+ptstr+'/ig  datm_atm_in > datm_atm_in_tmp')
    os.system('mv datm_atm_in_tmp datm_atm_in')
    os.system('sed -e s/CLM_USRDAT/1x1pt_'+options.site+'/ig  datm_atm_in > datm_atm_in_tmp')
    os.system('mv datm_atm_in_tmp datm_atm_in')
    os.system('mv clm1PT.CLM_USRDAT.stream.txt clm1PT.1x1pt_REFCASE.stream.txt')
    os.system('sed -e s/REFCASE/'+options.site+'/ig clm1PT.1x1pt_REFCASE.stream.txt > clm1PTstream.tmp')
    os.system('mv clm1PTstream.tmp clm1PT.'+ptstr+'_'+options.site+'.stream.txt')
    os.system('sed -e s/1x1pt/'+ptstr+'/ig clm1PT.'+ptstr+'_'+options.site+'.stream.txt > clm1PTstream.tmp')
    os.system('mv clm1PTstream.tmp clm1PT.'+ptstr+'_'+options.site+'.stream.txt')
    os.system('rm *REFCASE*')
    os.system('sed -e s/'+incasename+'/'+casename+'/ig  drv_in > drv_in_tmp')
    os.system('mv drv_in_tmp drv_in')
    os.system('sed -e s/REFCASE/'+options.site+'/ig  drv_in > drv_in_tmp')
    os.system('mv drv_in_tmp drv_in')
    
    #modify met stream file for correct years
    myinput  = open('clm1PT.'+ptstr+'_'+options.site+'.stream.txt')
    myoutput = open('clm1PT.'+ptstr+'_'+options.site+'.stream.txt.tmp','w')
    for s in myinput:
        if (s[0:22] == '            2000-01.nc'):
            for y in range(startyear,endyear+1):
                for m in range(1,13):
                    if (m < 10):
                        myoutput.write('            '+str(y)+'-0'+str(m)+'.nc\n')
                    else:
                        myoutput.write('            '+str(y)+'-'+str(m)+'.nc\n')
        elif (s[0:17] == '            2000-'):
            continue  #do nothing
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv clm1PT.'+ptstr+'_'+options.site+'.stream.txt.tmp clm1PT.'+ptstr+'_'+options.site+'.stream.txt')

    #modify presearo stream file to change to 1850-2000 file
    myinput  = open('presaero.stream.txt')
    myoutput = open('presaero.stream.txt.tmp','w')
    for s in myinput:
        if (s[0:22] == '            aerosoldep'):
            myoutput.write('            aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc\n')
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv presaero.stream.txt.tmp presaero.stream.txt')

    #modify datm_atm_in for correct years
    myinput  = open('datm_atm_in')
    myoutput = open('datm_atm_in_tmp','w')
    for s in myinput:
        if (s[0:10] == '   streams'):
            if (compset == 'I20TRCLM45CN'):
                myoutput.write("   streams        = 'clm1PT."+ptstr+"_"+options.site+".stream"+ \
                                   ".txt "+str(options.align_year)+" "+str(startyear)+" "+str(endyear)+" ',\n")  
            else:
                myoutput.write("   streams        = 'clm1PT."+ptstr+"_"+options.site+ \
                                   ".stream.txt 1 "+str(startyear)+" "+str(endyear)+" ',\n")
        elif (s[0:40] == "                    'presaero.stream.txt'"):
            if (compset != 'I20TRCLM45CN'):
                myoutput.write("                       'presaero.stream.txt 1 1 1'\n")
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv datm_atm_in_tmp datm_atm_in')

    #modify component .nml files
    nmlfiles=['atm','cpl','glc','ice','lnd','ocn']
    for mynml in nmlfiles:
        outfile = open(mynml+'_modelio.nml','w')
        outfile.write('&modelio\n')
        outfile.write('   diri    = "'+os.path.abspath('../..')+'/'+incasename+'/'+ \
                          mynml+'   "\n')
        outfile.write('   diro    = "./"\n') #'+os.path.abspath('.')+'   "\n')
        outfile.write('   logfile = "'+mynml+'.log   "\n')
        outfile.write('/\n')
        outfile.close()
          
    #make drv_in namelist modifications (run length for final spin/tranisent case)
    myinput  = open('drv_in')
    myoutput = open('drv_in_tmp','w')
    for s in myinput:
        if (s[0:8] == '  stop_n'):
            myoutput.write("  stop_n         = "+str(options.run_n)+'\n')                               
        elif (s[0:11] == '  restart_n'):
            myoutput.write("  restart_n      = "+str(options.run_n)+'\n')
        elif (s[0:10] == '  stop_ymd'):
            myoutput.write("  stop_ymd       = -999\n")
        elif (s[0:13] == '  restart_ymd'):
            myoutput.write("  restart_ymd    = -999\n")
        elif (s[0:12] == '  atm_cpl_dt'):
            myoutput.write("  atm_cpl_dt     = "+str(float(options.tstep)*3600)+'\n')
        elif (s[0:12] == '  lnd_cpl_dt'):
            myoutput.write("  lnd_cpl_dt     = "+str(float(options.tstep)*3600)+'\n')
        elif (s[0:12] == '  ice_cpl_dt'):
            myoutput.write("  atm_cpl_dt     = "+str(float(options.tstep)*3600)+'\n')
        elif (s[0:11] == '  start_ymd'):
            if (options.exit_spinup):
                myoutput.write("  start_ymd      = "+finidat_yst+'0101\n')
            else:
                myoutput.write(s)
        else:
            myoutput.write(s)
    myinput.close()
    myoutput.close()
    os.system('mv drv_in_tmp drv_in')
            
    #write a basic PBS script
    output = open(casename+'.run','w')
    output.write('#PBS -S /bin/bash\n')
    output.write('#PBS -V\n')
    output.write('#PBS -m ae\n')
    output.write('#PBS -N '+casename+'\n')
    output.write('#PBS -q '+options.queue+'\n')
    output.write("#PBS -l nodes="+str((int(options.np)-1)/8+1)+ \
                     ":ppn="+str(min(int(options.np),8))+"\n")  
    output.write('#PBS -l walltime=48:00:00\n')
    output.write("cd "+csmdir+'/run/'+casename+"/run\n")
    if (options.np == 1):
        output.write("./ccsm.exe > ccsm_log.txt\n")
    else:
        output.write("mpirun -np "+options.np+" --hostfile $PBS_NODEFILE ./ccsm.exe\n")
    output.close()

#---------------------------end of refcase ------------------------------------------------

#----- copy rpointers and restart files to current run directory prior to run model
if (options.finidat_case != ''):
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/' + \
              options.finidat_case+'.*'+finidat_yst+'* ' + rundir)
    os.system('cp -f '+runroot+'/'+options.finidat_case+'/run/'+ \
              'rpointer.* '+rundir)

#----- submit job if requested
if (options.no_submit == False):
    os.chdir(casedir)
    if (options.machine == "darwin-gnu"):
        os.system("./"+casename+".run")
    else:   
        os.system("qsub ./"+casename+".run")


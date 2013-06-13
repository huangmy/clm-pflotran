#!/usr/bin/env python
import os, sys, csv, time, math
from optparse import OptionParser
import Scientific.IO.NetCDF
from Scientific.IO import NetCDF
import numpy

print('\n')
print('Makepointdata.py version 0.2')
print('Utillity to create point-level data from 0.5x0.5 gridded datasets')
print('NOTE - Python libraries and tools needed: Scientific, numpy, nco: ncks')
print('Contact: ricciutodm@ornl.gov')

parser = OptionParser()
parser.add_option("--compset", dest="compset", default='I1850CLM45CN', \
                  help = "component set to use (required)")
parser.add_option("--site", dest="site", default='', \
                  help = '6-character FLUXNET code to run (required)')
parser.add_option("--sitegroup", dest="sitegroup", default="AmeriFlux", \
                  help = "site group to use (default AmeriFlux)")
parser.add_option("--regional", action="store_true", \
                   dest="regional", default=False, \
                   help="Flag for regional run (2x2 or greater)")
parser.add_option("--csmdir", dest="csmdir", default='..', \
                  help = "base CESM directory (default = ../)")
parser.add_option("--ccsm_input", dest="ccsm_input", \
                  default='../../ccsm_inputdata', \
                  help = "input data directory for CESM (required)")
parser.add_option("--grid_input", dest="grid_input", \
                  default='/ugrid', \
                  help = "input data directory for creating point grid/surface data)")
parser.add_option("--metdir", dest="metdir", default="none", \
                  help = 'subdirectory for met data forcing')
parser.add_option("--makemetdata", dest="makemet", default=False, \
		  help = 'Generate meteorology', action="store_true")
parser.add_option("--soilgrid", dest="soilgrid", default=False, \
                  help = 'Use gridded soil data', action="store_true")

(options, args) = parser.parse_args()


csmdir=os.path.abspath(options.csmdir)
options.ccsm_input = os.path.abspath(options.ccsm_input)

#------------------- get site information ----------------------------------
sitedatadir = os.path.abspath('./PTCLM_files/PTCLM_sitedata')
os.chdir(sitedatadir)
AFdatareader = csv.reader(open(options.sitegroup+'_sitedata.txt',"rb"))
for row in AFdatareader:
    if row[0] == options.site:
        lon=float(row[3])
        if (lon < 0):
            lon=360.0+float(row[3]) 
        lat=float(row[4])
        startyear=int(row[6])
        endyear=int(row[7])
        alignyear = int(row[8])
        if options.regional == True:
            numxpts=int(row[9])
            numypts=int(row[10])
            resx=float(row[11])
            resy=float(row[12])
        else:
            numxpts=1
            numypts=1
            resx=0.1      #longitudinal resolution (degrees) 
            resy=0.1      #latitudinal resolution (degrees)
        if (options.makemet):
            print(" Making meteorological data for site")
            metcmd = 'python '+csmdir+'/scripts/makemetdata.py' \
                          +' --site '+options.site+' --lat '+row[4]+' --lon '+ \
                          row[3]+' --ccsm_input '+options.ccsm_input+ \
                          ' --startyear '+row[6]+' --endyear '+row[7]+' --numxpts '+ \
                          str(numxpts)+' --numypts '+str(numypts)
            if (options.metdir != 'none'):
                metcmd = metcmd + ' --metdir '+options.metdir
            os.system(metcmd)
        else:
            print('Met data not requested.  Model will not run if data do not exist')

#get corresponding 0.5x0.5 degree grid cells
if (lon < 0):
    xgrid = 720+int(round(lon*2))
else:
    xgrid = int(round(lon*2))
ygrid = int(round(lat*2)+180)

#---------------------Create domain data --------------------------------------------------
print(' ------ Creating domain data ------')
domainfile_orig = options.ccsm_input+options.grid_input \
    +'/domain.360x720_ORCHIDEE0to360.100409.nc'
print('origin: '+ domainfile_orig)
domainfile_new = options.ccsm_input+'/share/domains/domain.clm/' \
    +'domain.lnd.'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'_navy.nc'
print('new: '+domainfile_new)
if (os.path.isfile(domainfile_new)):
    print('Warning:  Removing existing domain file: '+domainfile_new)
    os.system('rm -rf '+domainfile_new)
os.system('ncks -d ni,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d nj,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+domainfile_orig+' '+domainfile_new)
domainfile_new_nc = NetCDF.NetCDFFile(domainfile_new, "a")
frac = domainfile_new_nc.variables['frac']
frac_vals = frac.getValue()
mask = domainfile_new_nc.variables['mask']
mask_vals = mask.getValue()
xc = domainfile_new_nc.variables['xc']
xc_vals = xc.getValue()
yc = domainfile_new_nc.variables['yc']
yc_vals = yc.getValue()
xv = domainfile_new_nc.variables['xv']
xv_vals = xv.getValue()
xv.assignValue(xv_vals)
yv = domainfile_new_nc.variables['yv']
yv_vals = yv.getValue()
area = domainfile_new_nc.variables['area']
area_vals = area.getValue()
for i in range(0,numxpts):
    for j in range(0,numypts):
        frac_vals[j][i] = 1.0
        mask_vals[j][i] = 1
        xc_vals[j][i] = lon+i*resx
        yc_vals[j][i] = lat+j*resy
        xv_vals[j][i][0] = lon-resx/2+i*resx
        xv_vals[j][i][1] = lon+resx/2+i*resx
        xv_vals[j][i][2] = lon-resx/2+i*resx
        xv_vals[j][i][3] = lon+resx/2+i*resx
        yv_vals[j][i][0] = lat-resy/2+j*resy
        yv_vals[j][i][1] = lat-resy/2+j*resy
        yv_vals[j][i][2] = lat+resy/2+j*resy
        yv_vals[j][i][3] = lat+resy/2+j*resy
        area_vals[j][i] = resx*resy*math.pi/180*math.pi/180
frac.assignValue(frac_vals)
mask.assignValue(mask_vals)
xc.assignValue(xc_vals)
yc.assignValue(yc_vals)
xv.assignValue(xv_vals)
yv.assignValue(yv_vals)
area.assignValue(area_vals)
domainfile_new_nc.close()

os.system("cp -f "+domainfile_new+" "+options.ccsm_input+"/atm/datm7/domain.clm/")

#-------------------- create surface data ----------------------------------
print(' ------ Creating surface data ------')
if (lon < 0):
    xgrid = 720+int(round(lon*2))
else:
    xgrid = int(round(lon*2))
ygrid = int(round(lat*2)+180)

mysimyr=1850
if (options.compset == 'ICLM45CN'):
    mysimyr=2000
surffile_orig = options.ccsm_input+options.grid_input  \
    +'/surfdata_360x720cru_simyr1850_c130415.nc'
print('origin: '+ surffile_orig)
surffile_new = options.ccsm_input+'/lnd/clm2/surfdata/' \
    +'surfdata_'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'_simyr'+str(mysimyr)+'.nc'
print('new: '+ surffile_new)
if (os.path.isfile(surffile_new)):
    print('Warning:  Removing existing surface file ' + surffile_new)
    os.system('rm -rf '+surffile_new)
os.system('ncks -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+surffile_orig+' '+surffile_new)
    
surffile_new_nc = NetCDF.NetCDFFile(surffile_new, "a")
landfrac_pft = surffile_new_nc.variables['LANDFRAC_PFT']
landfrac_pft_vals = landfrac_pft.getValue()
pftdata_mask = surffile_new_nc.variables['PFTDATA_MASK']
pftdata_mask_vals = pftdata_mask.getValue()
longxy = surffile_new_nc.variables['LONGXY']
longxy_vals = longxy.getValue()
latixy = surffile_new_nc.variables['LATIXY']
latixy_vals = latixy.getValue()
area = surffile_new_nc.variables['AREA']
area_vals = area.getValue()
pct_wetland = surffile_new_nc.variables['PCT_WETLAND']
pct_wetland_vals = pct_wetland.getValue()
pct_lake = surffile_new_nc.variables['PCT_LAKE']
pct_lake_vals = pct_lake.getValue()
pct_glacier = surffile_new_nc.variables['PCT_GLACIER']
pct_glacier_vals = pct_glacier.getValue()
pct_urban = surffile_new_nc.variables['PCT_URBAN']
pct_urban_vals = pct_urban.getValue()

#input from site-specific information
soil_color = surffile_new_nc.variables['SOIL_COLOR']
soil_color_vals = soil_color.getValue()
pct_sand = surffile_new_nc.variables['PCT_SAND']
pct_sand_vals = pct_sand.getValue()
pct_clay = surffile_new_nc.variables['PCT_CLAY']
pct_clay_vals = pct_clay.getValue()
organic = surffile_new_nc.variables['ORGANIC']
organic_vals = organic.getValue()
fmax = surffile_new_nc.variables['FMAX']
fmax_vals = fmax.getValue()
pct_pft = surffile_new_nc.variables['PCT_PFT']
pct_pft_vals = pct_pft.getValue()
monthly_lai = surffile_new_nc.variables['MONTHLY_LAI']
monthly_lai_vals = monthly_lai.getValue()
monthly_sai = surffile_new_nc.variables['MONTHLY_SAI']
monthly_sai_vals = monthly_sai.getValue()
monthly_height_top = surffile_new_nc.variables['MONTHLY_HEIGHT_TOP']
monthly_height_top_vals = monthly_height_top.getValue()
monthly_height_bot = surffile_new_nc.variables['MONTHLY_HEIGHT_BOT']
monthly_height_bot_vals = monthly_height_bot.getValue()

npft = 17

#read file for site-specific PFT information
AFdatareader = csv.reader(open(options.sitegroup+'_pftdata.txt','rb'))
mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
for row in AFdatareader:
    #print(row[0], row[1], options.site)
    if row[0] == options.site:
        for thispft in range(0,5):
            mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])

#read file for site-specific soil information
if (options.soilgrid == False):
    AFdatareader = csv.reader(open(options.sitegroup+'_soildata.txt','rb'))
    for row in AFdatareader:
        if row[0] == options.site:
            mypct_sand = row[4]
            mypct_clay = row[5]

for i in range(0,numxpts):
    for j in range(0,numypts):
        landfrac_pft_vals[j][i] = 1.0
        pftdata_mask_vals[j][i] = 1
        longxy_vals[j][i] = lon+i*resx
        latixy_vals[j][i] = lat+j*resy
        area_vals[j][i] = 111.2*resy*111.321*math.cos((lon+i*resx)*math.pi/180)*resx
        pct_wetland_vals[j][i] = 0.0
        pct_lake_vals[j][i]    = 0.0
        pct_glacier_vals[j][i] = 0.0
        pct_urban_vals[j][i]   = 0.0
        soil_color_vals[j][i] = soil_color_vals[0][0]
        fmax_vals[j][i] = fmax_vals[0][0]
        for k in range(0,10):
            if (options.soilgrid == False):
                organic_vals[k][j][i] = organic_vals[k][0][0]
                pct_sand_vals[k][j][i]= mypct_sand
                pct_clay_vals[k][j][i]= mypct_clay
        for p in range(0,npft):
            pct_pft_vals[p][j][i] = mypft_frac[p]
            #print p, mypft_frac[p]
            for t in range(0,12):
                monthly_lai_vals[t][p][j][i] = monthly_lai[t][p][0][0]
                monthly_sai_vals[t][p][j][i] = monthly_sai[t][p][0][0]
                monthly_height_top_vals[t][p][j][i] = monthly_height_top[t][p][0][0]
                monthly_height_bot_vals[t][p][j][i] = monthly_height_bot[t][p][0][0]

landfrac_pft.assignValue(landfrac_pft_vals)
pftdata_mask.assignValue(pftdata_mask_vals)
longxy.assignValue(longxy_vals)
latixy.assignValue(latixy_vals)
area.assignValue(area_vals)
pct_wetland.assignValue(pct_wetland_vals)
pct_lake.assignValue(pct_lake_vals)
pct_glacier.assignValue(pct_glacier_vals)
pct_urban.assignValue(pct_urban_vals)
soil_color.assignValue(soil_color_vals)
fmax.assignValue(fmax_vals)
organic.assignValue(organic_vals)
pct_sand.assignValue(pct_sand_vals)
pct_clay.assignValue(pct_clay_vals)
pct_pft.assignValue(pct_pft_vals)
monthly_height_top.assignValue(monthly_height_top_vals)
monthly_height_bot.assignValue(monthly_height_bot_vals)
surffile_new_nc.close()


#-------------------- create pftdyn surface data ----------------------------------

if (options.compset == 'I20TRCLM45CN'):

    print('Creating dynpft data')

    pftdyn_orig = options.ccsm_input+options.grid_input  \
        +'/surfdata.pftdyn_0.5x0.5_simyr1850-2010.nc'
    pftdyn_new = options.ccsm_input+'/lnd/clm2/surfdata/' \
        +'surfdata.pftdyn_'+str(numxpts)+'x'+str(numypts)+'pt_'+options.site+'.nc'
    if (os.path.isfile(pftdyn_new)):
        print('Warning:  Removing existing pftdyn file')
        os.system('rm -rf '+pftdyn_new)
        os.system('ncks -d lsmlon,'+str(xgrid)+','+str(xgrid+numxpts-1)+' -d lsmlat,'+str(ygrid)+ \
          ','+str(ygrid+numypts-1)+' '+pftdyn_orig+' '+pftdyn_new)
        
    pftdyn_new_nc = NetCDF.NetCDFFile(pftdyn_new, "a")
    landfrac_pft = pftdyn_new_nc.variables['LANDFRAC_PFT']
    landfrac_pft_vals = landfrac_pft.getValue()
    pftdata_mask = pftdyn_new_nc.variables['PFTDATA_MASK']
    pftdata_mask_vals = pftdata_mask.getValue()
    longxy = pftdyn_new_nc.variables['LONGXY']
    longxy_vals = longxy.getValue()
    latixy = pftdyn_new_nc.variables['LATIXY']
    latixy_vals = latixy.getValue()
    area = pftdyn_new_nc.variables['AREA']
    area_vals = area.getValue()
    pct_wetland = pftdyn_new_nc.variables['PCT_WETLAND']
    pct_wetland_vals = pct_wetland.getValue()
    pct_lake = pftdyn_new_nc.variables['PCT_LAKE']
    pct_lake_vals = pct_lake.getValue()
    pct_glacier = pftdyn_new_nc.variables['PCT_GLACIER']
    pct_glacier_vals = pct_glacier.getValue()
    pct_urban = pftdyn_new_nc.variables['PCT_URBAN']
    pct_urban_vals = pct_urban.getValue()
    pct_pft = pftdyn_new_nc.variables['PCT_PFT']
    pct_pft_vals = pct_pft.getValue()
    grazing = pftdyn_new_nc.variables['GRAZING']
    grazing_vals = grazing.getValue()
    harvest_sh1 = pftdyn_new_nc.variables['HARVEST_SH1']
    harvest_sh1_vals = harvest_sh1.getValue()
    harvest_sh2 = pftdyn_new_nc.variables['HARVEST_SH2']
    harvest_sh2_vals = harvest_sh2.getValue()
    harvest_sh3 = pftdyn_new_nc.variables['HARVEST_SH3']
    harvest_sh3_vals = harvest_sh3.getValue()
    harvest_vh1 = pftdyn_new_nc.variables['HARVEST_VH1']
    harvest_vh1_vals = harvest_vh1.getValue()
    harvest_vh2 = pftdyn_new_nc.variables['HARVEST_VH2']
    harvest_vh2_vals = harvest_vh2.getValue()

    npft = 17

    #read file for site-specific PFT information
    AFdatareader = csv.reader(open(options.sitegroup+'_pftdata.txt','rb'))
    mypft_frac=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for row in AFdatareader:
        #print(row[0], row[1], options.site)
        if row[0] == options.site:
            for thispft in range(0,5):
                mypft_frac[int(row[2+2*thispft])]=float(row[1+2*thispft])

    if (os.path.exists(options.site+'_dynpftdata.txt')):
        dynexist = True
        DYdatareader = csv.reader(open(options.site+'_dynpftdata.txt','rb'))
        dim = (19,200)
        pftdata = numpy.zeros(dim)
        for row in DYdatareader:
            if row[0] == '1850':
                nrows=1
                for i in range(0,19):
                    pftdata[i][0] = float(row[i])
            elif row[0] != 'trans_year':
                nrows += 1
                for i in range(0,19):
                    pftdata[i][nrows-1] = float(row[i])
    else:
        dynexist = False
        print('Warning:  Dynamic pft file for site '+options.site+' does not exist')
        print('Using constant 1850 values')

    for i in range(0,numxpts):
        for j in range(0,numypts):
            landfrac_pft_vals[j][i] = 1.0
            pftdata_mask_vals[j][i] = 1
            longxy_vals[j][i] = lon+i*resx
            latixy_vals[j][i] = lat+j*resy
            area_vals[j][i] = 111.2*resy*111.321*math.cos((lon+i*resx)*math.pi/180)*resx
            pct_wetland_vals[j][i] = 0.0
            pct_lake_vals[j][i]    = 0.0
            pct_glacier_vals[j][i] = 0.0
            pct_urban_vals[j][i]   = 0.0
            thisrow = 0
            for t in range(0,161):     
                if (dynexist):
                    for p in range(0,npft):
                        pct_pft_vals[t][p][j][i] = 0.
                    harvest_thisyear = False
                    if pftdata[0][thisrow+1] == 1850+t:
                        thisrow = thisrow+1
                        harvest_thisyear = True
                    if (t == 0 or pftdata[16][thisrow] == 1):
                        harvest_thisyear = True
                    for k in range(0,5):
                        pct_pft_vals[t][int(pftdata[k*2+2][thisrow])][j][i] = \
                            pftdata[k*2+1][thisrow]
                        grazing_vals[t][j][i] = pftdata[17][thisrow]
                        if (harvest_thisyear):
                            harvest_sh1_vals[t][j][i] = pftdata[13][thisrow]
                            harvest_sh2_vals[t][j][i] = pftdata[14][thisrow]
                            harvest_sh3_vals[t][j][i] = pftdata[15][thisrow]
                            harvest_vh1_vals[t][j][i] = pftdata[11][thisrow]
                            harvest_vh2_vals[t][j][i] = pftdata[12][thisrow]
                        else:
                            harvest_sh1_vals[t][j][i] = 0.
                            harvest_sh2_vals[t][j][i] = 0.
                            harvest_sh3_vals[t][j][i] = 0.
                            harvest_vh1_vals[t][j][i] = 0.
                            harvest_vh2_vals[t][j][i] = 0.
                else:
                    for p in range(0,npft):
                        pct_pft_vals[t][p][j][i] = mypft_frac[p]
                        grazing_vals[t][j][i] = 0.
                        harvest_sh1_vals[t][j][i] = 0.
                        harvest_sh2_vals[t][j][i] = 0.
                        harvest_sh3_vals[t][j][i] = 0.
                        harvest_vh1_vals[t][j][i] = 0.
                        harvest_vh2_vals[t][j][i] = 0.

    landfrac_pft.assignValue(landfrac_pft_vals)
    pftdata_mask.assignValue(pftdata_mask_vals)
    longxy.assignValue(longxy_vals)
    latixy.assignValue(latixy_vals)
    area.assignValue(area_vals)
    pct_wetland.assignValue(pct_wetland_vals)
    pct_lake.assignValue(pct_lake_vals)
    pct_glacier.assignValue(pct_glacier_vals)
    pct_urban.assignValue(pct_urban_vals)
    pct_pft.assignValue(pct_pft_vals)
    grazing.assignValue(grazing_vals)
    harvest_sh1.assignValue(harvest_sh1_vals)
    harvest_sh2.assignValue(harvest_sh2_vals)
    harvest_sh3.assignValue(harvest_sh3_vals)
    harvest_vh1.assignValue(harvest_vh1_vals)
    harvest_vh2.assignValue(harvest_vh2_vals)

print "\n makepointdata.py successfully done! \n"
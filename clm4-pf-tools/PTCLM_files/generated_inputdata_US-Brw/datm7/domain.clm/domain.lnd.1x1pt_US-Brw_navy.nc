CDF       
      nj        ni        nv              Conventions       NCAR-CSM:CF-1.0    title         CCSM domain data:      user_comment      PStandard CCSM3.1/4.0 domain specification file created from CLM inputdata files:   source        $from CLM fraction and griddata files   history      Fri Jun  7 23:04:06 2013: ncks -d ni,407,407 -d nj,322,322 /Users/f9y/clm4_5_inputdata/ugrid/0.5x0.5data/domain.360x720_ORCHIDEE0to360.100409.nc /Users/f9y/clm4_5_inputdata/share/domains/domain.clm/domain.lnd.1x1pt_US-Brw_navy.nc

04/09/10 15:31:15 slevis:be1105en.ucar.ed       mkdatadomain_version      �$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/branch_tags/ccsm4_0_rel_tags/ccsm4_0_rel_03_clm3_7_10/models/lnd/clm/tools/mkdatadomain/addglobal.F90 $    SVN_url       �$HeadURL: https://svn-ccsm-models.cgd.ucar.edu/clm2/branch_tags/ccsm4_0_rel_tags/ccsm4_0_rel_03_clm3_7_10/models/lnd/clm/tools/mkdatadomain/addglobal.F90 $    mkdatadomain_version_Id       4$Id: addglobal.F90 13984 2009-01-20 05:54:15Z erik $   source_code       4$Id: addglobal.F90 13984 2009-01-20 05:54:15Z erik $   Land_Grid_Dataset         griddata_0360x0720.nc      Land_Fraction_Dataset         fracdata_0360x0720_ORCHIDEE.nc     NCO       4.3.1            area                   	long_name         $area of grid cell in radians squared   
coordinate        xc yc      units         radians2        
   frac                   	long_name         $fraction of grid cell that is active   
coordinate        xc yc      units         unitless   filter1       =error if frac> 1.0+eps or frac < 0.0-eps; eps = 0.1000000E-11      filter2       Jlimit frac to [fminval,fmaxval]; fminval= 0.1000000E-02 fmaxval=  1.000000          
   mask                   	long_name         land domain mask   
coordinate        xc yc      note      unitless   comment       70=ocean and 1=land, 0 indicates that cell is not active         
   xc                     	long_name         longitude of grid cell center      units         degrees_east   bounds        xv          
   xv                        	long_name         longitude of grid cell vertices    units         degrees_east         
    yc                     	long_name         latitude of grid cell center   units         degrees_north      bounds        yv          
@   yv                        	long_name         latitude of grid cell vertices     units         degrees_north            
H>ɍ���[�?�         @ik���-�@ij_��F@im�:)�z@ij_��F@im�:)�z@Q���
=q@Q�p��
>@Q�p��
>@Q��
=p�@Q��
=p�
CLM-PFLOTRAN Regression test suite
==================================

The regression test suite attempts to do some basic testing of clm-pflotran, including:

1. building the executable
2. ensuring all test problems run to completion without generating standard cesm error messages
3. basic comparison to a baseline



Running the tests
-----------------

1. Download the three repos: (i)PFLOTRAN repo; (ii)CLM repo; and (iii) Repo containing inputdata for regression tests


    ::
    
        export RTEST_DIR=$PWD
        hg clone https://bitbucket.org/clm_pflotran/pflotran-clm-trunk
        hg clone https://bitbucket.org/clm_pflotran/clm-pflotran-trunk
        hg clone https://bitbucket.org/clm_pflotran/clm-pflotran-data-trunk-testing
    
2. Test pflotran regression tests

    ::
    
        cd ${RTEST_DIR}/pflotran-clm-trunk/src/pflotran
        make test
    
3. Build libpflotran.a

    ::
    
        cd ${RTEST_DIR}/pflotran-clm-trunk/src/clm-pflotran
        ./link_files.sh; 
        make libpflotran.a
    
4. Create a configuration file for your machine and fill in ALL the fields:


    ::
    
        mkdir ~/.cesm
        hostname -s
        cat > ~/.cesm/clm-pflotran-machines.cfg <<EOF
        [XXX_MY_MACHINE_NAME_XXX]
        machine = userdefined
        compiler = gnu
        mpi_vendor = 
        max_np = 
        os = 
        GMAKE =
        CC =
        CXX =
        FC =
        MPICC =
        MPICXX =
        MPIFC =
        MPIEXEC =
        FFLAGS =
        SLIBS =
        NETCDF_PATH =
        BLAS_FLAGS =
        ldflags = 
        EOF


    Note: Replace `XXX_MY_MACHINE_NAME_XXX` with the output of `hostname -s`.


5. Update to 'clm_pflotran' branch

    ::
    
        cd $RTEST_DIR/clm-pflotran-trunk
        hg update clm_pflotran
    
6. Create a local.cfg file

    ::
    
        cd ${RTEST_DIR}/clm-pflotran-trunk/clm4-pf-tools/regression_tests
        cat > local.cfg <<EOF
        [petsc]
        petsc_dir = $PETSC_DIR
        petsc_arch = $PETSC_ARCH
        
        [pflotran]
        pflotran_dir = ${RTEST_DIR}/pflotran-clm-trunk
        
        [data]
        data_dir = ${RTEST_DIR}/clm-pflotran-data-trunk-testing
        EOF
        
7. Run the tests

    ::
    
        make test


Additional information:
- Regression test suite creates new CESM cases in $RTEST_DIR/clm-pflotran-trunk/test_cases

- By default, only the $RTEST_DIR/clm-pflotran-trunk/test_cases/common-executable case is build. All other tests reuse the executable: $RTEST_DIR/clm-pflotran-trunk/test_cases/common-executable/bld/cesm.exe

- For each CESM regression case, a corresponding shell script is created in $RTEST_DIR/clm-pflotran-trunk/test_cases/.

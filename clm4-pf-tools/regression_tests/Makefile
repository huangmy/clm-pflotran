# -*- mode: makefile; -*-

# Makefile for manually running clm-pflotran tests suite

PYTHON=python2.7
ifdef python
    PYTHON=$(python)
endif

MACHINES=machines.cfg
ifdef machines
    MACHINES=$(machines)
endif

LOCAL_CFG=local.cfg
ifdef local
    LOCAL_CFG=$(local)
endif

OPTS=
ifdef debug
    OPTS:=$(OPTS) --backtrace
endif

CESM_ROOT=../../

STANDARD_TESTS= \
    tests/sgrid-1x1x10-subsurface-th-noice-may.cfg \
    tests/ugrid-1x1x10-subsurface-th-noice-may.cfg \
    tests/ugrid-1x1x10-subsurface-th-noice-dec.cfg \
    tests/ugrid-1x1x10-subsurface-richards-may.cfg \
    tests/ugrid-13x26x10-subsurface-th-noice-dec-NGEE_SiteB-np-1.cfg \
    tests/ugrid-13x26x10-surface-subsurface-th-noice-dec-NGEE_SiteB.cfg
#    tests/ugrid-13x26x10-subsurface-th-noice-dec-NGEE_SiteB-np-4.cfg \
#    tests/ugrid-1x1x10-subsurface-th-ice-may.cfg \

COMMON_EXE= \
    tests/common-executable.cfg 

TEST_CMD=$(PYTHON) ./clm-pflotran-tests.py $(OPTS) \
    --cesm-root-dir $(CESM_ROOT) \
    --local-config-file $(LOCAL_CFG) \
    --machines-config-file $(MACHINES) \

test : FORCE
	$(TEST_CMD) --test-config-files $(COMMON_EXE) $(STANDARD_TESTS)

# assume the common executable has already been built
fast : FORCE
	$(TEST_CMD) --test-config-files $(STANDARD_TESTS)

tests/common-executable.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/sgrid-1x1x10-subsurface-th-noice-may.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-1x1x10-subsurface-th-noice-may.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-1x1x10-subsurface-th-ice-may.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-1x1x10-subsurface-th-noice-dec.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-1x1x10-subsurface-richards-may.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-13x26x10-subsurface-th-noice-dec-NGEE_SiteB-np-4.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-13x26x10-subsurface-th-noice-dec-NGEE_SiteB-np-1.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

tests/ugrid-13x26x10-surface-subsurface-th-noice-dec-NGEE_SiteB.cfg : FORCE
	$(TEST_CMD) --test-config-files $@

FORCE :


clean :
	-rm -f *.testlog *~

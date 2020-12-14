

PLTFMOPT=Makefile.$(shell uname)
include $(PLTFMOPT)

#                  +------- link to where EPICS is on this machine
#                  |        see below
EPICS_ROOT="curEPICS"

# at ANL we use
#rm curEPICS; ln -s /global/devel6_sl62/base/base-3.14.12.1 curEPICS
#or maybe
#rm curEPICS; ln -s /global/devel6_sl62/base/R3.14.10 curEPICS

# at gretina @ MSU we use
#rm curEPICS; ln -s /global/devel/base/R3.14.11 curEPICS

EPICS_LIBSUBDIR=linux-x86
EPICS_INCSUBDIR=Linux

#------------------------------------------ end modifications ---

EPICSLIBDIR=$(EPICS_ROOT)/lib
EPICSARCHLIBDIR=$(EPICSLIBDIR)/$(EPICS_LIBSUBDIR)
EPICSINCDIR=$(EPICS_ROOT)/include
EPICSOSINCDIR=$(EPICS_ROOT)/include/os/$(EPICS_INCSUBDIR)

CFLAGS=-I$(EPICSINCDIR) -I$(EPICSOSINCDIR) -I.

EPICSFLAGS=-L$(EPICSARCHLIBDIR) -lCom -Wl,"-rpath=$(EPICSARCHLIBDIR)"

#-------------------

all:	GEBMerge EventBuilder EventBuilder_trace

#------------------- GEBMegre

time_stamp.o:	time_stamp.c
		$(cc)  -c  $^

GEBMerge.c:	gdecomp.h GTMerge.h GTMerge_readnew.h 

GEBMerge:	GEBMerge.c time_stamp.o
		$(cc) -o $@ -lm -lz $^

#------------------- EventBuilder, the GEBSort

GTPrintCC.o:	GTPrint.c
		$(CC) -c -o $@  $^

GEBSort.cxx:	GTMerge.h GEBSort.h gdecomp.h GEBClient.h gretTapClient.h

GEBSort_nogeb.o:	GEBSort.cxx 
		$(CC) -DHAVE_VXWORKS=0  $(ROOTINC)  $(CFLAGS) -c GEBSort.cxx -o GEBSort_nogeb.o

time_stampCC.o:	time_stamp.c
		$(CC)  -c  -o time_stampCC.o time_stamp.c  

GEBHeaderCC.o:	GEBHeader.c
		$(cc) -c -o $@ $^

jta.o:	jta.c
		$(CC) -c -o $@  $^

bin_XA.o:	bin_XA.c
		$(CC) -c -o $@ $^

#rename GEBSort_nogeb to be EventBuilder
EventBuilder:	GTPrintCC.o GEBSort_nogeb.o time_stampCC.o GEBHeaderCC.o jta.o bin_XA.o 
		$(CC) $^ -DHAVE_VXWORKS=0 $(ROOTLIB)  -o $@ 


#------------------------------

bin_XA_trace.o:	bin_XA_trace.c
		$(CC) -c -o $@ $^

#rename GEBSort_nogeb to be EventBuilder
EventBuilder_trace:	GTPrintCC.o GEBSort_nogeb.o time_stampCC.o GEBHeaderCC.o jta.o bin_XA_trace.o 
		$(CC) $^ -DHAVE_VXWORKS=0 $(ROOTLIB)  -o $@ 


#------------------------------


clean:		
		rm -f *.o
		rm -f core.*
		rm -f *%
		rm -f a.out
		rm -f *.gch
		rm EventBuilder GEBMerge

#------------------------------


#GEBHeader.c:	GTMerge.h gdecomp.h

#GTPrint.c:	GTMerge.h gdecomp.h GEBSort.h



#str_decompCC.o:	str_decomp.c
#		$(CC)  -c  -o $@ $^

#utilsCC.o:	utils.c
#		$(CC)  -c  -o $@ $^

#util.o:		utils.c
#		$(cc)  -c  -o $@ $^

#str_decomp.c:	ctk.h

#GEBClient.o:	GEBClient.c 
#		$(CC) $(EPICSFLAGS) $(CFLAGS) -c -o $@  $^

#GEBClient.c:	GEBClient.h

#gretTapClient.o:	gretTapClient.c 
#			$(CC) $(EPICSFLAGS) $(CFLAGS) -c -o $@  $^

#gretTapClient.c:	GEBClient.h gretTapClient.h


#-------------------


#printEvent.c:	ctk.h

#-------------------

#GEBFilter:	GEBFilter.c GF_veto_cube.o get_a_seed.o
#		$(cc) -o $@ -lm -lz $^
#		$@  GEBFilter.chat GTDATA/merged.gtd_000 filtertest.gtd

#GEBFilter.c:	ctk.h gdecomp.h GEBFilter.h


#-------------------

#listTS:		listTS.c gdecomp.h
#		$(cc) -o $@ $^

#spe_funCC.o:	spe_fun.c
#		$(CC) -c -o spe_funCC.o  $^

#G4toMode2.c:	gdecomp.h

#-------------------


#debug:		
#		gdb -d ./ G4toMode2_EU `ls core.*`
#		gdb -d ./ GEBSort_nogeb `ls core.*`
#		gdb -d ./ GEBMerge `ls core.*`
#		gdb -d ./ GEBbunch `ls core.*`
#		gdb -d ./ AG_mode3_to_GT `ls core.*`
#		gdb -d ./ dgs_ecal2 `ls core.*`



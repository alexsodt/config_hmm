#CC = /usr/bin/g++
#F77 = /usr/bin/gfortran

CC = mpic++
F77 = mpif77

LDFLAGS = 

COPT = -g 

.C.o: 
	$(CC) $(GSLINC) \
		 -o $*.o -c $(COPT)  -std=gnu++11 $*.C

.c.o:
	gcc  $(QHULL)	$(COPT) $(CDEFN)	-c	$*.c

.f.o:
	$(F77)  -c $*.f -o $*.o

EXECS = cluster extractDimers extractDimerHMM printFullDimers HMMAssignment sortindex discretizedFourier.exe fourierExtract.exe make_dbg_data EncodeHMMLipidAllAtom curvFromModeRange multistateStats.exe cholHMM.exe estimatedFdR.exe
 
all: $(EXECS)

BLASOBJS = \
BLAS/caxpy.o  \
BLAS/ccopy.o  \
BLAS/cdotc.o  \
BLAS/cdotu.o  \
BLAS/cgbmv.o  \
BLAS/cgemm.o  \
BLAS/cgemv.o  \
BLAS/cgerc.o  \
BLAS/cgeru.o  \
BLAS/chbmv.o  \
BLAS/chemm.o  \
BLAS/chemv.o  \
BLAS/cher.o  \
BLAS/cher2.o  \
BLAS/cher2k.o  \
BLAS/cherk.o  \
BLAS/chpmv.o  \
BLAS/chpr.o  \
BLAS/chpr2.o  \
BLAS/crotg.o  \
BLAS/cscal.o  \
BLAS/csrot.o  \
BLAS/csscal.o  \
BLAS/cswap.o  \
BLAS/csymm.o  \
BLAS/csyr2k.o  \
BLAS/csyrk.o  \
BLAS/ctbmv.o  \
BLAS/ctbsv.o  \
BLAS/ctpmv.o  \
BLAS/ctpsv.o  \
BLAS/ctrmm.o  \
BLAS/ctrmv.o  \
BLAS/ctrsm.o  \
BLAS/ctrsv.o  \
BLAS/dasum.o  \
BLAS/daxpy.o  \
BLAS/dcabs1.o  \
BLAS/dcopy.o  \
BLAS/ddot.o  \
BLAS/dgbmv.o  \
BLAS/dgemm.o  \
BLAS/dgemv.o  \
BLAS/dger.o  \
BLAS/dnrm2.o  \
BLAS/drot.o  \
BLAS/drotg.o  \
BLAS/drotm.o  \
BLAS/drotmg.o  \
BLAS/dsbmv.o  \
BLAS/dscal.o  \
BLAS/dsdot.o  \
BLAS/dspmv.o  \
BLAS/dspr.o  \
BLAS/dspr2.o  \
BLAS/dswap.o  \
BLAS/dsymm.o  \
BLAS/dsymv.o  \
BLAS/dsyr.o  \
BLAS/dsyr2.o  \
BLAS/dsyr2k.o  \
BLAS/dsyrk.o  \
BLAS/dtbmv.o  \
BLAS/dtbsv.o  \
BLAS/dtpmv.o  \
BLAS/dtpsv.o  \
BLAS/dtrmm.o  \
BLAS/dtrmv.o  \
BLAS/dtrsm.o  \
BLAS/dtrsv.o  \
BLAS/dzasum.o  \
BLAS/dznrm2.o  \
BLAS/icamax.o  \
BLAS/idamax.o  \
BLAS/isamax.o  \
BLAS/izamax.o  \
BLAS/lsame.o  \
BLAS/sasum.o  \
BLAS/saxpy.o  \
BLAS/scabs1.o  \
BLAS/scasum.o  \
BLAS/scnrm2.o  \
BLAS/scopy.o  \
BLAS/sdot.o  \
BLAS/sdsdot.o  \
BLAS/sgbmv.o  \
BLAS/sgemm.o  \
BLAS/sgemv.o  \
BLAS/sger.o  \
BLAS/snrm2.o  \
BLAS/srot.o  \
BLAS/srotg.o  \
BLAS/srotm.o  \
BLAS/srotmg.o  \
BLAS/ssbmv.o  \
BLAS/sscal.o  \
BLAS/sspmv.o  \
BLAS/sspr.o  \
BLAS/sspr2.o  \
BLAS/sswap.o  \
BLAS/ssymm.o  \
BLAS/ssymv.o  \
BLAS/ssyr.o  \
BLAS/ssyr2.o  \
BLAS/ssyr2k.o  \
BLAS/ssyrk.o  \
BLAS/stbmv.o  \
BLAS/stbsv.o  \
BLAS/stpmv.o  \
BLAS/stpsv.o  \
BLAS/strmm.o  \
BLAS/strmv.o  \
BLAS/strsm.o  \
BLAS/strsv.o  \
BLAS/xerbla.o  \
BLAS/zaxpy.o  \
BLAS/zcopy.o  \
BLAS/zdotc.o  \
BLAS/zdotu.o  \
BLAS/zdrot.o  \
BLAS/zdscal.o  \
BLAS/zgbmv.o  \
BLAS/zgemm.o  \
BLAS/zgemv.o  \
BLAS/zgerc.o  \
BLAS/zgeru.o  \
BLAS/zhbmv.o  \
BLAS/zhemm.o  \
BLAS/zhemv.o  \
BLAS/zher.o  \
BLAS/zher2.o  \
BLAS/zher2k.o  \
BLAS/zherk.o  \
BLAS/zhpmv.o  \
BLAS/zhpr.o  \
BLAS/zhpr2.o  \
BLAS/zrotg.o  \
BLAS/zscal.o  \
BLAS/zswap.o  \
BLAS/zsymm.o  \
BLAS/zsyr2k.o  \
BLAS/zsyrk.o  \
BLAS/ztbmv.o  \
BLAS/ztbsv.o  \
BLAS/ztpmv.o  \
BLAS/ztpsv.o  \
BLAS/ztrmm.o  \
BLAS/ztrmv.o  \
BLAS/ztrsm.o  \
BLAS/ztrsv.o  

LAPACKOBJS = \
lapack/double/dgebak.o \
lapack/double/dgebal.o \
lapack/double/dgeev.o \
lapack/double/dgehd2.o \
lapack/double/dgehrd.o \
lapack/double/dgesv.o \
lapack/double/dgetf2.o \
lapack/double/dgetrf.o \
lapack/double/dgetri.o \
lapack/double/dgetrs.o \
lapack/double/dhseqr.o \
lapack/double/disnan.o \
lapack/double/dlabad.o \
lapack/double/dlacpy.o \
lapack/double/dladiv.o \
lapack/double/dlae2.o \
lapack/double/dlaev2.o \
lapack/double/dlaexc.o \
lapack/double/dlahqr.o \
lapack/double/dlahr2.o \
lapack/double/dlaisnan.o \
lapack/double/dlaln2.o \
lapack/double/dlange.o \
lapack/double/dlanst.o \
lapack/double/dlansy.o \
lapack/double/dlanv2.o \
lapack/double/dlapy2.o \
lapack/double/dlaqr0.o \
lapack/double/dlaqr1.o \
lapack/double/dlaqr2.o \
lapack/double/dlaqr3.o \
lapack/double/dlaqr4.o \
lapack/double/dlaqr5.o \
lapack/double/dlarfb.o \
lapack/double/dlarf.o \
lapack/double/dlarfg.o \
lapack/double/dlarft.o \
lapack/double/dlarfx.o \
lapack/double/dlartg.o \
lapack/double/dlascl.o \
lapack/double/dlaset.o \
lapack/double/dlasr.o \
lapack/double/dlasrt.o \
lapack/double/dlassq.o \
lapack/double/dlaswp.o \
lapack/double/dlasy2.o \
lapack/double/dlasyf.o \
lapack/double/dlatrd.o \
lapack/double/dorg2l.o \
lapack/double/dorg2r.o \
lapack/double/dorghr.o \
lapack/double/dorgql.o \
lapack/double/dorgqr.o \
lapack/double/dorgtr.o \
lapack/double/dorm2r.o \
lapack/double/dormhr.o \
lapack/double/dormqr.o \
lapack/double/dsteqr.o \
lapack/double/dsterf.o \
lapack/double/dsyev.o \
lapack/double/dsytd2.o \
lapack/double/dsytf2.o \
lapack/double/dsytrd.o \
lapack/double/dsytrf.o \
lapack/double/dsytri.o \
lapack/double/dtrevc.o \
lapack/double/dtrexc.o \
lapack/double/dtrti2.o \
lapack/double/dtrtri.o \
lapack/util/dlamch.o \
lapack/util/ieeeck.o \
lapack/util/iladlc.o \
lapack/util/iladlr.o \
lapack/util/ilaenv.o \
lapack/util/iparmq.o \
lapack/util/lsame.o \
lapack/util/xerbla.o

clean: 
	rm  *.o $(LAPACKOBJS) $(BLASOBJS) $(EXECS)

blas.a: $(BLASOBJS)
	ar -r blas.a $(BLASOBJS)

lapack.a: $(LAPACKOBJS)
	ar -r lapack.a $(LAPACKOBJS)

extractDimers: extractDimers.o util.o dcd.o alignSet.o pdb.o geometry.o comparison.o lapack.a blas.a proc_definition.o
	$(CC) $(LDFLAGS)  -o extractDimers extractDimers.o util.o dcd.o alignSet.o pdb.o geometry.o comparison.o proc_definition.o lapack.a blas.a -lboost_regex -lstdc++ -lgfortran

extractDimerHMM: extractDimerHMM.o util.o dcd.o alignSet.o pdb.o geometry.o lapack.a blas.a medoidLR.o comparison.o proc_definition.o
	$(CC) $(LDFLAGS)  -o extractDimerHMM extractDimerHMM.o util.o dcd.o alignSet.o pdb.o geometry.o medoidLR.o comparison.o proc_definition.o -lboost_regex lapack.a blas.a -lstdc++ -lgfortran

printFullDimers: printFullDimers.o util.o dcd.o alignSet.o pdb.o geometry.o lapack.a blas.a comparison.o proc_definition.o
	$(CC) $(LDFLAGS)  -o printFullDimers printFullDimers.o util.o dcd.o alignSet.o pdb.o geometry.o proc_definition.o lapack.a blas.a -lboost_regex -lstdc++ -lgfortran

HMMAssignment: HMMAssignment.o util.o dcd.o alignSet.o pdb.o geometry.o lapack.a blas.a medoidLR.o comparison.o proc_definition.o
	$(CC) $(LDFLAGS)  -o HMMAssignment HMMAssignment.o util.o dcd.o alignSet.o pdb.o geometry.o medoidLR.o comparison.o proc_definition.o lapack.a blas.a -lboost_regex -lstdc++ -lgfortran

cluster: cluster.o medoidLR.o util.o dcd.o alignSet.o pdb.o geometry.o lapack.a blas.a comparison.o proc_definition.o
	$(CC) -o cluster cluster.o medoidLR.o util.o dcd.o alignSet.o pdb.o geometry.o comparison.o lapack.a proc_definition.o blas.a -lboost_regex -lstdc++ -lgfortran

sortindex: sortindex.o
	$(CC) -o sortindex sortindex.o

discretizedFourier.exe: discretizedFourier.o util.o   alignSet.o dcd.o pdb.o blas.a lapack.a
	$(CC) -o discretizedFourier.exe discretizedFourier.o util.o    alignSet.o dcd.o pdb.o  $(LDFLAGS) $(GSLLIB) lapack.a blas.a -lstdc++ -lgfortran

fourierExtract.exe: fourierExtract.o util.o  alignSet.o dcd.o pdb.o blas.a lapack.a proc_definition.o
	$(CC) -o fourierExtract.exe fourierExtract.o util.o     alignSet.o dcd.o pdb.o  $(LDFLAGS) $(GSLLIB) proc_definition.o lapack.a blas.a -lboost_regex -lstdc++ -lgfortran

make_dbg_data: make_dbg_data.o util.o pdb.o dcd.o alignSet.o lapack.a blas.a
	$(CC) -o make_dbg_data make_dbg_data.o util.o pdb.o dcd.o alignSet.o lapack.a blas.a -lgfortran

EncodeHMMLipidAllAtom: EncodeHMMLipidAllAtom.o util.o dcd.o alignSet.o pdb.o lapack.a blas.a 
	$(CC) $(LDFLAGS) -o EncodeHMMLipidAllAtom EncodeHMMLipidAllAtom.o util.o dcd.o alignSet.o pdb.o lapack.a blas.a -lgfortran 

curvFromModeRange: curvFromModeRange.o util.o
	$(CC) -o curvFromModeRange curvFromModeRange.o util.o

multistateStats.exe: multistateStats.o util.o  alignSet.o dcd.o pdb.o blas.a lapack.a proc_definition.o
	$(CC) -o multistateStats.exe multistateStats.o util.o     alignSet.o dcd.o pdb.o  $(LDFLAGS) $(GSLLIB) proc_definition.o lapack.a blas.a -lboost_regex -lstdc++ -lgfortran

cholHMM.exe: cholHMM.o util.o
	$(CC) -o cholHMM.exe cholHMM.o util.o

estimatedFdR.exe: estimatedFdR.o util.o
	$(CC) -o estimatedFdR.exe estimatedFdR.o util.o

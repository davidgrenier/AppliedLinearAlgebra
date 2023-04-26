WARN = -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wundef -Werror -Wno-unused-variable
EXTRA = -Wold-style-cast -Wctor-dtor-privacy -Wstrict-overflow=5 
WARNPLUS = $(WARN) $(EXTRA) -Wlogical-op -Wstrict-null-sentinel -Wsign-conversion# -Wnoexcept
BOOST = -lboost_iostreams -lboost_system -lboost_filesystem
GNUPLOT = -isystem/home/david/projects/gnuplot-iostream
EIGEN = -isystem/usr/local/include/eigen3 #-lopenblas -fopenmp # compiler benchmark below exclude BLAS
GPP = g++
# Eigen householderPerf(g++), Naive householderPerf(g++), Compile time clang(g++)
# OPTI = -O0                                  # O: 360.824(369.177)ms, C: 540.425(517.698)ms, Co: 3.055(4.418)s
# OPTI = -O0 -mavx2 -mfma -ffp-contract=fast  # O: 186.659(173.879)ms, C: 283.883(283.372)ms, Co: 3.352(4.837)s
# OPTI = -O2							      # O: 6.715(6.586)ms, C: 22.034(22.527)ms, Co: 4.429(4.903)s
# OPTI = -O3								  # O: 6.554(6.586)ms, C: 25.107(14.345)ms, Co: 4.463(5.254)s
OPTI = -O1 -mavx2 -mfma -ffp-contract=fast  # O: 5.878(5.561)ms, C: 15.512(15.414)ms, Co: 4.393(4.800)s
# OPTI = -O2 -mavx2 -mfma -ffp-contract=fast  # O: 5.100(5.585)ms, C: 15.386(16.338)ms, Co: 4.903(5.569)s
# OPTI = -O3 -mavx2 -mfma -ffp-contract=fast  # O: 5.191(5.289)ms, C: 16.162(13.212)ms, Co: 4.985(6.021)s
CLANG = c++ $(EXTRA)
COMPILE = $(CLANG) -std=c++20 $(OPTI) $(WARN) $(EIGEN) $(GNUPLOT) $(BOOST)

scratch: scratch.cpp
	$(COMPILE) $^ -o $@

leastSquarePoly: leastSquarePoly.cpp
	$(COMPILE) $^ -o $@

householderPerf: householderPerf.cpp
	$(COMPILE) $^ -o $@

householderTest: householderTest.cpp
	$(COMPILE) $^ -o $@

deterioration: deterioration.cpp
	$(COMPILE) $^ -lX11 -o $@

legenApprox: legenApprox.cpp
	$(COMPILE) $^ -o $@

legendreQR: legendreQR.cpp
	$(COMPILE) $^ -o $@

patho2: patho2.cpp
	$(COMPILE) $^ -o $@

pathological: pathological.cpp
	$(COMPILE) $^ -o $@

gramSchmidt: gramSchmidt.cpp
	$(COMPILE) $^ -o $@

eigenDecHermit: eigenDecHermit.cpp
	$(COMPILE) $^ -o $@

svdPolarDecomp: svdPolarDecomp.cpp
	$(COMPILE) $^ -o $@

singular2x2: singular2x2.cpp
	$(COMPILE) $^ -o $@

hermitSvd: hermitSvd.cpp
	$(COMPILE) $^ -o $@

svd: svd.cpp
	$(COMPILE) $^ -o $@

frob: frob.cpp
	$(COMPILE) $^ -o $@

vandermonde: vandermonde.cpp
	$(COMPILE) $^ -o $@

norms: norms.cpp
	$(COMPILE) $^ -o $@

perturb: perturb.cpp
	$(COMPILE) $^ -o $@

lesson1: lesson1.cpp
	$(COMPILE) $^ -o $@

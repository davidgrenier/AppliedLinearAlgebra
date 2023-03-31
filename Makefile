WARN = -pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization -Wformat=2 -Winit-self -Wmissing-include-dirs -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wundef -Werror -Wno-unused-variable
EXTRA = -Wold-style-cast -Wsign-conversion -Wctor-dtor-privacy -Wstrict-overflow=5 
WARNPLUS = $(WARN) $(EXTRA) -Wlogical-op -Wstrict-null-sentinel# -Wnoexcept
BOOST = -lboost_iostreams -lboost_system -lboost_filesystem
GNUPLOT = -isystem/home/david/projects/gnuplot-iostream
EIGEN = -isystem/usr/local/include/eigen3 -lopenblas -fopenmp
GPP = g++ $(WARN) -std=c++20 -O3
CWARN = -Wall
CLANG = c++ $(WARN) $(EXTRA) -std=c++20 -O3
OPTI = $(CLANG) $(EIGEN) $(GNUPLOT) $(BOOST)

scratch: scratch.cpp
	$(OPTI) $^ -o $@

eigenDecHermit: eigenDecHermit.cpp
	$(OPTI) $^ -o $@

svdPolarDecomp: svdPolarDecomp.cpp
	$(OPTI) $^ -o $@

singular2x2: singular2x2.cpp
	$(OPTI) $^ -o $@

hermitSvd: hermitSvd.cpp
	$(OPTI) $^ -o $@

svd: svd.cpp
	$(OPTI) $^ -o $@

frob: frob.cpp
	$(OPTI) $^ -o $@

vandermonde: vandermonde.cpp
	$(OPTI) $^ -o $@

norms: norms.cpp
	$(OPTI) $^ -o $@

perturb: perturb.cpp
	$(OPTI) $^ -o $@

lesson1: lesson1.cpp
	$(OPTI) $^ -o $@
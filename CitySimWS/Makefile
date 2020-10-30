# ===========================================
# Project : CitySimWebS
# Autor : Sauvain Renaud
# Date : 20.03.2011
# ===========================================

TARGET = CitySimWebS
SRCDIR = src/

SOURCES = $(TARGET).cpp soapC.cpp soapServer.cpp 
OUTDIR = out/
FLAGS = -fPIC -lGLU -pipe -Wall -fopenmp -w -DWITH_GZIP  -DTIXML_USE_STL

INCLUDE = -I./src -I.. -I../NR -I../Sky -I../Tinyxml -I../VFCLibrary
LIBRARY = -L./src -L.. -L../NR -L../Sky -L../Tinyxml -L../VFCLibrary -lCitySim
GSOAPFLAGS = src/stdsoap2.cpp
OBJECTS = $(patsubst %.cpp,$(SRCDIR)%.cpp,$(SOURCES))

# == default ===========================================
all: 
	$(CXX) $(INCLUDE) -ggdb -lz -o $(OUTDIR)$(TARGET) $(OBJECTS) $(GSOAPFLAGS) $(FLAGS) $(LIBRARY)

# == clean ==========================================
clean:
	#rm -Rf $(OBJDIR) $(OUTDIR)

# == mkdir ==========================================
"%/":
	@[ -z "$@" -o -d "$@" ] || mkdir -p $@ && chmod g+w $@

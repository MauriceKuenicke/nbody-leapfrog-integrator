#
# Makefile for nbody:
#
#    make                   : compile program
#    make install           : install program in bin directory
#    make clean             : remove temporary files
#    make cleaner           : remove everything but source
#

# source and bin directory
SRCDIR     := src
BINDIR     := bin

# name of the executable
EXECUTABLE := nbody

all:
	cd $(SRCDIR); $(MAKE) $(EXECUTABLE)


install: all
	cp $(SRCDIR)/$(EXECUTABLE) $(BINDIR)

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*~

cleaner: clean
	rm -f $(BINDIR)/$(EXECUTABLE) $(SRCDIR)/$(EXECUTABLE)


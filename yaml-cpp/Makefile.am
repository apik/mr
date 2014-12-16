## Process this file with automake to produce Makefile.in .
#
# Author: Paul_Kunz@slac.stanford.edu
#

# The following is set; otherwise it would follow GNU conventions.
# AUTOMAKE_OPTIONS = foreign

# The convenience library to be built.
lib_LIBRARIES = libYAMLpp.a

libYAMLpp_a_SOURCES =				\
	src/binary.cpp	        src/emit.cpp	        src/emitterstate.cpp    src/memory.cpp       \
	src/node_data.cpp       src/ostream_wrapper.cpp src/regex.cpp           src/scantag.cpp      \
	src/singledocparser.cpp	src/convert.cpp         src/emitfromevents.cpp  src/emitterutils.cpp \
	src/nodebuilder.cpp     src/nodeevents.cpp      src/parse.cpp	        src/scanner.cpp	     \
	src/scantoken.cpp       src/stream.cpp          src/directives.cpp      src/emitter.cpp	     \
	src/exp.cpp		  src/node.cpp	       src/null.cpp	   src/parser.cpp	    src/scanscalar.cpp	src/simplekey.cpp  src/tag.cpp


libYAMLpp_a_CPPFLAGS =                          \
	-I$(top_srcdir) -Iinclude


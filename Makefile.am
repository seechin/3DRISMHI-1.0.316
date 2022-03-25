bin_PROGRAMS = gensolvent gmxtop2solute heatmap ts4sdump rismhi3d

gensolvent_SOURCES = src/gensolvent.cpp
gmxtop2solute_SOURCES = src/gmxtop2solute.cpp
heatmap_SOURCES = src/heatmap.cpp
ts4sdump_SOURCES = src/ts4sdump.cpp
rismhi3d_SOURCES = src/rismhi3d.cpp

LDADD = $(EXTRA_LDFLAGS) $(LDFLAGS)
rismhi3d_CXXFLAGS = $(CFLAGS) $(EXTRA_CFLAGS)
gensolvent_CXXFLAGS = $(CFLAGS) $(EXTRA_CFLAGS)

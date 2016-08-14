#g++ SimpleXyRd.cpp `ncxx4-config  --libs` -Wall -o SimpleXyRd
#g++ findBadBC.cpp `ncxx4-config  --libs` -Wall -o findBadBC

#gcc simple_xy_rd.c -Wall `nc-config --cflags --libs` -o simple_xy_rd
gcc findBadBC.c -Wall -g `nc-config --cflags --libs` -o findBadBC

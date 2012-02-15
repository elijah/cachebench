set logscale x 2
set grid
set term postscript
set title "Memory Hierarchy Performance"
set xlabel "Vector Length"
set ylabel "MB/Sec"
plot "llcbench_cache_read.dat" title 'read' with linespoints \
, "llcbench_cache_write.dat" title 'write' with linespoints \
, "llcbench_cache_rmw.dat" title 'rmw' with linespoints \
, "llcbench_cache_handread.dat" title 'handread' with linespoints \
, "llcbench_cache_handwrite.dat" title 'handwrite' with linespoints \
, "llcbench_cache_handrmw.dat" title 'handrmw' with linespoints \
, "llcbench_cache_memset.dat" title 'memset' with linespoints \
, "llcbench_cache_memcpy.dat" title 'memcpy' with linespoints \

from stability import *

rs = get_inner_at_target(200,4,1e-10,1e-3)

print ("P(RS[{},{},{}]) = {:1.3e}".format(rs.n,rs.k,rs.d,rs.P_result()))

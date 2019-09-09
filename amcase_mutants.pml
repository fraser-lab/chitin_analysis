bg white
viewport 1500,1500
fetch 3rm4
color grey, 3rm4
remove resn HOH
remove chain B
color red, resn 3RM
color mendelevium, resi 239+364
color tv_orange, resi 246
show sticks, resi 239+364+246 and not name N+C+O
show sticks, resi 136+138+140 and not name N+C+O
color aquamarine, resi 136+138+140 and not name N+C+O
set_view (\
     0.534795165,   -0.734789670,    0.417223930,\
     0.550100207,    0.677555561,    0.488162488,\
    -0.641388476,   -0.031554054,    0.766564131,\
    -0.000081837,    0.000013820, -183.666030884,\
     0.435390949,    1.959021568,  -16.750881195,\
  -512.875793457,  880.209716797,  -20.000000000 )

set ambient_occlusion_mode, 1
set shininess, 10
set reflect, 0.24
set cartoon_side_chain_helper, 1
